"""Substituent prefix names for multiplicative names (P-29.4): simple groups
come from a table; anything else is named by attaching the substituent to a
benzene carrier (with the principal group of the parent para to it, so
senior-class groups inside the prefix are demoted correctly) and stripping
the carrier from the resulting name.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, halogen_prefix, halogen_substituents
from ._multiplicative_groups import KETONE_LIKE

SUFFIX_CARRIERS = {
    "carboxylic_acid": ("C(=O)O", "benzoic acid"),
    "sulfonic_acid": ("S(=O)(=O)O", "benzene-1-sulfonic acid"),
    "amide": ("C(N)=O", "benzamide"),
    "hydrazide": ("C(=O)NN", "benzohydrazide"),
    "amidine": ("C(N)=N", "benzene-1-carboximidamide"),
    "sulfonamide": ("S(N)(=O)=O", "benzenesulfonamide"),
    "nitrile": ("C#N", "benzonitrile"),
    "aldehyde": ("C=O", "benzaldehyde"),
    "alcohol": ("O", "phenol"),
    "thiol": ("S", "benzenethiol"),
    "amine": ("N", "aniline"),
}

SIMPLE_PREFIXES = {
    "carboxylic_acid": "carboxy",
    "sulfonic_acid": "sulfo",
    "sulfinic_acid": "sulfino",
    "amide": "carbamoyl",
    "hydrazide": "hydrazinecarbonyl",
    "amidine": "carbamimidoyl",
    "sulfonamide": "sulfamoyl",
    "nitrile": "cyano",
    "aldehyde": "formyl",
    "ketone": "oxo",
    "thione": "sulfanylidene",
    "selone": "selanylidene",
    "tellone": "tellanylidene",
    "alcohol": "hydroxy",
    "thiol": "sulfanyl",
    "amine": "amino",
}

_RENAMED = {"phenylmethyl": "benzyl"}

_state = {"depth": 0}


def hook_suspended():
    return _state["depth"] > 0


def labels_in(atoms):
    """{atom: label entry} of the isotopically modified atoms among `atoms` while a labelled structure is named."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    return {a: e for a, e in context["labels"].items() if a in atoms} if context else {}


def label_marked(mol):
    """`mol` with each nuclide label written as an isotope number of its atom, so that canonical ranks and fragment
    SMILES tell structures that differ only in their modification apart."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context or not context["labels"]:
        return mol
    marked = Chem.RWMol(mol)
    for idx, entry in context["labels"].items():
        mass = int("".join(ch for ch in entry["skeleton"] or "" if ch.isdigit()) or 0)
        marked.GetAtomWithIdx(idx).SetIsotope(mass + 1000 * entry["H"].get("2H", 0) + 100000 * entry["H"].get("3H", 0))
    return marked


def consume_labels(atoms):
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if context:
        context["consumed"].update(a for a in atoms if a in context["labels"])


def subtree(mol, root, from_atom):
    seen = {root}
    stack = [root]
    while stack:
        atom = stack.pop()
        for n in mol.GetAtomWithIdx(atom).GetNeighbors():
            idx = n.GetIdx()
            if idx != from_atom and idx not in seen:
                seen.add(idx)
                stack.append(idx)
    return seen


def _copy_fragment(mol, atoms):
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        rw.RemoveAtom(idx)
    mapping = {old: new for new, old in enumerate(sorted(atoms))}
    return rw.GetMol(), mapping


def _carrier_smiles(mol, root, from_atom, suffix_group):
    atoms = subtree(mol, root, from_atom)
    fragment, mapping = _copy_fragment(mol, atoms)
    benzene = Chem.MolFromSmiles("c1ccccc1")
    combined = Chem.RWMol(Chem.CombineMols(benzene, fragment))
    offset = benzene.GetNumAtoms()
    combined.AddBond(0, offset + mapping[root], Chem.BondType.SINGLE)
    if suffix_group is not None:
        group_mol = Chem.MolFromSmiles(SUFFIX_CARRIERS[suffix_group][0])
        group_offset = combined.GetNumAtoms()
        combined = Chem.RWMol(Chem.CombineMols(combined, group_mol))
        combined.AddBond(3, group_offset, Chem.BondType.SINGLE)
    result = combined.GetMol()
    Chem.SanitizeMol(result)
    return Chem.MolToSmiles(result)


def _split_enclosed(text):
    pairs = {"(": ")", "[": "]", "{": "}"}
    if text and text[0] in pairs:
        depth = 0
        for i, ch in enumerate(text):
            if ch in pairs:
                depth += 1
            elif ch in pairs.values():
                depth -= 1
                if depth == 0:
                    return (text[1:i], True) if i == len(text) - 1 else (text, False)
    return text, False


def probe_name(smiles, name_function=None):
    """Name of an auxiliary probe molecule, with the multiplicative hook
    suspended so a probe containing identical rings is named substitutively."""
    if name_function is None:
        from .core import smiles_to_iupac as name_function
    _state["depth"] += 1
    try:
        return name_function(smiles)
    except UnsupportedStructure:
        raise
    except Exception as error:
        raise UnsupportedStructure(f"auxiliary probe could not be named: {error}") from error
    finally:
        _state["depth"] -= 1


def _has_ring(mol, root, from_atom):
    """A substituent holding a ring would compete with the benzene carrier for the parent of the probe name (P-44.1.2.2)."""
    return any(mol.GetAtomWithIdx(a).IsInRing() for a in subtree(mol, root, from_atom))


def _cites_stereo_in(atoms):
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    return bool(context) and (
        any(a in atoms for a in context["atoms"]) or any(a in atoms and b in atoms for a, b in context["bonds"])
    )


def cited_linker_substituents(mol, blocked, linker_atoms):
    """The substituents on `linker_atoms` outside `blocked` cite their own stereo elements (P-93.4.1); a no-op when
    none of them carries one."""
    from contextlib import nullcontext

    from ._common import adjacency, specified_stereo_elements, stereo_element_atoms
    from ._substituents import cited_branch_stereo

    graph = adjacency(mol)
    roots = [n for a in linker_atoms for n in graph[a] if n not in blocked]
    inside, stack = set(), list(roots)
    while stack:
        a = stack.pop()
        if a not in inside:
            inside.add(a)
            stack.extend(n for n in graph[a] if n not in blocked and n not in inside)
    stereo = specified_stereo_elements(mol)
    if not stereo or not inside & set(stereo_element_atoms(mol, stereo)):
        return nullcontext()
    return cited_branch_stereo(mol, graph, blocked, roots)


def prefix_name(mol, root, from_atom, suffix_group=None, name_function=None, groups=None):
    """(name, is_compound) of the substituent rooted at `root` (hanging off
    `from_atom`), named as a prefix in the context of a parent whose
    principal group class is `suffix_group` (None: plain ring parent)."""
    atom = mol.GetAtomWithIdx(root)
    if labels_in(subtree(mol, root, from_atom)) or _cites_stereo_in(subtree(mol, root, from_atom)):
        from ._common import adjacency
        from ._substituents import name_branch

        return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol, unsaturated=True)
    if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() == 1:
        return halogen_prefix(atom), False
    if not atom.IsInRing() and mol.GetBondBetweenAtoms(root, from_atom).GetBondTypeAsDouble() != 1.0:
        from ._common import adjacency
        from ._substituents import name_branch

        return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol, unsaturated=True)
    from ._common import CATIONIC_AMINES
    from ._polyfunctional import AMINIUM

    if suffix_group == "amine" and (AMINIUM.get() is True or CATIONIC_AMINES.get()):
        from ._common import adjacency
        from ._substituents import name_branch

        return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol, unsaturated=True)
    if not atom.IsInRing():
        from ._common import adjacency
        from ._hetero_prefixes import MONONUCLEAR_HYDRIDES, hetero_branch_name

        if atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES:
            result = hetero_branch_name(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol)
            if result is not None:
                return result
    if atom.IsInRing():
        from ._common import adjacency
        from ._substituents import name_branch

        return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol)
    if suffix_group in KETONE_LIKE:
        suffix_group = None
    if suffix_group is not None and suffix_group not in SUFFIX_CARRIERS:
        raise UnsupportedStructure(f"substituents of a {suffix_group} parent are not supported yet")
    if suffix_group is None and _has_ring(mol, root, from_atom):
        from ._common import adjacency
        from ._substituents import name_branch

        try:
            return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol)
        except UnsupportedStructure:
            pass
    name = probe_name(_carrier_smiles(mol, root, from_atom, suffix_group), name_function)
    if suffix_group is None:
        tail = "benzene"
        head = ""
    else:
        tail = SUFFIX_CARRIERS[suffix_group][1]
        head = "4-"
    if not name.startswith(head) or not name.endswith(tail) or len(name) <= len(head) + len(tail):
        if suffix_group is None:
            from ._common import adjacency
            from ._substituents import name_branch

            return name_branch(adjacency(mol), root, from_atom, halogen_substituents(mol), frozenset(), mol=mol, unsaturated=True)
        raise UnsupportedStructure(f"substituent prefix could not be extracted from {name!r}")
    raw = name[len(head) : len(name) - len(tail)]
    inner, enclosed = _split_enclosed(raw)
    if inner in _RENAMED:
        return _RENAMED[inner], False
    return inner, enclosed
