"""Substituent prefix names for multiplicative names (P-29.4): simple groups
come from a table; anything else is named by attaching the substituent to a
benzene carrier (with the principal group of the parent para to it, so
senior-class groups inside the prefix are demoted correctly) and stripping
the carrier from the resulting name.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure

SUFFIX_CARRIERS = {
    "carboxylic_acid": ("C(=O)O", "benzoic acid"),
    "sulfonic_acid": ("S(=O)(=O)O", "benzenesulfonic acid"),
    "amide": ("C(N)=O", "benzamide"),
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
    "sulfonamide": "sulfamoyl",
    "nitrile": "cyano",
    "aldehyde": "formyl",
    "ketone": "oxo",
    "alcohol": "hydroxy",
    "thiol": "sulfanyl",
    "amine": "amino",
}

_RENAMED = {"phenylmethyl": "benzyl"}

_state = {"depth": 0}


def hook_suspended():
    return _state["depth"] > 0


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


def prefix_name(mol, root, from_atom, suffix_group=None, name_function=None, groups=None):
    """(name, is_compound) of the substituent rooted at `root` (hanging off
    `from_atom`), named as a prefix in the context of a parent whose
    principal group class is `suffix_group` (None: plain ring parent)."""
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() == 1:
        return HALOGEN_PREFIXES[atom.GetAtomicNum()], False
    if not atom.IsInRing():
        from ._common import adjacency
        from ._hetero_prefixes import MONONUCLEAR_HYDRIDES, hetero_branch_name

        if atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES:
            result = hetero_branch_name(adjacency(mol), root, from_atom, {}, frozenset(), mol)
            if result is not None:
                return result
    if atom.IsInRing():
        from ._multiplicative_groups import classify
        from ._multiplicative_ring import ring_substituent_name

        ring_atoms = next(set(r) for r in mol.GetRingInfo().AtomRings() if root in r)
        result = ring_substituent_name(
            mol, ring_atoms, root, from_atom, groups if groups is not None else classify(mol), suffix_group, name_function
        )
        if result is not None:
            return result
    if suffix_group is not None and suffix_group not in SUFFIX_CARRIERS:
        raise UnsupportedStructure(f"substituents of a {suffix_group} parent are not supported yet")
    name = probe_name(_carrier_smiles(mol, root, from_atom, suffix_group), name_function)
    if suffix_group is None:
        tail = "benzene"
        head = ""
    else:
        tail = SUFFIX_CARRIERS[suffix_group][1]
        head = "4-"
    if not name.startswith(head) or not name.endswith(tail) or len(name) <= len(head) + len(tail):
        raise UnsupportedStructure(f"substituent prefix could not be extracted from {name!r}")
    raw = name[len(head) : len(name) - len(tail)]
    inner, enclosed = _split_enclosed(raw)
    if inner in _RENAMED:
        return _RENAMED[inner], False
    return inner, enclosed
