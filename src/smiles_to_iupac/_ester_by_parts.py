"""Ester names built from their two parts (P-65.6.3.2.1): the alkyl/aryl group
cited as a substituent group, then the anion of the acid part. The molecule is
split at the O-alkyl bond; the acid fragment is named as a free acid and its
'-ic acid' ending becomes '-ate', so every acid shape the engine can name also
works as an ester acyl part.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._numerals import multiplying_prefix
from ._substituents import BRANCH_STEREO, name_branch

_ESTER = Chem.MolFromSmarts("[CX3;!R](=O)[OX2;!R][#6]")
_FREE_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
_HYDROGEN_WORDS = {1: "hydrogen", 2: "dihydrogen", 3: "trihydrogen"}


def _carbonyl_with_heteroatom(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    has_double_o = any(
        n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )
    return has_double_o and any(
        n.GetAtomicNum() in (7, 8, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in atom.GetNeighbors()
    )


def _insert_before_ending(anion, text):
    for ending in ("oate", "ate"):
        if anion.endswith(ending):
            return anion[: -len(ending)] + text + ending
    raise UnsupportedStructure("the ester ending is not delimited")


def _anion_name(acid_name):
    if not acid_name.endswith(" acid") or " " in acid_name[: -len(" acid")].strip():
        raise UnsupportedStructure("the acid part of this ester is not named as a plain acid")
    stem = acid_name[: -len(" acid")]
    if not stem.endswith("ic"):
        raise UnsupportedStructure("the acid part of this ester has an unexpected acid name ending")
    return stem[:-2] + "ate"


def _branch_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def name_ester_by_parts(mol) -> str:
    from ._isotope_labels import split_isotopes

    split = split_isotopes(mol)
    if split is None:
        return _name_ester_parts(mol, {})
    clean, labels, _ = split
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms()):
        raise UnsupportedStructure("stereodescriptors of an isotopically modified ester are not supported yet (P-82.4)")
    return _name_ester_parts(clean, labels)


def _restore_labels(editable, clean, labels, removed):
    """Put the nuclide labels of the atoms kept in the acid part back on `editable`."""
    from ._isotope_labels import HYDROGEN_ISOTOPES

    for idx, entry in labels.items():
        if idx in removed:
            continue
        atom = editable.GetAtomWithIdx(idx)
        if entry["skeleton"]:
            atom.SetIsotope(int("".join(ch for ch in entry["skeleton"] if ch.isdigit())))
        total = sum(entry["H"].values())
        if total:
            atom.SetNumExplicitHs(clean.GetAtomWithIdx(idx).GetTotalNumHs() - total)
            atom.SetNoImplicit(True)
        for nuclide, count in entry["H"].items():
            for _ in range(count):
                hydrogen = Chem.Atom(1)
                hydrogen.SetIsotope(HYDROGEN_ISOTOPES[nuclide])
                editable.AddBond(idx, editable.AddAtom(hydrogen), Chem.BondType.SINGLE)


def _name_ester_parts(mol, labels) -> str:
    from ._substituents import ISOTOPE_LABELS
    from .core import smiles_to_iupac

    matches = mol.GetSubstructMatches(_ESTER)
    if not matches:
        raise UnsupportedStructure("an acyclic ester group is required for part-wise ester naming")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    arms = []
    for acyl_carbon, _, ester_oxygen, alkyl_carbon in matches:
        atoms = _branch_atoms(graph, alkyl_carbon, ester_oxygen)
        if acyl_carbon in atoms:
            raise UnsupportedStructure("a ring-closing ester (lactone) is not named part-wise")
        if any(_carbonyl_with_heteroatom(mol, a) for a in atoms):
            raise UnsupportedStructure("an acid or ester group inside the alkyl part outranks this ester")
        if any(mol.GetAtomWithIdx(a).GetIsotope() or mol.GetAtomWithIdx(a).GetNumRadicalElectrons() for a in atoms):
            raise UnsupportedStructure("isotopes and radicals in the alkyl part are not supported")
        arms.append(atoms)
    removed = set().union(*arms)
    if sum(len(a) for a in arms) != len(removed):
        raise UnsupportedStructure("the alkyl parts share atoms, so these esters are of a polyol, not a polyacid")

    ester_labels = {}
    free_hydrogens = []
    labels = dict(labels)
    if labels:
        bridging = {ester_oxygen: i for i, (_, _, ester_oxygen, _) in enumerate(matches)}
        carbonyl = {
            m.GetIdx()
            for acyl_carbon, _, _, _ in matches
            for m in mol.GetAtomWithIdx(acyl_carbon).GetNeighbors()
            if m.GetAtomicNum() == 8 and m.GetIdx() not in bridging
        }
        for idx in [i for i in labels if i in bridging or i in carbonyl]:
            entry = labels[idx]
            if idx in carbonyl and not entry["skeleton"] and sum(entry["H"].values()) == 1 and mol.GetAtomWithIdx(idx).GetTotalNumHs() == 1:
                free_hydrogens.extend(entry["H"])
                del labels[idx]
                continue
            if entry["H"] or not entry["skeleton"]:
                raise UnsupportedStructure("this isotopic modification of an ester oxygen is not supported yet (P-82.6.4)")
            ester_labels[idx] = (entry["skeleton"], bridging.get(idx))
            del labels[idx]
        if len({a for _, _, _, a in matches}) != len(matches):
            raise UnsupportedStructure("this isotopic modification of an ester oxygen is not supported yet (P-82.6.4)")

    editable = Chem.RWMol(mol)
    for _, _, ester_oxygen, _ in matches:
        oxygen = editable.GetAtomWithIdx(ester_oxygen)
        oxygen.SetNumExplicitHs(1)
        oxygen.SetNoImplicit(True)
    _restore_labels(editable, mol, labels, removed)
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    if len(Chem.GetMolFrags(acid)) != 1:
        raise UnsupportedStructure("the acid parts of these esters are separate groups")
    Chem.SanitizeMol(acid)
    anion = _anion_name(smiles_to_iupac(Chem.MolToSmiles(acid)))

    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    context = {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in removed and a.HasProp("_CIPCode")},
        "bonds": {
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
            for b in probe.GetBonds()
            if b.GetBeginAtomIdx() in removed and b.GetEndAtomIdx() in removed and b.HasProp("_CIPCode")
        },
        "used": set(),
    }
    named = {}
    arm_label = {i: nuclide for nuclide, i in ester_labels.values() if i is not None}
    token = BRANCH_STEREO.set(context)
    isotope_context = {"labels": {a: e for a, e in labels.items() if a in removed}, "consumed": set(), "mol": mol}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if isotope_context["labels"] else None)
    try:
        for arm, (acyl_carbon, _, ester_oxygen, alkyl_carbon) in enumerate(matches):
            name, compound = name_branch(graph, alkyl_carbon, ester_oxygen, halogens, mol=mol)
            locant = (f"{arm_label[arm][:-1]}O" if arm in arm_label else "O") if ester_labels else ""
            entry = named.setdefault((name, locant), [0, compound, []])
            entry[0] += 1
    finally:
        BRANCH_STEREO.reset(token)
        ISOTOPE_LABELS.reset(isotope_token)
    if set(isotope_context["labels"]) - isotope_context["consumed"]:
        raise UnsupportedStructure("an isotopically modified atom of the alkyl part is not cited by any supported name")
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of the alkyl part is not cited by any supported name")
    if ester_labels:
        counts = {}
        for nuclide, _ in ester_labels.values():
            counts[nuclide] = counts.get(nuclide, 0) + 1
        text = "(" + ",".join(f"{n}{c}" for n, c in sorted(counts.items())) + ")"
        anion = text + anion if anion == "carbonate" else _insert_before_ending(anion, text)
    grouped = {}
    for (name, locant), (count, compound, _) in named.items():
        grouped.setdefault(name, []).append((locant, count, compound))
    parts = []
    for name in sorted(grouped, key=alpha_sort_key):
        for locant, count, compound in sorted(grouped[name]):
            text = name if count == 1 else (
                f"{multiplying_prefix(count, compound=compound)}({name})" if compound else f"{multiplying_prefix(count, compound=compound)}{name}"
            )
            if locant:
                text = f"{','.join([locant] * count)}-{text}"
            parts.append(text)
    free = len(acid.GetSubstructMatches(_FREE_ACID)) - len(matches)
    if free > 0:
        if free not in _HYDROGEN_WORDS:
            raise UnsupportedStructure("too many free acid groups beside the esters")
        word = _HYDROGEN_WORDS[free]
        if free_hydrogens:
            if free != 1 or len(free_hydrogens) != 1:
                raise UnsupportedStructure("isotopic modification of several acid hydrogens is not supported yet")
            word = f"({free_hydrogens[0]})hydrogen"
        parts.append(word)
    elif free_hydrogens:
        raise UnsupportedStructure("an isotopically modified hydroxy hydrogen of this ester is not supported yet")
    return " ".join(parts + [anion])
