"""Carbamate esters R-O-C(=O)-NR'R'' as the parent of a polyfunctional molecule (P-65.2.1.1, P-65.6.3.2.1, P-41).

An ester outranks an amide, a nitrile, a ketone, an alcohol and every other group junior to it, so the rest of the molecule
is cited as substituent groups: the alkyl group of the ester before the word, the nitrogen substituents in front of
'carbamate' with no locants, because the nitrogen is the only substitutable position: '2-hydroxypropyl
(2-aminoethyl)carbamate (PIN)', 'tert-butyl benzyl[3-(pyridin-3-yl)propyl]carbamate'.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._multiplicative_text import enclose
from ._substituents import BRANCH_STEREO, format_mononuclear_prefixes, name_branch

_CARBAMATE = Chem.MolFromSmarts("[NX3;!R]-[CX3;!R](=O)-[OX2;!R]-[#6]")
_SENIOR_TO_CARBAMATE = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3;!$(C(N)(=O)O)](=O)[OX2H1]",
        "[CX3;!$(C(N)(=O)O)](=O)[OX2][#6]",
        "[CX3](=O)[OX2][CX3]=O",
        "[CX3](=O)[F,Cl,Br,I]",
        "[#16X4](=O)(=O)[OX2]",
        "[#15](=O)[OX2]",
        "[CX3](=O)([NX3])[OX2][#6].[CX3](=O)([NX3])[OX2][#6]",
    )
]


def _side(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _arm_name(mol, graph, root, parent, halogens, aromatic, codes):
    atoms = _side(graph, root, parent)
    context = {
        "atoms": {i: c for i, c in codes[0].items() if i in atoms},
        "bonds": {b: c for b, c in codes[1].items() if b[0] in atoms and b[1] in atoms},
        "used": set(),
    }
    token = BRANCH_STEREO.set(context)
    try:
        named = name_branch(graph, root, parent, halogens, aromatic, mol=mol)
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of a carbamate substituent is not cited by its name")
    return named, atoms


def name_carbamate_ester(mol) -> str:
    """The name of a molecule whose senior group is one acyclic carbamate ester, else raises."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("charged, isotopic and multi-fragment structures are not named as carbamate esters here")
    matches = mol.GetSubstructMatches(_CARBAMATE)
    if len(matches) != 1:
        raise UnsupportedStructure("exactly one acyclic carbamate ester is required")
    if any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_CARBAMATE):
        raise UnsupportedStructure("a group senior to the carbamate ester is present")
    nitrogen, carbon, _, ester_oxygen, alkyl = matches[0]
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    codes = (
        {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.HasProp("_CIPCode")},
        {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode") for b in probe.GetBonds() if b.HasProp("_CIPCode")},
    )
    alkyl_name, alkyl_atoms = _arm_name(mol, graph, alkyl, ester_oxygen, halogens, aromatic, codes)
    if carbon in alkyl_atoms or nitrogen in alkyl_atoms:
        raise UnsupportedStructure("a cyclic carbamate is not named as an ester")
    entries, covered = [], {nitrogen, carbon, ester_oxygen} | alkyl_atoms
    carbonyl_oxygen = next(n for n in graph[carbon] if n not in (nitrogen, ester_oxygen))
    covered.add(carbonyl_oxygen)
    for root in graph[nitrogen]:
        if root == carbon:
            continue
        (name, compound), atoms = _arm_name(mol, graph, root, nitrogen, halogens, aromatic, codes)
        if atoms & covered:
            raise UnsupportedStructure("a ring closes through the carbamate nitrogen")
        covered |= atoms
        entries.append((name, compound))
    if len(covered) != mol.GetNumAtoms():
        raise UnsupportedStructure("an atom of the molecule lies outside the carbamate and its substituents")
    alkyl_text = enclose(alkyl_name[0]) if alkyl_name[1] and " " in alkyl_name[0] else alkyl_name[0]
    if not entries:
        return f"{alkyl_text} carbamate"
    if len(entries) == 1 and entries[0][1]:
        return f"{alkyl_text} {enclose(entries[0][0])}carbamate"
    return f"{alkyl_text} {format_mononuclear_prefixes(entries)}carbamate"
