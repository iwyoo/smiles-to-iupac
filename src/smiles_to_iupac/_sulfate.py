"""Esters of sulfuric, sulfurous, selenic, selenous, telluric and tellurous acid, R groups identical or mixed, with
the partial esters cited with 'hydrogen' and their salts with the cation in front (P-67.1.3.2, P-72.2.2.2.1). One
implementation for the -ic and -ous acids of the three elements; the word assembly is shared with `_phosphate.py`."""

from rdkit import Chem

from ._acid_derivatives import protonated, without_hydrogen_word
from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import centre_nuclides
from ._phosphate import check_ester_organyl, cited_ester_groups, format_ester_entries
from ._salt import _MONOATOMIC_CATION_NAMES

_NOT_ESTER_ROOTS = {7, 8}
_ANION_WORDS = {16: ("sulfite", "sulfate"), 34: ("selenite", "selenate"), 52: ("tellurite", "tellurate")}


def _centre_atoms(mol, oxo):
    """Chalcogen atoms shaped like a fully or partially esterified acid with `oxo` double-bonded terminal oxygens and
    two further single-bonded oxygens: an ester oxygen, a hydroxyl or one deprotonated oxygen, with at least one ester."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ANION_WORDS or atom.GetDegree() != oxo + 2 or atom.GetFormalCharge() != 0 or atom.IsInRing():
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) != oxo + 2:
            continue
        orders = {o.GetIdx(): mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() for o in oxygens}
        double_os = [o for o in oxygens if orders[o.GetIdx()] == 2.0]
        single_os = [o for o in oxygens if orders[o.GetIdx()] == 1.0]
        if len(double_os) != oxo or len(single_os) != 2:
            continue
        if any(o.GetDegree() != 1 or o.GetFormalCharge() != 0 for o in double_os):
            continue
        charged_os = [o for o in single_os if o.GetFormalCharge() == -1]
        neutral_os = [o for o in single_os if o.GetFormalCharge() == 0]
        if len(charged_os) + len(neutral_os) != 2 or len(charged_os) > 1 or any(o.GetDegree() != 1 for o in charged_os):
            continue
        ester_os = [o for o in neutral_os if o.GetDegree() == 2]
        if not ester_os or any(
            n.GetAtomicNum() in _NOT_ESTER_ROOTS for o in ester_os for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx()
        ):
            continue
        matches.append(atom)
    return matches


def _name(mol, oxo) -> str:
    centres = _centre_atoms(mol, oxo)
    if len(centres) != 1:
        raise UnsupportedStructure("more than one group of this acid is not supported yet")
    (centre,) = centres

    group_oxygens = {n.GetIdx() for n in centre.GetNeighbors() if n.GetAtomicNum() == 8}
    charged_oxygens = {idx for idx in group_oxygens if mol.GetAtomWithIdx(idx).GetFormalCharge() == -1}

    frags = Chem.GetMolFrags(mol)
    anion_frag = next(frag for frag in frags if centre.GetIdx() in frag)
    other_frags = [frag for frag in frags if frag is not anion_frag]

    if charged_oxygens and not other_frags:
        return without_hydrogen_word(_name(protonated(mol, charged_oxygens), oxo))

    cation_name = None
    if charged_oxygens:
        if len(other_frags) != 1 or len(other_frags[0]) != 1:
            raise UnsupportedStructure("a salt with other than one monoatomic counter-ion is not supported yet")
        (cation_idx,) = other_frags[0]
        cation_atom = mol.GetAtomWithIdx(cation_idx)
        cation_name = _MONOATOMIC_CATION_NAMES.get((cation_atom.GetSymbol(), cation_atom.GetFormalCharge()))
        if cation_name is None or cation_atom.GetFormalCharge() != 1:
            raise UnsupportedStructure("only a single +1 monoatomic cation is supported for a partial-ester salt yet")
    elif other_frags:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    check_ester_organyl(mol, group_oxygens | {centre.GetIdx()}, bool(charged_oxygens), atoms=anion_frag)

    double_os = {
        n.GetIdx()
        for n in centre.GetNeighbors()
        if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(centre.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    }
    single_oxygens = [idx for idx in group_oxygens if idx not in double_os]
    ester_oxygens = [idx for idx in single_oxygens if mol.GetAtomWithIdx(idx).GetDegree() == 2]
    has_hydroxyl = len(single_oxygens) - len(ester_oxygens) - len(charged_oxygens) > 0

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != centre.GetIdx()]
        roots.append((root, oxygen_idx))

    entries = cited_ester_groups(mol, graph, roots, halogens, aromatic_atoms)
    ester_words = format_ester_entries(entries)
    nuclide = centre_nuclides(centre.GetIdx(), oxygens=True)
    word = _ANION_WORDS[centre.GetAtomicNum()][oxo == 2]
    anion_name = f"{ester_words} hydrogen {nuclide}{word}" if has_hydroxyl else f"{ester_words} {nuclide}{word}"
    return f"{cation_name} {anion_name}" if cation_name else anion_name


def has_sulfate_shape(mol) -> bool:
    return bool(_centre_atoms(mol, 2))


def has_sulfite_shape(mol) -> bool:
    return bool(_centre_atoms(mol, 1))


def name_sulfate(mol) -> str:
    return _name(mol, 2)


def name_sulfite(mol) -> str:
    return _name(mol, 1)
