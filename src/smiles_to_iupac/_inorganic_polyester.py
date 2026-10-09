"""One polyvalent alcohol component esterified by several sulfate, sulfite, phosphate or phosphite groups, each left
with its remaining acid hydroxyl groups (P-65.6.3.2.2, P-67.1.3.2): the component is cited once as a 'diyl' group before
the multiplied anion, 'ethane-1,2-diyl bis(hydrogen sulfate)'."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._diester_anions import _INORGANIC_ANIONS
from ._diester_ring_diyl import name_diester_ring_diyl

_CENTRES = (15, 16)


def _centre_match(mol, atom):
    """(centre, doubly bonded oxygen or the centre itself, ester oxygen, alcohol carbon) of an acid centre esterified
    once, or None."""
    if atom.GetAtomicNum() not in _CENTRES or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetTotalNumHs():
        return None
    esters, doubled = [], []
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != 8 or n.GetFormalCharge():
            return None
        order = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
        if order == 2.0 and n.GetDegree() == 1:
            doubled.append(n)
        elif order == 1.0 and n.GetDegree() == 2:
            esters.append(n)
        elif not (order == 1.0 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1):
            return None
    if len(esters) != 1 or (atom.GetAtomicNum(), len(doubled)) not in _INORGANIC_ANIONS:
        return None
    if atom.GetDegree() != {(16, 2): 4, (16, 1): 3, (15, 1): 4, (15, 0): 3}[(atom.GetAtomicNum(), len(doubled))]:
        return None
    carbon = next(n for n in esters[0].GetNeighbors() if n.GetIdx() != atom.GetIdx())
    if carbon.GetAtomicNum() != 6 or carbon.GetIsAromatic():
        return None
    return atom, doubled[0] if doubled else atom, esters[0], carbon


def inorganic_ester_matches(mol):
    """The esterified acid centres when every functional group of the molecule is one of them, else []."""
    if len(Chem.GetMolFrags(mol)) > 1:
        return []
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in _CENTRES:
            match = _centre_match(mol, atom)
            if match is None:
                return []
            matches.append(match)
    if len(matches) < 2:
        return []
    owned = {n.GetIdx() for m in matches for n in m[0].GetNeighbors()} | {m[0].GetIdx() for m in matches}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in owned:
            continue
        if atom.GetFormalCharge() or atom.GetIsotope() or atom.GetAtomicNum() not in (6, 8, *HALOGEN_PREFIXES):
            return []
        if atom.GetAtomicNum() == 8 and (atom.GetDegree() > 2 or atom.IsInRing()):
            return []
        if any(b.GetBondTypeAsDouble() == 2.0 and (b.GetBeginAtom().GetAtomicNum() != 6 or b.GetEndAtom().GetAtomicNum() != 6) for b in atom.GetBonds()):
            return []
    return matches


def has_inorganic_polyester_shape(mol) -> bool:
    return bool(inorganic_ester_matches(mol))


def name_inorganic_polyester(mol) -> str:
    matches = inorganic_ester_matches(mol)
    if not matches:
        raise UnsupportedStructure("several singly esterified sulfate, sulfite, phosphate or phosphite groups are required")
    return name_diester_ring_diyl(mol, matches)
