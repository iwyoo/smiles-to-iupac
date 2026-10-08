"""Mancude heteromonocycles with ring atoms of nonstandard bonding number (P-14.1.3, P-22.2.3, P-31.1.3.2): 1λ4,3-thiazine.
The Hantzsch-Widman name of the ring is kept, the λ symbol follows the locant of its atom, and the heteroatoms take the
lowest locants as a set, then in the order O, S, Se, Te, N, P, then the higher bonding number first (P-44.4.1.3.2)."""

from rdkit import Chem

from ._common import adjacency, ring_cycle
from ._ring_lambda_heterone import _SENIORITY, _stem

_LAMBDA_STANDARD = {"S": 2, "Se": 2, "Te": 2, "P": 3, "Cl": 1, "Br": 1, "I": 1}


def _match(mol):
    """(ring order, {ring atom: element}, {ring atom: bonding number}) when every ring atom of a substituent-free
    monocycle but the λ atoms takes part in exactly one ring double bond."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    info = mol.GetRingInfo()
    if info.NumRings() != 1 or not 5 <= len(info.AtomRings()[0]) <= 10 or mol.GetNumAtoms() != len(info.AtomRings()[0]):
        return None
    ring = list(info.AtomRings()[0])
    hetero, lam = {}, {}
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetNumRadicalElectrons():
            return None
        doubles = sum(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())
        if any(b.GetBondTypeAsDouble() == 3.0 for b in atom.GetBonds()):
            return None
        symbol = atom.GetSymbol()
        if symbol == "C" or symbol == "N":
            if doubles != 1 or (symbol == "N" and atom.GetTotalNumHs()) or (symbol == "C" and atom.GetTotalNumHs() > 1):
                return None
            if symbol == "N":
                hetero[idx] = "N"
        elif symbol in _LAMBDA_STANDARD:
            valence = atom.GetTotalValence()
            if doubles == 0 and valence == _LAMBDA_STANDARD[symbol]:
                hetero[idx] = symbol
            elif doubles in (1, 2) and valence != _LAMBDA_STANDARD[symbol]:
                hetero[idx] = symbol
                lam[idx] = valence
            else:
                return None
        elif symbol == "O" and not doubles:
            hetero[idx] = "O"
        else:
            return None
    if not lam:
        return None
    return ring_cycle(adjacency(mol), ring), hetero, lam


def has_lambda_ring_shape(mol) -> bool:
    return _match(mol) is not None


def name_lambda_ring(mol) -> str:
    order, hetero, lam = _match(mol)
    size = len(order)
    best = None
    for base in (order, order[::-1]):
        for start in range(size):
            candidate = base[start:] + base[:start]
            position = {a: i + 1 for i, a in enumerate(candidate)}
            locants = sorted(position[a] for a in hetero)
            seniority = tuple(
                position[a] for element in _SENIORITY for a in sorted((x for x, e in hetero.items() if e == element), key=position.get)
            )
            lam_locants = sorted(position[a] for a in lam)
            key = (locants, seniority, lam_locants, [-lam[a] for a in sorted(lam, key=position.get)])
            if best is None or key < best[0]:
                best = (key, position)
    position = best[1]
    cited = []
    for element in _SENIORITY:
        for a in sorted((x for x, e in hetero.items() if e == element), key=position.get):
            cited.append(f"{position[a]}λ{lam[a]}" if a in lam else str(position[a]))
    stem = _stem([hetero[a] for a in order if a in hetero], size, True)
    return f"{','.join(cited)}-{stem}"
