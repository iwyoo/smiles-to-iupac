"""P-62.2.3, P-62.2.1.1.1: the prefix 'anilino' is retained and substitution on it is allowed, so an amino nitrogen
carrying a benzene ring is named like the N-substituted aniline with the ending 'anilino' ('3-(N-methylanilino)phenol',
'4-chloroanilino')."""

from ._common import ring_cycle
from ._substituents import name_branch, substituents_for_ring


def _benzene_ring_of(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
        return None
    info = mol.GetRingInfo()
    for ring in info.AtomRings():
        if carbon in ring:
            if len(ring) != 6 or any(info.NumAtomRings(i) != 1 for i in ring):
                return None
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                return None
            return ring
    return None


def anilino_prefix(graph, mol, nitrogen, others, halogens, aromatic_atoms):
    """(name, True) for -N(R)-Ar with Ar a lone benzene ring, else None; the ring with the most substituents is the anilino
    ring (P-45.1) and any other nitrogen substituent is an N-prefix."""
    from ._amine import _aniline_candidate_key

    best = None
    for carbon in others:
        ring = _benzene_ring_of(mol, carbon)
        if ring is None:
            continue
        n_names = [name_branch(graph, other, nitrogen, halogens, aromatic_atoms, mol=mol) for other in others if other != carbon]
        ring_order = ring_cycle(graph, list(ring))
        for start in range(len(ring_order)):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                substituents = substituents_for_ring(graph, candidate, halogens, {nitrogen}, mol=mol)
                count = sum(len(v) for v in substituents.values())
                key = (-count, *_aniline_candidate_key(position_of[carbon], substituents, n_names))
                if best is None or key < best:
                    best = key
    if best is None:
        return None
    name = best[-1]
    return name[: -len("aniline")] + "anilino", True
