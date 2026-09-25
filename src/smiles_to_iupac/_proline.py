"""Naming of proline (pyrrolidine-2-carboxylic acid) using its P-103
retained name and L/D stereodescriptor, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained name 'proline'.
- P-103.1.2: "The atoms in proline are numbered as in pyrrolidine, the
  nitrogen atom being numbered '1', and the carbon atom bonded to the
  carboxy group is numbered '2'." Unlike every other common amino acid,
  proline's own amino group is the ring nitrogen itself (a secondary
  amine), not an exocyclic primary amine -- `_carboxylic_acid_amine.py`
  explicitly excludes rings and requires a primary amine, so it can never
  reach this shape; a separate module is needed, not an extension of that
  one (confirmed empirically during scoping, #1056 M2 step 1 -- both
  L-/D-proline currently misroute to `_ketone.py`'s hetero-ring-ketone
  path, an unrelated module that claims any saturated single-heteroatom
  5/6/7-membered ring purely by shape).
- P-103.1.3.1: 'L' corresponds to CIP 'S' at the alpha-carbon (proline is
  not the cysteine exception), confirmed against real PubChem L-/D-
  proline structures during scoping.

Scope, deliberately narrow: a plain pyrrolidine ring (5-membered, one ring
nitrogen bearing a single ring-hydrogen and no substituent) whose C-2 (the
ring carbon adjacent to nitrogen) bears an exocyclic -COOH and nothing
else, with every other ring carbon a plain -CH2-. Any substituent on the
ring beyond the C-2 carboxy group (e.g. hydroxyproline) is out of scope,
deferred to a follow-up.
"""

from ._carboxylic_acid_amine import _find_carboxylic_acid_carbon
from ._common import adjacency, specified_stereocenters

_ALPHA_TO_LD = {"S": "L", "R": "D"}


def _match(mol):
    """alpha_carbon_idx if `mol` is plain proline's exact shape, else
    None."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring = ring_info.AtomRings()[0]
    if len(ring) != 5 or any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return None
    ring_set = set(ring)

    nitrogens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    if len(nitrogens) != 1:
        return None
    (ring_n,) = nitrogens
    if any(mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in ring if a != ring_n):
        return None

    n_atom = mol.GetAtomWithIdx(ring_n)
    if n_atom.GetFormalCharge() != 0 or n_atom.GetDegree() != 2 or n_atom.GetTotalNumHs() != 1:
        return None

    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return None
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    acid_carbon_idx = acid_carbon.GetIdx()
    if acid_carbon_idx in ring_set:
        return None

    graph = adjacency(mol)
    alpha_candidates = [n for n in graph[acid_carbon_idx] if n in ring_set]
    if len(alpha_candidates) != 1:
        return None
    (alpha_carbon,) = alpha_candidates
    if ring_n not in graph[alpha_carbon]:
        return None

    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if alpha_atom.GetFormalCharge() != 0 or alpha_atom.GetDegree() != 3 or alpha_atom.GetTotalNumHs() != 1:
        return None

    for ring_carbon in ring_set - {ring_n, alpha_carbon}:
        atom = mol.GetAtomWithIdx(ring_carbon)
        if atom.GetFormalCharge() != 0 or atom.GetDegree() != 2 or atom.GetTotalNumHs() != 2:
            return None

    excluded = ring_set | {acid_carbon_idx, carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    if mol.GetNumAtoms() != len(excluded):
        return None
    return alpha_carbon


def has_proline_shape(mol) -> bool:
    return _match(mol) is not None


def name_proline(mol) -> str:
    alpha_carbon = _match(mol)
    stereo = specified_stereocenters(mol)
    if not stereo:
        return "proline"
    labels = {atom_idx: label for atom_idx, label in stereo}
    if alpha_carbon not in labels:
        return "proline"
    return f"{_ALPHA_TO_LD[labels[alpha_carbon]]}-proline"
