"""Naming of histidine (2-amino-3-(1H-imidazol-4-yl)propanoic acid) using
its P-103 retained name and L/D stereodescriptor, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained name 'histidine', formula
  2-amino-3-(1H-imidazol-4-yl)propanoic acid.
- P-103.1.3.1: 'L' corresponds to CIP 'S' at the alpha-carbon (histidine
  is not the cysteine exception), confirmed against real PubChem L-/D-
  histidine structures' CIP labels during scoping (#1056 M2 step 2).
- P-103.1.2's special pi/tau imidazole-ring numbering only matters when
  citing a *substituent* on the ring itself (e.g. a methylated histidine
  derivative) -- plain histidine's own retained name needs no internal
  ring locant at all, so that numbering system is out of scope here.

Unlike proline (`_proline.py`), histidine's amino/carboxyl groups are NOT
on the ring -- the imidazole ring is an ordinary side-chain substituent,
structurally much closer to phenylalanine's shape than proline's.
`_carboxylic_acid_amine.py` already has a working "aromatic-ring-on-a-
chain" path, but it only recognizes a *plain benzene* ring
(`is_plain_benzene_ring`); an imidazole ring falls through that check
into the module's generic ring rejection instead (confirmed empirically
during scoping -- misleadingly worded as "a -COOH and/or -NH2 group on/in
a ring," even though neither group is actually on the ring here).

Scope, deliberately narrow: the alpha-amino-acid backbone (acyclic,
`_carboxylic_acid_amine.py`'s own shared detection) with a side chain of
exactly one -CH2- bonded to a plain 5-membered aromatic ring bearing
exactly 2 nitrogens (one NH, one pyridine-type) and 3 carbons, only one of
which carries the -CH2- branch (the other two are plain aromatic -CH=).
Which specific ring nitrogen sits next to the substituted carbon is not
checked -- the NH's position is a fast tautomeric/drawing choice, not a
distinct molecule, so either arrangement is accepted. Any substituent on
the ring itself (would need the special pi/tau numbering) is out of
scope, deferred to a follow-up.
"""

from ._carboxylic_acid_amine import _amine_on_a_different_carbon, _find_carboxylic_acid_carbon
from ._common import adjacency, find_primary_amines, specified_stereocenters

_ALPHA_TO_LD = {"S": "L", "R": "D"}


def _is_imidazole_ring_with_one_substituent(mol, graph, ring_root, coming_from):
    ring_info = mol.GetRingInfo()
    ring = next((r for r in ring_info.AtomRings() if ring_root in r), None)
    if ring is None or len(ring) != 5:
        return False
    if any(not mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return False
    ring_set = set(ring)

    nitrogens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    carbons = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    if len(nitrogens) != 2 or len(carbons) != 3:
        return False

    nh_count = sum(1 for n in nitrogens if mol.GetAtomWithIdx(n).GetTotalNumHs() == 1)
    pyridine_count = sum(1 for n in nitrogens if mol.GetAtomWithIdx(n).GetTotalNumHs() == 0)
    if nh_count != 1 or pyridine_count != 1:
        return False
    if any(mol.GetAtomWithIdx(n).GetDegree() != 2 for n in nitrogens):
        return False

    substituted = [c for c in carbons if [n for n in graph[c] if n not in ring_set]]
    if substituted != [ring_root]:
        return False
    plain = [c for c in carbons if c != ring_root]
    if any(mol.GetAtomWithIdx(c).GetDegree() != 2 or mol.GetAtomWithIdx(c).GetTotalNumHs() != 1 for c in plain):
        return False

    exo = [n for n in graph[ring_root] if n not in ring_set]
    if exo != [coming_from]:
        return False
    if mol.GetAtomWithIdx(ring_root).GetDegree() != 3:
        return False
    return True


def _match(mol):
    """alpha_carbon_idx if `mol` is plain histidine's exact shape, else
    None."""
    if mol.GetRingInfo().NumRings() != 1:
        return None
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return None
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    amines = find_primary_amines(mol)
    if len(amines) != 1:
        return None
    if not _amine_on_a_different_carbon(mol, acid_carbon.GetIdx(), amines):
        return None
    (amine_n,) = amines

    graph = adjacency(mol)
    acid_carbon_idx = acid_carbon.GetIdx()
    alpha_carbon = graph[amine_n][0]
    if acid_carbon_idx not in graph[alpha_carbon]:
        return None

    exclude = {amine_n, acid_carbon_idx}
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if alpha_atom.GetAtomicNum() != 6 or alpha_atom.GetFormalCharge() != 0 or alpha_atom.GetTotalNumHs() != 1:
        return None
    side_neighbors = [n for n in graph[alpha_carbon] if n not in exclude]
    if len(side_neighbors) != 1:
        return None
    (ch2,) = side_neighbors
    ch2_atom = mol.GetAtomWithIdx(ch2)
    if ch2_atom.GetAtomicNum() != 6 or ch2_atom.GetFormalCharge() != 0 or ch2_atom.GetTotalNumHs() != 2:
        return None
    ring_neighbors = [n for n in graph[ch2] if n != alpha_carbon]
    if len(ring_neighbors) != 1:
        return None
    (ring_root,) = ring_neighbors
    if not _is_imidazole_ring_with_one_substituent(mol, graph, ring_root, ch2):
        return None

    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    ring = next(r for r in mol.GetRingInfo().AtomRings() if ring_root in r)
    accounted = exclude | excluded_acid_oxygens | {alpha_carbon, ch2} | set(ring)
    if mol.GetNumAtoms() != len(accounted):
        return None
    return alpha_carbon


def has_histidine_shape(mol) -> bool:
    return _match(mol) is not None


def name_histidine(mol) -> str:
    alpha_carbon = _match(mol)
    stereo = specified_stereocenters(mol)
    if not stereo:
        return "histidine"
    labels = {atom_idx: label for atom_idx, label in stereo}
    if alpha_carbon not in labels:
        return "histidine"
    return f"{_ALPHA_TO_LD[labels[alpha_carbon]]}-histidine"
