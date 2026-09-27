"""Naming of a chalcogen (S/Se/Te) ring-oxide on an otherwise plain,
unsubstituted monocyclic ring (mancude or saturated), per the IUPAC 2013
Recommendations:

- P-62.5's functional-class "oxide" nomenclature isn't limited to acyclic
  amines (`_amine_oxide.py`) -- confirmed via `tmp/bluebook/P7.txt:3135`
  ("1,2-thiazole 1-oxide") and `tmp/bluebook/P6a.txt:944` ("1,2-thiazinane
  1-oxide"): a single chalcogen oxide on an otherwise standard-valence
  ring is named by appending "<locant>-oxide" to the ring's own plain
  name, no lambda-convention involved (contrast P-22.2.7, whose own
  worked examples are all bare parent hydrides needing extra indicated
  hydrogen, never a discrete oxide substituent).
- Built exactly like `_amine_oxide.py`: remove the oxide oxygen, re-
  sanitize, and name the reduced ring via `_hetero_monocyclic.py`'s
  existing canonical-structure lookup -- not a parallel implementation.
  The sole heteroatom is always locant 1 for these unsubstituted single-
  heteroatom rings (P-22.2.2.1.2), so the oxide's own locant is fixed.
- Confirmed against PubChem CID 9548690: `C1=CS(=O)C=C1` -> "thiophene
  1-oxide" exactly. Removing the oxide oxygen from that structure and
  re-sanitizing reproduces plain thiophene's own canonical SMILES,
  verifying the reduction is chemically sound (the ring re-aromatizes).
  The same reduction naturally extends to a saturated ring, confirmed
  against PubChem CID 534965: `O=S1CCCCC1` -> "thiane 1-oxide".

Out of scope: a ring with more than one heteroatom (needs the parent's
own real numbering to place a locant other than 1), a charge-separated
oxide depiction (pyridine N-oxide's own `[N+]`/`[O-]` shape is
`_amine_oxide.py`'s neighboring but distinct territory), any substituent
elsewhere on the ring, and any fused/bridged ring system.
"""

from rdkit import Chem

from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic

_CHALCOGENS = {16, 34, 52}


def _reduced_ring(mol, oxide_oxygen_idx):
    reduced = Chem.RWMol(mol)
    reduced.RemoveAtom(oxide_oxygen_idx)
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)
    return reduced_mol


def _find_ring_oxide(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = set(ring_info.AtomRings()[0])
    candidates = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetIdx() in ring_atoms
        and atom.GetAtomicNum() in _CHALCOGENS
        and atom.GetFormalCharge() == 0
        and atom.GetIsotope() == 0
    ]
    if len(candidates) != 1:
        return None
    (heteroatom,) = candidates
    oxide_neighbors = [
        n
        for n in heteroatom.GetNeighbors()
        if n.GetIdx() not in ring_atoms
        and n.GetAtomicNum() == 8
        and n.GetFormalCharge() == 0
        and n.GetIsotope() == 0
        and n.GetDegree() == 1
    ]
    if len(oxide_neighbors) != 1:
        return None
    (oxide_oxygen,) = oxide_neighbors
    if mol.GetBondBetweenAtoms(heteroatom.GetIdx(), oxide_oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
        return None
    return oxide_oxygen.GetIdx()


def has_hetero_ring_oxide_shape(mol) -> bool:
    oxide_oxygen_idx = _find_ring_oxide(mol)
    if oxide_oxygen_idx is None:
        return False
    return has_hetero_monocyclic_name(_reduced_ring(mol, oxide_oxygen_idx))


def name_hetero_ring_oxide(mol) -> str:
    oxide_oxygen_idx = _find_ring_oxide(mol)
    base_name = name_hetero_monocyclic(_reduced_ring(mol, oxide_oxygen_idx))
    return f"{base_name} 1-oxide"
