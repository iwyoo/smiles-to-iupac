"""Naming of a chalcogen (S/Se/Te) ring-oxide on an otherwise plain,
unsubstituted retained-name fused aromatic bicyclic ring -- the fused
analogue of `_hetero_ring_oxide.py`'s single-ring case (P-62.5's
functional-class "oxide" nomenclature), per the IUPAC 2013
Recommendations:

- Same construction as the monocyclic case: remove the oxide oxygen,
  re-sanitize, and name the reduced ring via `_heteroaromatic_fused.py`'s
  existing retained-name lookup -- not a parallel implementation.
- The oxide's own locant is the chalcogen atom's fixed ring position in
  the retained name's own systematic numbering, which for the two
  relevant retained names ('1-benzothiophene'/'2-benzothiophene' -- the
  only chalcogen-heteroatom entries in that table; quinoline/isoquinoline/
  indole have a nitrogen instead, and O itself can't carry a further oxide
  the way S/Se/Te can) happens to equal the leading disambiguating numeral
  already in the name itself: S is ring position 1 in 1-benzothiophene,
  position 2 in 2-benzothiophene, per each isomer's own Blue Book
  systematic numbering.
- Confirmed against PubChem CID 5383918: `C1=CC=C2C(=C1)C=CS2=O` ->
  "1-benzothiophene 1-oxide" exactly. Removing the oxide oxygen and
  re-sanitizing reproduces plain 1-benzothiophene's own canonical SMILES.

Out of scope: a second ring heteroatom, any substituent elsewhere on the
ring, a charge-separated oxide depiction, a polycyclic (3+ ring) system,
and any hypervalent-heterocycle base needing the lambda convention.
"""

from rdkit import Chem

from ._heteroaromatic_fused import has_retained_heteroaromatic_fused_name, name_retained_heteroaromatic_fused

_CHALCOGENS = {16, 34, 52}


def _reduced_ring(mol, oxide_oxygen_idx):
    reduced = Chem.RWMol(mol)
    reduced.RemoveAtom(oxide_oxygen_idx)
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)
    return reduced_mol


def _find_fused_ring_oxide(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 2:
        return None
    ring_atoms = {idx for ring in ring_info.AtomRings() for idx in ring}
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


def _match(mol):
    oxide_oxygen_idx = _find_fused_ring_oxide(mol)
    if oxide_oxygen_idx is None:
        return None
    reduced = _reduced_ring(mol, oxide_oxygen_idx)
    if not has_retained_heteroaromatic_fused_name(reduced):
        return None
    base_name = name_retained_heteroaromatic_fused(reduced)
    if not base_name[0].isdigit():
        return None
    locant = base_name.split("-", 1)[0]
    return base_name, locant


def has_fused_hetero_ring_oxide_shape(mol) -> bool:
    return _match(mol) is not None


def name_fused_hetero_ring_oxide(mol) -> str:
    base_name, locant = _match(mol)
    return f"{base_name} {locant}-oxide"
