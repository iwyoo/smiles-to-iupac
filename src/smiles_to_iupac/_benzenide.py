"""Naming of the benzenide anion (phenyl anion, P-72.2.2.1, Chapter P-7,
https://iupac.qmul.ac.uk/BlueBook/P7.html): a benzene ring with one ring
carbon carrying the sole formal charge -1 names as 'benzenide (PIN)' --
`tmp/bluebook/P7.txt` lines 908-911 and 946-947 give this exact worked
example directly. Unlike the cyclopentadienide anion (`_cyclopentadienide.py`),
this needs no locant citation at all: every ring position is equivalent by
symmetry, so the anion center's own numbering is never essential
(mirroring how a lone ring double bond's locant is never cited in
`_cyclic_unsaturated.py`, P-14.3.3).

Scope, deliberately narrow: exactly one 6-membered, all-carbon, monocyclic
aromatic ring, exactly one ring carbon at formal charge -1, no other
charge/isotope anywhere, and no atom outside the ring at all (no
substituents). Any other ring size, any substituent, or more than one
anionic center is out of scope for this pass.
"""

_RING_SIZE = 6


def _is_benzenide(mol) -> bool:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return False
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _RING_SIZE or mol.GetNumAtoms() != _RING_SIZE:
        return False
    for idx in ring_atoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return False
    if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
        return False

    anions = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetFormalCharge() == -1]
    if len(anions) != 1:
        return False
    others_neutral = all(
        mol.GetAtomWithIdx(idx).GetFormalCharge() == 0 for idx in ring_atoms if idx != anions[0]
    )
    if not others_neutral:
        return False
    return all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)


def has_benzenide_shape(mol) -> bool:
    return _is_benzenide(mol)


def name_benzenide(mol) -> str:
    return "benzenide"
