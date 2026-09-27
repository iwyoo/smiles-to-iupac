"""Naming of the cyclopentadienide anion (P-72.2.2.1, Chapter P-7,
https://iupac.qmul.ac.uk/BlueBook/P7.html): a 5-membered all-carbon ring
with one formal-charge -1 ring carbon and two ring C=C double bonds (the
neutral-parent-equivalent "cyclopenta-1,3-diene" shape) names as
'cyclopenta-2,4-dien-1-ide (PIN)' -- `tmp/bluebook/P7.txt` lines 916-919
and 946-949 give this exact worked example directly ('cyclopentadienide',
the P-76 delocalized-anion alternative name, is noted alongside but is
not the PIN). `_carbanide.py` is explicitly acyclic-only, so this
ring+diene shape needs its own module rather than a tweak to that one's
scope.

RDKit perceives the anion's 6 pi-electron ring as aromatic (bond order
1.5 throughout, e.g. `c1cc[cH-]c1`), not as the two explicit ring double
bonds the Blue Book's own name is built from -- so the diene locants are
read off a Kekulized copy of the ring, fixing the anion at locant 1
(P-72.2.2.1, mirroring how a suffix group anchors ring numbering
elsewhere in this project) and choosing ring direction for the lowest
ene-locant set via `_cyclic_unsaturated.py`'s own bond-locant helper.

Scope, deliberately narrow (mirroring the issue's own worked example):
exactly one 5-membered, all-carbon, monocyclic ring, exactly one ring
carbon at formal charge -1, no other charge/isotope anywhere, and no atom
outside the ring at all (no substituents). Any other ring size, any
substituent, or more than one anionic center is out of scope for this
pass.
"""

from rdkit import Chem

from ._common import adjacency, ring_cycle
from ._cyclic_unsaturated import _bond_locant


def _find_ring_anion(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 5 or mol.GetNumAtoms() != 5:
        return None
    for idx in ring_atoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            return None
    if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
        return None

    anions = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetFormalCharge() == -1]
    if len(anions) != 1:
        return None
    others_neutral = all(
        mol.GetAtomWithIdx(idx).GetFormalCharge() == 0 for idx in ring_atoms if idx != anions[0]
    )
    if not others_neutral:
        return None
    if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
        return None
    return anions[0], ring_atoms


def has_cyclopentadienide_shape(mol) -> bool:
    return _find_ring_anion(mol) is not None


def name_cyclopentadienide(mol) -> str:
    anion_idx, ring_atoms = _find_ring_anion(mol)

    kekulized = Chem.Mol(mol)
    Chem.Kekulize(kekulized, clearAromaticFlags=True)
    graph = adjacency(kekulized)
    ring_order = ring_cycle(graph, list(ring_atoms))
    start = ring_order.index(anion_idx)
    rotated = ring_order[start:] + ring_order[:start]

    double_bonds = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in kekulized.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
    ]

    best_locants = None
    for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
        locants = sorted(_bond_locant(candidate, bond) for bond in double_bonds)
        if best_locants is None or locants < best_locants:
            best_locants = locants

    loc_str = ",".join(str(loc) for loc in best_locants)
    return f"cyclopenta-{loc_str}-dien-1-ide"
