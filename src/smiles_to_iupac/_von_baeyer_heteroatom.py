"""Skeletal replacement ('a' prefix) naming of a saturated von Baeyer
bicyclic ring system (`_bicyclic.py`) containing exactly one ring
heteroatom, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-23.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  skeletal replacement ('a') nomenclature names a heteroatom in an
  otherwise-carbon von Baeyer ring system by citing a nondetachable
  prefix ('oxa', 'thia', 'aza', ...) with the locant of the replaced
  skeletal atom, directly in front of the 'bicyclo[x.y.z]alkane' parent
  (e.g. '7-oxabicyclo[2.2.1]heptane', '1-azabicyclo[2.2.2]octane' —
  quinuclidine). The bridge-length brackets and the ring's numbering
  (P-23.2.3) are unaffected by which skeletal atom is a heteroatom; only
  which of the numbering choices `_bicyclic.py`'s
  `iter_bicyclic_numberings` produces is preferred changes: the
  heteroatom must get the lowest locant available among them, ahead of
  substituent locants (mirrors the substituent-locant tie-break already
  used for a plain hydrocarbon bicyclic).
- Table 2.8's element seniority order governs which heteroatom is 'a'-1
  when several *different* kinds are present; moot here since this module
  only ever accepts exactly one.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one skeletal heteroatom (of any kind), or a heteroatom outside
  the ring core (e.g. a substituent -OH/-NH2 alongside a plain hydrocarbon
  bicyclic — a different module's territory).
- Any heteroatom other than O, N, or S.
- Unsaturation, charged/isotopic atoms, or anything else
  `_bicyclic.name_bicycloalkane`'s validation already rejects for the
  all-carbon case.
"""

from ._bicyclic import _candidate_key, bicyclic_parent_name, iter_bicyclic_numberings
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, non_single_bonds
from ._cyclic import _substituents_for_ring

_HETEROATOM_PREFIXES = {8: "oxa", 7: "aza", 16: "thia"}
_ALLOWED_ATOMIC_NUMS = {6, *_HETEROATOM_PREFIXES, *HALOGEN_PREFIXES}


def _ring_heteroatoms(mol, core):
    bh1, bh2, bridges = core
    core_atoms = {bh1, bh2} | {atom for bridge in bridges for atom in bridge}
    return [atom for atom in core_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() in _HETEROATOM_PREFIXES]


def has_single_ring_heteroatom_shape(mol, core) -> bool:
    """True iff `core` (a `_bicyclic.find_bicyclic_core` result) has exactly
    one O/N/S skeletal ring atom, the shape this module accepts."""
    return len(_ring_heteroatoms(mol, core)) == 1


def name_von_baeyer_heteroatom(mol, core) -> str:
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a single skeletal O/N/S (P-23.2.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated bicyclic ring systems are not supported yet (see "
            "P-31.1.4, unsaturated von Baeyer ring systems)"
        )

    ring_heteroatoms = _ring_heteroatoms(mol, core)
    all_heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(ring_heteroatoms) != 1 or set(all_heteroatoms) != set(ring_heteroatoms):
        raise UnsupportedStructure(
            "exactly one skeletal ring heteroatom is supported (P-23.2.1's "
            "'a'-prefix ordering for two or more heteroatoms, and any "
            "heteroatom outside the ring skeleton, are out of scope here)"
        )
    (heteroatom_idx,) = ring_heteroatoms
    a_prefix = _HETEROATOM_PREFIXES[mol.GetAtomWithIdx(heteroatom_idx).GetAtomicNum()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    parent = bicyclic_parent_name(core)

    best_key = None
    for full_order in iter_bicyclic_numberings(core):
        heteroatom_locant = full_order.index(heteroatom_idx) + 1
        substituents = _substituents_for_ring(graph, full_order, halogens)
        key = _candidate_key(
            parent,
            substituents,
            heteroatom_locant=heteroatom_locant,
            nondetachable_prefix=f"{heteroatom_locant}-{a_prefix}",
        )
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]
