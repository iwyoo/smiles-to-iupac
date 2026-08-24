"""Skeletal replacement ('a' prefix) naming of a saturated monospiro ring
system (`_spiro.py`) containing exactly one ring heteroatom, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-24.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) and
  the companion spiro-nomenclature document
  (https://iupac.qmul.ac.uk/spiro/sp16p.html): "The numbering of the spiro
  hydrocarbon ring system is never modified by the introduction of
  heteroatoms, but low locants must be attributed to heteroatoms if there
  is a choice." I.e. the ring-size-based numbering priority
  (`_spiro.py`'s `iter_monospiro_numberings`) is unchanged by a
  heteroatom; the heteroatom's locant is minimized only among the
  numbering choices that rule already leaves open (direction of traversal
  within each ring; which of two equal-size rings is numbered first),
  ahead of substituent locants -- cross-checked against real named
  compounds: 6-oxaspiro[4.5]decane (heteroatom in the larger ring, not
  10-oxa; PubChem CID 12630235) and 1-oxaspiro[4.5]decane (heteroatom in
  the smaller ring, adjacent to the spiro atom, not 4-oxa; PubChem CID
  79743's 7,9-dione derivative confirms the same ring-1 numbering).
- The spiro atom itself must stay carbon: P-24.2.1 numbers it using all
  four of its ring bonds (two shared with each ring, no room for an
  exocyclic substituent or a lone pair), a valence no neutral O/N/S atom
  can supply -- so this module's single heteroatom is always at a
  non-spiro ring position.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one ring heteroatom (of any kind), a heteroatom at the spiro
  atom itself, or a heteroatom outside the two rings entirely.
- Any heteroatom other than O, N, or S.
- Combination with polyspiro (`_polyspiro.py`) or branched polyspiro
  systems, unsaturation, or anything else `_spiro.name_monospiro`'s
  validation already rejects for the all-carbon case.
"""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, non_single_bonds
from ._cyclic import _substituents_for_ring
from ._spiro import _candidate_key, iter_monospiro_numberings

_HETEROATOM_PREFIXES = {8: "oxa", 7: "aza", 16: "thia"}
_ALLOWED_ATOMIC_NUMS = {6, *_HETEROATOM_PREFIXES, *HALOGEN_PREFIXES}


def _ring_atoms(mol):
    atom_rings = mol.GetRingInfo().AtomRings()
    return set(atom_rings[0]) | set(atom_rings[1])


def has_single_ring_heteroatom_shape(mol, spiro_atom) -> bool:
    """True iff `mol` has exactly one O/N/S atom, and it's a ring atom
    other than `spiro_atom` itself (the shape this module accepts) --
    excludes a heteroatom that's an exocyclic substituent (e.g. a
    hydroxyl on an otherwise all-carbon spiro skeleton), which is a
    different module's territory."""
    if mol.GetAtomWithIdx(spiro_atom).GetAtomicNum() in _HETEROATOM_PREFIXES:
        return False
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(heteroatoms) != 1:
        return False
    return heteroatoms[0].GetIdx() in _ring_atoms(mol)


def name_spiro_heteroatom(mol, spiro_atom) -> str:
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a single ring O/N/S (P-24.2.1) and "
                "halogen substituents (P-35.2.1) are not supported yet"
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
            "unsaturated spiro ring systems are not supported yet (see "
            "P-31.1.5, unsaturated alicyclic spiro ring systems)"
        )

    heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(heteroatoms) != 1 or heteroatoms[0] not in _ring_atoms(mol):
        raise UnsupportedStructure(
            "exactly one *ring* heteroatom is supported (two or more ring "
            "heteroatoms, 'a'-prefix seniority ordering among different "
            "kinds, and an exocyclic heteroatom substituent like -OH, are "
            "out of scope here)"
        )
    (heteroatom_idx,) = heteroatoms
    if heteroatom_idx == spiro_atom:
        raise UnsupportedStructure(
            "a heteroatom at the spiro atom itself is not supported here "
            "(a neutral O/N/S can't supply the four ring bonds a spiro "
            "atom needs; see module docstring)"
        )
    a_prefix = _HETEROATOM_PREFIXES[mol.GetAtomWithIdx(heteroatom_idx).GetAtomicNum()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    for parent, full_order in iter_monospiro_numberings(mol, spiro_atom):
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
