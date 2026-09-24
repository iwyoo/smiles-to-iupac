"""Naming of the biphenyl ring assembly (two benzene rings joined by a single
bond, sharing no atom), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-28.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): two
  identical cyclic parent hydrides joined directly by a single bond are named
  as a ring assembly using the parent hydride's name preceded by 'bi'; two
  benzene rings joined this way have the retained name 'biphenyl'. The
  locants of the two ring-joining atoms are compulsorily cited as '1,1'' even
  when the compound is otherwise unsubstituted, giving the PIN
  '1,1'-biphenyl'.
- P-28.2.2 / P-14.3.2 (numbering): each ring is numbered independently,
  starting at its own point of attachment (locant 1), with one ring's
  locants left unprimed and the other's primed. Substituent locants are
  chosen to be as low as possible as a set, where an unprimed locant is
  considered lower than the same number primed -- this module reuses the
  same lowest-locant-set/citation-order tie-break machinery as every other
  ring module (`_cyclic.py`, `_aromatic.py`), just with string locants (e.g.
  "4", "4'") instead of plain integers so ordinary tuple/string comparison
  already encodes that "unprimed < primed at the same number" rule.
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the same
  way as in every other ring module.

Scope, deliberately narrow: only two *identical* benzene rings connected by
exactly one single (non-aromatic) bond, each bearing at most simple
substituents (halogens, alkyl). A different pair of rings or non-benzene
ring assemblies are out of scope and fall through to `UnsupportedStructure`
elsewhere in the dispatch chain. Three to six benzene rings in an
unbranched chain (terphenyl etc.) are `_ring_assembly_chain.py`'s job
instead -- a separate composite-locant numbering scheme (P-28.3), not a
generalization of this module's own primed-locant one.
"""

from ._common import (
    adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._substituents import format_substituent_prefixes, name_branch


def _bond_between(bond, atoms_a, atoms_b):
    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    return (x in atoms_a and y in atoms_b) or (x in atoms_b and y in atoms_a)


def find_ring_assembly_core(mol):
    """Return (ring0_atoms, ring1_atoms, attach0, attach1) if `mol` is
    exactly two disjoint 6-membered all-carbon aromatic rings joined by one
    single bond, else None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    for ring in atom_rings:
        if len(ring) != 6:
            return None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                return None

    ring0, ring1 = set(atom_rings[0]), set(atom_rings[1])
    if ring0 & ring1:
        return None

    connecting = [bond for bond in mol.GetBonds() if _bond_between(bond, ring0, ring1)]
    if len(connecting) != 1:
        return None
    bond = connecting[0]
    if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() != 1.0:
        return None

    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    attach0, attach1 = (x, y) if x in ring0 else (y, x)
    return atom_rings[0], atom_rings[1], attach0, attach1


def _numberings_from_attachment(graph, ring_atoms, attach, prime):
    cycle = ring_cycle(graph, list(ring_atoms))
    start = cycle.index(attach)
    rotated = cycle[start:] + cycle[:start]
    suffix = "'" if prime else ""
    for seq in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
        yield {atom: f"{position}{suffix}" for position, atom in enumerate(seq, start=1)}


def _candidate_key(locants, ring_atoms, graph, halogens, mol=None):
    substituents = {}
    for atom, position in locants.items():
        branch_roots = [n for n in graph[atom] if n not in ring_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]

    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if prefix:
        prefix += "-"
    name = prefix + "1,1'-biphenyl"
    return locant_set, citation_locants, name


def name_ring_assembly(mol, core) -> str:
    validate_atoms_and_bonds(mol)

    ring0_atoms, ring1_atoms, attach0, attach1 = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(ring0_atoms) | set(ring1_atoms)

    best_key = None
    best_name = None
    for (ring_a, attach_a), (ring_b, attach_b) in (
        ((ring0_atoms, attach0), (ring1_atoms, attach1)),
        ((ring1_atoms, attach1), (ring0_atoms, attach0)),
    ):
        for locants_a in _numberings_from_attachment(graph, ring_a, attach_a, prime=False):
            for locants_b in _numberings_from_attachment(graph, ring_b, attach_b, prime=True):
                locants = {**locants_a, **locants_b}
                key = _candidate_key(locants, ring_atoms, graph, halogens, mol=mol)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]

    return best_name
