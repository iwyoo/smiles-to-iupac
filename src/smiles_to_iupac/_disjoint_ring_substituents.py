"""Naming of two disjoint (mutually unconnected) plain ring systems joined
only through an acyclic bridge -- e.g. dicyclohexylmethane, 1-cyclohexyl-
3-phenylpropane-shaped molecules -- per the IUPAC 2013 Recommendations:

- P-44.1.1 (Chapter P-4): a ring or ring system is always senior to a
  chain as the parent hydride, so one of the two rings is the parent and
  the entire bridge, with the other ring at its far end, becomes a single
  compound substituent prefix (P-29.4.1) -- reusing `name_branch`'s
  general P-46 chain-selection mechanism, extended (see
  `_substituents.py`'s `_longest_chains_from_root` docstring) to let a
  chain walk terminate at a separate plain ring instead of raising.
- Between an aromatic and a saturated ring, the aromatic ring is senior
  (P-44.1.1's own seniority-of-rings order); between two rings of the same
  class the choice doesn't change the resulting name for the molecules in
  scope here (each such case in this module's own tests is a symmetric
  pair of identical rings), so either is picked.
- Scope: exactly two disjoint plain rings (no shared atom, no direct
  ring-to-ring bond -- already claimed by `_ring_assembly.py` before this
  module runs), each bearing no substituent of its own besides the single
  bridge attachment -- a ring with any other substituent, or a fused/
  spiro/polycyclic ring, is out of scope and raises. Cross-checked against
  real PubChem structures: 'C1CCCCC1CCC1CCCCC1' (CID 76838,
  '2-cyclohexylethylcyclohexane'), 'C1CCCCC1CCCc1ccccc1' (CID 561990,
  '3-cyclohexylpropylbenzene'), 'c1ccccc1CCCc1ccccc1' (CID 14125,
  '3-phenylpropylbenzene'), 'C1CCCCC1CC1CCCCC1' (CID 76644, PubChem's own
  name 'cyclohexylmethylcyclohexane' omits the compound-substituent
  parentheses this project's own established convention always keeps,
  e.g. '(chloromethyl)benzene' -- P-16.3.3 calls for them, so this module
  follows that existing project convention over PubChem's own looser
  string).
"""

from ._common import UnsupportedStructure, adjacency
from ._numerals import alkane_name
from ._substituents import _simple_ring_substituent, name_branch


def find_disjoint_ring_pair_core(mol):
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    ring_a, ring_b = (set(r) for r in atom_rings)
    if ring_a & ring_b:
        return None
    return ring_a, ring_b


def _ring_attachment(mol, graph, ring):
    external = [(atom, n) for atom in ring for n in graph[atom] if n not in ring]
    if len(external) != 1:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic substituent, alongside a "
            "second disjoint ring elsewhere in the molecule, is not "
            "supported yet (see P-29.2)"
        )
    (attach_atom, bridge_atom), = external
    if mol.GetBondBetweenAtoms(attach_atom, bridge_atom).GetBondTypeAsDouble() != 1.0:
        # e.g. a ring-assembly-ylidene shape (two rings joined by a C=C),
        # already claimed elsewhere in core.py if it fits that shape --
        # `name_branch`'s own chain walk has no bond-order awareness for
        # this first, ring-external bond either (see its docstring), so it
        # must be checked here rather than silently dropped.
        raise UnsupportedStructure(
            "a non-single bond joining a ring to the rest of the molecule "
            "is not supported in this disjoint-ring-pair shape"
        )
    return attach_atom, bridge_atom


def name_disjoint_ring_pair(mol, core) -> str:
    ring_a, ring_b = core
    graph = adjacency(mol)
    aromatic_atoms = frozenset(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())

    attach_a, bridge_a = _ring_attachment(mol, graph, ring_a)
    attach_b, bridge_b = _ring_attachment(mol, graph, ring_b)
    shape_a = _simple_ring_substituent(graph, attach_a, bridge_a, aromatic_atoms, mol=mol)
    shape_b = _simple_ring_substituent(graph, attach_b, bridge_b, aromatic_atoms, mol=mol)
    if shape_a is None or shape_b is None:
        raise UnsupportedStructure(
            "a fused, spiro, or otherwise non-simple ring, alongside a "
            "second disjoint ring elsewhere in the molecule, is not "
            "supported yet (see P-23/P-24/P-25)"
        )
    size_a, aromatic_a, _ = shape_a
    size_b, aromatic_b, _ = shape_b

    if aromatic_b and not aromatic_a:
        attach_a, bridge_a, size_a, aromatic_a = attach_b, bridge_b, size_b, aromatic_b

    name, is_compound = name_branch(graph, bridge_a, attach_a, {}, aromatic_atoms, mol=mol)
    display_name = f"({name})" if is_compound else name
    parent = "benzene" if aromatic_a else "cyclo" + alkane_name(size_a)
    return f"{display_name}{parent}"
