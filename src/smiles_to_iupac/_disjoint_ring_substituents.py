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
- The parent ring is the senior one by P-44.2.1 (heteroatoms, then more skeletal
  atoms), then the aromatic one (P-44.4.1.1); for two otherwise equal rings the lower
  parent attachment locant wins, then the alphanumerically earlier name --
  the Blue Book gives no further criterion for that tie. A bridge may
  carry multiple bonds and may join the parent ring through a double
  bond, e.g. '(cyclohexylidenemethyl)benzene' (P-44.4.1.1's own example).
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

from ._common import UnsupportedStructure, adjacency, alpha_sort_key
from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
from ._multiplicative_ring import monocycle_spec, numberings
from ._multiplicative_text import enclose
from ._ring_system_seniority import ring_seniority_key
from ._substituents import name_branch


def find_disjoint_ring_pair_core(mol):
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    ring_a, ring_b = (set(r) for r in atom_rings)
    if ring_a & ring_b:
        return None
    # P-44.1.2.1: a bridge atom of Si, P, B, ... outranks carbon, so a ring is not the parent
    if any(a.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not a.IsInRing() for a in mol.GetAtoms()):
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
    return attach_atom, bridge_atom


def _is_plain_monocycle(graph, ring, attach, bridge):
    members = set(ring)
    return all(
        {n for n in graph[a] if n not in members} == ({bridge} if a == attach else set()) for a in ring
    ) and all(sum(1 for n in graph[a] if n in members) == 2 for a in ring)


def name_disjoint_ring_pair(mol, core) -> str:
    ring_a, ring_b = core
    graph = adjacency(mol)
    aromatic_atoms = frozenset(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())

    rings = []
    for ring in (ring_a, ring_b):
        attach, bridge = _ring_attachment(mol, graph, ring)
        if not _is_plain_monocycle(graph, ring, attach, bridge):
            raise UnsupportedStructure(
                "a fused, spiro, or otherwise non-simple ring, alongside a "
                "second disjoint ring elsewhere in the molecule, is not "
                "supported yet (see P-23/P-24/P-25)"
            )
        spec = monocycle_spec(mol, ring)
        if spec is None:
            raise UnsupportedStructure(
                "an unsaturated or otherwise unsupported monocycle, alongside a "
                "second disjoint ring elsewhere in the molecule, is not "
                "supported yet (see P-31.1.4, P-44.2)"
            )
        rings.append((attach, bridge, spec))

    def rank(spec):
        return ring_seniority_key(mol, set(spec.cycle))

    best = min(rank(spec) for _, _, spec in rings)
    names = []
    for attach, bridge, spec in rings:
        if rank(spec) != best:
            continue
        name, is_compound = name_branch(graph, bridge, attach, {}, aromatic_atoms, mol=mol, unsaturated=True)
        display_name = enclose(name) if is_compound else name
        if spec.hetero is None:
            names.append((0, f"{display_name}{spec.parent}"))
            continue
        locant = min(numberings(spec), key=lambda loc: loc[attach])[attach]
        joiner = "-" if spec.parent[0].isdigit() or spec.parent[0] == "Δ" else ""
        names.append((locant, f"{locant}-{display_name}{joiner}{spec.parent}"))
    return min(names, key=lambda item: (item[0], alpha_sort_key(item[1])))[1]
