"""Naming of partially unsaturated monocyclic all-carbon hydrocarbons
(cycloalkenes, cycloalkadienes, ...) per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-31.1.4.2 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  a ring double bond is denoted the same way as in an acyclic chain --
  changing 'ane' to 'ene' -- with locants as low as possible given to the
  double bonds as a set; numbering starts from any ring atom (there is no
  fixed starting point as in a chain), so this is always achievable at
  locant '1' for a single double bond.
- P-31.1.1.2 (Chapter P-3): two or more double bonds take a multiplying
  prefix ('di', 'tri', ...) before 'ene', with the euphonic 'a' inserted
  after the parent stem when a locant-bearing multiplying prefix follows
  (e.g. 'cyclohexa-1,3-diene (PIN)'), exactly as for the acyclic case (see
  `_unsaturated.py`).
- P-14.3.3 (Chapter P-1): a locant is cited only when essential. For a
  *single* ring double bond, the ring's free choice of starting atom and
  direction always allows it to be numbered '1' regardless of any
  substituents present, so that locant is never essential and is always
  omitted (e.g. 'cyclohexene', '3-methylcyclohexene' -- never
  'cyclohex-1-ene' or '3-methylcyclohex-1-ene'). This differs from the
  acyclic rule (P-14.3.4.2(d)), where omission is restricted to short,
  unsubstituted chains, because an acyclic chain's numbering start is not
  fully free. With two or more double bonds the locant set is no longer
  redundant (different relative spacings of the double bonds are
  genuinely different structures), so it is always cited.
- P-14.4(e) (Chapter P-1): among numbering choices that are otherwise
  tied, lower locants go to the double bonds as a set before substituent
  locants (mirrored from `_unsaturated.py`'s acyclic tie-break).
- P-29.4 / P-46 (Chapter P-2, P-4) and P-35.2.1 (Chapter P-3): ring
  substituents (simple/compound alkyl groups, halogens) are named exactly
  as in `_cyclic.py`.

Scope, deliberately narrow: a single monocyclic, all-carbon ring bearing
one or more carbon-carbon ring double bonds (no triple bonds, no
exocyclic unsaturation, no heteroatoms, no fused/bridged/spiro
combination, no aromaticity -- benzene itself is `_aromatic.py`'s). Any
molecule not exactly matching this shape falls through to another module
in `core.py` (most commonly `_cyclic.py`, whose existing "unsaturated
rings are not supported yet" message still applies to out-of-scope cases
like a ring triple bond).
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch


def _ring_cycle(graph, ring_atoms):
    """Order a simple ring's atoms into a cyclic sequence by walking its
    bonds (works regardless of bond order -- adjacency doesn't encode it)."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def find_cyclic_unsaturated_core(mol):
    """Return the ring's atom indices if `mol` is a single monocyclic,
    all-carbon, non-aromatic ring carrying at least one ring C=C double
    bond and no unsaturation anywhere else, else None (deliberately
    permissive about the invalid neighboring shapes this module itself
    rejects, e.g. a ring triple bond -- returning None there lets it fall
    through to `_cyclic.py`'s existing rejection instead, see module
    docstring)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    ring_set = set(ring_atoms)
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None

    has_ring_double_bond = False
    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if order != 2.0 or a not in ring_set or b not in ring_set:
            return None
        has_ring_double_bond = True
    if not has_ring_double_bond:
        return None
    return list(ring_atoms)


def _double_bonds(mol):
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
    ]


def _ring_bond_locants(ring_order, double_bonds):
    n = len(ring_order)
    position = {atom: i + 1 for i, atom in enumerate(ring_order)}
    locants = []
    for a, b in double_bonds:
        pa, pb = position[a], position[b]
        locants.append(n if {pa, pb} == {1, n} else min(pa, pb))
    return sorted(locants)


def _substituents_for_ring(graph, ring_order, halogens):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(ring_size, ene_locants, grouped):
    parent_stem = "cyclo" + alkane_name(ring_size)[:-3]
    prefix = format_substituent_prefixes(grouped)
    count = len(ene_locants)
    if count == 1:
        # P-14.3.3: always achievable at locant '1', so never cited (see
        # module docstring).
        return prefix + parent_stem + "ene"
    word = numerical_term(count) + "ene"
    loc_str = ",".join(str(loc) for loc in ene_locants)
    return prefix + parent_stem + "a-" + loc_str + "-" + word


def _candidate_key(ring_size, ene_locants, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(ring_size, ene_locants, grouped)
    return ene_locant_set, locant_set, citation_locants, name


def name_cyclic_unsaturated(mol, ring_atoms) -> str:
    validate_atoms_and_bonds(mol)
    ring_set = set(ring_atoms)

    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if order != 2.0:
            raise UnsupportedStructure(
                "only carbon-carbon double bonds are supported in a "
                "partially unsaturated ring (see P-31.1.4.2); a triple "
                "bond is not supported yet"
            )
        if a not in ring_set or b not in ring_set:
            raise UnsupportedStructure(
                "a double bond outside the ring (an exocyclic alkylidene "
                "substituent) is not supported yet"
            )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    double_bonds = _double_bonds(mol)
    ring_order = _ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            ene_locants = _ring_bond_locants(candidate, double_bonds)
            substituents = _substituents_for_ring(graph, candidate, halogens)
            key = _candidate_key(ring_size, ene_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name
