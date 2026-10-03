"""Naming of partially unsaturated monocyclic all-carbon hydrocarbons
(cycloalkenes, cycloalkadienes, cycloalkynes, mixed cycloalkenynes, ...)
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-31.1.3.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  a ring double or triple bond is denoted the same way as in an acyclic
  chain -- changing 'ane' to 'ene'/'yne' -- with one multiple bond always
  allocated locant '1' (numbering starts from any ring atom, there being
  no fixed starting point as in a chain). With both ene and yne bonds
  present, locants as low as possible go to the full set first, even if
  this gives 'yne' a lower locant than 'ene'; if a choice remains, low
  locants go to the double bonds specifically -- the same two-step rule
  as P-31.1.1.1 for chains, and directly worked in this section's own
  example '(1Z)-cyclopentadec-1-en-4-yne (PIN)'.
- P-31.1.1.2 (Chapter P-3): two or more multiple bonds of one kind take a
  multiplying prefix ('di', 'tri', ...) before 'ene'/'yne', with the
  euphonic 'a' inserted after the parent stem when a locant-bearing
  multiplying prefix follows (e.g. 'cyclohexa-1,3-diene (PIN)'), exactly
  as for the acyclic case (see `_unsaturated.py`, which this module
  mirrors for the ene/yne suffix-construction logic).
- P-14.3.3 (Chapter P-1): a locant is cited only when essential. For a
  ring bearing exactly one multiple bond total, the ring's free choice of
  starting atom and direction always allows it to be numbered '1'
  regardless of any substituents present, so that locant is never
  essential and is always omitted (e.g. 'cyclohexene', 'cyclooctyne',
  '3-methylcyclohexene' -- never 'cyclohex-1-ene' or 'cyclooct-1-yne').
  This differs from the acyclic rule (P-14.3.4.2(d)), where omission is
  restricted to short, unsubstituted chains, because an acyclic chain's
  numbering start is not fully free. With two or more multiple bonds the
  locant set is no longer redundant (different relative spacings are
  genuinely different structures), so it is always cited.
- P-14.4(e) (Chapter P-1): among numbering choices that are otherwise
  tied, lower locants go to the multiple bonds as a set, then to the
  double bonds specifically, before substituent locants (mirrored from
  `_unsaturated.py`'s acyclic tie-break).
- P-29.4 / P-46 (Chapter P-2, P-4) and P-35.2.1 (Chapter P-3): ring
  substituents (simple/compound alkyl groups, halogens) are named exactly
  as in `_cyclic.py`.
- P-91.2.2 (Chapter P-9): a specified ring C=C double bond gets a bare
  `(Z)-`/`(E)-` prefix, no locant (redundant for the same reason the
  bond's own locant is already omitted above), but only when the ring is
  large enough for the geometry to be a genuine stereogenic unit --
  confirmed directly from the primary source text: a 3-7-membered ring's
  double bond is always 'Z' (the only physically realizable form) and the
  descriptor is correspondingly always omitted, while an 8+-membered
  ring's requires it. `_common.specified_double_bond_stereo` already
  delegates this exact size distinction to RDKit's own stereo perception
  (`Chem.FindPotentialStereo` reports no element at all for the small-ring
  case regardless of input markers), so no ring-size branching is
  hand-coded here. A ring with two or more double bonds needing
  locant-bearing citation is a separate, unverified case, out of scope.

Scope: a single monocyclic, all-carbon ring bearing one or more
carbon-carbon ring double and/or triple bonds, plus any exocyclic
substituents `name_branch` can name -- including ylidene/ylidyne and
enyl/ynyl groups (P-29.2, P-32.1.1). No heteroatoms, no fused/bridged/spiro
combination, no aromaticity (benzene itself is `_aromatic.py`'s).
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    ring_cycle,
    specified_double_bond_stereo,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import format_substituent_prefixes, ring_and_branch_stereo_display, substituents_for_ring

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_VALID_ORDERS = (_ENE_ORDER, _YNE_ORDER)


def find_cyclic_unsaturated_core(mol):
    """Return the ring's atom indices if `mol` is a single monocyclic,
    all-carbon, non-aromatic ring carrying at least one ring C=C or C#C
    bond and no unsaturation anywhere else, else None (deliberately
    permissive about the invalid neighboring shapes this module itself
    rejects, e.g. an exocyclic double bond -- returning None there lets it
    fall through to `_cyclic.py`'s existing rejection instead, see module
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

    has_ring_multi_bond = False
    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if order not in _VALID_ORDERS:
            return None
        if a in ring_set and b in ring_set:
            has_ring_multi_bond = True
    if not has_ring_multi_bond:
        return None
    return list(ring_atoms)


def _multi_bonds(mol):
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() in _VALID_ORDERS
    ]


def _double_bonds(mol):
    return [(a, b) for a, b, order in _multi_bonds(mol) if order == _ENE_ORDER]


def _bond_locant(ring_order, bond_atoms):
    n = len(ring_order)
    position = {atom: i + 1 for i, atom in enumerate(ring_order)}
    pa, pb = position[bond_atoms[0]], position[bond_atoms[1]]
    return n if {pa, pb} == {1, n} else min(pa, pb)


def _ring_bond_locants(ring_order, double_bonds):
    return sorted(_bond_locant(ring_order, bond) for bond in double_bonds)


def _ring_multi_bond_locants(ring_order, bonds):
    """(ene_locants, yne_locants), both sorted, for every ring multiple
    bond under this ring ordering."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _bond_locant(ring_order, (a, b))
        (ene if order == _ENE_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)




def _unsaturation_suffix(ene_locants, yne_locants):
    """Locant-and-suffix string (e.g. '1,3-dien-5-yne') plus whether the
    parent stem needs its euphonic trailing 'a' -- mirrors
    `_unsaturated.py`'s helper of the same name; see P-31.1.1.1/
    P-31.1.3.1 for the shared 'ene' precedes 'yne' rule."""
    ene_count, yne_count = len(ene_locants), len(yne_locants)
    ene_word = multiplied_word(ene_count, "ene")
    yne_word = multiplied_word(yne_count, "yne")

    if ene_count and yne_count:
        ene_part = ene_word[:-1] if yne_word[0] in "aeiouy" else ene_word
        ene_loc_str = ",".join(str(loc) for loc in ene_locants)
        yne_loc_str = ",".join(str(loc) for loc in yne_locants)
        body = f"{ene_loc_str}-{ene_part}-{yne_loc_str}-{yne_word}"
    elif ene_count:
        body = f"{','.join(str(loc) for loc in ene_locants)}-{ene_word}"
    else:
        body = f"{','.join(str(loc) for loc in yne_locants)}-{yne_word}"

    needs_stem_a = (ene_count >= 2) if ene_count else (yne_count >= 2)
    return body, needs_stem_a


def _name_from_substituents(ring_size, ene_locants, yne_locants, grouped):
    parent_stem = "cyclo" + alkane_name(ring_size)[:-3]
    prefix = format_substituent_prefixes(grouped)
    total_count = len(ene_locants) + len(yne_locants)
    if total_count == 1:
        # P-14.3.3: always achievable at locant '1', so never cited (see
        # module docstring).
        return prefix + parent_stem + ("ene" if ene_locants else "yne")
    body, needs_stem_a = _unsaturation_suffix(ene_locants, yne_locants)
    return prefix + parent_stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(ring_size, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(ring_size, ene_locants, yne_locants, grouped)
    return combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def name_cyclic_unsaturated_yl(mol, ring_atoms, root) -> str:
    """Name a monocyclic all-carbon partially unsaturated ring, bearing no
    exocyclic substituent other than the free valence itself, as a '-yl'
    substituent group (P-29.2) with that free valence fixed at position 1
    -- e.g. 'cyclohex-3-en-1-yl'. Unlike the parent-hydride name
    (`name_cyclic_unsaturated`), a single ring double bond's locant is
    never redundant here: fixing the free valence at '1' removes the free
    choice of numbering start that otherwise always lets a lone double
    bond land on '1', so its locant must always be cited."""
    graph = adjacency(mol)
    double_bonds = _double_bonds(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    start = ring_order.index(root)
    rotated = ring_order[start:] + ring_order[:start]

    best_locants = None
    for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
        ene_locants = _ring_bond_locants(candidate, double_bonds)
        if best_locants is None or ene_locants < best_locants:
            best_locants = ene_locants

    parent_stem = "cyclo" + alkane_name(ring_size)[:-3]
    count = len(best_locants)
    if count == 1:
        return f"{parent_stem}-{best_locants[0]}-en-1-yl"
    word = numerical_term(count) + "en"
    loc_str = ",".join(str(loc) for loc in best_locants)
    return f"{parent_stem}a-{loc_str}-{word}-1-yl"


def name_cyclic_unsaturated(mol, ring_atoms) -> str:
    validate_atoms_and_bonds(mol)
    ring_set = set(ring_atoms)

    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if order not in _VALID_ORDERS:
            raise UnsupportedStructure(
                "only carbon-carbon double and triple bonds are "
                "supported in a partially unsaturated ring (see "
                "P-31.1.1.1)"
            )

    bonds = [(a, b, order) for a, b, order in _multi_bonds(mol) if a in ring_set and b in ring_set]

    # A specified tetrahedral stereocenter and a specified ring C=C E/Z
    # element are mutually exclusive shapes here (P-92 vs. P-93) --
    # `Chem.FindPotentialStereo`'s own element types decide which of
    # `_common`'s two dedicated checks applies, so neither one is called
    # on a molecule shaped for the other (each would otherwise reject the
    # other's shape outright as "beyond" its own kind, see #776).
    stereo_elements = Chem.FindPotentialStereo(mol)
    specified_elements = [e for e in stereo_elements if e.specified == Chem.StereoSpecified.Specified]
    tetrahedral_specified = specified_elements and all(
        e.type == Chem.StereoType.Atom_Tetrahedral for e in specified_elements
    )

    bond_stereo = None
    stereo = None
    if tetrahedral_specified:
        stereo = specified_stereocenters(mol)
    else:
        bond_stereo = specified_double_bond_stereo(mol)
        if bond_stereo is not None:
            # RDKit's own `Chem.FindPotentialStereo` already tells apart a
            # 3-7-membered ring's genuinely non-stereogenic double bond (no
            # element reported, `bond_stereo` stays None regardless of input
            # markers) from an 8+-membered ring's genuinely stereogenic one
            # (P-91.2.2) -- no ring-size branching needed here at all.
            if len(_multi_bonds(mol)) != len(bonds):
                raise UnsupportedStructure(
                    "a specified E/Z double bond alongside an exocyclic multiple bond is not supported yet (see P-93)"
                )
            if any(order == _YNE_ORDER for _, _, order in bonds):
                raise UnsupportedStructure(
                    "a specified double-bond E/Z stereo element combined with "
                    "a triple bond is not supported yet (see P-93)"
                )
            if len(bond_stereo) != 1 or len(bonds) != 1:
                raise UnsupportedStructure(
                    "a ring with more than one double bond, or one where not "
                    "every double bond is specified, is not supported yet -- "
                    "P-91.2.2's multi-bond citation rule needs separate "
                    "verification"
                )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    branch_stereo = None
    if stereo is not None and any(atom not in ring_set for atom, _ in stereo):
        branch_stereo = ring_and_branch_stereo_display(graph, ring_order, frozenset(), stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch, other than the "
                "ring's own attachment atom paired with it, is not "
                "supported yet (see P-92)"
            )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            ene_locants, yne_locants = _ring_multi_bond_locants(candidate, bonds)
            substituents = substituents_for_ring(graph, candidate, halogens, mol=mol, unsaturated=True)
            if branch_stereo is not None:
                branch_ring_atom, display, _ring_r_or_s = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            key = _candidate_key(ring_size, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if branch_stereo is not None:
        branch_ring_atom, _display, ring_r_or_s = branch_stereo
        if ring_r_or_s is not None:
            locant = best_position_of[branch_ring_atom]
            return f"({locant}{ring_r_or_s})-{best_name}"
        return best_name
    if stereo is not None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    if bond_stereo is not None:
        # P-91.2.2's own worked example ('(Z)-cyclooctene', '(E)-cyclooctene')
        # cites no locant -- redundant for the same reason the ring's own
        # double-bond locant is already omitted in the parent name
        # (P-14.3.3).
        ((_, code),) = bond_stereo
        return f"({code})-{best_name}"
    return best_name
