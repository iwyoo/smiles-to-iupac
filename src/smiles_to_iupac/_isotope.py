"""Naming of isotopically substituted compounds (the '(<nuclide>)' isotope
descriptor, P-82.2.1) for a mononuclear methane parent, an unbranched
multi-carbon alkane chain parent, and a monocyclic saturated all-carbon
ring parent, each bearing one or more isotopic modifications, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-82.2.1 (Chapter P-8, https://iupac.qmul.ac.uk/BlueBook/P8.html): the name
  of an isotopically substituted compound is formed by adding the nuclide
  symbol(s), enclosed in parentheses, preceded by any necessary locant(s),
  before the part of the name that is isotopically substituted. Confirmed
  worked examples: CH3-2H -> '(2H1)methane (PIN)'; (14C)methane (PIN);
  trichloro(12C)methane (PIN); dichloro(2H2)methane (PIN). The isotope's
  count is always cited as a right subscript to its symbol, even for a
  single atom ('2H1', never a bare '2H') -- "the number of atoms
  substituted is always specified as a right subscript ... even in case of
  monosubstitution".
- Methane (a mononuclear parent) never needs a locant, the same way
  `_acyclic.py`'s P-14.3.4.2(a) never cites a locant on a single-carbon
  parent (e.g. 'chloromethane').
- Multi-carbon chains (P-82.6 Locants, `tmp/bluebook/P8.txt` lines 510-579,
  and P-82.5.2's numbering-priority rule, same file lines ~465-490):
  locants ARE cited whenever needed to specify the modified position.
  Confirmed worked examples: '(2-14C)butane (PIN)' [not '(3-14C)butane' --
  P-82.5.2's own numbering-priority example, an unhalogenated, otherwise
  plain butane chain] and '1,1,1-trifluoro(2-2H1)ethane (PIN)' [a
  halogenated ethane chain]. P-82.5.2 states directly: "the starting point
  and the direction of numbering ... are chosen so as to give lowest
  locants to the modified atoms or groups considered together in one
  series" -- i.e. the isotope locant and any halogen substituent locants
  are minimized together as one combined series when picking chain
  numbering direction, exactly the way `_acyclic.py`/`_selenoate.py`
  already pick numbering direction for their own substituent locants.
- This module deliberately stays conservative about locant *omission* on a
  chain: only two shapes are confirmed by the worked examples above (an
  unhalogenated chain of 3+ carbons, and a halogenated 2-carbon chain whose
  halogen count already forces citation on its own, i.e. more than one
  halogen substituent instance). A 2-carbon chain with a single isotope
  modification and no halogens (e.g. plain deuterioethane) or with exactly
  one halogen substituent is NOT covered by any confirmed worked example
  here -- P-82.6.1.1's own text ("locants are omitted if no locants are
  necessary in unmodified names") suggests these might also omit the
  locant, but without a worked example to confirm the exact omitted form,
  this module raises `UnsupportedStructure` for them rather than guessing.
- Within a chain, deuterium may now be spread across more than one chain
  carbon (a locant *set*, mirroring the '(1,1,1,3,3-2H5)pentan-2-one'-style
  names cited by P-82.3, though this module still only handles a plain
  alkane parent, not a ketone), each position independently carrying 1-4
  deuterium atoms, always with its own locant repeated once per atom --
  the same citation style `_halogen_prefix_for_positions` already uses for
  multiple halogens at one position. A single labeled skeletal carbon
  (12C/13C/14C) is still the only carbon-isotope shape supported; more
  than one such carbon (its own locant set) and simultaneous carbon-isotope
  + deuterium modification on a chain (P-82.3) both stay out of scope for a
  follow-up task. Structure-verified via PubChem: `[2H]C([2H])C([2H])[2H]`
  (PubChem's own name: "1,1,2,2-tetradeuterioethane"),
  `[2H]C([2H])CC([2H])[2H]` (PubChem: "1,1,3,3-tetradeuteriopropane") --
  PubChem uses its own systematic 'deuterio' prefix style rather than the
  Blue Book's parenthesized nuclide descriptor, but the structures match
  exactly.
- P-82.3.1 (order of nuclide symbols): a mononuclear methane parent may now
  carry a skeletal carbon isotope (12C/13C/14C) and deuterium substitution
  at the same time -- "when isotopes of different elements are present ...
  their symbols are arranged in alphabetical order", confirmed by the
  worked example `methan(2H,18O)ol (PIN)` (a different element pair, but
  the same alphabetical-ordering rule applies to carbon+hydrogen: 'C'
  before 'H', e.g. '(13C,2H1)methane'). Structure-verified via PubChem:
  `[13CH3][2H]` -> CID 10949861, PubChem's own name
  "deuterio(113C)methane" (PubChem's systematic style, not a PIN, but
  confirms the combined structure is valid and distinct from either
  isotope alone).
- The same carbon-isotope + deuterium combination is now also supported on
  a multi-carbon chain (P-82.5.2, `tmp/bluebook/P8.txt` lines 466-490),
  confirmed by two worked examples: `(3-14C,2,2-2H2)butane (PIN)` [not
  `(2-14C,3,3-2H2)butane`] and `(2-14C,3-2H1)butane (PIN)` [not
  `(3-14C,2-2H1)butane`]. P-82.5.2 numbers the chain to give lowest
  locants to "the modified atoms or groups considered together in one
  series" -- i.e. the carbon-isotope locant and the deuterium locant(s)
  are minimized together as a single combined series (alongside any
  halogen locants), the same way this module already combines
  multi-position deuterium locants with halogen locants. Locants are
  always cited for this combination (never omitted) -- both worked
  examples cite them even though butane's own unmodified name needs none.
- A skeletal carbon isotope may now also be spread across more than one
  chain carbon (its own locant set), mirroring the deuterium locant-set
  rule above, as long as every labeled carbon shares the same nuclide
  (mixing e.g. 13C and 14C on the same chain stays out of scope). The
  isotope count is cited as a subscript on the nuclide symbol only when
  more than one carbon carries it (e.g. '13C2'), matching the
  single-carbon form's bare 'nC' (no '1' subscript). Structure-verified
  via PubChem: `[13CH3]CC[13CH3]` -> "(1,4-13C2)butane" (CID 13378975),
  `[13CH3][13CH2]CC` -> "(1,2-13C2)butane" (CID 90969531),
  `[13CH3]C[13CH2]C` -> "(1,3-13C2)butane" (CID 13378977). It combines
  with deuterium/halogen locants using the same combined-locant-set
  numbering rule as everything else in this module.

- P-82.6.1.1 ("locants are omitted if no locants are necessary in
  unmodified names") also covers a monocyclic saturated all-carbon ring
  (`_name_ring`, mirroring `_name_chain`'s combined-locant-set numbering
  but trying every rotation/direction of the ring, since a ring has no
  fixed terminus): unlike a 3+-carbon chain (whose interior vs. terminal
  positions are chemically distinct even before any modification, so a
  locant is still needed to say *which* one), every position of a plain,
  otherwise-unsubstituted monocyclic ring is equivalent to every other, so
  a *single* isotopically modified position needs no locant at all --
  e.g. `[2H]C1CCCCC1` -> `(2H1)cyclohexane` (no '1-' prefix). Two or more
  distinct modified ring positions, or any coexisting halogen, still cite
  locants the same way the chain case does (the modification breaks the
  ring's symmetry, so a locant becomes necessary to specify the pattern).
  Deuterium and a skeletal carbon isotope both modifying the *same* single
  ring position (with no halogen) is treated the same way as an
  unsubstituted single position for this omission rule (extending the
  mononuclear methane case's own precedent, P-82.3.1's alphabetical
  '(13C,2H1)methane' ordering, to a monocyclic ring) -- not itself
  confirmed by a registered structure, so this module raises
  `UnsupportedStructure` for that specific untested combination rather
  than guessing at its exact formatting.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any parent other than an unbranched methane/alkane chain or a
  monocyclic saturated all-carbon ring -- branched chains, fused/spiro/
  polycyclic rings, an aromatic ring, and any parent other than a plain
  hydrocarbon are not supported.
- Any isotope other than deuterium (2H) among hydrogen atoms (e.g. tritium).
- Any isotopically labeled halogen.
- Mixing different skeletal carbon isotope nuclides (e.g. 13C and 14C
  together) on the same chain/ring.
- A 2-carbon chain whose isotope-locant citation isn't settled by a
  confirmed worked example (see above) -- this only applies to a single
  deuterium-bearing position; two or more positions always cite locants.
- A ring carbon isotope and deuterium both modifying the exact same single
  ring position, with no halogen (see above) -- unconfirmed formatting.
- Any heteroatom (other than halogens), charge, or additional isotopic
  modification.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    longest_chains,
    lowest_locant_set,
    ring_cycle,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import format_substituent_prefixes

_CARBON_ISOTOPES = {12, 13, 14}
_MAX_DEUTERIUMS = 4


def has_isotope_shape(mol) -> bool:
    """True if the molecule contains any isotopically labeled atom at all
    (deuterium, a heavy carbon isotope, or otherwise), regardless of
    whether the rest of the molecule is in scope. Used by `core.py` to
    route here before every other branch, none of which recognize an
    explicit, isotopically labeled atom (RDKit represents 2H as its own
    atom, atomic number 1) at all."""
    return any(atom.GetIsotope() != 0 for atom in mol.GetAtoms())


def _validate_other_atoms(other_atoms, mol, bonded_carbon_idx):
    for atom in other_atoms:
        if atom.GetAtomicNum() not in {1, *HALOGEN_PREFIXES}:
            raise UnsupportedStructure(
                "this module only supports deuterium and halogen substituents; "
                "other atoms are not supported yet (P-82.2.1)"
            )
        bond = mol.GetBondBetweenAtoms(bonded_carbon_idx(atom), atom.GetIdx())
        if bond is None or bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("substituents must be singly bonded")
        if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetIsotope() != 0:
            raise UnsupportedStructure("an isotopically labeled halogen is not supported yet")


def _halogen_prefix_for_positions(halogen_positions):
    substituents: dict[int, list] = {}
    for position, atomic_num in halogen_positions:
        substituents.setdefault(position, []).append((HALOGEN_PREFIXES[atomic_num], False))
    grouped = group_substituents(substituents)
    return format_substituent_prefixes(grouped)


def _name_methane(mol, carbon) -> str:
    if carbon.GetFormalCharge() != 0 or carbon.GetIsAromatic():
        raise UnsupportedStructure("an unsupported carbon shape for methane")

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() != carbon.GetIdx()]
    _validate_other_atoms(other_atoms, mol, lambda atom: carbon.GetIdx())

    hydrogens = [atom for atom in other_atoms if atom.GetAtomicNum() == 1]
    non_deuteriums = [atom for atom in hydrogens if atom.GetIsotope() != 2]
    if non_deuteriums:
        raise UnsupportedStructure(
            "only deuterium (2H) is supported among hydrogen isotopes; "
            "other hydrogen isotopes (e.g. tritium) are not supported yet"
        )
    deuterium_count = len(hydrogens)

    carbon_isotope = carbon.GetIsotope()
    if carbon_isotope != 0 and carbon_isotope not in _CARBON_ISOTOPES:
        raise UnsupportedStructure(
            "only 12C/13C/14C skeletal carbon isotopes are supported (P-82.2.1)"
        )
    if deuterium_count > _MAX_DEUTERIUMS:
        raise UnsupportedStructure("methane cannot carry more than 4 deuterium atoms")

    if carbon_isotope != 0 and deuterium_count:
        isotope_descriptor = f"{carbon_isotope}C,2H{deuterium_count}"
    elif carbon_isotope != 0:
        isotope_descriptor = f"{carbon_isotope}C"
    elif deuterium_count:
        isotope_descriptor = f"2H{deuterium_count}"
    else:
        raise UnsupportedStructure("no isotopically labeled atom found")

    halogens = [atom for atom in other_atoms if atom.GetAtomicNum() in HALOGEN_PREFIXES]
    halogen_counts: dict[str, int] = {}
    for atom in halogens:
        prefix = HALOGEN_PREFIXES[atom.GetAtomicNum()]
        halogen_counts[prefix] = halogen_counts.get(prefix, 0) + 1
    halogen_prefix = "".join(
        prefix if count == 1 else numerical_term(count) + prefix
        for prefix, count in sorted(halogen_counts.items())
    )

    return f"{halogen_prefix}({isotope_descriptor})methane"


def _collect_substituents(mol, skeleton_set, shape_label):
    """{other_atom_idx -> its single skeleton-carbon neighbor's idx} after
    checking every non-skeleton atom is a single leaf bonded directly to a
    skeleton carbon -- shared by `_name_chain` (skeleton = the chain) and
    `_name_ring` (skeleton = the ring)."""
    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in skeleton_set]
    bonded_carbon = {}
    for atom in other_atoms:
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(neighbors) != 1 or neighbors[0] not in skeleton_set:
            raise UnsupportedStructure(
                f"every substituent must be a single leaf atom bonded directly "
                f"to a {shape_label} carbon (unbranched chains/plain monocyclic "
                f"rings only)"
            )
        bonded_carbon[atom.GetIdx()] = neighbors[0]
    _validate_other_atoms(other_atoms, mol, lambda atom: bonded_carbon[atom.GetIdx()])
    return other_atoms, bonded_carbon


def _collect_isotopes(mol, skeleton, other_atoms, bonded_carbon, shape_label):
    """(deuterium_counts, carbon_isotope_positions, carbon_isotope_symbol,
    has_carbon_isotope, has_deuterium, deuterium_count, multi_position_deuterium)
    after validating every skeleton carbon's own shape and every isotope's
    kind/count -- shared by `_name_chain`/`_name_ring`."""
    for carbon_idx in skeleton:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        if carbon.GetFormalCharge() != 0 or carbon.GetIsAromatic():
            raise UnsupportedStructure(f"an unsupported carbon shape for this {shape_label}")

    hydrogens = [atom for atom in other_atoms if atom.GetAtomicNum() == 1]
    non_deuteriums = [atom for atom in hydrogens if atom.GetIsotope() != 2]
    if non_deuteriums:
        raise UnsupportedStructure(
            "only deuterium (2H) is supported among hydrogen isotopes; "
            "other hydrogen isotopes (e.g. tritium) are not supported yet"
        )
    deuterium_counts: dict[int, int] = {}
    for atom in hydrogens:
        carbon_idx = bonded_carbon[atom.GetIdx()]
        deuterium_counts[carbon_idx] = deuterium_counts.get(carbon_idx, 0) + 1
    for count in deuterium_counts.values():
        if count > _MAX_DEUTERIUMS:
            raise UnsupportedStructure(f"a single {shape_label} position cannot carry more than 4 deuterium atoms")
    deuterium_count = len(hydrogens)

    carbon_isotope_positions = [c for c in skeleton if mol.GetAtomWithIdx(c).GetIsotope() != 0]
    carbon_isotopes = {mol.GetAtomWithIdx(c).GetIsotope() for c in carbon_isotope_positions}
    for isotope in carbon_isotopes:
        if isotope not in _CARBON_ISOTOPES:
            raise UnsupportedStructure(
                "only 12C/13C/14C skeletal carbon isotopes are supported (P-82.2.1)"
            )
    if len(carbon_isotopes) > 1:
        raise UnsupportedStructure(
            f"mixing different skeletal carbon isotope nuclides (e.g. 13C "
            f"and 14C together) is not supported yet on the same {shape_label}"
        )
    has_carbon_isotope = bool(carbon_isotope_positions)
    has_deuterium = deuterium_count > 0
    if not has_carbon_isotope and not has_deuterium:
        raise UnsupportedStructure("no isotopically labeled atom found")

    multi_position_deuterium = len(deuterium_counts) > 1
    carbon_isotope_count = len(carbon_isotope_positions)
    carbon_isotope_nuclide = next(iter(carbon_isotopes)) if has_carbon_isotope else None
    carbon_isotope_symbol = (
        f"{carbon_isotope_nuclide}C" + (str(carbon_isotope_count) if carbon_isotope_count > 1 else "")
        if has_carbon_isotope
        else None
    )
    return (
        deuterium_counts,
        carbon_isotope_positions,
        carbon_isotope_symbol,
        has_carbon_isotope,
        has_deuterium,
        deuterium_count,
        multi_position_deuterium,
    )


def _pick_best_locants(candidates, bonded_carbon, deuterium_counts, carbon_isotope_positions, halogens):
    """Try every candidate atom-order list (a chain's two directions, or a
    ring's every rotation x 2 directions) and return the
    (carbon_isotope_locants, deuterium_locants, halogen_positions) whose
    combined locant set (P-82.5.2: isotope and halogen locants minimized
    together as one series) is lowest."""
    best_key = None
    best = None
    for candidate in candidates:
        halogen_positions = [
            (candidate.index(bonded_carbon[h.GetIdx()]) + 1, h.GetAtomicNum()) for h in halogens
        ]
        deuterium_locants = sorted(
            candidate.index(carbon_idx) + 1
            for carbon_idx, count in deuterium_counts.items()
            for _ in range(count)
        )
        carbon_isotope_locants = sorted(candidate.index(c) + 1 for c in carbon_isotope_positions)
        combined = lowest_locant_set(
            carbon_isotope_locants + deuterium_locants + [pos for pos, _ in halogen_positions]
        )
        candidate_result = (carbon_isotope_locants, deuterium_locants, halogen_positions)
        if best_key is None or combined < best_key:
            best_key = combined
            best = candidate_result
    return best


def _name_chain(mol, carbons) -> str:
    carbon_graph = carbon_adjacency(mol)
    chains = longest_chains(carbon_graph)
    chain = chains[0]
    if len(chain) != len(carbons):
        raise UnsupportedStructure(
            "a branched carbon skeleton is not supported yet (unbranched "
            "chains only)"
        )
    chain_set = set(chain)

    other_atoms, bonded_carbon = _collect_substituents(mol, chain_set, "chain")
    (
        deuterium_counts,
        carbon_isotope_positions,
        carbon_isotope_symbol,
        has_carbon_isotope,
        has_deuterium,
        deuterium_count,
        multi_position_deuterium,
    ) = _collect_isotopes(mol, chain, other_atoms, bonded_carbon, "chain")

    halogens = [atom for atom in other_atoms if atom.GetAtomicNum() in HALOGEN_PREFIXES]

    carbon_isotope_locants, deuterium_locants, halogen_positions = _pick_best_locants(
        (chain, list(reversed(chain))), bonded_carbon, deuterium_counts, carbon_isotope_positions, halogens
    )

    chain_length = len(chain)
    halogen_count = len(halogen_positions)

    descriptor_parts = []
    if has_carbon_isotope:
        locants_str = ",".join(str(loc) for loc in carbon_isotope_locants)
        descriptor_parts.append(f"{locants_str}-{carbon_isotope_symbol}")

    if has_deuterium:
        deuterium_symbol = f"2H{deuterium_count}"
        if has_carbon_isotope or multi_position_deuterium:
            cite_locant = True
        elif chain_length >= 3:
            cite_locant = True
        elif halogen_count == 0:
            raise UnsupportedStructure(
                "an unhalogenated 2-carbon chain with a single isotopic "
                "modification is not covered by a confirmed worked example yet "
                "(P-82.6.1.1's locant-omission rule for this shape is unconfirmed)"
            )
        elif halogen_count == 1:
            raise UnsupportedStructure(
                "a 2-carbon chain with exactly one halogen substituent alongside "
                "an isotopic modification is not covered by a confirmed worked "
                "example yet"
            )
        else:
            cite_locant = True

        descriptor_parts.append(
            f"{','.join(str(loc) for loc in deuterium_locants)}-{deuterium_symbol}"
            if cite_locant
            else deuterium_symbol
        )

    isotope_descriptor = ",".join(descriptor_parts)

    halogen_prefix = _halogen_prefix_for_positions(halogen_positions)

    return f"{halogen_prefix}({isotope_descriptor}){alkane_name(chain_length)}"


def _name_ring(mol, ring_atoms) -> str:
    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_set = set(ring_order)
    ring_size = len(ring_order)

    other_atoms, bonded_carbon = _collect_substituents(mol, ring_set, "ring")
    (
        deuterium_counts,
        carbon_isotope_positions,
        carbon_isotope_symbol,
        has_carbon_isotope,
        has_deuterium,
        deuterium_count,
        multi_position_deuterium,
    ) = _collect_isotopes(mol, ring_order, other_atoms, bonded_carbon, "ring")

    halogens = [atom for atom in other_atoms if atom.GetAtomicNum() in HALOGEN_PREFIXES]

    candidates = []
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        candidates.append(rotated)
        candidates.append(list(reversed(rotated)))
    carbon_isotope_locants, deuterium_locants, halogen_positions = _pick_best_locants(
        candidates, bonded_carbon, deuterium_counts, carbon_isotope_positions, halogens
    )
    halogen_count = len(halogen_positions)

    # P-82.6.1.1: a monocyclic ring's positions are all equivalent to one
    # another before any modification (unlike a 3+-carbon chain's own
    # chemically distinct interior/terminal positions), so a *single*
    # modified position needs no locant at all -- see module docstring.
    modified_positions = set(carbon_isotope_locants) | set(deuterium_locants)
    single_position_no_halogen = len(modified_positions) <= 1 and halogen_count == 0
    if single_position_no_halogen and has_carbon_isotope and has_deuterium:
        raise UnsupportedStructure(
            "a skeletal carbon isotope and deuterium both modifying the "
            "exact same single ring position, with no halogen, is not "
            "covered by a confirmed worked example yet (whether the "
            "locant is omitted here, extending the mononuclear methane "
            "precedent to a ring, is unconfirmed)"
        )

    descriptor_parts = []
    if has_carbon_isotope:
        cite_locant = not single_position_no_halogen
        descriptor_parts.append(
            f"{','.join(str(loc) for loc in carbon_isotope_locants)}-{carbon_isotope_symbol}"
            if cite_locant
            else carbon_isotope_symbol
        )

    if has_deuterium:
        deuterium_symbol = f"2H{deuterium_count}"
        cite_locant = not single_position_no_halogen
        descriptor_parts.append(
            f"{','.join(str(loc) for loc in deuterium_locants)}-{deuterium_symbol}"
            if cite_locant
            else deuterium_symbol
        )

    isotope_descriptor = ",".join(descriptor_parts)
    halogen_prefix = _halogen_prefix_for_positions(halogen_positions)

    return f"{halogen_prefix}({isotope_descriptor})cyclo{alkane_name(ring_size)}"


def name_isotope(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbons:
        raise UnsupportedStructure(
            "this module only supports a hydrocarbon parent hydride; "
            "other parent hydrides are not supported yet (P-82.2.1)"
        )

    if len(carbons) == 1:
        return _name_methane(mol, carbons[0])

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = ring_info.AtomRings()[0]
        ring_bond_indices = ring_info.BondRings()[0]
        is_plain_saturated_ring = set(ring_atoms) == {c.GetIdx() for c in carbons} and all(
            mol.GetBondWithIdx(b).GetBondTypeAsDouble() == 1.0 for b in ring_bond_indices
        )
        if is_plain_saturated_ring:
            return _name_ring(mol, ring_atoms)

    return _name_chain(mol, carbons)
