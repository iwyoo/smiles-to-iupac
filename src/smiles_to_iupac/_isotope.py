"""Naming of isotopically substituted compounds (the '(<nuclide>)' isotope
descriptor, P-82.2.1) for a mononuclear methane parent, and for an
unbranched multi-carbon alkane chain parent bearing a single isotopic
modification, per the IUPAC 2013 Recommendations ("the Blue Book"):

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
  + deuterium modification (P-82.3, confirmed distinct by
  '(2-14C,3-2H1)butane (PIN)') both stay out of scope for a follow-up task.
  Structure-verified via PubChem: `[2H]C([2H])C([2H])[2H]` (PubChem's own
  name: "1,1,2,2-tetradeuterioethane"), `[2H]C([2H])CC([2H])[2H]`
  (PubChem: "1,1,3,3-tetradeuteriopropane") -- PubChem uses its own
  systematic 'deuterio' prefix style rather than the Blue Book's
  parenthesized nuclide descriptor, but the structures match exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any parent other than an unbranched methane/alkane chain -- branched
  chains, rings, and any parent other than a plain hydrocarbon are not
  supported.
- Any isotope other than deuterium (2H) among hydrogen atoms (e.g. tritium).
- The carbon isotope and deuterium substitution occurring at the same time
  (P-82.3), or any isotopically labeled halogen.
- More than one distinct chain position bearing a skeletal carbon isotope.
- A 2-carbon chain whose isotope-locant citation isn't settled by a
  confirmed worked example (see above) -- this only applies to a single
  deuterium-bearing position; two or more positions always cite locants.
- Any heteroatom (other than halogens), charge, or additional isotopic
  modification.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
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
    if carbon_isotope != 0 and deuterium_count:
        raise UnsupportedStructure(
            "a skeletal carbon isotope combined with deuterium substitution "
            "is not supported yet (P-82.3)"
        )
    if deuterium_count > _MAX_DEUTERIUMS:
        raise UnsupportedStructure("methane cannot carry more than 4 deuterium atoms")

    if carbon_isotope != 0:
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

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in chain_set]
    bonded_carbon = {}
    for atom in other_atoms:
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(neighbors) != 1 or neighbors[0] not in chain_set:
            raise UnsupportedStructure(
                "every substituent must be a single leaf atom bonded directly "
                "to a chain carbon (unbranched chains only)"
            )
        bonded_carbon[atom.GetIdx()] = neighbors[0]
    _validate_other_atoms(other_atoms, mol, lambda atom: bonded_carbon[atom.GetIdx()])

    for carbon_idx in chain:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        if carbon.GetFormalCharge() != 0 or carbon.GetIsAromatic():
            raise UnsupportedStructure("an unsupported carbon shape for this chain")

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
            raise UnsupportedStructure("a single chain position cannot carry more than 4 deuterium atoms")
    deuterium_count = len(hydrogens)

    carbon_isotope_positions = [c for c in chain if mol.GetAtomWithIdx(c).GetIsotope() != 0]
    for c in carbon_isotope_positions:
        isotope = mol.GetAtomWithIdx(c).GetIsotope()
        if isotope not in _CARBON_ISOTOPES:
            raise UnsupportedStructure(
                "only 12C/13C/14C skeletal carbon isotopes are supported (P-82.2.1)"
            )
    if len(carbon_isotope_positions) > 1:
        raise UnsupportedStructure(
            "more than one isotopically labeled skeletal carbon (a locant "
            "set) is not supported yet"
        )
    if carbon_isotope_positions and deuterium_count:
        raise UnsupportedStructure(
            "a skeletal carbon isotope combined with deuterium substitution "
            "is not supported yet (P-82.3)"
        )

    multi_position_deuterium = len(deuterium_counts) > 1

    if carbon_isotope_positions:
        isotope_carbon = carbon_isotope_positions[0]
        isotope_descriptor_base = f"{mol.GetAtomWithIdx(isotope_carbon).GetIsotope()}C"
    elif deuterium_count and not multi_position_deuterium:
        (isotope_carbon,) = deuterium_counts
        isotope_descriptor_base = f"2H{deuterium_count}"
    elif not deuterium_count:
        raise UnsupportedStructure("no isotopically labeled atom found")

    halogens = [atom for atom in other_atoms if atom.GetAtomicNum() in HALOGEN_PREFIXES]

    best_key = None
    best = None
    for candidate in (chain, list(reversed(chain))):
        halogen_positions = [
            (candidate.index(bonded_carbon[h.GetIdx()]) + 1, h.GetAtomicNum()) for h in halogens
        ]
        if multi_position_deuterium:
            deuterium_locants = sorted(
                candidate.index(carbon_idx) + 1
                for carbon_idx, count in deuterium_counts.items()
                for _ in range(count)
            )
            combined = lowest_locant_set(deuterium_locants + [pos for pos, _ in halogen_positions])
            candidate_result = (deuterium_locants, halogen_positions)
        else:
            isotope_locant = candidate.index(isotope_carbon) + 1
            combined = lowest_locant_set([isotope_locant] + [pos for pos, _ in halogen_positions])
            candidate_result = (isotope_locant, halogen_positions)
        if best_key is None or combined < best_key:
            best_key = combined
            best = candidate_result

    chain_length = len(chain)

    if multi_position_deuterium:
        deuterium_locants, halogen_positions = best
        isotope_descriptor = f"{','.join(str(loc) for loc in deuterium_locants)}-2H{deuterium_count}"
    else:
        isotope_locant, halogen_positions = best
        halogen_count = len(halogen_positions)

        if chain_length >= 3:
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

        isotope_descriptor = (
            f"{isotope_locant}-{isotope_descriptor_base}" if cite_locant else isotope_descriptor_base
        )

    halogen_prefix = _halogen_prefix_for_positions(halogen_positions)

    return f"{halogen_prefix}({isotope_descriptor}){alkane_name(chain_length)}"


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
    return _name_chain(mol, carbons)
