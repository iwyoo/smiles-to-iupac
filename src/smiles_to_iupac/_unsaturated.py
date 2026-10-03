"""Naming of acyclic hydrocarbons containing one or more carbon-carbon double
and/or triple bonds (alkenes, alkynes, and mixed enynes) on the principal
chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-31.1.1.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): the
  presence of a double or triple bond in an otherwise saturated chain is
  denoted by changing the 'ane' ending to 'ene' or 'yne'. Locants as low as
  possible are given to multiple bonds as a set, even though this may at
  times give 'yne' endings lower locants than 'ene' endings; if a choice
  remains, preference for low locants is given to the double bonds. In
  names, 'ene' always precedes 'yne', with elision of the final letter 'e'
  in 'ene' when directly followed by 'yne' (e.g. 'pent-3-en-1-yne (PIN)',
  'pent-1-en-4-yne (PIN)').
- P-31.1.1.2 (Chapter P-3): the multiplying prefixes 'di', 'tri', etc. are
  placed before 'ene'/'yne' to indicate the number of multiple bonds of each
  kind. For euphonic reasons, when 'ene'/'yne' are preceded by a multiplying
  prefix and a locant, the letter 'a' is inserted after the parent stem,
  e.g. 'buta-1,3-diene (PIN)', 'nona-1,3,5,7-tetraene (PIN)'. There is no
  elision of a multiplying prefix's own final 'a' before 'ene'/'yne'
  ('tetraene', not 'tetrene').
- P-31.1.2.2.1 (Chapter P-3): systematic acyclic examples combining the two
  rules above, e.g. 'hexa-2,3-diene (PIN)', 'hexa-1,3-dien-5-yne (PIN)'.
- P-31.1.2.1 (Chapter P-3): 'acetylene' is the retained preferred IUPAC name
  for unsubstituted HC#CH (substitution of any kind is disallowed for this
  retained name); once substituted, chain-extended, or combined with a
  second multiple bond, the systematic '-yne' name is used instead.
- P-14.3.4.2(d) / P-14.3.3 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  the locant '1' is omitted only for unsubstituted, two- or three-carbon
  alkenes/alkynes carrying a single multiple bond ('ethene (PIN)',
  'acetylene (PIN)', 'propene (PIN)', 'propyne (PIN)'); with a substituent,
  a longer chain, or a second multiple bond, the chain is always long enough
  that more than one position is possible, so the locant(s) are always
  essential and this omission never applies.
- P-44.3.2 / P-44.4.1.1 / P-44.4.1.2 (Chapter P-4,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): chain length is chosen first
  (as in `_acyclic.py`), then the chain with the greater number of multiple
  bonds, then of double bonds. A multiple bond left off the chain is cited
  inside a substituent prefix -- 'ethenyl', 'prop-1-en-2-yl', 'methylidene',
  ... -- via `name_branch` (P-32.1.1, P-29.2), e.g. '3-methylidenepentane'.
- P-14.4(e) / P-44.4.1.10, P-44.4.1.10.1 (Chapter P-1 / P-4): numbering
  direction is chosen to give the lowest locants to the full set of multiple
  bonds (ene and yne together) ahead of substituent locants; if a choice
  still remains, lower locants go to the double ('ene') bonds specifically
  (same rule restated for the naming/citation format in P-31.1.1.1 above).
- P-45.2 (Chapter P-4) / P-14.4(f,g) (Chapter P-1): any remaining tie is
  broken exactly as in `_acyclic.py`: by substituent count, then lowest
  locant set, then lowest locants in citation order.
- P-29.4 / P-46 (Chapter P-2, P-4): substituents on the chain are named the
  same way as for alkanes/cycloalkanes, via `_substituents.py`, so simple
  *and* branched ("compound") substituents (e.g. isopropyl) are supported
  here too.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  are never skeletal atoms (P-44.3), so the principal chain is found over
  carbon-carbon connectivity only (`carbon_adjacency`, see `_common.py`)
  while substituent detection still uses the full atom graph.

Unsaturation in a ring is out of scope and raises `UnsupportedStructure`.

- P-91.3 / P-93 (Chapter P-9, https://iupac.qmul.ac.uk/BlueBook/P9.html):
  when a C=C double bond's
  geometry is specified in the input (`/`/`\`), a locanted "(nE)-"/"(nZ)-"
  prefix is added to the whole name, e.g. "(2E)-but-2-ene",
  "(2Z)-2-chlorobut-2-ene" -- the primary source's own worked example is
  "(2Z)-but-2-ene (PIN)" (P-91.3), and P-91.3 states plainly that "in
  preferred IUPAC names, stereodescriptors, preceded by a locant, must be
  cited to specify each stereogenic unit" for acyclic chains; the locant is
  never dropped just because a single double bond leaves nothing to
  disambiguate (that omission is reserved for specific ring systems --
  three- to seven-membered unsaturated alicyclics, von Baeyer, spiro, fused,
  and cyclophane systems -- per P-91.2.2, none of which apply here). CIP
  priority computation is delegated entirely to RDKit
  (`_common.specified_double_bond_stereo`), mirroring `_alcohol.py`'s R/S
  handling: a non-stereogenic double bond, or one left unspecified in the
  input, is not a new rejection case -- it's named exactly as before (no
  prefix).
- Generalizing the above the same way R/S was generalized from one
  stereocenter to many: two or more C=C double bonds, *all*
  specified, are cited together in one parenthesized group, ascending
  locant order, e.g. "(2E,4E)-hexa-2,4-diene", "(2Z,4E)-hexa-2,4-diene"
  (both confirmed via PubChem PUG REST). A specified double bond alongside a
  triple bond, a partially-specified set of double bonds (some with a
  slash marker, some left as a plain double bond), or a tetrahedral
  stereocenter is out of scope and raises `UnsupportedStructure`
  explicitly.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    bond_locant,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    specified_double_bond_stereo,
    substituent_locant_set_and_citation,
    unsaturation_suffix,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_VALID_ORDERS = (_ENE_ORDER, _YNE_ORDER)

def _unsaturation_suffix_from_citations(ene_citations, yne_citations):
    """Like `unsaturation_suffix` (`_common.py`), but each locant is a (primary_locant,
    display) pair (`_common.von_baeyer_unsaturation_citations`) instead of
    a plain integer -- P-31.1.4.2(1)'s compound-locant display (e.g.
    '1(7)') for a von Baeyer bicyclic/polycyclic parent, sorted by each
    bond's own primary (never-parenthesized) locant."""
    ene_citations = sorted(ene_citations)
    yne_citations = sorted(yne_citations)
    ene_count, yne_count = len(ene_citations), len(yne_citations)
    ene_word = multiplied_word(ene_count, "ene")
    yne_word = multiplied_word(yne_count, "yne")

    ene_str = ",".join(display for _, display in ene_citations)
    yne_str = ",".join(display for _, display in yne_citations)
    if ene_count and yne_count:
        body = f"{ene_str}-{ene_word[:-1]}-{yne_str}-{yne_word}"
    elif ene_count:
        body = f"{ene_str}-{ene_word}"
    else:
        body = f"{yne_str}-{yne_word}"

    needs_stem_a = (ene_count >= 2) if ene_count else (yne_count >= 2)
    return body, needs_stem_a


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    prefix = format_substituent_prefixes(grouped)
    if not ene_locants and not yne_locants:
        return prefix + alkane_name(chain_length)
    stem = alkane_name(chain_length)[:-3]
    single_bond = len(ene_locants) + len(yne_locants) == 1
    if single_bond and chain_length <= 3 and not prefix:
        # P-14.3.4.2(d): the locant is omittable only when unsubstituted, a
        # single multiple bond, and the chain is short enough that no other
        # position is possible.
        if chain_length == 2 and yne_locants:
            return "acetylene"
        return stem + ("ene" if ene_locants else "yne")
    if single_bond and chain_length == 2:
        suffix = "ene" if ene_locants else "yne"
        if len(grouped) == 1:
            name = next(iter(grouped))
            info = grouped[name]
            if len(info["locants"]) == 1 and not info["compound"]:
                # P-14.3.4.2(b): a homogeneous two-carbon chain bearing
                # exactly one substituent has only one possible structure
                # regardless of numbering direction, so both the multiple
                # bond's and the substituent's locants are omittable, e.g.
                # 'fluoroethyne (PIN)' for fluoroacetylene (P-31.1.2.1).
                return name + stem + suffix
        # A two-carbon chain has only one possible bond position (the
        # C1=C2/C1#C2 pair), so the bond's own locant is always omittable
        # here -- even with 2+ substituents, whose own locants (already
        # baked into `prefix`) are still needed to distinguish isomers
        # like "1,2-" from "1,1-".
        return prefix + stem + suffix
    body, needs_stem_a = unsaturation_suffix(ene_locants, yne_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-14.4(e)/P-44.4.1.10 then P-45.2.1-P-45.2.3,
    most-preferred first."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
    # The multiple bonds' locants (as a set, then 'ene' locants specifically)
    # outrank every substituent criterion; substituent count and lower
    # locants are preferred, so negate the count to sort every field in
    # ascending "most preferred first" order.
    return (
        (
            -len(combined_locant_set),
            -len(ene_locants),
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _bond_atoms(mol, bond_idx):
    bond = mol.GetBondWithIdx(bond_idx)
    return bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()


def _on_chain_bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for the multiple bonds lying on `chain`;
    the rest are cited inside substituent prefixes."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = bond_locant(chain, (a, b))
        if locant is not None:
            (ene if order == _ENE_ORDER else yne).append(locant)
    return ene, yne


def name_acyclic_unsaturated(mol) -> str:
    validate_atoms_and_bonds(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see smiles_to_iupac._cyclic)"
        )

    bonds = non_single_bonds(mol)
    if not bonds or any(order not in _VALID_ORDERS for _, _, order in bonds):
        raise UnsupportedStructure(
            "this module handles only carbon-carbon double and triple bonds; "
            "a bond order other than double or triple is not supported (see "
            "P-31.1.1.1)"
        )

    stereo = specified_double_bond_stereo(mol)
    if stereo is not None:
        if any(order == _YNE_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a specified double-bond E/Z stereo element combined with a "
                "triple bond is not supported yet (see P-93)"
            )
        ene_bond_count = sum(1 for _, _, order in bonds if order == _ENE_ORDER)
        if len(stereo) != ene_bond_count:
            raise UnsupportedStructure(
                "a non-stereogenic or unspecified double bond alongside one "
                "or more specified double-bond E/Z elements is not "
                "supported yet -- every double bond in the molecule must be "
                "specified (see P-93)"
            )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    best_key = None
    best_name = None
    best_candidate = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            ene_locants, yne_locants = _on_chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, mol=mol, unsaturated=True)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_candidate = key, name, candidate

    if stereo is not None:
        if any(bond_locant(best_candidate, _bond_atoms(mol, bond_idx)) is None for bond_idx, _ in stereo):
            raise UnsupportedStructure(
                "a specified E/Z double bond inside a substituent group is not supported yet (see P-93)"
            )
        # P-91.3: a locant always precedes each stereodescriptor, cited in
        # ascending locant order (a bare "(E)-"/"(Z)-" is only for the ring
        # systems P-91.2.2 lists, not acyclic chains).
        labels = sorted(
            (bond_locant(best_candidate, (mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx(),
                                            mol.GetBondWithIdx(bond_idx).GetEndAtomIdx())), code)
            for bond_idx, code in stereo
        )
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
