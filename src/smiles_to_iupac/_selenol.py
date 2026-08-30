"""Naming of selenols (the '-selenol' suffix, -SeH) on acyclic saturated or
unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.1.1, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'selenol' is the next chalcogen analogue in the same family as 'thiol'
  (`_thiol.py`) — 'ol'/'thiol'/'selenol'/'tellurol' are all the same
  substitutive-suffix construction, just with O/S/Se/Te respectively.
  This module mirrors `_thiol.py`'s own first-pass scope (single -SH,
  before that module's later dithiol/monocyclic extensions) exactly,
  swapping the chalcogen atomic number (34 instead of 16) and suffix word
  ('selenol' instead of 'thiol') and otherwise reusing the same
  chain-search/locant machinery from `_common.py` directly (this module is
  new enough to use the shared `longest_chains`/`bond_locant`/
  `bond_locants`/`group_substituents`/`multiplied_word` helpers rather
  than adding yet another private copy of them -- `_thiol.py` itself still
  has its own pre-`_common.py` copies, a known, separate cleanup debt not
  addressed here). Confirmed via PubChem PUG REST: CID 440764 (`C[SeH]`)
  -> "methaneselenol", CID 5252527 (`CC[SeH]`) -> "ethaneselenol" (the
  same P-14.3.4.2(b) two-carbon locant omission as 'ethanethiol'), CID
  71373846 (`CCC[SeH]`) -> "propane-1-selenol", CID 157368643
  (`CC(C)C[SeH]`) -> "2-methylpropane-1-selenol" (branched), CID 15821407
  (`C=CC[SeH]`) -> "prop-2-ene-1-selenol" (unsaturated chain).
- P-35.2.1 (Chapter P-3): a halogen substituent on the carbon chain is
  prefix-only and coexists freely with the -SeH suffix, the same
  mechanism already independently verified in `_thiol.py`/`_alcohol.py`/
  several other modules -- no PubChem-listed halogenated selenol was
  found to independently confirm this specific combination (every
  candidate SMILES tried came back as CID 0), so this is a reviewed
  (eyeballed), not independently verified, extension along an otherwise
  already-confirmed single axis.
- Two -SeH groups (a diselenol) mirrors `_thiol.py`'s dithiol extension --
  the locant/suffix machinery already generalizes over a list of
  selenol locants, so only the single-group guard needed lifting.
  Confirmed via PubChem PUG REST: `[SeH]CC[SeH]` ->
  "ethane-1,2-diselenol", `[SeH]CCC[SeH]` -> "propane-1,3-diselenol".
- Three or more -SeH groups (triselenol, tetraselenol, ...) lifts the
  group-count cap entirely, the same way `_thiol.py` itself generalizes
  to "dithiol/trithiol/..." without citing a specific 3+-group PubChem
  worked example of its own -- `_suffix_body`/`_name_from_substituents`
  below are already fully generalized over an arbitrary-length
  `se_locants` list (the same shared mechanism `_alcohol.py`'s own
  confirmed polyol support, e.g. glycerol, uses), so lifting the cap adds
  no new code path to verify. No triselenol compound (three -SeH groups
  on a single chain) was found registered in PubChem to independently
  confirm this specific chalcogen (unlike `_thiol.py`'s own trithiol
  case, which is in the same boat) -- this is a reviewed, not
  independently structure-verified, generalization along an
  already-confirmed axis, matching `_thiol.py`'s own precedent exactly.
- A single, otherwise-unsubstituted saturated monocyclic ring also works,
  mirroring `_thiol.py`'s own monocyclic support: confirmed via PubChem
  structure match, `C1CCCCC1[SeH]` -> "cyclohexaneselenol",
  `C1CCCC1[SeH]` -> "cyclopentaneselenol". The ring machinery below
  (`_ring_name_from_substituents`/`_ring_candidate_key`) is ported
  directly from `_thiol.py`'s own already-generalized-over-multiple-
  groups ring code, so a multi-selenol ring is supported too, the same
  way `_thiol.py` itself doesn't cite a specific multi-group ring PubChem
  example either -- trusting the shared, already-exercised locant
  machinery rather than re-deriving it.

Scope, deliberately narrow (mirrors `_thiol.py`'s own group-count- and
ring-generalized scope): one or more -SeH groups on an acyclic chain, or
on a single saturated carbon ring, with no other heteroatom (in
particular no -OH, -SH, or amine nitrogen) anywhere in the molecule.
Explicitly out of scope (raise `UnsupportedStructure`): polycyclic/spiro/
unsaturated rings, an -SeH on a substituent branch off an otherwise-
unsubstituted ring, a selenide (-Se- ether-analogue) or any other
selenium-oxidation-state group, a tellurol or other chalcogen atom, and
any oxygen or nitrogen atom at all.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_cycle,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_YNE_BOND_ORDER = 3.0
_SELENIUM = 34
_ALLOWED_ATOMIC_NUMS = {6, _SELENIUM, *HALOGEN_PREFIXES}


def has_selenol_shape(mol) -> bool:
    return any(atom.GetAtomicNum() == _SELENIUM for atom in mol.GetAtoms())


def _validate_and_collect_selenols(mol):
    selenols = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a selenol selenium (P-63.1.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting -OH, -SH, or amine nitrogen "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for selenols"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure("aromatic rings are out of scope for this module")
        elif atomic_num == _SELENIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a selenium bonded to more than one heavy atom (e.g. a "
                    "selenide) is out of scope; only an isolated selenol "
                    "(-SeH) is supported (P-63.1.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("a selenium double-bonded to carbon is not a selenol")
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "a -Se- atom that isn't a simple selenol (-SeH) is out "
                    "of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a selenol must be attached to a carbon atom")
            selenols.add(atom.GetIdx())
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
        )
    if not selenols:
        raise UnsupportedStructure(
            "no selenol (-SeH) group found; this module only handles selenols"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return selenols


def _reject_eneselenol_carbons(graph, selenols, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for se_idx in selenols:
        (carbon,) = graph[se_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a selenol on a carbon that is also part of a C=C/C#C bond "
                "is out of scope for this module"
            )


def _suffix_body(ene_locants, yne_locants, se_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(se_locants), multiplied_word(len(se_locants), "selenol")))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, se_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited.
        selenol_word = multiplied_word(len(se_locants), "selenol")
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + selenol_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(se_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (here, the sole -SeH) in total omits the locant, e.g.
        # 'ethaneselenol'.
        return alkane_name(2) + "selenol"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, se_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, se_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    se_locant_set = lowest_locant_set(se_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, se_locants, ene_locants, yne_locants, grouped)
    return (
        (
            se_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _se_locants(position_of, selenols, graph):
    locants = []
    for se in selenols:
        (carbon,) = graph[se]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, selenols):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in selenols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _substituents_for_ring(graph, ring_order, halogens, selenols):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in selenols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, se_locants, grouped):
    stem = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    selenol_word = multiplied_word(len(se_locants), "selenol")

    if total_subs == 0 and len(se_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexaneselenol'.
        return stem + selenol_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(se_locants))
    return f"{prefix}{stem}-{loc_str}-{selenol_word}"


def _ring_candidate_key(ring_size, se_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    se_locant_set = lowest_locant_set(se_locants)
    name = _ring_name_from_substituents(ring_size, se_locants, grouped)
    return se_locant_set, locant_set, citation_locants, name


def _name_cyclic_selenol(mol, selenols):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, selenols, graph)
            substituents = _substituents_for_ring(graph, candidate, halogens, selenols)
            key = _ring_candidate_key(ring_size, se_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_selenol(mol) -> str:
    selenols = _validate_and_collect_selenols(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, _YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_eneselenol_carbons(graph, selenols, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        raise UnsupportedStructure(
            "polycyclic/spiro selenols are not supported yet (this module "
            "only handles acyclic chains and a single saturated ring)"
        )
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_selenols = {s for s in selenols if next(iter(graph[s])) in ring_atoms}
        if ring_selenols != selenols:
            raise UnsupportedStructure(
                "a selenol on a substituent branch chain rather than the "
                "ring itself is not supported yet"
            )
        return _name_cyclic_selenol(mol, selenols)

    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _se_locants(position_of, selenols, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every selenol-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain, "
            "or an -SeH expressed as a 'selanyl' substituent prefix, is "
            "not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            se_locants = _se_locants(position_of, selenols, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, selenols)
            key, name = _candidate_key(chain_length, se_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
