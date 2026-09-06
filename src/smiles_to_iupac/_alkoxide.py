"""Naming of alkoxide anions (R-O(-), the '-olate' suffix), restricted to a
single such group on an acyclic carbon skeleton, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-72.2.2.2.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html): an
  anion formed by removing a hydron from the chalcogen atom of a hydroxy
  compound (a plain '-ol', not a carbonyl/anion tautomer pair like the
  carboxylate/thioate/selenoate modules) is named by adding 'ate' to the
  '-ol' suffix -- e.g. 'propan-2-ol' -> 'propan-2-olate (PIN)' (confirmed
  worked example; PubChem CID 3260420 structure match, and PubChem's own
  auto-generated name agrees exactly here, unlike some of this project's
  other isotope/anion modules).
- Retained names: "methoxide, ethoxide, propoxide, butoxide, tert-butoxide,
  phenoxide (but not isopropoxide), and aminoxide" are preferred IUPAC
  names for CH3-O(-), C2H5-O(-), (n-)C3H7-O(-), (n-)C4H9-O(-),
  (CH3)3C-O(-), C6H5-O(-), and H2N-O(-) respectively (P-72.2.2.2.2's own
  text, `tmp/bluebook/P7.txt` lines 1108-1111). This module implements the
  first five (methoxide/ethoxide/propoxide/butoxide/tert-butoxide) for
  their exact plain, unhalogenated, unsaturated, terminal-oxygen shape
  only -- isopropoxide is explicitly NOT the PIN for propan-2-ol's anion
  (the worked example 'propan-2-olate (PIN)' confirms this), so a
  non-terminal oxygen always falls through to the systematic '-olate'
  path. Phenoxide (-O(-) directly on a benzene ring, `_name_phenoxide`)
  is implemented as a separate aromatic-ring path -- P-63.8.1 confirms it
  "may be substituted in the same way as the corresponding alcohols"
  (`tmp/bluebook/P6.txt` lines 3091-3111, worked example "lithium
  phenoxide (PIN)"), mirroring `_alcohol.py`'s `_name_phenol` exactly.
  Aminoxide (a carbon-free H2N-O(-) shape) is not implemented.
- Structure-verified via PubChem: `CC[O-]` (CID 119440), `CCC[O-]` (CID
  12543515), `CC(C)[O-]` (CID 3260420) -- PubChem's own generated names use
  the systematic '-olate' form even for the retained-name cases (e.g.
  "ethanolate" for CC[O-]), which is expected (this project has repeatedly
  observed PubChem's generator not always picking the Blue Book's specific
  retained-name PIN) -- the retained names themselves come directly from
  the Blue Book's own text, not from PubChem.
- Mirrors `_carboxylate.py`/`_selenoate.py`'s acyclic, single-group,
  halogen+unsaturation-coexisting scope, but unlike those, the oxygen here
  sits as a substituent ON a chain carbon rather than being the chain's own
  terminus, so its locant follows the same citation rules as `_alcohol.py`'s
  '-ol' suffix (P-14.3.4.2(a) mononuclear-parent locant omission; otherwise
  the locant is always cited on a chain of 2+ carbons once the retained
  names' narrow, unsubstituted shape doesn't apply).
- A branched carbon skeleton is named via `_alcohol.py`'s own general
  substituent-branch machinery (`name_branch`), reused here unchanged --
  e.g. `CC(C)C[O-]` -> "2-methylpropan-1-olate" (PubChem structure/name
  match), the same skeleton as `_alcohol.py`'s own "2-methylpropan-1-ol"
  with 'ate' appended. A retained name (methoxide/ethoxide/propoxide/
  butoxide) only applies to a truly unbranched chain using every carbon in
  the molecule -- a branched skeleton whose longest chain happens to match
  a retained name's length (e.g. a 5-carbon branched skeleton with a
  4-atom longest chain) must fall through to the systematic path instead.

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring anywhere in the molecule other than a single plain benzene ring
  bearing the -O(-) itself or one chain substituent (acyclic-only,
  mirrors the other anion modules; see `_name_phenoxide`/
  `_name_phenyl_chain_alkoxide` for the two benzene-ring paths).
- More than one -O(-) group, or any oxygen that isn't the single alkoxide
  anion (an ether, a second alkoxide, a carbonyl).
- A carbon bearing the anionic oxygen that's also double-bonded to another
  chalcogen (O/S/Se/Te) -- that shape is a carboxylate/thioate/selenoate-
  style carbonyl-anion tautomer pair (P-72.2.2.2.1.1), handled by its own
  module, not this one's plain hydroxy-compound anion.
- Any charged or radical atom other than the single alkoxide oxygen's
  formal charge -1.
- Any other heteroatom (N, S, ...), aromatic ring, or isotopic modification.

P-91.3/P-92: a molecule with one
or more *specified* tetrahedral stereocenters -- every one on the
principal chain itself, no unspecified one alongside them, and no
C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
ascending locant order, same pattern as `_sulfonic_acid.py`/`_thiol.py`
(chain only, this module has no ring support). The alkoxide oxygen itself
is never a potential stereocenter (a monovalent, negatively-charged
terminal atom), confirmed via RDKit `FindPotentialStereo` on
`CC[C@@H](C)[O-]`.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_chain_attachment,
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_RETAINED_ALKOXIDES = {1: "methoxide", 2: "ethoxide", 3: "propoxide", 4: "butoxide"}
_CHALCOGENS = (8, 16, 34, 52)


def _find_alkoxide_oxygens(mol):
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        if atom.GetFormalCharge() != -1 or atom.GetDegree() != 1:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6:
            continue
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx())
        if bond.GetBondTypeAsDouble() != 1.0:
            continue
        # Exclude a carboxylate/thioate/selenoate-style carbon (a carbonyl-
        # anion tautomer pair, P-72.2.2.2.1.1) -- its own module already
        # handles that shape, and it isn't a plain hydroxy-compound anion
        # (P-72.2.2.2.2) at all.
        if any(
            n.GetAtomicNum() in _CHALCOGENS
            and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in carbon.GetNeighbors()
        ):
            continue
        matches.append(atom)
    return matches


def has_alkoxide_shape(mol) -> bool:
    return bool(_find_alkoxide_oxygens(mol))


def _find_alkoxide_group(mol):
    matches = _find_alkoxide_oxygens(mol)
    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one alkoxide (-O(-)) group is required; zero or "
            "multiple such groups are not supported yet (P-72.2.2.2.2)"
        )
    return matches[0]


def _suffix_body(ene_locants, yne_locants, o_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append(([o_locant], "olate"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    elide_stem = words[0][0] in "aeiouy"
    return "-".join(parts), elide_stem


def _name_from_substituents(chain_length, o_locant, ene_locants, yne_locants, grouped):
    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited, however many substituents there are.
        stem = alkane_name(1)[:-1]  # 'olate' starts with a vowel
        return format_substituent_prefixes(grouped, omit_locants=True) + stem + "olate"

    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, o_locant)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, o_locant, ene_locants, yne_locants, substituents):
    """Sort key: lowest locant to the principal characteristic group (the
    alkoxide oxygen) first, mirroring `_alcohol.py`'s own '-ol' priority,
    then ene/yne locants, then substituent-prefix locants (P-45.2)."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, o_locant, ene_locants, yne_locants, grouped)
    return (
        (
            (o_locant,),
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded_atoms):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _is_tert_butoxide(mol, oxygen_idx, bonds, halogen_atoms):
    if bonds or halogen_atoms:
        return False
    (oxygen_carbon,) = mol.GetAtomWithIdx(oxygen_idx).GetNeighbors()
    if oxygen_carbon.GetDegree() != 4:
        return False
    methyls = [n for n in oxygen_carbon.GetNeighbors() if n.GetIdx() != oxygen_idx]
    return len(methyls) == 3 and all(m.GetAtomicNum() == 6 and m.GetDegree() == 1 for m in methyls)


def _name_acyclic_alkoxide(mol, oxygen_idx, excluded_atoms, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch is out of scope), and the winning candidate's own
    locants are used to format a "(<locant><R/S>,...)-" prefix onto the
    final name. A genuine stereocenter always requires a branched skeleton
    (module docstring's retained names/tert-butoxide are all unbranched or
    symmetric), so `stereo` never fires alongside those fast paths."""
    halogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in HALOGEN_PREFIXES]
    if stereo is None and _is_tert_butoxide(mol, oxygen_idx, bonds, halogen_atoms):
        return "tert-butoxide"

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    carbon_graph = carbon_adjacency(mol)
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    (oxygen_carbon,) = [n.GetIdx() for n in mol.GetAtomWithIdx(oxygen_idx).GetNeighbors()]

    eligible = []
    for chain in chains:
        if oxygen_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            oxygen_carbon in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the alkoxide carbon (and/or multiple bonds) does not lie on a "
            "single longest carbon chain"
        )

    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if (
        stereo is None
        and chain_length == total_carbons
        and chain_length in _RETAINED_ALKOXIDES
        and not bonds
        and not halogen_atoms
    ):
        terminal_positions = {chains[0][0], chains[0][-1]}
        if oxygen_carbon in terminal_positions:
            return _RETAINED_ALKOXIDES[chain_length]

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            o_locant = position_of[oxygen_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded_atoms)
            key, name = _candidate_key(chain_length, o_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _validate_and_prepare_alkoxide(mol, aromatic_ring_atoms=frozenset()):
    """(oxygen, excluded_atoms, bonds, stereo) after checking the molecule
    fits this module's scope (see module docstring).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_alkoxide`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_alkoxide`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors `_alcohol.py`'s
    equivalent aromatic-exemption pattern."""
    oxygen = _find_alkoxide_group(mol)
    excluded_atoms = {oxygen.GetIdx()}
    stereo = specified_stereocenters(mol)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in excluded_atoms:
            continue
        if atomic_num not in {6, *HALOGEN_PREFIXES}:
            raise UnsupportedStructure(
                "heteroatoms other than the alkoxide's own oxygen "
                "(P-72.2.2.2.2) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "single alkoxide anion oxygen is not supported yet"
            )
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1)"
            )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_atoms
        and b[1] not in excluded_atoms
        and (b[0] not in aromatic_ring_atoms or b[1] not in aromatic_ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    return oxygen, excluded_atoms, bonds, stereo


def _phenoxide_name_from_substituents(grouped):
    # P-63.8.1/P-72.2.2.2.2: the retained name 'phenoxide' stands for the
    # whole ring+O(-) system and "may be substituted in the same way as
    # the corresponding alcohols" (`tmp/bluebook/P6.txt` lines 3091-3111)
    # -- mirrors `_alcohol.py`'s `_name_phenol`'s identical treatment of
    # 'phenol', so the O(-)'s own ring locant is never cited, only other
    # substituents'.
    if not grouped:
        return "phenoxide"
    return f"{format_substituent_prefixes(grouped)}phenoxide"


def _phenoxide_candidate_key(o_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _phenoxide_name_from_substituents(grouped)
    return o_locant, locant_set, citation_locants, name


def _substituents_for_ring(graph, ring_order, halogens, excluded):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_phenoxide(mol, ring_atoms, oxygen):
    """P-63.8.1/P-72.2.2.2.2: -O(-) attached directly to a benzene ring
    carbon -- e.g. 'phenoxide' (PubChem CID 998, structure match; PubChem's
    own auto-generated name), '4-methylphenoxide' (structure-verified
    against the corresponding phenol, CID 2879, per the Blue Book's own
    "substituted the same way as the corresponding alcohols" text -- see
    module docstring). Mirrors `_alcohol.py`'s `_name_phenol` exactly, with
    'phenoxide' as the retained parent name in place of 'phenol'."""
    if specified_stereocenters(mol):
        raise UnsupportedStructure("a specified stereocenter alongside phenoxide is not supported yet")

    graph = adjacency(mol)
    oxygen_idx = oxygen.GetIdx()
    (oxygen_carbon,) = [n.GetIdx() for n in oxygen.GetNeighbors()]
    halogens = halogen_substituents(mol)
    excluded = {oxygen_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            o_locant = position_of[oxygen_carbon]
            substituents = _substituents_for_ring(graph, candidate, halogens, excluded)
            key = _phenoxide_candidate_key(o_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_alkoxide(mol, ring_atoms):
    """Name an alkoxide whose -O(-) lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropan-1-olate. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_alcohol.py`'s `_name_phenyl_chain_alcohol`. Narrower than the
    acyclic path above: no chain unsaturation and no specified
    stereocenter -- the retained names (methoxide/.../tert-butoxide) never
    apply here since they require every carbon in the molecule to be part
    of the alkoxide's own chain, which a benzene-ring substituent always
    violates."""
    oxygen, excluded_atoms, bonds, stereo = _validate_and_prepare_alkoxide(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if stereo:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "alkoxide chain is not supported yet"
        )
    non_ring_unsaturation = [b for b in bonds if b[0] not in ring_atoms or b[1] not in ring_atoms]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "alkoxide chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain alkoxide is not supported yet"
        )
    (oxygen_carbon,) = [n.GetIdx() for n in oxygen.GetNeighbors()]
    chain, branches = longest_branched_chain_through(graph, oxygen_carbon, ring_atoms, excluded_atoms)
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    halogens = halogen_substituents(mol)
    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        o_locant = position_of[oxygen_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, o_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_alkoxide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            oxygen = _find_alkoxide_group(mol)
            (oxygen_carbon,) = [n.GetIdx() for n in oxygen.GetNeighbors()]
            if oxygen_carbon in ring_atoms:
                return _name_phenoxide(mol, ring_atoms, oxygen)
            return _name_phenyl_chain_alkoxide(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring is out of scope for this acyclic-only module"
        )

    oxygen, excluded_atoms, bonds, stereo = _validate_and_prepare_alkoxide(mol)
    return _name_acyclic_alkoxide(mol, oxygen.GetIdx(), excluded_atoms, bonds, stereo)
