"""Naming of ketones (the '-one' suffix, C=O with two carbon substituents) on
acyclic saturated or unsaturated carbon chains and on simple monocyclic
saturated rings, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-33.4, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'one' is the preselected suffix for a ketone carbonyl, cited in a combined
  chain/suffix name the same way 'ol'/'amine' are for -OH/-NH2
  (`_alcohol.py`/`_amine.py`), e.g. 'propan-2-one', 'cyclohexanone'. Table
  3.3 ranks 'one' senior to 'ol' and 'amine' but junior to the carboxylic
  acid/ester/amide/nitrile/aldehyde suffixes; a carbonyl carbon shaped like
  an aldehyde (only one carbon neighbor) or a carboxylic acid/ester/amide
  (a second oxygen on the same carbon) is rejected outright rather than
  silently named as if it were a ketone. A coexisting hydroxyl (-OH), being
  junior to 'one', is *not* rejected: it is cited as the 'hydroxy'
  substituent prefix instead (P-41), e.g. 'CC(=O)CCO' ->
  '4-hydroxybutan-2-one'.
- P-44.4.1.8 / P-45.2: suffix locants are minimized before 'ene'/'yne'
  locants, which are minimized before substituent-prefix locants — same
  ordering as `_alcohol.py`/`_amine.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  ketone suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- P-91.3/P-92 (`tasks/ketone-stereocenter-naming.md`): a molecule with one or
  more *specified* tetrahedral stereocenters -- every one on the principal
  chain/ring itself, no unspecified one alongside them, and no C=C/C#N
  double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix, ascending
  locant order, e.g. '(3R)-3-chloropentan-2-one',
  '(2R)-2-chlorocyclohexan-1-one', same pattern as `_carboxylic_acid.py`/
  `_aldehyde.py` (chain) and `_alcohol.py`'s `_name_cyclic_alcohol` (ring).

Unlike -OH/-NH2, a ketone carbonyl carbon always has exactly two carbon
neighbors and no hydrogens, so it can never be the sole substituent on a
mononuclear (P-14.3.4.2(a)) or homogeneous two-carbon (P-14.3.4.2(b)) chain
— those locant-omission special cases from `_alcohol.py` don't apply here
and are simply absent below; even the seemingly unambiguous 'propan-2-one'
(acetone) still cites its locant, confirmed against PubChem.

This module otherwise mirrors `_alcohol.py`'s scope restrictions:

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated carbonyl oxygen or a
  singly-bonded hydroxyl (ethers, and any other oxygen shape).
- A carbonyl carbon with fewer than two carbon neighbors (aldehyde) or an
  aromatic carbonyl carbon (aryl ketone) — a separate module's territory.
- Any other heteroatom (N, S, ...) — this module only resolves the 'one'
  vs. 'ol' seniority competition (Table 3.3) between a ketone carbonyl and a
  coexisting hydroxyl; only C, halogen, and ketone/hydroxyl oxygen atoms are
  accepted at all.
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- -one on a von Baeyer polycyclic or spiro skeleton — deferred, same as
  `_alcohol.py`.

`tasks/hetero-monocyclic-ketone-naming.md`: a narrow extra path handles one
or more ketone carbonyls on an otherwise unsubstituted, saturated,
single-heteroatom (N/O/S) monocyclic ring of size 5-7 -- reusing
`_hetero_monocyclic.py`'s retained-name table (piperidine/pyrrolidine/
oxane/thiane/... for the bare stem) with the ring heteroatom always fixed
at locant 1 (only its own H, if any, may be present -- no other exocyclic
substituent), and the ketone locant set minimized over the two possible
numbering directions. PubChem-confirmed: `O=C1CCNCC1` ->
'piperidin-4-one' (CID 33721), `O=C1CCCNC1` -> 'piperidin-3-one' (CID
33722), `O=C1CCOCC1` -> 'oxan-4-one' (CID 121599), `O=C1CCSCC1` ->
'thian-4-one' (CID 66173), `O=C1CNC(=O)N1` (hydantoin) ->
'imidazolidine-2,4-dione' (CID 10006, a 1,3-diazole ring -- out of scope
here, single-heteroatom only). This path is deliberately narrow: 2+
heteroatom rings (e.g. piperazine), Se/Te heteroatoms, mancude (aromatic)
rings, 3-/4-membered rings, and any substituent other than the ring
heteroatom's own indicated hydrogen and the ketone carbonyl(s) themselves
(alkyl, halogen, hydroxyl, ...) are all out of scope -- these fall through
to this module's existing carbocyclic-only validation, which raises its
own (more general) error.

`tasks/hetero-ring-ketone-lactam-routing.md`: a ketone carbonyl directly
bonded to the ring heteroatom (a lactam, e.g. piperidin-2-one) fits this
same shape and this module already names it correctly -- the only real
blocker was `core.py`'s dispatch order, since that shape also looks
amide-shaped (a carbonyl plus a singly-bonded N with 0-2 carbon
substituents) to the earlier `has_amide_shape` check, which would
otherwise send it to `_amide.py`'s unconditional "any ring is a lactam,
out of scope" guard. `core.py` now checks this module's
`has_hetero_ring_ketone_shape` first. PubChem-confirmed: `O=C1CCCCN1` ->
'piperidin-2-one' (CID 12665), `O=C1CCCN1` -> 'pyrrolidin-2-one' (CID
12025), `O=C1CCCCCN1` -> 'azepan-2-one' (CID 7768). A symmetric,
unsubstituted cyclic imide (e.g. succinimide, `O=C1CCC(=O)N1`) fits the
same shape too (both ring carbonyls bonded to the same N) and is named
here as a plain ring dione ('pyrrolidine-2,5-dione', CID 11439) rather
than via `_imide.py`'s acyclic-only "N-acyl amide" construction
(P-66.6.3), which was never meant to cover the cyclic case.

`tasks/two-heteroatom-1-4-saturated-ring-naming.md`: the same hetero-ring
path extends to the 1,4-related two-heteroatom 6-membered saturated rings
that have their own retained/systematic name
(`_hetero_monocyclic.py`'s `saturated_two_heteroatom_1_4_ring_name`) --
morpholine (N+O), piperazine (N+N), thiomorpholine (N+S). The
higher-priority heteroatom (O or S over N, per P-22.2.1's element
seniority -- confirmed via PubChem's own '4-methylmorpholine'/
'4-methylthiomorpholine', both citing the ring N as locant 4, i.e. O/S
always wins locant 1) is fixed at locant 1 when the two elements differ;
for piperazine's two identical nitrogens, either one may be locant 1, so
both are tried alongside both directions. The other heteroatom always
lands at locant 4 regardless of direction (a 6-ring's antipodal position
is 3 steps either way), so only the ketone locant set varies between
candidates -- same minimization as the single-heteroatom path.
PubChem-confirmed: `O=C1COCCN1` -> 'morpholin-3-one' (CID 66953),
`O=C1CNCCN1` -> 'piperazin-2-one' (CID 231360), `O=C1CNC(=O)CN1` ->
'piperazine-2,5-dione' (CID 7817), `O=C1CSCCN1` -> 'thiomorpholin-3-one'
(CID 88402).

`tasks/dioxane-naming.md`: the same identical-element tie already handled
for piperazine's two nitrogens also covers O+O (1,4-dioxane) -- no new
branch needed, just adding `frozenset(("O", "O"))` to
`_TWO_HETERO_RING_ELEMENT_PAIRS` and to
`_hetero_monocyclic.py`'s own name table. PubChem-confirmed:
`O=C1COCCO1` -> '1,4-dioxan-2-one' (CID 18233), `O=C1COC(=O)CO1` ->
'1,4-dioxane-2,5-dione' (CID 65432).

`tasks/oxathiane-dithiane-naming.md`: extending to O+S (1,4-oxathiane) and
S+S (1,4-dithiane) surfaced a latent bug in `_TWO_HETERO_PRIORITY` -- it
had treated O and S as equal-priority (both outranking N only), which
never mattered while O and S never shared a pair, but silently picked
whichever heteroatom happened to come first in ring-atom order once O+S
became possible. Fixed to the strict P-22.2.1 order O > S > N (Table 2.8),
confirmed via PubChem's own '1,4-oxathian-3-one' (O at locant 1, S at
locant 4). S+S ties through the existing identical-element branch, same as
O+O and piperazine. PubChem-confirmed: `O=C1COCCS1` -> '1,4-oxathian-3-one'
(CID 15238324), `O=C1CSCCS1` -> '1,4-dithian-2-one' (CID 542724),
`O=C1CSC(=O)CS1` -> '1,4-dithiane-2,5-dione' (CID 319007).

Unlike -OH/-NH2, a ketone carbon can never itself also be a C=C/C#C alkene
carbon (its two remaining bonds, after the C=O double bond, are already
committed to its two required carbon substituents — a ketone carbon with a
third, double-bonded C=C neighbor would have five bonds), so there is no
enol-analogous "enone" case to scope out here: an alpha,beta-unsaturated
ketone like 'CC(=O)C=CC' (the double bond adjacent to, but not on, the
carbonyl carbon) is a perfectly nameable 'pent-3-en-2-one' and is supported
via the same suffix-locant-priority mechanism as `_alcohol.py`'s
'pent-4-en-1-ol'.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    specified_stereocenters,
)
from ._hetero_monocyclic import saturated_ring_name, saturated_two_heteroatom_1_4_ring_name
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}
_HETERO_RING_ELEMENTS = {7: "N", 8: "O", 16: "S"}
_HETERO_RING_SIZES = (5, 6, 7)
_TWO_HETERO_RING_ELEMENT_PAIRS = {
    frozenset(("N", "O")),
    frozenset(("N", "N")),
    frozenset(("N", "S")),
    frozenset(("O", "O")),
    frozenset(("O", "S")),
    frozenset(("S", "S")),
}
_TWO_HETERO_RING_SIZE = 6
# P-22.2.1 element seniority for locant 1 (Table 2.8's replacement-nomenclature
# order O > S > N, confirmed via PubChem's own 'oxathian-3-one' citing O at
# locant 1 and S at locant 4): the strictly lower value always wins, so O
# outranks S which outranks N.
_TWO_HETERO_PRIORITY = {"O": 0, "S": 1, "N": 2}


def _validate_and_collect_ketones(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (ketones, hydroxyls): the set of carbonyl-oxygen atom indices,
    and the set of any coexisting hydroxyl-oxygen atom indices. A hydroxyl is
    junior to 'one' in Table 3.3's suffix seniority order, so it is cited as
    the 'hydroxy' substituent prefix instead of competing for the suffix
    (P-41)."""
    ketones = set()
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a ketone carbonyl oxygen (P-33.4) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; only an isolated ketone "
                    "carbonyl or hydroxyl is supported (Table 3.3, P-33.4)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a ketone/hydroxyl oxygen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() == 1.0:
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        "an oxygen that isn't a carbonyl (=O) or hydroxyl "
                        "(-OH) is out of scope for this module"
                    )
                hydroxyls.add(atom.GetIdx())
                continue
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring (an aryl ketone) is out "
                    "of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a carbonyl carbon with fewer than two carbon neighbors "
                    "(an aldehyde or terminal carbonyl) is a more senior "
                    "characteristic group than a plain ketone (Table 3.3), "
                    "which this module does not attempt to disambiguate"
                )
            ketones.add(atom.GetIdx())
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if not ketones:
        raise UnsupportedStructure("no ketone (C=O) group found; this module only handles ketones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return ketones, hydroxyls


def _suffix_body(ene_locants, yne_locants, one_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'one' endings
    (e.g. '4-en-1-one')."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(one_locants), multiplied_word(len(one_locants), "one")))

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


def _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, one_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    one_locant_set = lowest_locant_set(one_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped)
    return (
        (
            one_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _one_locants(position_of, ketones, graph):
    locants = []
    for o in ketones:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, ketones):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in ketones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_ketone(mol, ketones, hydroxyls, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`/`_aldehyde.py`'s identical treatment),
    and the winning candidate's own locants are used to format a
    "(<locant><R/S>,...)-" prefix onto the name, ascending locant order
    (P-91.3)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _one_locants(position_of, ketones, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _one_locants({a: i + 1 for i, a in enumerate(c)}, ketones, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every ketone-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more C=O groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, ketones)
            key, name = _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _ring_cycle(graph, ring_atoms):
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _substituents_for_ring(graph, ring_order, halogens, ketones):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in ketones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, one_locants, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    one_word = multiplied_word(len(one_locants), "one")
    elide = one_word[0] in "aeiouy"
    stem = parent[:-1] if elide else parent

    if total_subs == 0 and len(one_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanone'.
        return stem + one_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(one_locants))
    return f"{prefix}{stem}-{loc_str}-{one_word}"


def _ring_candidate_key(ring_size, one_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    one_locant_set = lowest_locant_set(one_locants)
    name = _ring_name_from_substituents(ring_size, one_locants, grouped)
    return one_locant_set, locant_set, citation_locants, name


def _name_cyclic_ketone(mol, ketones, hydroxyls, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must lie on
    the ring itself (P-92: a stereocenter on a substituent branch is out
    of scope, mirroring `_alcohol.py`'s `_name_cyclic_alcohol`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the ring "
            "itself is not supported yet (see P-92)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            if one_locants is None:
                raise UnsupportedStructure(
                    "a ketone not on the ring itself (e.g. on a substituent "
                    "branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, ketones)
            key = _ring_candidate_key(ring_size, one_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _hetero_ring_heteroatom(mol):
    """The ring-atom index of the sole non-carbon heteroatom in this
    molecule's one ring, if it looks like a saturated single-heteroatom
    monocyclic ketone shape (N/O/S in a 5-, 6-, or 7-membered ring) --
    None if it doesn't match that shape at all, in which case the caller
    falls through to the existing carbocyclic-only path (which raises its
    own, more general error for whatever doesn't fit)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) not in _HETERO_RING_SIZES:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    (heteroatom,) = heteroatoms
    if mol.GetAtomWithIdx(heteroatom).GetAtomicNum() not in _HETERO_RING_ELEMENTS:
        return None
    return heteroatom


def _validate_and_collect_hetero_ring_ketone(mol, heteroatoms):
    """(ring_order, ketones, elements_by_atom) for a hetero-ring ketone --
    validates that every ring heteroatom (`heteroatoms`, 1 or 2 atom
    indices) carries no substituent beyond its own indicated hydrogen,
    every other ring atom is a plain CH2 or an unsubstituted ketone
    carbonyl carbon, and the ring itself is fully saturated (P-22.2.1's
    plain retained-name ring shape)."""
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    graph = adjacency(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_set = set(ring_atoms)

    elements_by_atom = {}
    for heteroatom in heteroatoms:
        hetero_atom_obj = mol.GetAtomWithIdx(heteroatom)
        element = _HETERO_RING_ELEMENTS[hetero_atom_obj.GetAtomicNum()]
        elements_by_atom[heteroatom] = element
        if hetero_atom_obj.GetFormalCharge() != 0 or hetero_atom_obj.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        expected_h = 1 if element == "N" else 0
        if hetero_atom_obj.GetDegree() != 2 or hetero_atom_obj.GetTotalNumHs() != expected_h:
            raise UnsupportedStructure(
                "a ring heteroatom bearing a substituent is out of scope for "
                "this module's hetero-ring ketone path (see "
                "tasks/hetero-monocyclic-ketone-naming.md)"
            )

    ketones = set()
    for atom_idx in ring_atoms:
        if atom_idx in elements_by_atom:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        exo = [n for n in graph[atom_idx] if n not in ring_set]
        if not exo:
            if atom.GetTotalNumHs() != 2:
                raise UnsupportedStructure(
                    "a substituent other than a ketone carbonyl is out of "
                    "scope for this module's hetero-ring ketone path"
                )
            continue
        if len(exo) != 1:
            raise UnsupportedStructure(
                "a substituent other than a ketone carbonyl is out of "
                "scope for this module's hetero-ring ketone path"
            )
        (other,) = exo
        other_atom = mol.GetAtomWithIdx(other)
        bond = mol.GetBondBetweenAtoms(atom_idx, other)
        if other_atom.GetAtomicNum() != 8 or other_atom.GetDegree() != 1 or bond.GetBondTypeAsDouble() != 2.0:
            raise UnsupportedStructure(
                "a substituent other than a ketone carbonyl is out of "
                "scope for this module's hetero-ring ketone path"
            )
        ketones.add(other)
    if not ketones:
        raise UnsupportedStructure("no ketone (C=O) group found; this module only handles ketones")

    for i in range(len(ring_order)):
        a, b = ring_order[i], ring_order[(i + 1) % len(ring_order)]
        if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "an unsaturated hetero ring is out of scope for this "
                "module's hetero-ring ketone path"
            )

    return ring_order, ketones, elements_by_atom


def _best_one_locants(graph, ring_order, ketones, starts):
    """Lowest ketone locant set over every (start, direction) candidate in
    `starts` (each start rotated to position 1, tried both directions) --
    shared by the single- and two-heteroatom hetero-ring ketone paths."""
    best_locants = None
    for start in starts:
        rotated_start = ring_order.index(start)
        rotated = ring_order[rotated_start:] + ring_order[:rotated_start]
        for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            if one_locants is None:
                raise UnsupportedStructure(
                    "a ketone not on the ring itself (e.g. on a substituent "
                    "branch) is not supported yet"
                )
            locant_set = lowest_locant_set(one_locants)
            if best_locants is None or locant_set < best_locants:
                best_locants = locant_set
    return best_locants


def _hetero_ring_ketone_name(stem, best_locants):
    one_word = multiplied_word(len(best_locants), "one")
    elide = one_word[0] in "aeiouy"
    base = stem[:-1] if elide else stem
    loc_str = ",".join(str(loc) for loc in best_locants)
    return f"{base}-{loc_str}-{one_word}"


def _name_hetero_cyclic_ketone(mol, heteroatom):
    graph = adjacency(mol)
    ring_order, ketones, elements_by_atom = _validate_and_collect_hetero_ring_ketone(mol, {heteroatom})
    ring_size = len(ring_order)
    element = elements_by_atom[heteroatom]
    stem = saturated_ring_name(element, ring_size)
    if stem is None:
        raise UnsupportedStructure(
            f"no retained/Hantzsch-Widman name for a {ring_size}-membered "
            f"{element}-heteroatom saturated ring (P-22.2.1)"
        )
    best_locants = _best_one_locants(graph, ring_order, ketones, [heteroatom])
    return _hetero_ring_ketone_name(stem, best_locants)


def _hetero_ring_two_heteroatoms(mol):
    """(het1, het2) ring-atom indices for a saturated, 6-membered,
    1,4-related two-heteroatom ketone shape (morpholine/piperazine/
    thiomorpholine's element pairs only -- see
    `_TWO_HETERO_RING_ELEMENT_PAIRS`), or None if it doesn't match that
    shape at all."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _TWO_HETERO_RING_SIZE:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _TWO_HETERO_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    ring_order = _ring_cycle(graph, list(ring_atoms))
    if abs(ring_order.index(het1) - ring_order.index(het2)) != _TWO_HETERO_RING_SIZE // 2:
        # Not the 1,4 (antipodal) relationship morpholine/piperazine/
        # thiomorpholine need -- e.g. a 1,2- or 1,3-diheteroatom ring,
        # which has no retained name and is out of scope.
        return None
    return het1, het2


def _name_two_hetero_cyclic_ketone(mol, het1, het2):
    graph = adjacency(mol)
    ring_order, ketones, elements_by_atom = _validate_and_collect_hetero_ring_ketone(mol, {het1, het2})
    stem = saturated_two_heteroatom_1_4_ring_name((elements_by_atom[het1], elements_by_atom[het2]))
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this two-heteroatom saturated ring "
            "(P-22.2.1)"
        )
    if elements_by_atom[het1] == elements_by_atom[het2]:
        starts = [het1, het2]
    else:
        starts = [het1 if _TWO_HETERO_PRIORITY[elements_by_atom[het1]] < _TWO_HETERO_PRIORITY[elements_by_atom[het2]] else het2]
    best_locants = _best_one_locants(graph, ring_order, ketones, starts)
    return _hetero_ring_ketone_name(stem, best_locants)


def has_hetero_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the narrow hetero-ring-ketone shape
    (single heteroatom, see `_hetero_ring_heteroatom`, or the 1,4
    two-heteroatom morpholine/piperazine/thiomorpholine/1,4-dioxane
    shape, see `_hetero_ring_two_heteroatoms`) -- used by `core.py` to
    route ahead of both `has_ester_shape` and `has_amide_shape`, since a
    ketone directly bonded to a ring oxygen/sulfur (a lactone, e.g.
    oxan-2-one/1,4-dioxan-2-one) or ring nitrogen (a lactam, e.g.
    piperidin-2-one/morpholin-3-one) would otherwise look ester- or
    amide-shaped to those checks and get rejected by `_ester.py`'s
    acyclic-only construction or `_amide.py`'s ring-always-out-of-scope
    guard before ever reaching this module."""
    return _hetero_ring_heteroatom(mol) is not None or _hetero_ring_two_heteroatoms(mol) is not None


def name_ketone(mol) -> str:
    hetero_atom = _hetero_ring_heteroatom(mol)
    if hetero_atom is not None:
        return _name_hetero_cyclic_ketone(mol, hetero_atom)
    two_heteroatoms = _hetero_ring_two_heteroatoms(mol)
    if two_heteroatoms is not None:
        return _name_two_hetero_cyclic_ketone(mol, *two_heteroatoms)
    ketones, hydroxyls = _validate_and_collect_ketones(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the ketone oxygen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in ketones and b[1] not in ketones]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior ketone/aldehyde "
                "form and is out of scope for this module (P-31.1.4.2.4)"
            )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_ketone(mol, ketones, hydroxyls, bonds, stereo)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_ketone(mol, ketones, hydroxyls, stereo)
    raise UnsupportedStructure(
        "polycyclic and spiro ketones are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
