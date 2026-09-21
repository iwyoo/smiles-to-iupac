"""Naming of alcohols (the '-ol' suffix, -OH) on acyclic saturated or
unsaturated carbon chains, on simple monocyclic saturated rings, and on a
chain hanging off an otherwise-plain monocyclic ring (saturated, or
carrying ring C=C unsaturation of its own), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-44.1.1 ring-vs-chain competition (P-44.1.2.2): when the ring itself
  bears no -OH at all, there is no actual competition for which structure
  is the senior parent -- the chain trivially captures the maximum number
  of principal characteristic groups (all of them), so it's the parent and
  the ring is cited as a plain "cyclo..." substituent prefix (P-29.3.3),
  e.g. cyclohexylmethanol. When the ring itself carries one or more C=C
  double bonds (and still no -OH of its own), it's instead cited via
  `_cyclic_unsaturated.name_cyclic_unsaturated_yl`, with the attachment
  fixed at locant 1 (P-29.2) -- e.g. '(cyclohex-3-en-1-yl)methanol'; this
  path still requires the ring to have no other exocyclic branch (see
  `_common.ring_chain_attachment`), so a ring bearing both the -OH chain and
  another substituent of its own stays unsupported regardless of
  saturation. When the ring's own -OH count is at least the
  chain substituent's -OH count, the ring is always the senior parent:
  either P-44.1.1 settles it outright (the ring captures strictly more of
  the principal characteristic group), or the two counts tie and
  P-44.1.2.2 resolves the tie in the ring's favor (no chain-length
  comparison, unlike the 1993 recommendations). Since the winner is fixed
  for this whole "ring count >= chain count" shape, this module doesn't
  build or compare a real chain-parent candidate name; it simply names the
  ring as parent and cites the chain (with *all* of its own -OH's, if more
  than one) as a "(hydroxy...alkyl)"/"(dihydroxy...alkyl)" substituent
  prefix, reusing the `{oxygen_idx: "hydroxy"}` trick already used by
  `_carboxylic_acid.py`/`_amide.py`/`_aldehyde.py`/`_ketone.py` -- mapping
  multiple chain hydroxyls this way lets the shared substituent-grouping
  machinery multiply the "hydroxy" prefix exactly as it already does for
  repeated halogens, no new logic needed. The reverse case (chain -OH
  count greater than the ring's, regardless of how many -OH's the
  outcompeted ring has of its own) is supported too: the chain becomes the
  senior parent and the ring is cited as a substituent via
  `_substituents.name_branch`'s `_ring_substituent_with_hydroxyls`, which
  cites the ring's own hydroxyls together with a "di"/"tri" multiplying
  prefix (e.g. '1-(4-hydroxycyclohexyl)ethane-1,2-diol', PubChem CID
  21395558) -- completing this module's whole ring-vs-chain competition.
- P-92 (Chapter P-9, https://iupac.qmul.ac.uk/BlueBook/P9.html): an
  acyclic (chain) alcohol whose molecule has one or more stereo elements
  overall -- every one a specified tetrahedral stereocenter located on the
  principal chain itself, no unspecified one, and no C=C/C#N double-bond
  E/Z stereo anywhere -- gets a "(<locant><R/S>)-" prefix, e.g.
  '(2R)-butan-2-ol', '(3R)-pent-1-en-3-ol' (both Blue Book worked
  examples). Two or more stereocenters are cited together in one
  parenthesized group, ascending locant order, comma-separated (P-91.3),
  e.g. '(2R,3R)-3-chlorobutan-2-ol' (PubChem CID 12575191) -- CIP priority
  computation itself is delegated entirely to RDKit
  (`_common.specified_stereocenters`); this module only formats the
  resulting label(s) using the chain locants already computed for the
  winning numbering (P-92 doesn't get its own say in *which* numbering
  wins -- it's purely descriptive once the chain/locants are otherwise
  fixed). The same mechanism extends to a plain monocyclic ring whose
  -OH's are all on the ring itself (no exocyclic hydroxyl chain): every
  specified stereocenter must lie on the ring itself, e.g.
  '(1S,2S)-2-methylcyclohexan-1-ol' (PubChem CID 642632; PubChem's own
  redundant relative "trans-"/"cis-" prefix is dropped, matching this
  project's existing acyclic convention of citing R/S alone). An acyclic
  chain's specified tetrahedral stereocenter(s) may also coexist with
  specified C=C double-bond E/Z element(s), cited together in the same
  ascending-locant group (`_common.specified_stereo_elements`), e.g.
  '(2Z,5R,7E)-nona-2,7-dien-5-ol' (a Blue Book worked example, P-91.3).
  Any stereo element beyond this (a specified element mixed with an
  unspecified one, a stereocenter or E/Z double bond on a substituent
  branch rather than the chain itself, an E/Z double bond coexisting with
  a ring, or any stereocenter at all on a polycyclic/spiro skeleton or
  alongside a ring-vs-chain hydroxyl comparison) raises
  `UnsupportedStructure` explicitly.

- P-33.2.1, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'ol' is the preselected suffix for -OH, ranked 14th (out of 17) in Table
  3.3's seniority-for-citation-as-suffix order. Table 3.3 also lists several
  suffixes senior to 'ol' that this module must not be fooled into
  misreading as a plain alcohol: carboxylic acid/oic acid ('-CO-OH'), amide,
  nitrile, and in particular 'al' ('-CHO') and 'one' ('=O') — any oxygen
  double-bonded to a carbon (aldehyde, ketone, or the carbonyl half of a
  carboxylic acid) is rejected outright rather than silently named as if it
  were an alcohol.
- P-30 Introduction (Chapter P-3): "a characteristic group systematically
  introduced as a suffix attached to a parent hydride, for example ...
  ethanol" is the Blue Book's own example of this construction: the parent
  hydride name ('ethane') with its final 'e' elided before the vowel-initial
  suffix 'ol' -> 'ethanol'.
- P-44.1.1 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  senior parent structure has the *maximum number* of principal
  characteristic groups (here, -OH) — this is applied before chain length
  (P-44.3.2) is even considered, e.g. 'pentane-1,4-diol' is preferred over a
  longer chain that would only capture one of two -OH groups on its own
  chain, the second -OH being pushed into a branch as a 'hydroxy' prefix.
  This module does NOT implement that full generality: it only accepts a
  candidate chain that is simultaneously (a) one of the longest carbon
  chains (P-44.3.2, as in `_acyclic.py`/`_unsaturated.py`) and (b) carries
  every -OH-bearing carbon in the molecule (mirroring `_unsaturated.py`'s
  existing, analogous restriction that every multiple bond must lie on a
  candidate longest chain). A molecule where the true PIN would require a
  shorter chain to capture more -OH groups, or would push a leftover -OH
  into a 'hydroxy' substituent prefix, is out of scope and raises
  `UnsupportedStructure`.
- P-44.4.1 (Chapter P-4): once the parent chain/ring is otherwise fixed, its
  *numbering* is chosen by criteria (a)-(l) in order; in particular
  criterion (h) (P-44.4.1.8, "the lower locant for an attached group
  expressed as a suffix") outranks criterion (j) (P-44.4.1.10, the 'ene'/
  'yne' locants) which in turn outranks substituent-prefix locants (P-45.2,
  applied after all of P-44.4.1 — see `_acyclic.py`/`_unsaturated.py`).
  Concretely: the -OH locant set is minimized first; only if that leaves a
  choice are 'ene'/'yne' locants minimized; only if that still leaves a
  choice are substituent-prefix locants minimized. Example confirmed against
  this exact ordering: 'pent-4-en-1-ol' (not 'pent-1-en-5-ol' /
  'pent-1-en-4-ol'-numbered-from-the-alkene-end) — the -OH gets locant 1
  even though this gives the double bond the *higher* available locant (4
  rather than 1).
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-alcohol name (e.g. 'pent-4-en-1-ol')
  reuses the same mechanics as `_unsaturated.py` (ending replaces 'ane'
  entirely, multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene'
  before 'yne', euphonic stem 'a' before a multiplied ending).
- P-14.3.4.2(a)/(b) (Chapter P-1, as already applied to halogens in
  `_acyclic.py`): a substituent prefix's own locant is omitted from a
  mononuclear (one-carbon) parent regardless of how many
  substituents/suffixes it carries ('methanol'); the -OH suffix's own
  locant is additionally omitted from a saturated two-carbon chain
  regardless of any other substituent present or where it sits, since
  P-44.4.1.8 always forces the suffix locant to be C1 there ('ethanol';
  'ethane-1,2-diol' still needs locants because it has two -OH's, not
  one; a substituent prefix's own locant, unlike the suffix's, is *not*
  omitted at two carbons since its position is still real distinguishing
  information -- '2-chloroethanol', not 'chloroethanol'). Both forms
  confirmed directly via PubChem PUG REST.
- P-35.2.1 (Chapter P-3): halogen substituents are prefix-only and coexist
  freely with the -OH suffix (they never compete for suffix status), reusing
  `halogen_substituents`/`format_substituent_prefixes` unchanged.
- P-29.3.3 (Chapter P-2): on an acyclic chain, a simple alkoxy ether (-O-R,
  R a plain unbranched saturated alkyl group, e.g. '-OCH3') coexists freely
  with the -OH suffix too, cited as an ordinary 'methoxy'/'ethoxy'/...
  substituent prefix (also covers what P-66.6.5.2 calls a 'hemiacetal',
  RR'C(OH)(OR'') -- there is no distinct hemiacetal nomenclature, it's just
  this same alcohol-plus-alkoxy-ether combination). Confirmed via PubChem
  PUG REST: CID 3015637 (`CC(O)OC`) -> "1-methoxyethanol", CID 8109
  (`OCCCOCC`) -> "3-ethoxypropan-1-ol", CID 8107 (`OCCCCOC`) ->
  "4-methoxybutan-1-ol", CID 12486323 (`OCC(OC)COC`) ->
  "2,3-dimethoxypropan-1-ol" (more than one such ether is fine, each gets
  its own prefix). On a two-carbon chain this can collide with the
  documented P-14.3.4.2(a)/(b) locant-citation edge case above (e.g.
  PubChem's own '2-methoxyethanol' omits the -OH locant); this module
  keeps its existing policy of always citing the locant once any
  substituent is present, an accepted, reviewed divergence rather than a
  PubChem-confirmed one for that specific chain length, same as the
  halogen/ring-substituent cases already handled that way.

- P-31.1.3: a monocyclic ring bearing a hydroxyl and exactly one C=C ring
  double bond -- e.g. 'cyclohex-2-en-1-ol', 'cyclohex-3-en-1-ol', both
  confirmed via PubChem PUG REST. The hydroxyl's own locant is never
  omittable here even as the ring's sole substituent (a competing ring
  double bond means '1' must still be cited, unlike the fully saturated
  'cyclohexanol' case above); the numbering direction is then chosen to
  minimize the double bond's own locant -- mirrors `_ketone.py`'s
  identical extension (same underlying `_cyclic_unsaturated.py` pattern).
  Deliberately narrow: a ring triple bond, and any other substituent
  alongside the ring double bond, are both still explicitly rejected
  pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that is not an isolated, singly-bonded -OH with exactly one H,
  or a simple alkoxy ether as described above (any C=O — aldehyde, ketone,
  or the carbonyl of a carboxylic acid/amide/ester; a branched, cyclic, or
  unsaturated alkoxy R; an ether attached to a ring rather than an acyclic
  chain).
- Two or more *different* characteristic-group types (e.g. an alcohol and an
  amine) — not applicable here since only C, halogen, -OH-shaped oxygen,
  and simple-alkoxy-ether oxygen atoms are accepted at all; any other
  heteroatom (N, S, ...) is rejected.
- -OH on an aromatic ring (phenol-type) combined with a *second* hydroxyl
  elsewhere (a second ring -OH, or a separate chain -OH), or an alkoxy
  ether alongside it — a single ring -OH on an otherwise-plain benzene
  ring is handled by `_name_phenol` (P-63.1.1), a separate follow-up from
  this narrower multi-group combination.
- -OH directly on a heteroaromatic ring (pyridine/pyrrole/furan/thiophene)
  — unlike benzene's free numbering, a heteroaromatic ring's numbering
  must fix the heteroatom at locant 1, so `_name_phenol`'s retained-name
  construction does not apply there; only the chain-substituent shape
  (`_name_phenyl_chain_alcohol`, e.g. `pyridin-3-ylmethanol`) extends to
  a heteroaromatic ring (P-616 M2 step 8), mirroring `_thiol.py`'s/
  `_selenol.py`'s/`_tellurol.py`'s/`_ketone.py`'s identical
  generalization.
- -OH on a von Baeyer polycyclic (bicyclic through pentacyclic) or spiro
  skeleton — deferred; those modules' internal numbering would need real
  integration work to prioritize a suffix locant correctly.
- -OH on a carbon that is also part of a C=C/C#C bond (an enol) — this is a
  further, deliberate narrowing beyond what full generality would allow (the
  Blue Book's numbering machinery *would* handle this case, combining suffix
  and 'ene'/'yne' locant priority the same way as any other case), but the
  keto-enol-tautomer-adjacent naming nuances were not independently verified
  here, so it is scoped out rather than guessed at.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    heteroaromatic_monocycle_name,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereo_elements,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    two_separate_rings_with_plain_aromatic_substituent,
)
from ._bicyclic import (
    _candidate_key as _bicyclic_candidate_key,
    bicyclic_parent_name,
    find_bicyclic_core,
    iter_bicyclic_numberings,
)
from ._cyclic_unsaturated import name_cyclic_unsaturated_yl
from ._numerals import alkyl_name
from ._polycyclic import (
    _candidate_key as _polycyclic_candidate_key,
    find_polycyclic_core,
    iter_polycyclic_candidates,
)
from ._spiro import (
    _candidate_key as _spiro_candidate_key,
    find_monospiro_atom,
    iter_monospiro_numberings,
)
from ._substituents import (
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    ring_branch_stereo_display,
    substituents_for_chain,
    substituents_for_ring,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _plain_terminal_alkyl_length(mol, graph, start, coming_from):
    """Walk a branch outward, like `linear_branch`, but only succeeds if
    every atom on the way (including `start`) is a plain carbon -- returns
    None as soon as a heteroatom, a fork, or a cycle is hit. Used to tell
    a simple terminal alkyl group (the 'R' of an alkoxy substituent, e.g.
    '-OCH3') apart from a direction that actually leads back into the rest
    of the molecule (e.g. toward a coexisting -OH)."""
    length = 0
    previous, current, visited = coming_from, start, set()
    while True:
        if current in visited:
            return None
        visited.add(current)
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None
        length += 1
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) == 0:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]


def _ether_oxygens(mol, graph):
    """Simple ethers (-O-, single-bonded to two carbons) where exactly one
    side is a plain unbranched terminal alkyl chain (P-29.3.3's 'alkoxy'
    substituent, e.g. 'methoxy', 'ethoxy') -- {ether_o_idx: alkoxy_name}.
    A degree-2 oxygen where neither or both sides qualify (an ambiguous or
    unsupported shape, e.g. a symmetrical ether) is left out here and
    falls through to the generic heteroatom rejection below."""
    ethers = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 2:
            continue
        if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
            continue
        neighbors = atom.GetNeighbors()
        if any(n.GetAtomicNum() != 6 for n in neighbors):
            continue
        n1, n2 = neighbors
        len1 = _plain_terminal_alkyl_length(mol, graph, n1.GetIdx(), atom.GetIdx())
        len2 = _plain_terminal_alkyl_length(mol, graph, n2.GetIdx(), atom.GetIdx())
        if (len1 is None) == (len2 is None):
            continue
        length = len1 if len1 is not None else len2
        ethers[atom.GetIdx()] = _alkoxy_name(length)
    return ethers


_CONTRACTED_ALKOXY_NAMES = {1: "methoxy", 2: "ethoxy", 3: "propoxy", 4: "butoxy"}


def _alkoxy_name(length):
    """P-29.3.3: the alkoxy prefix for C1-C4 is a contracted retained form
    ('methoxy', not 'methyloxy'); C5 and up uses the full alkyl name plus
    'oxy' ('pentyloxy')."""
    if length in _CONTRACTED_ALKOXY_NAMES:
        return _CONTRACTED_ALKOXY_NAMES[length]
    return alkyl_name(length) + "oxy"


def _validate_and_collect_hydroxyls(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (hydroxyls, ethers): the set of hydroxyl-oxygen atom
    indices, and a {ether_o_idx: alkoxy_name} dict for any simple alkoxy
    ether (see `_ether_oxygens`).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring or heteroaromatic monocycle (pyridine/furan/thiophene/pyrrole)
    with exactly one exocyclic attachment -- exempted here wholesale from
    the per-atomic-number checks below (already independently verified by
    that shape check itself) so `name_alcohol`'s benzene-ring-substituent
    path (see `_name_phenyl_chain_alcohol`) and its `two_separate_rings_
    with_plain_aromatic_substituent` shape (#624, mirroring `_ketone.py`'s
    #622) can both reuse this same validation for the rest of the
    molecule. Empty by default, so every other caller's behavior is
    unchanged."""
    graph = adjacency(mol)
    ethers = _ether_oxygens(mol, graph)
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in aromatic_ring_atoms:
            if atomic_num == 6:
                has_carbon = True
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a hydroxyl oxygen (P-33.2.1), a "
                "simple alkoxy ether (P-29.3.3), and halogen substituents "
                "(P-35.2.1) are not supported yet"
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
            if atom.GetIdx() in ethers:
                continue
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope unless it is a simple alkoxy "
                    "ether (P-29.3.3) with a plain unbranched terminal "
                    "alkyl on exactly one side"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "an oxygen double-bonded to carbon indicates a more "
                    "senior characteristic group (aldehyde, ketone, or a "
                    "carboxylic acid) than a plain alcohol, which this "
                    "module does not attempt to disambiguate (Table 3.3's "
                    "suffix seniority order, P-33.2.1)"
                )
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "an -O- atom that isn't a simple hydroxyl (-OH) is out "
                    "of scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a hydroxyl must be attached to a carbon atom")
            hydroxyls.add(atom.GetIdx())
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
    if not hydroxyls:
        raise UnsupportedStructure("no hydroxyl (-OH) group found; this module only handles alcohols")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return hydroxyls, ethers


def _reject_enol_carbons(graph, hydroxyls, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o_idx in hydroxyls:
        (carbon,) = graph[o_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol-type structure) is out of scope for this module"
            )


def _name_from_substituents(chain_length, oh_locants, ene_locants, yne_locants, grouped):
    # Only a truly mononuclear parent (a single carbon) makes every
    # substituent's own locant trivially '1' -- at chain_length == 2 the
    # suffix's own locant is separately omittable (see
    # `name_from_substituents`/`should_omit_mononuclear_locants`), but a
    # substituent prefix can still sit at either carbon, so its locant
    # stays meaningful and must not be omitted, e.g. '2-aminoethanethiol'
    # (PubChem-verified) keeps its '2-' even though the thiol locant drops.
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return prefix + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(len(oh_locants), "ol"), oh_locants
    )


def _candidate_key(chain_length, oh_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    oh_locant_set = lowest_locant_set(oh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, oh_locants, ene_locants, yne_locants, grouped)
    return (
        (
            oh_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _oh_locants(position_of, hydroxyls, graph):
    locants = []
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants




def _best_acyclic_alcohol_candidate(
    mol, hydroxyls, bonds, stereo=None, ethers=None, extra_names=None, required_atoms=frozenset()
):
    """(best_name, best_position_of) -- the winning chain numbering
    (P-44.4.1.8's own suffix-locant-first tie-break) and its
    fully-formatted name, factored out of `_name_acyclic_alcohol` (which
    just adds the stereo-descriptor prefix on top) so `_isotope_alcohol.py`
    can reuse the identical numbering decision: P-82.5.1 states the
    presence of isotopic nuclides is "considered last" among numbering
    criteria, so an isotope label never changes which candidate wins here,
    and that module only needs `best_position_of` to independently compute
    its own isotope descriptor's locant.

    `ethers`: optional {ether_o_idx: alkoxy_name} (see `_ether_oxygens`)
    -- merged into `halogens` so `name_branch` resolves each ether oxygen
    directly to its alkoxy prefix name instead of recursing into it (its
    far-side terminal alkyl is already invisible to `carbon_adjacency`,
    so it never competes for the principal chain).

    `extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection.
    `required_atoms`: additional carbon atoms (e.g. every demoted amine's
    own carbon neighbor) that a candidate chain must also carry -- both
    empty/None by default so existing callers are unaffected."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(ethers or {}), **(extra_names or {})}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [idx for kind, idx, _ in stereo if kind == "atom"] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _oh_locants(position_of, hydroxyls, graph) is None:
            continue
        if not required_atoms <= set(chain):
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in chain for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _oh_locants({a: i + 1 for i, a in enumerate(c)}, hydroxyls, graph) is not None
            and required_atoms <= set(c)
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every hydroxyl-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more -OH groups (P-44.1.1), or an -OH expressed as a "
            "'hydroxy' substituent prefix, is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locants = _oh_locants(position_of, hydroxyls, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = substituents_for_chain(graph, candidate, halogens, hydroxyls, mol=mol)
            key, name = _candidate_key(chain_length, oh_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    return best_name, best_position_of


def _name_acyclic_alcohol(
    mol, hydroxyls, bonds, stereo=None, ethers=None, extra_names=None, required_atoms=frozenset()
):
    """`stereo`: None, or a list of ("atom"/"bond", idx, "R"/"S"/"E"/"Z")
    from `specified_stereo_elements` -- if given, only chain candidates
    that include *every* tetrahedral stereocenter are eligible (P-92: a
    stereocenter on a substituent branch rather than the principal chain
    is out of scope, see module docstring; a double-bond E/Z element's
    atoms are already required to lie on the chain via `bonds`, so no
    separate check is needed for those), and the winning candidate's own
    locants for each element are used to format a
    "(<locant><R/S/E/Z>,...)-" prefix onto the name, ascending locant
    order (P-91.3, including when both kinds coexist). See
    `_best_acyclic_alcohol_candidate` for `ethers`/`extra_names`/
    `required_atoms`."""
    best_name, best_position_of = _best_acyclic_alcohol_candidate(
        mol, hydroxyls, bonds, stereo, ethers, extra_names, required_atoms
    )
    if stereo is not None:
        labels = []
        for kind, idx, code in stereo:
            if kind == "atom":
                locant = best_position_of[idx]
            else:
                bond = mol.GetBondWithIdx(idx)
                locant = min(best_position_of[bond.GetBeginAtomIdx()], best_position_of[bond.GetEndAtomIdx()])
            labels.append((locant, code))
        labels.sort()
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _ring_name_from_substituents(ring_size, oh_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size, ene_locants, yne_locants, prefix, total_subs, multiplied_word(len(oh_locants), "ol"), oh_locants
    )


def _ring_candidate_key(ring_size, oh_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    oh_locant_set = lowest_locant_set(oh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, oh_locants, ene_locants, yne_locants, grouped)
    return oh_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _phenol_name_from_substituents(grouped):
    # Like the mancude-ring 'benzenesulfonic acid'/'benzoic acid' cases
    # (`_sulfonic_acid.py`/`_carboxylic_acid.py`), the retained name
    # 'phenol' stands for the whole ring+OH system, so the ring's own
    # numbering is always free to start at the -OH carbon and that
    # locant is never cited, even when other substituents need theirs --
    # e.g. '4-methylphenol' (PubChem CID 2879), not '4-methylphenol-1-ol'.
    if not grouped:
        return "phenol"
    return f"{format_substituent_prefixes(grouped)}phenol"


def _phenol_candidate_key(oh_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _phenol_name_from_substituents(grouped)
    return oh_locant, locant_set, citation_locants, name


def _name_phenol(mol, ring_atoms):
    """P-63.1.1: -OH attached directly to a benzene ring carbon -- e.g.
    'phenol' (PubChem CID 996), '4-methylphenol' (CID 2879), '2-
    chlorophenol' (CID 7245). Mirrors `_sulfonic_acid.py`'s
    `_name_benzenesulfonic_acid` exactly, with the retained name 'phenol'
    standing in for 'benzenesulfonic acid'. Narrower than the general
    ring case: exactly one hydroxyl (directly on the ring, no coexisting
    alkoxy ether), and no specified stereocenter -- more than one ring
    hydroxyl (resorcinol-style) is a separate follow-up."""
    hydroxyls, ethers = _validate_and_collect_hydroxyls(mol, aromatic_ring_atoms=ring_atoms)
    if ethers:
        raise UnsupportedStructure("an alkoxy ether alongside phenol is not supported yet")
    if len(hydroxyls) != 1:
        raise UnsupportedStructure(
            "more than one hydroxyl directly on a benzene ring is not "
            "supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure("a specified stereocenter alongside phenol is not supported yet")

    (oh_oxygen,) = hydroxyls
    graph = adjacency(mol)
    (oh_carbon,) = graph[oh_oxygen]
    halogens = halogen_substituents(mol)
    excluded = {oh_oxygen}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locant = position_of[oh_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _phenol_candidate_key(oh_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _ring_branch_stereo_display(graph, ring_order, hydroxyls, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    return ring_branch_stereo_display(graph, ring_order, hydroxyls, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms)


def _name_cyclic_alcohol(mol, hydroxyls, stereo=None, bonds=(), ring_atoms=None, aromatic_atoms=frozenset()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must
    normally lie on the ring itself (P-92: a stereocenter on a
    substituent branch is out of scope, mirroring
    `_name_acyclic_alcohol`'s identical chain-only restriction), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3) -- same mechanism as the acyclic case, since
    P-92 doesn't affect which numbering wins. The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_aromatic.py`'s own single-
    branch-stereocenter case): exactly one stereocenter on the ring's
    sole substituent branch instead embeds a bracketed descriptor into
    that substituent's own name, in place of the usual ring-locant
    prefix.
    `bonds`: ring C=C double bonds (P-31.1.3), empty by default -- see
    `_ketone.py`'s identical `bonds` parameter for the shared reasoning
    (hydroxyl locant fixed first, then minimized ene locant).

    `ring_atoms`/`aromatic_atoms`: when the hydroxyl-bearing ring reaches
    this function as one half of a `two_separate_rings_with_plain_
    aromatic_substituent` shape (#624, mirroring `_ketone.py`'s #622), the
    caller passes the hydroxyl-bearing ring's own atoms explicitly
    (RDKit's SSSR would otherwise list either of the two disjoint rings
    first) along with the other, aromatic ring's atoms, cited as a plain
    substituent (phenyl or a heteroaromatic monocycle) via `name_branch`
    the same way an ordinary alkyl ring substituent already is.
    `ring_atoms=None` (default) preserves the original single-ring
    dispatch unchanged."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    if ring_atoms is None:
        ring_atoms = list(mol.GetRingInfo().AtomRings()[0])
    else:
        ring_atoms = list(ring_atoms)
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(
            graph, ring_order, hydroxyls, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms
        )
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(
        substituents_for_ring(graph, ring_order, halogens, hydroxyls, mol=mol, aromatic_atoms=aromatic_atoms).values()
    ):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "hydroxyl is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locants = _oh_locants(position_of, hydroxyls, graph)
            if oh_locants is None:
                raise UnsupportedStructure(
                    "a hydroxyl not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = substituents_for_ring(
                graph, candidate, halogens, hydroxyls, mol=mol, aromatic_atoms=aromatic_atoms
            )
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, oh_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _von_baeyer_core_atoms(bicyclic_core, polycyclic_core):
    if bicyclic_core is not None:
        bh1, bh2, bridges = bicyclic_core
        return {bh1, bh2} | {atom for bridge in bridges for atom in bridge}
    branch_atoms, bridges = polycyclic_core
    return set(branch_atoms) | {atom for _, _, path in bridges for atom in path}


def _name_von_baeyer_alcohol(mol, hydroxyls, stereo, bonds, bicyclic_core, polycyclic_core, ring_count):
    """P-23.2.1's bicyclic/polycyclic ('bicyclo[x.y.z]alkane'/
    'tricyclo[...]alkane'/...) numbering (P-23.2.3-P-23.2.6), extended
    with a single -OH suffix the same way `_von_baeyer_heteroatom.py`
    extends it with a skeletal replacement heteroatom: iterate every von
    Baeyer-valid numbering (`iter_bicyclic_numberings`/
    `iter_polycyclic_candidates`), compute the -OH carbon's own locant per
    candidate, and rank it ahead of substituent locants via
    `_bicyclic_candidate_key`/`_polycyclic_candidate_key`'s new
    `suffix_locant` parameter. Unlike a chain's own principal
    characteristic-group carbon (always the fixed C1, its locant never
    cited, P-14.3.3), a von Baeyer parent's atom 1 is always a bridgehead
    by P-23.2.3's own numbering rule regardless of where the -OH sits, so
    the -OH's own locant is never structurally forced and is always
    explicitly cited, e.g. 'bicyclo[2.2.1]heptan-2-ol' (confirmed against
    PubChem's own IUPACName for CID 19809, 'OC1CC2CCC1C2')."""
    if len(hydroxyls) != 1:
        raise UnsupportedStructure(
            "more than one hydroxyl on a von Baeyer bicyclic/polycyclic "
            "ring system is not supported yet"
        )
    if bonds:
        raise UnsupportedStructure(
            "an unsaturated von Baeyer bicyclic/polycyclic ring system is "
            "not supported yet (see P-31.1.4)"
        )
    if stereo is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a von Baeyer bicyclic/"
            "polycyclic alcohol is not supported yet (see P-92)"
        )

    (oh_oxygen,) = hydroxyls
    graph = adjacency(mol)
    (oh_carbon,) = graph[oh_oxygen]
    core_atoms = _von_baeyer_core_atoms(bicyclic_core, polycyclic_core)
    if oh_carbon not in core_atoms:
        raise UnsupportedStructure(
            "a hydroxyl not on the bicyclic/polycyclic ring system itself "
            "(e.g. on a substituent branch) is not supported yet"
        )

    halogens = halogen_substituents(mol)
    best_key = None
    if bicyclic_core is not None:
        base_parent = bicyclic_parent_name(bicyclic_core)
        for full_order in iter_bicyclic_numberings(bicyclic_core):
            oh_locant = full_order.index(oh_carbon) + 1
            substituents = substituents_for_ring(graph, full_order, halogens, hydroxyls, mol=mol)
            parent = base_parent[:-1] + f"-{oh_locant}-ol"
            key = _bicyclic_candidate_key(parent, substituents, suffix_locant=oh_locant)
            if best_key is None or key < best_key:
                best_key = key
        return best_key[-1]

    for full_order, parent, outer_key in iter_polycyclic_candidates(polycyclic_core, ring_count):
        oh_locant = full_order.index(oh_carbon) + 1
        substituents = substituents_for_ring(graph, full_order, halogens, hydroxyls, mol=mol)
        suffixed_parent = parent[:-1] + f"-{oh_locant}-ol"
        key = outer_key + _polycyclic_candidate_key(suffixed_parent, substituents, suffix_locant=oh_locant)
        if best_key is None or key < best_key:
            best_key = key
    if best_key is None:
        raise UnsupportedStructure(
            "this polycyclic topology is not supported yet (disjoint ring "
            "systems joined only by an acyclic linker are out of scope; "
            "see _polycyclic.py's name_polycycloalkane for the analogous "
            "non-alcohol guard)"
        )
    return best_key[-1]


def _name_monospiro_alcohol(mol, hydroxyls, stereo, spiro_atom):
    """P-24.2.1's monospiro numbering (`_spiro.py`), extended with a
    single -OH suffix the same way `_name_von_baeyer_alcohol` extends the
    von Baeyer bicyclic/polycyclic numbering above: iterate every
    P-24.2.1-valid numbering (`iter_monospiro_numberings`), compute the
    -OH carbon's own locant per candidate, and rank it ahead of
    substituent locants via `_spiro_candidate_key`'s `suffix_locant`
    parameter. Confirmed against a real registered structure (PubChem CID
    90762054, 'OC1CCCC2(C1)CCCCC2')."""
    if len(hydroxyls) != 1:
        raise UnsupportedStructure(
            "more than one hydroxyl on a monospiro ring system is not "
            "supported yet"
        )
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "an unsaturated monospiro ring system is not supported yet "
            "(see P-31.1.5)"
        )
    if stereo is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a monospiro alcohol is "
            "not supported yet (see P-92)"
        )

    (oh_oxygen,) = hydroxyls
    graph = adjacency(mol)
    (oh_carbon,) = graph[oh_oxygen]
    ring_atoms = {atom for ring in mol.GetRingInfo().AtomRings() for atom in ring}
    if oh_carbon not in ring_atoms:
        raise UnsupportedStructure(
            "a hydroxyl not on the monospiro ring system itself (e.g. on "
            "a substituent branch) is not supported yet"
        )

    halogens = halogen_substituents(mol)
    best_key = None
    for parent, full_order in iter_monospiro_numberings(mol, spiro_atom):
        oh_locant = full_order.index(oh_carbon) + 1
        substituents = substituents_for_ring(graph, full_order, halogens, hydroxyls, mol=mol)
        suffixed_parent = parent[:-1] + f"-{oh_locant}-ol"
        key = _spiro_candidate_key(suffixed_parent, substituents, suffix_locant=oh_locant)
        if best_key is None or key < best_key:
            best_key = key
    return best_key[-1]


def _name_ring_substituent_chain_alcohol(mol, hydroxyls):
    """Name an alcohol whose -OH lies entirely on a single unbranched chain
    hanging off one atom of an otherwise-plain monocyclic ring (the ring
    itself bears no -OH) -- e.g. cyclohexylmethanol. The ring is cited as a
    "cyclo..." substituent prefix (P-29.3.3) on the chain, which is the
    parent hydride (see module docstring). The ring may carry ring C=C
    unsaturation of its own (P-31.1.3 + P-29.2), named via
    `name_cyclic_unsaturated_yl` with the attachment fixed at locant 1 --
    e.g. '(cyclohex-3-en-1-yl)methanol'."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, hydroxyls)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    anchor_oxygen = next(iter(hydroxyls))
    (anchor_carbon,) = graph[anchor_oxygen]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, hydroxyls, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a hydroxyl outside the single unbranched chain hanging "
                "off the ring is not supported yet"
            )
    # `ring_atom` is named specially below (`ring_name`, aware of ring
    # unsaturation via `name_cyclic_unsaturated_yl`) -- drop it here so
    # the generic per-position loop doesn't also re-name it plainly.
    branches_by_atom = {
        chain[position - 1]: [r for r in roots if r != ring_atom]
        for position, roots in branches.items()
    }
    branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}

    ring_has_double_bond = any(
        bond.GetBondTypeAsDouble() == 2.0
        and bond.GetBeginAtomIdx() in ring_atoms
        and bond.GetEndAtomIdx() in ring_atoms
        for bond in mol.GetBonds()
    )
    if ring_has_double_bond:
        ring_name = name_cyclic_unsaturated_yl(mol, ring_atoms, ring_atom)
        ring_is_compound = True
    else:
        ring_name = "cyclo" + alkyl_name(len(ring_atoms))
        ring_is_compound = False
    chain_length = len(chain)

    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        oh_locants = _oh_locants(position_of, hydroxyls, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, ring_is_compound))
        key, name = _candidate_key(chain_length, oh_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_hydroxy_chain_alcohol(mol, hydroxyls):
    """Name an alcohol where the ring itself bears at least as many -OH's
    as the chain hanging off exactly one ring atom does (P-44.1.1: the
    candidate with the greater count of the principal characteristic group
    -OH is senior; P-44.1.2.2 resolves an exact tie in the ring's favor)
    -- e.g. 2-(hydroxymethyl)cyclohexan-1-ol,
    4-(hydroxymethyl)cyclohexane-1,2-diol,
    4-(1,2-dihydroxyethyl)cyclohexane-1,2-diol. The ring is the parent; the
    chain is cited as a '(hydroxy...alkyl)' substituent prefix, reusing the
    {oxygen_idx: "hydroxy"} trick already used by
    `_carboxylic_acid.py`/`_amide.py`/`_aldehyde.py`/`_ketone.py` -- mapping
    every chain hydroxyl this way lets `name_branch` group and multiply the
    "hydroxy" prefix exactly as it already does for repeated halogens.

    The chain hanging off the ring may itself be branched (P-44.3.2: the
    chain-wins case below uses `longest_branched_chain_through`, the same
    general branch-absorbing search `_name_ring_substituent_chain_alcohol`
    already uses for its own, ring-has-no-competing-OH case) -- no
    PubChem-registered example of a *carbon*-branched competing chain was
    found (most real compounds with this shape are unbranched, e.g.
    PubChem CID 82849893's '3-(2-hydroxycyclohexyl)propane-1,2-diol'), but
    the rule is the identical P-44.3.2 search already verified for the
    unbranched and ring-substituent-only cases."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, hydroxyls)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    ring_hydroxyls = {o for o in hydroxyls if next(iter(graph[o])) in ring_atoms}
    chain_hydroxyls = hydroxyls - ring_hydroxyls

    if len(ring_hydroxyls) < len(chain_hydroxyls):
        # P-44.1.1: the chain captures strictly more -OH's, so it's the
        # senior parent and the ring (with its own one or more -OH's) is
        # cited as a substituent instead -- mirrors
        # `_name_ring_substituent_chain_alcohol` exactly, substituting
        # the ring's own name_branch-computed name for the plain
        # "cyclo..." one that function uses.
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{o: "hydroxy" for o in ring_hydroxyls}}, mol=mol
        )
        chain, branches = longest_branched_chain_through(graph, chain_root, ring_atoms, hydroxyls, halogens=halogen_substituents(mol))
        branches_by_atom = {
            chain[position - 1]: [r for r in roots if r != ring_atom]
            for position, roots in branches.items()
        }
        branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locants = _oh_locants(position_of, chain_hydroxyls, graph)
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            substituents.setdefault(position_of[chain_root], []).append((ring_name, ring_is_compound))
            key, name = _candidate_key(chain_length, oh_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{o: "hydroxy" for o in chain_hydroxyls}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            oh_locants = _oh_locants(position_of, ring_hydroxyls, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, oh_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_alcohol(mol, ring_atoms):
    """Name an alcohol whose -OH lies entirely on a single unbranched chain
    hanging off one atom of a benzene or heteroaromatic-monocycle ring --
    e.g. 2-phenylethanol, pyridin-3-ylmethanol (P-616 M2 step 8). The ring
    is cited as a 'phenyl' (or, if the ring's other atoms each carry a
    single halogen, e.g. '4-chlorophenyl') or heteroaromatic substituent
    prefix (via `name_branch`'s aromatic-ring recognition) on the chain,
    which is the parent hydride, mirroring `_carboxylic_acid.py`'s
    `_name_phenyl_chain_carboxylic_acid` and this module's own
    `_name_ring_substituent_chain_alcohol` for a plain saturated ring.
    Narrower than either: exactly one -OH, no coexisting alkoxy ether, no
    chain unsaturation, no specified stereocenter, and no non-halogen ring
    substituent alongside the chain -- each is a separate follow-up rather
    than being combined with this first slice."""
    hydroxyls, ethers = _validate_and_collect_hydroxyls(mol, aromatic_ring_atoms=ring_atoms)
    if ethers:
        raise UnsupportedStructure(
            "an alkoxy ether alongside a benzene-ring-substituent alcohol "
            "chain is not supported yet"
        )
    if len(hydroxyls) != 1:
        raise UnsupportedStructure(
            "more than one hydroxyl alongside a benzene-ring substituent "
            "is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "alcohol chain is not supported yet"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "alcohol chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain alcohol is not "
            "supported yet"
        )
    ring_atom, chain_root = attachment
    if chain_root in hydroxyls:
        raise UnsupportedStructure(
            "a hydroxyl directly on the benzene ring (phenol-type) "
            "combined with a second hydroxyl elsewhere is not supported "
            "yet -- a single ring hydroxyl alone is handled by "
            "`_name_phenol`, out of scope for this chain-parent module"
        )
    anchor_oxygen = next(iter(hydroxyls))
    (anchor_carbon,) = graph[anchor_oxygen]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, hydroxyls, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a hydroxyl outside the single unbranched chain hanging "
                "off the benzene ring is not supported yet"
            )
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        oh_locants = _oh_locants(position_of, hydroxyls, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, oh_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_alcohol(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            hydroxyls, ethers = _validate_and_collect_hydroxyls(mol, aromatic_ring_atoms=ring_atoms)
            graph = adjacency(mol)
            ring_hydroxyls = {o for o in hydroxyls if next(iter(graph[o])) in ring_atoms}
            if is_benzene and not ethers and hydroxyls == ring_hydroxyls and len(hydroxyls) == 1:
                return _name_phenol(mol, ring_atoms)
            if is_heteroaromatic and ring_hydroxyls:
                raise UnsupportedStructure(
                    "a hydroxyl directly on a heteroaromatic ring is not "
                    "supported yet (a heteroaromatic ring's numbering must "
                    "fix the heteroatom at locant 1, unlike benzene's free "
                    "numbering, so 'phenol'-style naming does not apply "
                    "here)"
                )
            return _name_phenyl_chain_alcohol(mol, ring_atoms)

    # Two separate simple monocycles joined by one direct bond, one a plain
    # benzo/heteroaromatic ring with no substituent of its own (P-25 M1
    # step 3, #624) -- e.g. 2-phenylcyclohexan-1-ol -- reuses the existing
    # single-ring `_name_cyclic_alcohol` dispatch below with the aromatic
    # ring's atoms passed through as an exemption/substituent, the same
    # generalization #622 made for `_ketone.py`. Narrower than that: only
    # the "every hydroxyl is ring-borne" split is handled below (matching
    # this shape's only verified real structures) -- a hydroxyl entirely
    # on a chain hanging off the non-aromatic ring, with the ring itself
    # bearing no -OH, would need the ring cited as a substituent on that
    # chain, but the ring would then carry *two* exocyclic branches (the
    # aromatic ring and the chain) and would need its own name_branch-
    # computed compound name instead of a plain "cyclo..." prefix --
    # deferred as a separate, more involved follow-up.
    aromatic_shape = None
    if ring_info.NumRings() == 2:
        aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
    aromatic_atoms = aromatic_shape[1] if aromatic_shape is not None else frozenset()

    hydroxyls, ethers = _validate_and_collect_hydroxyls(mol, aromatic_ring_atoms=aromatic_atoms)
    graph = adjacency(mol)
    # Exclude a bond entirely inside the aromatic-substituent ring above --
    # its aromatic bond order (1.5) isn't a real chain/ring 'ene'/'yne'
    # bond either, and is already accounted for by naming that ring via
    # `name_branch` instead (mirrors `_ketone.py`'s identical exclusion).
    all_non_single = [
        b for b in non_single_bonds(mol) if not (b[0] in aromatic_atoms and b[1] in aromatic_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enol_carbons(graph, hydroxyls, bonds)

    num_rings = ring_info.NumRings()
    if num_rings == 0:
        # A ring is never involved past this point, so a specified C=C
        # double-bond E/Z element may coexist with a specified tetrahedral
        # stereocenter (P-91.3) -- every other branch below keeps using
        # `specified_stereocenters`, which still rejects that combination
        # (rings never have both: unsaturated rings are rejected outright
        # just below).
        stereo = specified_stereo_elements(mol)
        return _name_acyclic_alcohol(mol, hydroxyls, bonds, stereo, ethers)

    stereo = specified_stereocenters(mol)
    if ethers:
        raise UnsupportedStructure(
            "an alkoxy ether coexisting with a cyclic alcohol structure "
            "is not supported yet"
        )
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_hydroxyls = {o for o in hydroxyls if next(iter(graph[o])) in ring_atoms}
        chain_hydroxyls = hydroxyls - ring_hydroxyls
        if bonds:
            # A ring C=C double bond coexisting with a suffix -OH is
            # supported for the "chain is parent, ring has no -OH of its
            # own" shape (P-31.1.3 + P-29.2, see
            # `_name_ring_substituent_chain_alcohol`) and for the "-OH is
            # on the ring itself, no separate chain -OH" shape (P-31.1.3,
            # see `_name_cyclic_alcohol`'s own `bonds` parameter, mirroring
            # `_ketone.py`'s identical extension); a ring double bond
            # alongside *both* a ring -OH and a separate chain -OH, or
            # unsaturation reaching outside the ring (an exocyclic double
            # bond, or a triple bond), stays unsupported.
            ring_only_double_bonds = all(
                order == _ENE_ORDER and a in ring_atoms and b in ring_atoms
                for a, b, order in bonds
            )
            if not ring_only_double_bonds or (ring_hydroxyls and chain_hydroxyls):
                raise UnsupportedStructure(
                    "unsaturated rings are not supported yet (see P-31.1.3, "
                    "cycloalkenes and cycloalkynes)"
                )
        if not ring_hydroxyls:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_alcohol(mol, hydroxyls)
        if not chain_hydroxyls:
            return _name_cyclic_alcohol(mol, hydroxyls, stereo, bonds)
        if stereo is not None:
            raise UnsupportedStructure(
                "a stereocenter alongside a ring-vs-chain hydroxyl "
                "comparison is not supported yet (see P-92)"
            )
        return _name_ring_with_hydroxy_chain_alcohol(mol, hydroxyls)
    if num_rings == 2 and aromatic_shape is not None:
        ring_atoms, _, _, _ = aromatic_shape
        chain_hydroxyls = hydroxyls - {o for o in hydroxyls if next(iter(graph[o])) in ring_atoms}
        if chain_hydroxyls:
            raise UnsupportedStructure(
                "a hydroxyl on a chain hanging off the ring, with the ring "
                "itself bearing no hydroxyl of its own, alongside this "
                "two-ring aromatic-substituent shape is not supported yet"
            )
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_alcohol(
            mol, hydroxyls, stereo, bonds, ring_atoms=ring_atoms, aromatic_atoms=aromatic_atoms
        )

    # A von Baeyer bicyclic or polycyclic (ring_count>=3) skeleton (P-23,
    # see `_bicyclic.py`/`_polycyclic.py`) tries this dedicated numbering
    # extension before the generic spiro/fused rejection below -- a spiro
    # atom's degree-4 shape (shared by both its rings) never matches
    # `find_bicyclic_core`/`find_polycyclic_core`'s own degree-2/3 core
    # requirement, so this never misfires on a genuine spiro/fused system.
    bicyclic_core = find_bicyclic_core(mol)
    polycyclic_core = None
    von_baeyer_ring_count = None
    if bicyclic_core is None:
        for candidate_ring_count in (3, 4, 5, 6):
            polycyclic_core = find_polycyclic_core(mol, candidate_ring_count)
            if polycyclic_core is not None:
                von_baeyer_ring_count = candidate_ring_count
                break
    if bicyclic_core is not None or polycyclic_core is not None:
        return _name_von_baeyer_alcohol(
            mol, hydroxyls, stereo, bonds, bicyclic_core, polycyclic_core, von_baeyer_ring_count
        )

    # A monospiro skeleton (P-24.2.1, `_spiro.py`) also tries its own
    # dedicated numbering extension before the generic fused-ring
    # rejection below -- `find_monospiro_atom` only matches the exact
    # two-rings-sharing-one-atom shape, so this never misfires on a
    # fused (non-spiro) polycyclic system.
    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return _name_monospiro_alcohol(mol, hydroxyls, stereo, spiro_atom)

    if stereo is not None:
        raise UnsupportedStructure(
            "a stereocenter on a polycyclic/spiro skeleton is not "
            "supported yet (see P-92)"
        )
    raise UnsupportedStructure(
        "polycyclic and fused-ring alcohols are not supported yet (P-23/"
        "P-25 numbering integration with a suffix group is future work)"
    )
