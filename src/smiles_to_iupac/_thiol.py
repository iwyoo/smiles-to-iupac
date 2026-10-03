"""Naming of thiols (the '-thiol' suffix, -SH) on acyclic saturated or
unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.1.1 / P-65.1.1, Table 3.3 (Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'thiol' is the preselected
  suffix for -SH, the sulfur analogue of 'ol' — structurally the same
  substitutive-suffix construction, just with S instead of O, and (per
  Table 3.3's overall seniority order: acids > esters > amides > nitriles >
  aldehydes > ketones > alcohols > thiols > amines) ranked directly below
  alcohols and above amines. Unlike 'ol', 'thiol' begins with a consonant,
  so the parent hydride's final 'e' is never elided before it (P-16.3.3):
  'methane' + 'thiol' -> 'methanethiol', not 'methanthiol'.
- P-44.1.1 / P-44.4.1 / P-45.2 (Chapter P-4): same principal-chain and
  numbering machinery as `_alcohol.py` (maximum number of -SH groups first,
  then -SH locants ahead of ene/yne locants ahead of substituent-prefix
  locants) — this module deliberately mirrors that module's structure so
  the two stay easy to compare.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_alcohol.py` apply here too (mononuclear parent, or a saturated
  two-carbon chain, regardless of other substituents), e.g.
  'ethanethiol'.
- P-35.2.1 (Chapter P-3): halogen substituents are prefix-only and coexist
  freely with the -SH suffix, same as in `_alcohol.py`.
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters -- every one on the
  principal chain/ring itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, e.g. '(2R)-butane-2-thiol',
  '(1S,2R)-2-methylcyclohexane-1-thiol', same pattern as `_alcohol.py`/
  `_ketone.py`/`_sulfonic_acid.py` (chain and ring both). Unlike
  `_sulfinic_acid.py`'s sulfur, a thiol's -SH sulfur is monovalent (bonded
  only to its one carbon and one H) and can never itself be a
  stereocenter, so this support is unconditional, not chemistry-limited.

Scope, deliberately narrow (first pass at this functional group, mirroring
how `_amide.py`/`_nitrile.py`/etc. each started in isolation before any
cross-suffix seniority work): one or more -SH groups (P-63.1.1's
dithiol/trithiol/... multiplication) on an acyclic chain or on a single
saturated carbon ring, mirroring `_alcohol.py`'s polyol and monocyclic
support, with no other heteroatom (in particular no -OH or amine
nitrogen) anywhere in the molecule -- Table 3.3's alcohol/thiol/amine
seniority coexistence is future work, tracked as a separate roadmap item,
same as the analogous `multi-carbonyl-seniority.md` split for aldehyde/
ketone.

P-31.1.3: a monocyclic ring bearing a thiol and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-thiol', 'cyclohex-3-ene-1-thiol',
both confirmed via PubChem PUG REST (note 'thiol' starts with a
consonant, so 'ene' never elides here, unlike the ketone/alcohol
'-en-1-one'/'-en-1-ol' pattern -- 'pent-4-ene-1-thiol' already confirms
the same rule on the acyclic path). The thiol's own locant is never
omittable here even as the ring's sole substituent, mirroring
`_ketone.py`/`_alcohol.py`'s identical extension. Deliberately narrow: a
ring triple bond, and any other substituent alongside the ring double
bond, are both still explicitly rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`):
more than one thiol, or a specified stereocenter, on a von Baeyer
polycyclic or spiro skeleton (a single -SH on such a skeleton is
supported, P-23/P-24 numbering integration via `_polycyclic_suffix.py`
with `elide_e=False`, since 'thiol' begins with a consonant); ring
unsaturation reaching outside the ring or a ring triple bond, an -SH on a
substituent branch off an otherwise-unsubstituted *saturated* ring, a
sulfide (-S- ether-analogue) or any other sulfur-oxidation-state group
(sulfonic acid, etc.), and any oxygen or nitrogen atom at all. Two
*aromatic*-ring cases:
`_name_phenyl_chain_thiol` names one or more -SH groups on a chain
hanging off a plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-thiol', '3-phenylpropane-1,2-dithiol'), mirroring
`_alcohol.py`/`_sulfonic_acid.py`'s identical benzene-ring-substituent
path -- narrower than the acyclic path: no chain unsaturation, no
specified stereocenter. `_name_benzenethiol` names a single -SH directly
on a benzene ring carbon (with or without other ring substituents), e.g.
'benzenethiol' (PubChem CID 7969), '2-methylbenzenethiol' (CID 8712),
mirroring `_name_cyclic_thiol`'s ring-numbering search with the retained
name 'benzene' as stem -- the -SH's own locant is never cited here,
unlike the cycloalkane case; two or more -SH groups directly on the ring
(a dithiophenol-style structure) remain out of scope.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    all_chains,
    bond_locant,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    heteroaromatic_monocycle_name,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    most_multiple_bonds,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ordered_chain,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_chain_attachments_with_halogens,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    two_separate_rings_with_plain_aromatic_substituent,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
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
_ALLOWED_ATOMIC_NUMS = {6, 16, *HALOGEN_PREFIXES}


def _validate_and_collect_thiols(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring or heteroaromatic monocycle (pyridine/furan/thiophene/
    pyrrole) with exactly one exocyclic attachment -- exempted here
    wholesale from the per-atomic-number checks below (already
    independently verified by that shape check itself) so `name_thiol`'s
    benzene-ring-substituent path (see `_name_phenyl_chain_thiol`) and its
    `two_separate_rings_with_plain_aromatic_substituent` shape (#628,
    mirroring `_ketone.py`'s #622/`_alcohol.py`'s #624) can both reuse
    this same validation for the rest of the molecule. Empty by default,
    so every other caller's behavior is unchanged."""
    thiols = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in aromatic_ring_atoms:
            if atomic_num == 6:
                has_carbon = True
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a thiol sulfur (P-63.1.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting -OH or amine nitrogen needs "
                "Table 3.3 seniority-coexistence handling not yet "
                "implemented for thiols"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num == 16:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a sulfur bonded to more than one heavy atom (e.g. a "
                    "sulfide) is out of scope; only an isolated thiol "
                    "(-SH) is supported (Table 3.3, P-63.1.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "a sulfur double-bonded to carbon is not a thiol"
                )
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "an -S- atom that isn't a simple thiol (-SH) is out of "
                    "scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a thiol must be attached to a carbon atom")
            thiols.add(atom.GetIdx())
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
    if not thiols:
        raise UnsupportedStructure("no thiol (-SH) group found; this module only handles thiols")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return thiols


def _reject_enethiol_carbons(graph, thiols, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for s_idx in thiols:
        (carbon,) = graph[s_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a thiol on a carbon that is also part of a C=C/C#C bond "
                "is out of scope for this module"
            )


def _name_from_substituents(chain_length, sh_locants, ene_locants, yne_locants, grouped):
    # Only chain_length == 1 omits a substituent prefix's own locant too
    # (see `_alcohol.py`'s equivalent comment) -- at chain_length == 2 a
    # substituent can still sit at either carbon, e.g. '2-aminoethanethiol'
    # (PubChem-verified) keeps its '2-' even though the thiol locant drops.
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return prefix + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(len(sh_locants), "thiol"), sh_locants
    )


def _candidate_key(chain_length, sh_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    sh_locant_set = lowest_locant_set(sh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, sh_locants, ene_locants, yne_locants, grouped)
    return (
        (
            sh_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _sh_locants(position_of, thiols, graph):
    locants = []
    for s in thiols:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants

def _ring_name_from_substituents(ring_size, sh_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size, ene_locants, yne_locants, prefix, total_subs, multiplied_word(len(sh_locants), "thiol"), sh_locants
    )


def _ring_candidate_key(ring_size, sh_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    sh_locant_set = lowest_locant_set(sh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, sh_locants, ene_locants, yne_locants, grouped)
    return sh_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, thiols, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    return ring_branch_stereo_display(graph, ring_order, thiols, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms)


def _name_cyclic_thiol(mol, thiols, stereo=None, bonds=(), ring_atoms=None, aromatic_atoms=frozenset()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_alcohol.py`'s `_name_cyclic_alcohol`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_alcohol.py`/`_ketone.py`'s
    own case): exactly one stereocenter on the ring's sole substituent
    branch instead embeds a bracketed descriptor into that substituent's
    own name, in place of the usual ring-locant prefix.

    `ring_atoms`/`aromatic_atoms`: when the thiol-bearing ring reaches
    this function as one half of a `two_separate_rings_with_plain_
    aromatic_substituent` shape (#628, mirroring `_ketone.py`'s #622/
    `_alcohol.py`'s #624), the caller passes the thiol-bearing ring's own
    atoms explicitly (RDKit's SSSR would otherwise list either of the two
    disjoint rings first) along with the other, aromatic ring's atoms,
    cited as a plain substituent (phenyl or a heteroaromatic monocycle)
    via `name_branch` the same way an ordinary alkyl ring substituent
    already is. `ring_atoms=None` (default) preserves the original
    single-ring dispatch unchanged."""
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
            graph, ring_order, thiols, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms
        )
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(
        substituents_for_ring(graph, ring_order, halogens, thiols, mol=mol, aromatic_atoms=aromatic_atoms).values()
    ):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "thiol is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, thiols, graph)
            substituents = substituents_for_ring(
                graph, candidate, halogens, thiols, mol=mol, aromatic_atoms=aromatic_atoms
            )
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, sh_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzenethiol_name_from_substituents(sh_locants, grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -SH carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # e.g. '2-methylbenzenethiol' (PubChem CID 8712), not
    # '2-methylbenzene-1-thiol' -- mirrors `_sulfonic_acid.py`'s
    # identical 'benzenesulfonic acid' treatment.
    thiol_word = multiplied_word(len(sh_locants), "thiol")
    if not grouped:
        return "benzene" + thiol_word
    return f"{format_substituent_prefixes(grouped)}benzene{thiol_word}"


def _benzenethiol_candidate_key(sh_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    sh_locant_set = lowest_locant_set(sh_locants)
    name = _benzenethiol_name_from_substituents(sh_locants, grouped)
    return sh_locant_set, locant_set, citation_locants, name


def _name_benzenethiol(mol, ring_atoms, exempt_atoms=None):
    """P-63.1.1: -SH attached directly to a benzene ring carbon -- e.g.
    'benzenethiol' (PubChem CID 7969), '2-methylbenzenethiol' (CID 8712).
    Mirrors `_name_cyclic_thiol`'s ring-numbering search, with the
    retained name 'benzene' as stem in place of 'cyclo' + alkane_name.
    Narrower than that saturated-ring case: only a single -SH directly on
    the ring is verified here (two or more direct ring thiols, a
    dithiophenol-style structure, remain out of scope pending a
    PubChem-confirmed example)."""
    thiols = _validate_and_collect_thiols(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
    if len(thiols) != 1:
        raise UnsupportedStructure(
            "more than one thiol directly on the benzene ring is not "
            "supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzenethiol is not "
            "supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, thiols, graph)
            substituents = substituents_for_ring(graph, candidate, halogens, thiols, mol=mol)
            key = _benzenethiol_candidate_key(sh_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_thiol(mol, ring_atoms):
    """Name one or more -SH groups lying entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropane-1-thiol,
    3-phenylpropane-1,2-dithiol (PubChem CID 154223538). The ring is
    cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_alcohol.py`'s `_name_phenyl_chain_alcohol`; multiple -SH
    locants are handled by the same `_sh_locants`/`_candidate_key`
    machinery already used by the acyclic-chain path above. Narrower
    than the acyclic path above: no chain unsaturation and no specified
    stereocenter -- each is a separate follow-up (see
    tasks/phenyl-substituent-on-thiol-chain.md's scope note)."""
    thiols = _validate_and_collect_thiols(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "thiol chain is not supported yet"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "thiol chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_chain_attachments_with_halogens(graph, rings, set(), halogens)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain thiol is not "
            "supported yet"
        )
    if any(chain_root in thiols for _, chain_root in attachment):
        raise UnsupportedStructure(
            "a thiol directly on the benzene ring (thiophenol-type) uses "
            "a separate construction, out of scope for this chain-parent "
            "module"
        )
    anchor_sulfur = next(iter(thiols))
    (anchor_carbon,) = graph[anchor_sulfur]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, thiols, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for s in thiols:
        (carbon,) = graph[s]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a thiol outside the single unbranched chain hanging off "
                "the benzene ring is not supported yet"
            )
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        sh_locants = _sh_locants(position_of, thiols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, sh_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_thiol(mol, thiols):
    """Name one or more -SH groups lying entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no thiol) -- e.g. cyclohexylmethanethiol.
    The ring is cited as a "cyclo..." substituent prefix (P-29.3.3) on the
    chain, which is the parent hydride, mirroring `_name_phenyl_chain_
    thiol` above and `_alcohol.py`'s `_name_ring_substituent_chain_
    alcohol`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, thiols)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    anchor_sulfur = next(iter(thiols))
    (anchor_carbon,) = graph[anchor_sulfur]
    chain, branches = longest_branched_chain_through(graph, anchor_carbon, ring_atoms, thiols, halogens=halogen_substituents(mol))
    chain_set = set(chain)
    for s in thiols:
        (carbon,) = graph[s]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a thiol outside the single branched chain hanging off "
                "the ring is not supported yet"
            )
    branches_by_atom = {
        chain[position - 1]: [r for r in roots if r != ring_atom]
        for position, roots in branches.items()
    }
    branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}

    ring_name = "cyclo" + alkyl_name(len(ring_atoms))
    chain_length = len(chain)

    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        sh_locants = _sh_locants(position_of, thiols, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, sh_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_thiol_chain_thiol(mol, thiols):
    """Name a thiol compound where the ring itself bears at least as many
    -SH's as a single unbranched chain hanging off exactly one ring atom
    does (P-44.1.1: the candidate with the greater count of the principal
    characteristic group thiol is senior; P-44.1.2.2 resolves an exact
    tie in the ring's favor) -- e.g. 2-(2-sulfanylethyl)cyclohexane-1-
    thiol, 1-(sulfanylmethyl)cyclohexane-1,2-dithiol. The ring is the
    parent; the chain is cited as a '(sulfanyl...alkyl)' substituent
    prefix, mirroring `_alcohol.py`'s `_name_ring_with_hydroxy_chain_
    alcohol` and `_amine.py`'s `_name_ring_with_amine_chain_amine`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, thiols)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, thiols)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_thiols = {s for s in thiols if next(iter(graph[s])) in chain_set}
    ring_thiols = thiols - chain_thiols
    if len(ring_thiols) < len(chain_thiols):
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{s: "sulfanyl" for s in ring_thiols}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, chain_thiols, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, sh_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{s: "sulfanyl" for s in chain_thiols}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, ring_thiols, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, sh_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_von_baeyer_or_spiro_thiol(mol, thiols, stereo, bonds):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single -SH suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` (P-16.3.3: 'thiol' begins with a consonant, so the
    parent stem's final 'e' is kept, e.g. 'bicyclo[2.2.1]heptane-2-thiol',
    not '...heptan-2-thiol' -- PubChem CID 13487780). Mirrors
    `_amine.py`'s/`_ketone.py`'s own bicyclic/polycyclic-before-spiro
    dispatch order and restrictions: exactly one thiol on the ring system
    itself; ring unsaturation composes on the bicyclic/polycyclic branch
    (matches `_ketone.py`'s bicyclic/polycyclic-vs-monospiro split) but not the monospiro
    branch, which has no base mechanism yet (P-31.1.5)."""
    if len(thiols) != 1:
        raise UnsupportedStructure(
            "more than one thiol on a von Baeyer bicyclic/polycyclic or "
            "monospiro ring system is not supported yet"
        )
    (sulfur,) = thiols
    graph = adjacency(mol)
    (thiol_carbon,) = graph[sulfur]

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
        return name_von_baeyer_suffix(
            mol,
            thiol_carbon,
            thiols,
            "thiol",
            "thiol",
            bicyclic_core,
            polycyclic_core,
            von_baeyer_ring_count,
            elide_e=False,
            stereo=stereo,
            bonds=bonds,
        )

    if bonds:
        raise UnsupportedStructure(
            "an unsaturated monospiro ring system is not supported yet "
            "(see P-31.1.5)"
        )

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro_suffix(
            mol, thiol_carbon, thiols, "thiol", "thiol", spiro_atom, elide_e=False, stereo=stereo
        )

    raise UnsupportedStructure(
        "polycyclic and fused-ring thiols are not supported yet (P-23/"
        "P-25 numbering integration with a suffix group is future work)"
    )


def name_thiol(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        anchors = list(_validate_and_collect_thiols(mol, aromatic_ring_atoms=union))
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzenethiol(mol, host, union)
        return _name_phenyl_chain_thiol(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            if is_benzene:
                ring_thiols = _validate_and_collect_thiols(mol, aromatic_ring_atoms=ring_atoms)
                if len(ring_thiols) == 1:
                    (only_s,) = ring_thiols
                    (only_s_carbon,) = adjacency(mol)[only_s]
                    if only_s_carbon in ring_atoms:
                        return _name_benzenethiol(mol, ring_atoms)
            return _name_phenyl_chain_thiol(mol, ring_atoms)

    # Two separate simple monocycles joined by one direct bond, one a plain
    # benzo/heteroaromatic ring with no substituent of its own (P-25 M2
    # step 1, #628) -- e.g. 2-phenylcyclohexane-1-thiol -- reuses the
    # existing single-ring `_name_cyclic_thiol` dispatch below with the
    # aromatic ring's atoms passed through as an exemption/substituent,
    # the same generalization #622/#624 made for `_ketone.py`/
    # `_alcohol.py`. Narrower than that: only the "every thiol is
    # ring-borne" split is handled below (matching this shape's only
    # verified real structures), mirroring #624's identical scoping
    # decision -- a thiol entirely on a chain hanging off the
    # non-aromatic ring, with the ring itself bearing none, is deferred as
    # a separate, more involved follow-up.
    aromatic_shape = None
    if ring_info.NumRings() == 2:
        aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
    aromatic_atoms = aromatic_shape[1] if aromatic_shape is not None else frozenset()

    thiols = _validate_and_collect_thiols(mol, aromatic_ring_atoms=aromatic_atoms)
    stereo = specified_stereocenters(mol)
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
    _reject_enethiol_carbons(graph, thiols, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 2 and aromatic_shape is not None:
        ring_atoms, _, _, _ = aromatic_shape
        chain_thiols = thiols - {s for s in thiols if next(iter(graph[s])) in ring_atoms}
        if chain_thiols:
            raise UnsupportedStructure(
                "a thiol on a chain hanging off the ring, with the ring "
                "itself bearing no thiol of its own, alongside this "
                "two-ring aromatic-substituent shape is not supported yet"
            )
        ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
        return _name_cyclic_thiol(
            mol, thiols, stereo, ring_bonds, ring_atoms=ring_atoms, aromatic_atoms=aromatic_atoms
        )
    if num_rings > 1:
        return _name_von_baeyer_or_spiro_thiol(mol, thiols, stereo, bonds)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
        if any(order == _YNE_ORDER for _, _, order in ring_bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a thiol is not "
                "supported yet -- only a ring double bond is in scope for "
                "this first pass (see P-31.1.3)"
            )
        ring_thiols = {s for s in thiols if next(iter(graph[s])) in ring_atoms}
        if not ring_thiols:
            return _name_acyclic_thiol(mol, thiols, bonds, stereo)
        if not bonds and ring_thiols and ring_thiols != thiols:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain thiol "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_thiol_chain_thiol(mol, thiols)
        if ring_thiols != thiols:
            raise UnsupportedStructure(
                "a thiol on a substituent branch chain rather than the "
                "ring itself is not supported yet"
            )
        return _name_cyclic_thiol(mol, thiols, stereo, ring_bonds)

    return _name_acyclic_thiol(mol, thiols, bonds, stereo)


def _name_acyclic_thiol(
    mol, thiols, bonds, stereo=None, extra_names=None, required_atoms=frozenset(), carbon_graph=None
):
    """`extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection. `None`
    keeps the original halogens-only behavior unchanged. `required_atoms`:
    additional carbon atoms (e.g. every demoted amine's own carbon
    neighbor) that a candidate chain must also carry -- empty by default
    so existing callers are unaffected. `carbon_graph`: the carbon-only
    graph to search for the principal chain -- defaults to
    `carbon_adjacency(mol)` (unchanged behavior); `_ether_thiol.py` passes
    one with a coexisting ether's alkoxy-branch component already removed,
    since an ether oxygen (not itself a carbon) would otherwise leave that
    branch as a separate component that could wrongly outrank the real
    thiol-bearing chain in `longest_chains`' global-diameter search."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    chains = all_chains(carbon_graph if carbon_graph is not None else carbon_adjacency(mol))
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _sh_locants(position_of, thiols, graph) is None:
            continue
        chain_set = set(chain)
        if not required_atoms <= chain_set:
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _sh_locants({a: i + 1 for i, a in enumerate(c)}, thiols, graph) is not None
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every thiol-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -SH groups, or an -SH expressed as a "
            "'sulfanyl' substituent prefix, is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, thiols, graph)
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, thiols, mol=mol)
            key, name = _candidate_key(chain_length, sh_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def has_thiol_shape(mol) -> bool:
    """Excludes a thiophene ring's own aromatic sulfur: otherwise it would
    satisfy this check ahead of `has_selenol_shape` in `core.py`'s
    chalcogen-chain dispatch (P-616 M4)."""
    return any(atom.GetAtomicNum() == 16 and not atom.GetIsAromatic() for atom in mol.GetAtoms())
