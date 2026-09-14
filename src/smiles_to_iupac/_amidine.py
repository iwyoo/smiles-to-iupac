"""Naming of primary amidines (the '-imidamide' suffix, terminal
-C(=NH)NH2) on acyclic saturated or unsaturated carbon chains, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-66.4.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'imidamide' is the preferred suffix for -C(=NH)NH2 (functional
  replacement of the amide suffix's =O by =NH). Confirmed via PubChem PUG
  REST: CID 67171 (`CC(=N)N`) -> "ethanimidamide", CID 4362807
  (`CCCCCC(=N)N`) -> "hexanimidamide", CID 68047 (`C(=N)N`) ->
  "methanimidamide" -- all matching the Blue Book's own worked examples
  directly ('ethanimidamide (PIN)', 'methanimidamide (PIN)',
  P6a.pdf 1874-1901).
- Like an amide carbon (`_amide.py`), an amidine carbon is always a chain
  terminus: after its imino nitrogen (=NH) and amino nitrogen (-NH2), it
  has room for at most one more substituent, which must be another chain
  carbon (or nothing, for methanimidamide, HC(=NH)NH2). So the parent
  chain's numbering is never a locant-minimization choice for the
  amidine group itself: whichever chain end carries the amidine carbon
  simply becomes C1, and its suffix locant is never cited (P-14.3.3),
  mirroring `_amide.py` exactly.
- P-31.0/P-31.1.1.1-.2: construction of the 'ene'/'yne' portion of a
  combined unsaturated-amidine name reuses the same mechanics as
  `_amide.py` (ending replaces 'ane' entirely, multiplying prefixes
  'di'/'tri' for >=2 bonds of a kind, 'ene' before 'yne', euphonic stem
  'a' before a multiplied ending); 'imidamide' is appended directly onto
  the last such segment. Confirmed via PubChem PUG REST: CID 4056230
  (`C=CCC(=N)N`) -> "but-3-enimidamide".
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the 'imidamide' suffix, reusing `halogen_substituents`/
  `format_substituent_prefixes` unchanged. Confirmed via PubChem PUG
  REST: CID 35602 (`ClCC(=N)N`) -> "2-chloroethanimidamide".
- P-92 stereocenters: unlike
  `_amide.py`'s carbonyl (always non-stereogenic), this module's own
  C=NH imine bond is *always* flagged by RDKit's
  `Chem.FindPotentialStereo` as an unspecified potential Bond_Double
  stereo element, regardless of substituents -- confirmed on
  `CC(=N)N`/`CCC(=N)N`/`C(=N)N` alike. That means any specified chain
  tetrahedral stereocenter always coexists with this unspecified imine
  bond, and `_common.specified_stereocenters` correctly rejects the
  combination as partially specified (P-92/P-93) rather than silently
  dropping either one -- same conclusion as `_sulfinic_acid.py`'s
  sulfur. This module only ever explicitly rejects a specified
  stereocenter rather than attempting to cite one.

- P-66.4.1.4.1 (Chapter P-6a): the amino (-NH2) nitrogen may carry 0-2 and
  the imino (=NH) nitrogen 0-1 plain, unbranched, unsubstituted, saturated
  alkyl substituents (mirroring `_hydrazide.py`'s identical N-/N'-alkyl
  support) or a single plain, unsubstituted benzene ring (mirroring
  `_amide.py`'s/`_hydrazide.py`'s own N-phenyl pattern), cited with 'N'/
  'N'' locants -- confirmed via PubChem PUG REST and the primary text's
  own worked examples: 'N-phenylbenzenecarboximidamide (PIN)',
  'N'-methyl-N,N-diphenylbenzenecarboximidamide (PIN)',
  'N'-ethyl-N-methylbenzenecarboximidamide (PIN)' (P-66.4.1.4.1's own
  confirmed worked examples: "The locant N refers to the amino group and
  N' refers to the imino group"). Identically-named substituents (on the
  same or different nitrogen) share one multiplying prefix, e.g.
  'N,N-diphenyl...'.

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_amide.py`'s own first-pass scope:
- A branched, unsaturated, or otherwise-substituted N-/N'-alkyl
  substituent, or an alkyl-and-phenyl combination sharing the same
  nitrogen alongside a chain benzene-ring substituent elsewhere - each a
  separate follow-up.
- An amidine on/in a ring (e.g. 'cyclohexanecarboximidamide') - a
  separate module's territory.
- More than one amidine group in the same molecule.
- Any oxygen anywhere (functional-group coexistence with amidine is
  entirely unverified in this first pass).
- A carbon with other than exactly one imino nitrogen and one amino
  nitrogen, or more than one carbon neighbor (a ketone-shaped carbon is
  not a valid amidine carbon).
- An aromatic carbon or ring anywhere in the molecule.
- Any other heteroatom (O, S, ...).
"""

from rdkit import Chem

from ._common import (
    elides_before,
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
    linear_branch,
    longest_branched_chain,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix
from ._substituents import (
    alpha_sort_key,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
)

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def has_amidine_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded nitrogen and a
    singly-bonded nitrogen at the same time -- an amidine-family
    (-C(=N-)N<) pattern -- regardless of each nitrogen's own substitution
    or whether the rest of the molecule is in scope. Used by `core.py` to
    route ahead of the imine dispatch, since an amidine carbon would
    otherwise look imine-shaped to that module (both have a C=N double
    bond), and that module's own validation doesn't reject a second,
    singly-bonded nitrogen the way it should for this shape -- unlike
    `_amide.py`'s aldehyde fallback, which does reject a stray nitrogen
    outright. Deliberately not restricted to the primary -C(=NH)NH2
    shape this module actually names: this project's `_validate_and_
    collect_amidine` below is what raises a precise `UnsupportedStructure`
    for an N-substituted amidine, geminal diamidine (e.g. guanidine), or
    amidoxime instead of letting any of them fall through to a silent
    imine-module misname."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        has_imino_like = any(
            n.GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
        has_amino_like = any(
            n.GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in atom.GetNeighbors()
        )
        if has_imino_like and has_amino_like:
            return True
    return False


def _amidine_cores(mol):
    """Every -C(=N-)(N<) amidine-shaped carbon in `mol`, regardless of
    each nitrogen's own substitution or whether the rest of the molecule
    is in scope: (carbon, imino_nitrogen, amino_nitrogen, imino_subs,
    amino_subs). `imino_subs`/`amino_subs` are that nitrogen's own
    substituent atom indices (besides the amidine carbon itself) -- 0-1
    for the imino =N- (room for only one more bond besides the double
    bond to the amidine carbon), 0-2 for the amino -N< (room for up to
    two, P-66.4.1.4.1's confirmed 'N,N-diphenyl...' worked example)."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        imino_candidates = []
        amino_candidates = []
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() != 7:
                continue
            order = mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
            if order == 2.0:
                imino_candidates.append(n)
            elif order == 1.0:
                amino_candidates.append(n)
        if len(imino_candidates) != 1 or len(amino_candidates) != 1:
            continue
        (imino_n,) = imino_candidates
        (amino_n,) = amino_candidates
        imino_subs = tuple(x.GetIdx() for x in imino_n.GetNeighbors() if x.GetIdx() != atom.GetIdx())
        amino_subs = tuple(x.GetIdx() for x in amino_n.GetNeighbors() if x.GetIdx() != atom.GetIdx())
        cores.append((atom.GetIdx(), imino_n.GetIdx(), amino_n.GetIdx(), imino_subs, amino_subs))
    return cores


def _validate_and_collect_amidine(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amidine_carbon, imino_nitrogen, amino_nitrogen,
    imino_subs, amino_subs): the single -C(=N-)(N<) carbon/imino-
    nitrogen/amino-nitrogen atom indices, and each nitrogen's own 0-1/0-2
    substituent atom indices (P-66.4.1.1: 'N' is the amino locant, 'N''
    the imino one).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_amidine`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_amidine`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors `_amide.py`'s
    `_validate_and_collect_amide`."""
    cores = _amidine_cores(mol)
    if len(cores) != 1:
        raise UnsupportedStructure(
            "no primary amidine (-C(=N-)(N<)) group found, or more than "
            "one, which this module does not attempt to disambiguate"
        )
    (amidine_carbon, imino_nitrogen, amino_nitrogen, imino_subs, amino_subs) = cores[0]
    if len(imino_subs) > 1:
        raise UnsupportedStructure(
            "the amidine's imino (=N-) nitrogen has room for at most one "
            "substituent"
        )
    if len(amino_subs) > 2:
        raise UnsupportedStructure(
            "the amidine's amino (-N<) nitrogen has room for at most two "
            "substituents"
        )
    for n_idx, subs in ((imino_nitrogen, imino_subs), (amino_nitrogen, amino_subs)):
        for sub_idx in subs:
            bond = mol.GetBondBetweenAtoms(n_idx, sub_idx)
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "an amidine nitrogen substituent must be singly bonded"
                )
            if mol.GetAtomWithIdx(sub_idx).GetAtomicNum() != 6:
                raise UnsupportedStructure(
                    "an amidine nitrogen substituent other than a "
                    "carbon-based group is out of scope for this module"
                )

    known_nitrogen_atoms = {imino_nitrogen, amino_nitrogen}
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a primary amidine's two nitrogens "
                "(P-66.4.1.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 7:
            if atom.GetIdx() not in known_nitrogen_atoms:
                raise UnsupportedStructure(
                    "a nitrogen other than the amidine's own imino/amino "
                    "pair needs Table 3.3 seniority handling not yet "
                    "implemented here"
                )
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

    amidine_carbon_atom = mol.GetAtomWithIdx(amidine_carbon)
    carbon_neighbors = [n for n in amidine_carbon_atom.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "an amidine carbon with more than one carbon neighbor is "
            "not a valid terminal amidine carbon (a ketone-shaped "
            "carbon is out of scope for this module)"
        )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amidine_carbon, imino_nitrogen, amino_nitrogen, imino_subs, amino_subs


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'imidamide'
    ending (e.g. '3-enimidamide'); the amidine group's own locant is never
    cited (P-14.3.3, see module docstring)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["imidamide"]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    if segments:
        locant_parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words[:-1])
        ]
        body = "-".join(locant_parts) + words[-1]
    else:
        body = words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The amidine
    group's own locant isn't part of this key: candidates are
    pre-filtered so the amidine carbon always sits at C1 (see
    `_name_acyclic_amidine`)."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
    return (
        (
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _collect_n_alkyl(full_carbon_graph, mol, amino_subs, imino_subs):
    """Validate amino's (0-2) and imino's (0-1) alkyl substituents -- each
    must be a plain, unbranched, unsubstituted, saturated alkyl chain,
    mirroring `_hydrazide.py`'s identical `_collect_n_alkyl` -- and return
    (entries, n_alkyl_atoms): `entries` is a list of ("N"/"N'", name)
    pairs (P-66.4.1.1's N/N' convention -- 'N' is the amino locant, 'N''
    the imino one) ready for `_format_n_prefix`, and `n_alkyl_atoms` is
    the full set of atom indices spanned by every N-substituent, to
    exclude from the principal-chain search below."""
    entries = []
    n_alkyl_atoms = set()
    for locant_label, roots in (("N", amino_subs), ("N'", imino_subs)):
        for root in roots:
            length = linear_branch(full_carbon_graph, root, None)
            if length is None:
                raise UnsupportedStructure("a branched N-substituent is not supported yet")
            atoms = set()
            previous, current = None, root
            while current is not None:
                atoms.add(current)
                neighbors = [n for n in full_carbon_graph[current] if n != previous]
                previous, current = current, (neighbors[0] if neighbors else None)
            if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
                raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
            # `full_carbon_graph` only sees carbon-carbon bonds, so a
            # halogen (or any other heteroatom) hanging off the chain is
            # invisible to the walk above and would otherwise pass through
            # silently, misnaming e.g. a -CH2CF3 N-substituent as plain
            # 'ethyl' with the three fluorines dropped entirely (found via
            # real-data testing: 'CCC(=N)N(CC)CC(F)(F)F' was misnamed
            # 'N,N-diethylpropanimidamide'). The only legitimate neighbor
            # of a chain atom outside `atoms` is the amidine nitrogen
            # itself, at `root`.
            for atom_idx in atoms:
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                    if neighbor.GetAtomicNum() == 1 or neighbor.GetIdx() in atoms:
                        continue
                    if neighbor.GetAtomicNum() == 7:
                        continue
                    raise UnsupportedStructure(
                        "a substituted N-substituent (e.g. bearing a "
                        "halogen) is not supported yet; only a plain, "
                        "unsubstituted alkyl N-substituent is in scope"
                    )
            entries.append((locant_label, alkyl_name(length)))
            n_alkyl_atoms |= atoms
    return entries, n_alkyl_atoms


def _format_n_prefix(entries):
    """entries: [("N"/"N'", name), ...] -> the assembled "N-"/"N'-" prefix
    string, grouping identically-named substituents (whether on the same
    nitrogen or split across both) under one shared multiplying prefix,
    e.g. [("N", "phenyl"), ("N", "phenyl")] -> "N,N-diphenyl" (P-66.4.1.4.1's
    confirmed 'N,N-diphenylbenzenecarboximidamide' worked example). '' if
    entries is empty. Mirrors `_hydrazide.py`'s identical helper."""
    if not entries:
        return ""
    grouped = {}
    for locant, name in entries:
        grouped.setdefault(name, []).append(locant)
    parts = []
    for name in sorted(grouped, key=alpha_sort_key):
        locants = sorted(grouped[name])
        multiplier = multiplying_prefix(len(locants)) if len(locants) > 1 else ""
        parts.append(f"{','.join(locants)}-{multiplier}{name}")
    return "-".join(parts)


def _name_acyclic_amidine(
    mol,
    amidine_carbon,
    excluded,
    bonds,
    amino_subs=(),
    imino_subs=(),
    extra_n_entries=(),
    extra_excluded_carbons=frozenset(),
):
    """`extra_n_entries`/`extra_excluded_carbons`: an already-resolved
    ("N"/"N'", name) pair (e.g. a phenyl N-substituent, see
    `_name_amidine_with_n_phenyl`) and the atom indices it spans, merged
    in alongside the ordinary alkyl entries below -- both empty by
    default, so every other caller's behavior is unchanged."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    full_carbon_graph = carbon_adjacency(mol)
    n_entries, n_alkyl_atoms = _collect_n_alkyl(full_carbon_graph, mol, amino_subs, imino_subs)
    n_entries = n_entries + list(extra_n_entries)
    n_prefix = _format_n_prefix(n_entries)
    n_alkyl_atoms = n_alkyl_atoms | set(extra_excluded_carbons)

    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_alkyl_atoms}
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if amidine_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the amidine-bearing carbon (and/or a multiple bond) does not "
            "lie on a single longest carbon chain; a shorter principal "
            "chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != amidine_carbon:
                # The amidine carbon must sit at C1 (see module docstring);
                # a direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name

    if n_prefix:
        separator = "-" if best_name[0].isdigit() else ""
        best_name = f"{n_prefix}{separator}{best_name}"
    return best_name


def _name_phenyl_chain_amidine(mol, ring_atoms):
    """Name a primary amidine whose -C(=NH)NH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropanimidamide. The ring
    is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_amide.py`'s `_name_phenyl_chain_amide`. Narrower than the
    acyclic path above: no chain unsaturation and no specified
    stereocenter."""
    amidine_carbon, imino_nitrogen, amino_nitrogen, imino_subs, amino_subs = _validate_and_collect_amidine(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if imino_subs or amino_subs:
        raise UnsupportedStructure(
            "an N-/N'-alkyl-substituted amidine alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol) is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside this amidine's own "
            "always-unspecified C=NH imine bond is not supported yet "
            "(see P-92/P-93, module docstring)"
        )

    excluded = {imino_nitrogen, amino_nitrogen}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "amidine chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain amidine is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain(graph, amidine_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an amidine directly attached to the benzene ring (the "
            "'benzamidine'-style naming) uses a separate construction, out "
            "of scope for this acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def _name_amidine_with_n_phenyl(mol, ring_atoms):
    """Name an amidine whose amino (N) or imino (N') nitrogen carries a
    plain, unsubstituted benzene ring as a direct substituent -- e.g.
    'N'-phenylethanimidamide' for `CC(=Nc1ccccc1)N` (ring on the
    double-bonded/imino nitrogen; PubChem PUG REST IUPACName match, and
    its own stored structure round-trips to the same tautomer, so this
    verification is trustworthy). The other tautomer, `CC(=N)Nc1ccccc1`
    (ring on the singly-bonded/amino nitrogen instead -- a genuinely
    different explicit structure, not merely a resonance form), is
    verified directly from the primary text instead: P-66.4.1.4.1's own
    confirmed worked example 'N-phenylbenzenecarboximidamide (PIN)' is
    exactly this shape. PubChem cannot cross-check that specific tautomer
    independently here -- querying it for `CC(=N)Nc1ccccc1` silently
    returns the *other* tautomer's stored structure/name instead
    (`CC(=NC1=CC=CC=C1)N`, i.e. it renormalizes the amidine tautomer
    before naming), so its answer for that query is not evidence about
    the structure actually asked about. A mixed case with two different
    substituents doesn't have this ambiguity and round-trips cleanly,
    e.g. `CC(=NC)Nc1ccccc1` -> 'N'-methyl-N-phenylethanimidamide'
    (PubChem PUG REST IUPACName match, structure round-trip confirmed).
    Mirrors `_hydrazide.py`'s `_name_hydrazide_with_n_phenyl`: the ring is
    just another N-substituent name, merged into
    `_name_acyclic_amidine`'s existing N-/N'-prefix assembly. The amino
    nitrogen may carry a second, different substituent alongside the
    phenyl (P-66.4.1.4.1's confirmed 'N,N-diphenyl...' pattern allows up
    to two on N); the imino nitrogen has room for only the one."""
    amidine_carbon, imino_nitrogen, amino_nitrogen, imino_subs, amino_subs = _validate_and_collect_amidine(
        mol, aromatic_ring_atoms=ring_atoms
    )
    stereo_check = specified_stereocenters(mol)
    if stereo_check is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside this amidine's own "
            "always-unspecified C=NH imine bond is not supported yet "
            "(see P-92/P-93, module docstring)"
        )
    if set(imino_subs) & ring_atoms:
        if len(imino_subs) != 1:
            raise UnsupportedStructure(
                "the amidine's imino (=N-) nitrogen has room for at most "
                "one substituent"
            )
        extra_n_entries = (("N'", "phenyl"),)
        imino_subs, amino_subs = (), amino_subs
    else:
        if not (set(amino_subs) & ring_atoms):
            raise UnsupportedStructure(
                "a benzene ring elsewhere in the molecule is out of scope "
                "for this module"
            )
        extra_n_entries = (("N", "phenyl"),)
        remaining = tuple(s for s in amino_subs if s not in ring_atoms)
        imino_subs, amino_subs = imino_subs, remaining

    excluded = {imino_nitrogen, amino_nitrogen}
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and (b[0] not in ring_atoms or b[1] not in ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    return _name_acyclic_amidine(
        mol,
        amidine_carbon,
        excluded,
        bonds,
        amino_subs=amino_subs,
        imino_subs=imino_subs,
        extra_n_entries=extra_n_entries,
        extra_excluded_carbons=ring_atoms,
    )


def name_amidine(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            graph = adjacency(mol)
            attachment = ring_chain_attachment(graph, ring_atoms, set())
            if attachment is not None and mol.GetAtomWithIdx(attachment[1]).GetAtomicNum() == 7:
                return _name_amidine_with_n_phenyl(mol, ring_atoms)
            return _name_phenyl_chain_amidine(mol, ring_atoms)
    amidine_carbon, imino_nitrogen, amino_nitrogen, imino_subs, amino_subs = _validate_and_collect_amidine(mol)
    if specified_stereocenters(mol) is not None:
        # The amidine's own C=NH imine bond is always flagged by RDKit's
        # `Chem.FindPotentialStereo` as an unspecified potential
        # Bond_Double stereo element, regardless of substituents (module
        # docstring) -- so any specified chain stereocenter here always
        # coexists with that unspecified imine bond, and
        # `specified_stereocenters` correctly rejects the combination
        # (P-92/P-93) instead of the silent drop this project's
        # stereodescriptor safety net exists to fix.
        raise UnsupportedStructure(
            "a specified stereocenter alongside this amidine's own "
            "always-unspecified C=NH imine bond is not supported yet "
            "(see P-92/P-93, module docstring)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amidine on/in a ring is out of scope for this acyclic-only "
            "module"
        )

    excluded = {imino_nitrogen, amino_nitrogen}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    return _name_acyclic_amidine(mol, amidine_carbon, excluded, bonds, amino_subs=amino_subs, imino_subs=imino_subs)
