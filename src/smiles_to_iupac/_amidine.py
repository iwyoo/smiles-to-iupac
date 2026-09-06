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

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_amide.py`'s own first-pass scope:
- An amidine nitrogen (either the imino =NH or the amino -NH2) with any
  substituent other than its hydrogen(s) (N-substituted amidine) -
  deferred entirely; only a primary amidine is supported here.
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
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ordered_chain,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    specified_stereocenters,
)
from ._numerals import alkane_name
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


def _validate_and_collect_amidine(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amidine_carbon, imino_nitrogen, amino_nitrogen): the
    single -C(=NH)NH2 carbon/imino-nitrogen/amino-nitrogen atom indices.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_amidine`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_amidine`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors `_amide.py`'s
    `_validate_and_collect_amide`."""
    has_carbon = False
    amidine_carbons = set()
    imino_by_carbon = {}
    amino_by_carbon = {}
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
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "an amidine nitrogen must be attached to exactly one "
                    "carbon atom"
                )
            (carbon,) = carbon_neighbors
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx())
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 2.0 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
                imino_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
                continue
            if bond_order == 1.0 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2:
                amino_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
                continue
            raise UnsupportedStructure(
                "an amidine nitrogen with any substituent other than its "
                "hydrogen(s) (N-substituted amidine) is out of scope for "
                "this module"
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

    for carbon_idx in set(imino_by_carbon) | set(amino_by_carbon):
        iminos = imino_by_carbon.get(carbon_idx, [])
        aminos = amino_by_carbon.get(carbon_idx, [])
        if len(iminos) != 1 or len(aminos) != 1:
            raise UnsupportedStructure(
                "a nitrogen pattern that isn't exactly one imino (=NH) and "
                "one amino (-NH2) nitrogen on the same carbon is not a "
                "plain primary amidine, which this module does not attempt "
                "to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "an amidine carbon with more than one carbon neighbor is "
                "not a valid terminal amidine carbon (a ketone-shaped "
                "carbon is out of scope for this module)"
            )
        amidine_carbons.add(carbon_idx)

    if not amidine_carbons:
        raise UnsupportedStructure(
            "no primary amidine (-C(=NH)NH2) group found; this module only "
            "handles primary amidines"
        )
    if len(amidine_carbons) > 1:
        raise UnsupportedStructure(
            "more than one amidine group is out of scope for this module"
        )
    (amidine_carbon,) = amidine_carbons
    (imino_nitrogen,) = imino_by_carbon[amidine_carbon]
    (amino_nitrogen,) = amino_by_carbon[amidine_carbon]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amidine_carbon, imino_nitrogen, amino_nitrogen


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
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
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
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_amidine(mol, amidine_carbon, excluded, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
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
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
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
    amidine_carbon, imino_nitrogen, amino_nitrogen = _validate_and_collect_amidine(
        mol, aromatic_ring_atoms=ring_atoms
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
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, excluded)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside an "
            "amidine is not supported yet"
        )
    if chain[-1] != amidine_carbon:
        raise UnsupportedStructure(
            "the amidine carbon must be the chain's far terminus from the "
            "benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an amidine directly attached to the benzene ring (the "
            "'benzamidine'-style naming) uses a separate construction, out "
            "of scope for this acyclic-chain-parent module"
        )

    ordered = list(reversed(chain))
    chain_length = len(ordered)
    position_of = {atom: i + 1 for i, atom in enumerate(ordered)}
    substituents = {
        position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def name_amidine(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_amidine(mol, ring_atoms)
    amidine_carbon, imino_nitrogen, amino_nitrogen = _validate_and_collect_amidine(mol)
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

    return _name_acyclic_amidine(mol, amidine_carbon, excluded, bonds)
