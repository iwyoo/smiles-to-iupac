"""Naming of sulfonic acids (the '-sulfonic acid' suffix, -SO3H) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.3.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'sulfonic acid' is the preselected suffix for -S(=O)(=O)-OH. Table 3.3's
  overall seniority order clusters every acid class (carboxylic, sulfonic,
  sulfinic, ...) ahead of esters/amides/nitriles/aldehydes/ketones/alcohols/
  thiols/amines -- sulfonic acid is senior to all of those, but *junior* to
  carboxylic acid within the acid cluster itself (verified: carboxylic acid
  is listed first among the 'Suffix' acid classes). This module doesn't
  implement that acid-vs-acid seniority yet (see scope note below).
- Unlike 'ol'/'thiol'/'amine' (a single suffix word replacing/following the
  parent hydride's final 'e'), 'sulfonic acid' is cited as the parent
  hydride name followed directly by the two-word suffix with no elision --
  'methane' + 'sulfonic acid' -> 'methanesulfonic acid' -- since 'sulfonic'
  begins with a consonant, matching `_thiol.py`'s reasoning for 'thiol'.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_thiol.py`/`_alcohol.py` apply (mononuclear parent, or a saturated
  two-carbon chain, regardless of other substituents), e.g.
  'ethanesulfonic acid'.
- P-44.4.1.8 / P-45.2: the -SO3H locant is minimized before ene/yne locants,
  which are minimized before substituent-prefix locants -- same ordering
  as every other suffix module here.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -SO3H suffix.
- P-91.3/P-92: a molecule
  with one or more *specified* tetrahedral stereocenters -- every one on
  the principal chain/ring itself, no unspecified one alongside them, and
  no C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-"
  prefix, ascending locant order, e.g. '(2R)-butane-2-sulfonic acid',
  '(1S,2S)-2-chlorocyclohexane-1-sulfonic acid', same pattern as
  `_ketone.py` (chain and ring both).

Scope, deliberately narrow (first pass at this functional group, mirroring
how `_thiol.py` started): a single -SO3H on an acyclic chain or on a
single saturated carbon ring (monocyclic, mirroring `_thiol.py`'s own
monocyclic support), with no other heteroatom anywhere in the molecule
except the sulfonic acid group's own three oxygens -- acid-vs-acid
seniority (vs. a coexisting carboxylic acid) and any other Table 3.3
seniority coexistence is future work, tracked under
`multi-carbonyl-seniority.md`.

P-31.1.3: a monocyclic ring bearing a sulfonic acid and exactly one C=C
ring double bond -- e.g. 'cyclohex-2-ene-1-sulfonic acid',
'cyclohex-3-ene-1-sulfonic acid', both confirmed via PubChem PUG REST.
The sulfonic acid's own locant is never omittable here even as the
ring's sole substituent, mirroring `_ketone.py`/`_thiol.py`'s identical
extension. Deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): a specified
stereocenter on a von Baeyer polycyclic or spiro skeleton (a single
sulfonic acid on such a skeleton is supported, P-23/P-24 numbering
integration via `_polycyclic_suffix.py` with `elide_e=False`, mirroring
`_thiol.py`); unsaturation reaching outside the ring or a ring triple
bond, an -SO3H on a substituent branch off an otherwise-unsubstituted
*saturated* ring, two or more -SO3H groups, and a sulfonic acid on a
carbon that is also part of a C=C/C#C bond. Two *aromatic*-ring cases:
`_name_phenyl_chain_sulfonic_acid` names a -SO3H chain hanging off a
single plain, unsubstituted benzene ring (e.g.
'3-phenylpropane-1-sulfonic acid'), mirroring `_ketone.py`/
`_aldehyde.py`'s identical benzene-ring-substituent path -- narrower than
the acyclic path: no chain unsaturation and no specified stereocenter.
`_name_benzenesulfonic_acid` names -SO3H directly on a benzene ring
carbon (with or without other ring substituents), e.g.
'benzenesulfonic acid' (PubChem CID 7371), '2-methylbenzenesulfonic
acid' (CID 6925), mirroring `_name_cyclic_sulfonic_acid`'s ring-
numbering search with the retained name 'benzene' as stem -- the
-SO3H's own locant is never cited here, unlike the cycloalkane case.

`_name_phenyl_chain_sulfonic_acid`'s chain-substituent shape also
extends to a simple heteroaromatic monocycle (pyridine/furan/thiophene/
pyrrole, P-29.3.4.1), e.g. 'pyridin-3-ylmethanesulfonic acid'. A
heteroaromatic ring with -SO3H bonded directly to it (unlike benzene)
stays out of scope: benzene's own numbering is free to start at the
-SO3H carbon regardless, but a heteroaromatic ring's numbering must fix
the heteroatom at locant 1 and search for the -SO3H's own lowest locant
relative to it -- ring-locant-search machinery not built here yet.
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
    name_from_substituents,
    non_single_bonds,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_branch_attachments,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    substituents_for_ring,
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    ring_branch_stereo_display,
    substituents_for_chain,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfonic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonic acid group: bonded to exactly one
    carbon, two double-bonded (terminal) oxygens, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        hydroxyl_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_sulfonic_acid_shape(mol) -> bool:
    return bool(_sulfonic_sulfur_atoms(mol))


def _validate_and_collect_sulfonic_acids(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene or heteroaromatic-monocycle ring with exactly one exocyclic
    attachment -- exempted from the aromatic-atom rejection below so
    `name_sulfonic_acid`'s ring-substituent path (see
    `_name_phenyl_chain_sulfonic_acid`) can reuse this same validation for
    the rest of the molecule. Empty by default, so every other caller's
    behavior is unchanged."""
    sulfur_atoms = _sulfonic_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfonic acid (-SO3H) group found; this module only "
            "handles sulfonic acids"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfonic acid group is out of scope for this "
            "module"
        )
    sulfonic_atom_idxs = set()
    for s in sulfur_atoms:
        sulfonic_atom_idxs.add(s.GetIdx())
        sulfonic_atom_idxs.update(n.GetIdx() for n in s.GetNeighbors() if n.GetAtomicNum() == 8)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in sulfonic_atom_idxs and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfonic acid group (P-65.3.1), "
                "a heteroaromatic ring's own heteroatom (P-29.3.4.1), and "
                "halogen substituents (P-35.2.1) are not supported "
                "yet -- in particular a coexisting carboxylic acid or "
                "other characteristic group needs acid-vs-acid Table 3.3 "
                "seniority handling not yet implemented here"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (sulfur,) = sulfur_atoms
    (carbon,) = (n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    return sulfur.GetIdx(), carbon.GetIdx()


def _reject_enesulfonic_carbon(graph, so3h_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so3h_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfonic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
        )


def _name_from_substituents(chain_length, so3h_locant, ene_locants, yne_locants, grouped):
    # Only chain_length == 1 omits a substituent prefix's own locant too
    # (see `_alcohol.py`'s equivalent comment).
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    return prefix + name_from_substituents(chain_length, ene_locants, yne_locants, "sulfonic acid", [so3h_locant], substituted=bool(grouped))


def _candidate_key(chain_length, so3h_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so3h_locant, ene_locants, yne_locants, grouped)
    return (
        (
            so3h_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )

def _ring_name_from_substituents(ring_size, so3h_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size, ene_locants, yne_locants, prefix, total_subs, "sulfonic acid", [so3h_locant]
    )


def _ring_candidate_key(ring_size, so3h_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, so3h_locant, ene_locants, yne_locants, grouped)
    return so3h_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=None):
    return ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)


def _name_cyclic_sulfonic_acid(mol, sulfur_idx, so3h_carbon, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_ketone.py`'s `_name_cyclic_ketone`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_alcohol.py`/`_ketone.py`/
    `_thiol.py`'s own case): exactly one stereocenter on the ring's sole
    substituent branch instead embeds a bracketed descriptor into that
    substituent's own name, in place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(substituents_for_ring(graph, ring_order, halogens, excluded, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "sulfonic acid is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so3h_locant = position_of[so3h_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, so3h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _benzenesulfonic_acid_name_from_substituents(grouped):
    # Unlike the cycloalkane case, the mancude ring's own numbering is
    # always free to start at the -SO3H carbon (P-14.3.3-style), so its
    # locant is never cited even when other substituents need theirs,
    # e.g. '2-methylbenzenesulfonic acid' (PubChem CID 6925), not
    # '2-methylbenzene-1-sulfonic acid' -- mirrors `_carboxylic_acid.py`'s
    # identical 'benzoic acid' treatment.
    if not grouped:
        return "benzenesulfonic acid"
    return f"{format_substituent_prefixes(grouped)}benzenesulfonic acid"


def _benzenesulfonic_acid_candidate_key(so3h_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _benzenesulfonic_acid_name_from_substituents(grouped)
    return so3h_locant, locant_set, citation_locants, name


def _name_benzenesulfonic_acid(mol, ring_atoms, exempt_atoms=None):
    """P-65.3.1: -SO3H attached directly to a benzene ring carbon -- e.g.
    'benzenesulfonic acid' (PubChem CID 7371), '2-methylbenzenesulfonic
    acid' (CID 6925), '4-methylbenzenesulfonic acid' (CID 6101). Mirrors
    `_name_cyclic_sulfonic_acid`'s ring-numbering search exactly, with the
    aromatic retained name 'benzene' as stem in place of 'cyclo' +
    alkane_name; an aromatic ring has no ene/yne ring-bond locants of its
    own, so those are always empty here."""
    sulfur_idx, so3h_carbon = _validate_and_collect_sulfonic_acids(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside benzenesulfonic acid is "
            "not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so3h_locant = position_of[so3h_carbon]
            substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
            key = _benzenesulfonic_acid_candidate_key(so3h_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _name_phenyl_chain_sulfonic_acid(mol, ring_atoms):
    """Name a sulfonic acid whose -SO3H lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene or heteroaromatic-monocycle ring -- e.g.
    3-phenylpropane-1-sulfonic acid, pyridin-3-ylmethanesulfonic acid.
    The ring is cited as a 'phenyl'/heteroaromatic substituent prefix (via
    `name_branch`'s aromatic-ring recognition) on the chain, which is the
    parent hydride, mirroring `_ketone.py`'s `_name_phenyl_chain_ketone`.
    Narrower than the acyclic path above: no chain unsaturation and no
    specified stereocenter -- each is a separate follow-up."""
    sulfur_idx, so3h_carbon = _validate_and_collect_sulfonic_acids(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "sulfonic acid chain is not supported yet"
        )
    excluded = {sulfur_idx}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if sulfur_idx not in (b[0], b[1]) and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "sulfonic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain sulfonic acid is "
            "not supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, so3h_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so3h_locant = position_of[so3h_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, so3h_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_sulfonic_acid(mol, sulfur_idx, so3h_carbon):
    """Name a sulfonic acid whose -SO3H lies entirely on a single branched
    chain hanging off one atom of an otherwise-plain saturated monocyclic
    ring (the ring itself bears no sulfonic acid) -- e.g.
    cyclohexylmethanesulfonic acid. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_sulfonic_acid` above and
    `_ketone.py`'s `_name_ring_substituent_chain_ketone`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])
    excluded = {sulfur_idx}

    attachment = ring_chain_attachment(graph, ring_atoms, excluded)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment

    chain, branches = longest_branched_chain_through(graph, so3h_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
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
        so3h_locant = position_of[so3h_carbon]
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, so3h_locant, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_von_baeyer_or_spiro_sulfonic_acid(mol, sulfur_idx, so3h_carbon, bonds, stereo):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single sulfonic acid (-SO3H) suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` ('sulfonic acid' begins with a consonant, P-16.3.3,
    same as `_thiol.py`'s 'thiol' -- e.g. 'bicyclo[2.2.1]heptane-2-
    sulfonic acid', PubChem CID 14579878). Mirrors `_sulfinic_acid.py`'s
    own bicyclic/polycyclic-before-spiro dispatch order and restrictions:
    exactly one sulfonic acid on the ring system itself; ring
    unsaturation composes on the bicyclic/polycyclic branch (mirrors
    `_ketone.py`'s identical extension) but not the monospiro branch,
    which has no base mechanism yet (P-31.1.5)."""
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
            so3h_carbon,
            {sulfur_idx},
            "sulfonic acid",
            "sulfonic acid",
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
            mol, so3h_carbon, {sulfur_idx}, "sulfonic acid", "sulfonic acid", spiro_atom, elide_e=False, stereo=stereo
        )

    raise UnsupportedStructure(
        "polycyclic and fused-ring sulfonic acids are not supported yet "
        "(P-23/P-25 numbering integration with a suffix group is future "
        "work)"
    )


def name_sulfonic_acid(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        union = set().union(*aromatic_rings)
        sulfur_idx, _ = _validate_and_collect_sulfonic_acids(mol, aromatic_ring_atoms=union)
        anchors = [sulfur_idx]
        host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
        if host is not None:
            return _name_benzenesulfonic_acid(mol, host, union)
        return _name_phenyl_chain_sulfonic_acid(mol, union)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            sulfur_idx, so3h_carbon = _validate_and_collect_sulfonic_acids(mol, aromatic_ring_atoms=ring_atoms)
            if so3h_carbon in ring_atoms:
                if is_heteroaromatic:
                    # A heteroaromatic ring's numbering must fix the
                    # heteroatom at locant 1 and search for the -SO3H's own
                    # lowest locant relative to it -- unlike benzene, where
                    # any ring-numbering start is equally valid, so no
                    # locant is ever cited.
                    raise UnsupportedStructure(
                        "a sulfonic acid directly attached to a "
                        "heteroaromatic ring is not supported yet"
                    )
                return _name_benzenesulfonic_acid(mol, ring_atoms)
            return _name_phenyl_chain_sulfonic_acid(mol, ring_atoms)
    sulfur_idx, so3h_carbon = _validate_and_collect_sulfonic_acids(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and sulfur_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 2:
        # The two S=O double bonds are always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enesulfonic_carbon(graph, so3h_carbon, bonds)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 1:
        return _name_von_baeyer_or_spiro_sulfonic_acid(mol, sulfur_idx, so3h_carbon, bonds, stereo)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
        if any(order == _YNE_ORDER for _, _, order in ring_bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a sulfonic "
                "acid is not supported yet -- only a ring double bond is "
                "in scope for this first pass (see P-31.1.3)"
            )
        if so3h_carbon not in ring_atoms:
            if not bonds:
                if stereo is not None:
                    raise UnsupportedStructure(
                        "a stereocenter on a substituent branch rather "
                        "than the ring itself is not supported yet (see "
                        "P-92)"
                    )
                return _name_ring_substituent_chain_sulfonic_acid(mol, sulfur_idx, so3h_carbon)
            raise UnsupportedStructure(
                "a sulfonic acid on a substituent branch chain rather "
                "than the ring itself is not supported yet"
            )
        return _name_cyclic_sulfonic_acid(mol, sulfur_idx, so3h_carbon, stereo, ring_bonds)

    return _name_acyclic_sulfonic_acid(mol, sulfur_idx, so3h_carbon, bonds, stereo)


def _name_acyclic_sulfonic_acid(
    mol, sulfur_idx, so3h_carbon, bonds, stereo=None, extra_names=None, required_atoms=frozenset()
):
    """`extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted thiol's 'sulfanyl'), reused
    by `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection. `None`
    keeps the original halogens-only behavior unchanged. `required_atoms`:
    carbon atoms (e.g. every demoted thiol's own carbon) that a candidate
    chain must also carry, mirroring how a plain sulfonic acid always
    requires just its own carbon -- empty by default so existing callers
    are unaffected."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    excluded = {sulfur_idx}
    chains = all_chains(carbon_adjacency(mol))
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if so3h_carbon not in chain_set:
            continue
        if not required_atoms <= chain_set:
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            so3h_carbon in c
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the sulfonic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so3h_locant = position_of[so3h_carbon]
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, so3h_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
