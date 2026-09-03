"""Naming of nitriles (the '-nitrile' suffix, terminal -C#N) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-66.5, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  'nitrile' is the preselected suffix for -C#N, cited in a combined
  chain/suffix name the same way 'al'/'one'/'ol'/'amine' are for their own
  groups (`_aldehyde.py`/`_ketone.py`/`_alcohol.py`/`_amine.py`), e.g.
  'propanenitrile', 'pent-4-enenitrile'. Table 3.3 ranks 'nitrile' junior to
  'amide' and senior to 'al'; this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior or junior group, since
  only C, N (the nitrile nitrogen itself), and halogen atoms are accepted at
  all (any oxygen routes elsewhere in `core.py`, or is out of scope here).
- Unlike 'al'/'one'/'amine', 'nitrile' begins with a consonant, so the
  euphonic elision `_aldehyde.py` applies before a vowel-initial suffix
  never triggers here: the stem's final 'e' (from the parent alkane/alkene/
  alkyne name) is always kept, e.g. 'ethanenitrile', not 'ethannitrile';
  'pent-4-enenitrile', not 'pent-4-ennitrile'.
- Like -CHO (`_aldehyde.py`), a -C#N carbon is always a chain terminus (it
  has exactly one carbon neighbor, being otherwise saturated by the triple
  bond to N), so it is never a genuine locant choice: whichever end of the
  principal chain bears it is numbered C1 (P-44.4.1.8, suffix locants
  minimized first), and that locant is always '1' and never cited (P-14.3.3):
  'propanenitrile', not 'propane-1-nitrile'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-nitrile name reuses the same mechanics as
  `_unsaturated.py`/`_aldehyde.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'nitrile' suffix
  word is appended directly onto the last such segment, since it never
  carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'nitrile' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one nitrile group (a dinitrile, P-66.6.3) - future work.
- Any oxygen - routed to a different module by `core.py`, or out of scope
  entirely if this module is called directly on one.
- A nitrile nitrogen that isn't a plain, isolated -C#N (any degree other
  than 1, or a bond order other than 3.0 to its one carbon neighbor).
- A nitrile carbon with other than exactly one carbon neighbor: zero (bare
  HC#N) or two-or-more (not a valid nitrile shape) - neither is a terminal
  substitutive -C#N.
- An aromatic nitrile carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory, except for one narrow case: a
  nitrile's chain hanging off a single plain, unsubstituted benzene ring
  with no other substituent on the ring (`_name_phenyl_chain_nitrile`,
  e.g. '3-phenylpropanenitrile'), mirroring `_aldehyde.py`/`_amide.py`'s
  identical benzene-ring-substituent path. Narrower than the acyclic path:
  no chain unsaturation, no specified stereocenter, and no substituted
  benzene/naphthalene - each a separate follow-up.
- -C#N on a ring (the 'carbonitrile' suffix, P-66.5.1.2, a substituent-style
  name rather than this module's parent-chain suffix) - deferred entirely;
  only an acyclic terminal -C#N is supported here.
- Any other heteroatom (O, S, ...), or a nitrile carbon entangled with
  another heteroatom.

P-91.3/P-92: a molecule with one
or more *specified* tetrahedral stereocenters -- every one on the
principal chain itself, no unspecified one alongside them, and no
C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
ascending locant order, same pattern as `_amine.py`/`_thiol.py` (chain
only, this module has no ring support). The nitrile carbon itself (sp,
triple-bonded to nitrogen) is never a potential stereocenter, confirmed
via RDKit `FindPotentialStereo` on `CC[C@@H](C)C#N`.
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
    is_plain_benzene_ring,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ordered_chain,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_NITRILE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def has_nitrile_shape(mol) -> bool:
    """True iff `mol` has at least one nitrogen atom of degree 1, triple-
    bonded to a carbon (a plain -C#N pattern), regardless of whether the
    rest of the molecule is in scope. Used by `core.py` to route ahead of
    the amine dispatch within the no-oxygen nitrogen branch, since a
    nitrile nitrogen would otherwise look like an unhandled shape to
    `name_amine`'s own bond-order check."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        (neighbor,) = atom.GetNeighbors()
        if bond.GetBondTypeAsDouble() == _NITRILE_ORDER and neighbor.GetAtomicNum() == 6:
            return True
    return False


def _validate_and_collect_nitriles(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return the set of nitrile-nitrogen atom indices.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_nitrile`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_nitrile`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
    nitriles = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitrile nitrogen (P-66.5) and "
                "halogen substituents (P-35.2.1) are not supported yet"
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
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a nitrogen bonded to more than one heavy atom is not a "
                    "plain nitrile (-C#N) and is out of scope for this "
                    "module (P-66.5)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a nitrile nitrogen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != _NITRILE_ORDER:
                raise UnsupportedStructure(
                    "a nitrogen not triple-bonded to its one carbon neighbor "
                    "is not a nitrile and is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a nitrile carbon on an aromatic ring is out of scope "
                    "for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a nitrile carbon with other than exactly one carbon "
                    "neighbor is out of scope for this module (Table 3.3)"
                )
            nitriles.add(atom.GetIdx())
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
    if not nitriles:
        raise UnsupportedStructure("no nitrile (-C#N) group found; this module only handles nitriles")
    if len(nitriles) > 1:
        raise UnsupportedStructure("more than one nitrile group (a dinitrile) is not supported yet")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return nitriles


def _suffix_body(ene_locants, yne_locants, nitrile_count):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'nitrile'
    ending (e.g. '4-enenitrile'). Like `_aldehyde.py`'s '-al', 'nitrile'
    never carries its own locant (see module docstring), so it is glued
    directly onto the preceding word instead of getting a hyphenated locant
    segment of its own; unlike '-al', it never elides a preceding stem's
    trailing 'e' either, since 'nitrile' starts with a consonant."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    nitrile_word = multiplied_word(nitrile_count, "nitrile")
    words = [word for _, word in segments] + [nitrile_word]

    locanted_parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words[:-1])
    ]
    body = "-".join(locanted_parts) + words[-1] if locanted_parts else words[-1]
    return body


def _name_from_substituents(chain_length, nitrile_count, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, nitrile_count)
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The nitrile locant set still drives
    orientation choice even though it is never printed (see module
    docstring)."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    nitrile_locant_set = lowest_locant_set(nitrile_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(nitrile_locants), ene_locants, yne_locants, grouped)
    return (
        (
            nitrile_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _nitrile_locants(position_of, nitriles, graph):
    locants = []
    for n in nitriles:
        (carbon,) = graph[n]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, nitriles):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in nitriles]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_nitrile(mol, nitriles, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch is out of scope), and the winning candidate's own
    locants are used to format a "(<locant><R/S>,...)-" prefix onto the
    final name."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _nitrile_locants(position_of, nitriles, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in position_of for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the nitrile-bearing carbon (and/or multiple bond, and/or a "
            "stereocenter on a substituent branch, see P-92) does not lie "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing the -C#N group (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            nitrile_locants = _nitrile_locants(position_of, nitriles, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, nitriles)
            key, name = _candidate_key(chain_length, nitrile_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_nitrile(mol, ring_atoms):
    """Name a nitrile whose -C#N lies entirely on a single unbranched chain
    hanging off one atom of an otherwise-plain, unsubstituted benzene ring
    -- e.g. 3-phenylpropanenitrile. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than the
    acyclic path above: no chain unsaturation and no specified stereocenter
    -- each is a separate follow-up (see
    tasks/phenyl-substituent-on-nitrile-chain.md's scope note)."""
    nitriles = _validate_and_collect_nitriles(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "nitrile chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in nitriles and b[1] not in nitriles and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "nitrile chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain nitrile is not supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, nitriles)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "nitrile is not supported yet"
        )
    (nitrile_nitrogen,) = nitriles
    (nitrile_carbon,) = graph[nitrile_nitrogen]
    if chain[-1] != nitrile_carbon:
        raise UnsupportedStructure(
            "the nitrile carbon must be the chain's far terminus from the "
            "benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -C#N group directly attached to the benzene ring (the "
            "'carbonitrile' suffix, P-66.5.1.2) uses a separate naming "
            "construction, out of scope for this acyclic-chain-parent "
            "module"
        )

    ordered = list(reversed(chain))
    chain_length = len(ordered)
    position_of = {atom: i + 1 for i, atom in enumerate(ordered)}
    halogens = halogen_substituents(mol)
    substituents = {
        position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, 1, [], [], grouped)


def name_nitrile(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_nitrile(mol, ring_atoms)
    nitriles = _validate_and_collect_nitriles(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    # Exclude each C#N nitrile bond itself: `non_single_bonds` reports it as
    # order 3.0 same as a C#C, but it isn't a chain 'yne' bond (one endpoint
    # is the nitrile nitrogen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in nitriles and b[1] not in nitriles]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a nitrile on a ring (the 'carbonitrile' suffix, P-66.5.1.2) is "
            "out of scope for this module; only an acyclic terminal -C#N is "
            "supported"
        )
    return _name_acyclic_nitrile(mol, nitriles, bonds, stereo)
