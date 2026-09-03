"""Naming of tellones (the '-tellone' suffix, C=Te with two carbon
substituents), the tellurium analogue of ketones, on acyclic saturated or
unsaturated carbon chains and on simple monocyclic saturated rings, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-64.6.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Chalcogen analogues of ketones, pseudoketones and heterones are named
  by using the following suffixes and prefixes: ... =Te '-tellone' ..."
  -- this module mirrors `_thione.py`'s own '-thione' (C=S) suffix
  machinery exactly, just with tellurium instead of sulfur (a tellone
  carbon, like a ketone/thione carbon, always has two carbon neighbors
  and no hydrogens). Confirmed via PubChem PUG REST structure match:
  `CC(=[Te])C` -> "propane-2-tellone" (the same structure as the Blue
  Book's own worked example for the sulfur case, "propane-2-thione
  (PIN)").
- Two or more C=Te groups (a ditellone) mirrors `_thione.py`'s own
  dithione generalization -- the locant/suffix machinery (ported
  unchanged from `_thione.py`) already generalizes over a list of tellone
  locants.
- Like '-thione', '-tellone' begins with a consonant, so the parent
  hydride's terminal 'e' is never elided (P-16.3.3).

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_thione.py`/`_ketone.py` (chain and ring both). A tellone's C=Te
  carbon is double-bonded to tellurium exactly like a ketone's C=O carbon
  to oxygen (confirmed via RDKit's `Chem.FindPotentialStereo` to never
  itself be a potential stereocenter), so this support is unconditional.

Scope, deliberately narrow (mirrors `_thione.py`'s own scope exactly):
one or more C=Te groups on an acyclic chain or a single saturated
monocyclic ring, with no other heteroatom (in particular no ketone C=O,
thione C=S, hydroxyl -OH, or any other chalcogen) anywhere in the
molecule. Explicitly out of scope (raise `UnsupportedStructure`): any
oxygen or sulfur at all, a tellone carbon with fewer than two carbon
neighbors (a telluroaldehyde, a different suffix), an aromatic tellone
carbon, polycyclic/spiro/unsaturated rings, and a tellone on a
substituent branch off an otherwise-unsubstituted *saturated* ring. One
narrow *aromatic*-ring exception: `_name_phenyl_chain_tellone` names a
single tellone whose chain hangs off a plain, unsubstituted benzene ring
(e.g. '1-phenylpropane-2-tellone'), mirroring `_thione.py`'s identical
benzene-ring-substituent path -- narrower than the acyclic path: exactly
one tellone, no chain unsaturation, and no specified stereocenter, and a
tellone carbon directly attached to the ring (an aryl tellone) stays out
of scope for this module.
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
    ring_cycle,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_TELLURIUM = 52
_ALLOWED_ATOMIC_NUMS = {6, _TELLURIUM, *HALOGEN_PREFIXES}


def has_tellone_shape(mol) -> bool:
    """True iff `mol` has at least one tellurium double-bonded to a carbon
    (a tellone, C=Te) -- unlike a plain -TeH/-Te- tellurium (`_tellurol.py`'s/
    `_telluride.py`'s own loose "any tellurium atom" checks), this is precise
    enough to route correctly regardless of check order."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _TELLURIUM:
            continue
        if atom.GetDegree() == 1:
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() == 2.0:
                return True
    return False


def _validate_and_collect_tellones(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_tellone`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_tellone`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged."""
    tellones = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a tellone tellurium (P-64.6.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting ketone C=O or hydroxyl -OH "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for tellones"
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
        elif atomic_num == _TELLURIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a tellurium bonded to more than one heavy atom (e.g. a "
                    "telluroether/telluride) is out of scope; only an isolated "
                    "tellone is supported (P-64.6.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a tellone tellurium must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "a tellurium that isn't a tellone (C=Te) is out of scope for "
                    "this module (e.g. a tellurol, -TeH)"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a tellone on an aromatic ring is out of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a tellone carbon with fewer than two carbon neighbors "
                    "(a telluroaldehyde) is a different suffix, which "
                    "this module does not attempt to disambiguate"
                )
            tellones.add(atom.GetIdx())
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
    if not tellones:
        raise UnsupportedStructure("no tellone (C=Te) group found; this module only handles tellones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return tellones


def _suffix_body(ene_locants, yne_locants, tellone_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(tellone_locants), multiplied_word(len(tellone_locants), "tellone")))

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


def _name_from_substituents(chain_length, tellone_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, tellone_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, tellone_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    tellone_locant_set = lowest_locant_set(tellone_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, tellone_locants, ene_locants, yne_locants, grouped)
    return (
        (
            tellone_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _tellone_locants(position_of, tellones, graph):
    locants = []
    for s in tellones:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, tellones):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in tellones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_tellone(mol, tellones, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_thione.py`'s identical treatment), and the winning
    candidate's own locants are used to format a "(<locant><R/S>,...)-"
    prefix onto the name, ascending locant order (P-91.3)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _tellone_locants(position_of, tellones, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _tellone_locants({a: i + 1 for i, a in enumerate(c)}, tellones, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every tellone-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more C=Te groups is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            tellone_locants = _tellone_locants(position_of, tellones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, tellones)
            key, name = _candidate_key(chain_length, tellone_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _substituents_for_ring(graph, ring_order, halogens, tellones):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in tellones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, tellone_locants, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    tellone_word = multiplied_word(len(tellone_locants), "tellone")
    elide = tellone_word[0] in "aeiouy"
    stem = parent[:-1] if elide else parent

    if total_subs == 0 and len(tellone_locants) == 1:
        # P-14.3.3: the sole substituent on an otherwise unsubstituted ring
        # has no locant to distinguish, e.g. 'cyclohexanetellone'.
        return stem + tellone_word

    prefix = format_substituent_prefixes(grouped)
    loc_str = ",".join(str(loc) for loc in sorted(tellone_locants))
    return f"{prefix}{stem}-{loc_str}-{tellone_word}"


def _ring_candidate_key(ring_size, tellone_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    tellone_locant_set = lowest_locant_set(tellone_locants)
    name = _ring_name_from_substituents(ring_size, tellone_locants, grouped)
    return tellone_locant_set, locant_set, citation_locants, name


def _name_cyclic_tellone(mol, tellones, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must lie on
    the ring itself (P-92: a stereocenter on a substituent branch is out
    of scope, mirroring `_thione.py`'s `_name_cyclic_thione`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
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
            tellone_locants = _tellone_locants(position_of, tellones, graph)
            if tellone_locants is None:
                raise UnsupportedStructure(
                    "a tellone not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, tellones)
            key = _ring_candidate_key(ring_size, tellone_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_tellone(mol, ring_atoms):
    """Name a tellone whose C=Te lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 1-phenylpropane-2-tellone. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_thione.py`'s `_name_phenyl_chain_thione`. Narrower than the acyclic
    path above: exactly one tellone, no chain unsaturation, and no
    specified stereocenter -- each is a separate follow-up rather than
    being combined with the ring case in this first slice. A tellone
    carbon directly attached to the ring (an aryl tellone) stays out of
    scope, same as the module's existing acyclic-carbonyl check above."""
    tellones = _validate_and_collect_tellones(mol, aromatic_ring_atoms=ring_atoms)
    if len(tellones) != 1:
        raise UnsupportedStructure(
            "more than one tellone group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "tellone chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in tellones and b[1] not in tellones and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "tellone chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain tellone is not supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, tellones)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "tellone is not supported yet"
        )
    (tellone_tellurium,) = tellones
    (tellone_carbon,) = graph[tellone_tellurium]
    if tellone_carbon == chain_root:
        raise UnsupportedStructure(
            "a tellone carbon directly attached to the benzene ring (an "
            "aryl tellone) is out of scope for this module (see the "
            "separate aromatic-ring module)"
        )
    if tellone_carbon not in chain:
        raise UnsupportedStructure(
            "the tellone carbon must lie on the chain hanging off the "
            "benzene ring for this benzene-substituent path"
        )

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        tellone_locants = _tellone_locants(position_of, tellones, graph)
        substituents = {
            position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
        }
        key, name = _candidate_key(chain_length, tellone_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_tellone(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_tellone(mol, ring_atoms)
    tellones = _validate_and_collect_tellones(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in tellones and b[1] not in tellones]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_tellone(mol, tellones, bonds, stereo)
    if num_rings == 1:
        if bonds:
            raise UnsupportedStructure(
                "unsaturated rings are not supported yet (see P-31.1.3, "
                "cycloalkenes and cycloalkynes)"
            )
        return _name_cyclic_tellone(mol, tellones, stereo)
    raise UnsupportedStructure(
        "polycyclic and spiro tellones are not supported yet"
    )
