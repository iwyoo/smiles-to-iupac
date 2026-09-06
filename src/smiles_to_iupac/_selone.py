"""Naming of selones (the '-selone' suffix, C=Se with two carbon
substituents), the selenium analogue of ketones, on acyclic saturated or
unsaturated carbon chains and on simple monocyclic saturated rings, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-64.6.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Chalcogen analogues of ketones, pseudoketones and heterones are named
  by using the following suffixes and prefixes: ... =Se '-selone' ..." --
  this module mirrors `_thione.py`'s own '-thione' (C=S) suffix machinery
  exactly, just with selenium instead of sulfur (a selone carbon, like a
  ketone/thione carbon, always has two carbon neighbors and no
  hydrogens). Confirmed via PubChem PUG REST structure match:
  `CC(=[Se])C` -> "propane-2-selone" (the same structure as the Blue
  Book's own worked example for the sulfur case, "propane-2-thione
  (PIN)"), `CC(=[Te])C`'s own sibling case in `_tellone.py` confirms the
  pattern generalizes across chalcogens. The Blue Book's own worked
  example for this exact suffix is "hexane-3-selone (PIN)"
  (`tmp/bluebook/P6.txt`), matching the structure of `CCC(=[Se])CC`.
- Two or more C=Se groups (a diselone) mirrors `_thione.py`'s own
  dithione generalization -- the locant/suffix machinery (ported
  unchanged from `_thione.py`) already generalizes over a list of selone
  locants.
- Like '-thione', '-selone' begins with a consonant, so the parent
  hydride's terminal 'e' is never elided (P-16.3.3).

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_thione.py`/`_ketone.py` (chain and ring both). A selone's C=Se carbon
  is double-bonded to selenium exactly like a ketone's C=O carbon to
  oxygen (confirmed via RDKit's `Chem.FindPotentialStereo` to never
  itself be a potential stereocenter), so this support is unconditional.

Scope, deliberately narrow (mirrors `_thione.py`'s own scope exactly):
one or more C=Se groups on an acyclic chain or a single saturated
monocyclic ring, with no other heteroatom (in particular no ketone C=O,
thione C=S, hydroxyl -OH, or any other chalcogen) anywhere in the
molecule.

P-31.1.3: a monocyclic ring bearing a selone and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-selone', 'cyclohex-3-ene-1-selone',
both confirmed via PubChem PUG REST. Mirrors `_thione.py`'s identical
extension; deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): any oxygen or
sulfur at all, a selone carbon with fewer than two carbon neighbors (a
selenoaldehyde, a different suffix), an aromatic selone carbon,
polycyclic/spiro rings, unsaturation reaching outside the ring or a ring
triple bond, and a selone on a substituent
branch off an otherwise-unsubstituted *saturated* ring. One narrow
*aromatic*-ring exception: `_name_phenyl_chain_selone` names a single
selone whose chain hangs off a plain, unsubstituted benzene ring (e.g.
'1-phenylpropane-2-selone'), mirroring `_thione.py`'s identical benzene-
ring-substituent path -- narrower than the acyclic path: exactly one
selone, no chain unsaturation, and no specified stereocenter, and a
selone carbon directly attached to the ring (an aryl selone) stays out
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
from ._substituents import alpha_sort_key, branch_atom_locant, format_substituent_prefixes, name_branch

_SELENIUM = 34
_ALLOWED_ATOMIC_NUMS = {6, _SELENIUM, *HALOGEN_PREFIXES}


def has_selone_shape(mol) -> bool:
    """True iff `mol` has at least one selenium double-bonded to a carbon
    (a selone, C=Se) -- unlike a plain -SeH/-Se- selenium (`_selenol.py`'s/
    `_selenide.py`'s own loose "any selenium atom" checks), this is precise
    enough to route correctly regardless of check order."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SELENIUM:
            continue
        if atom.GetDegree() == 1:
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() == 2.0:
                return True
    return False


def _validate_and_collect_selones(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_selone`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_selone`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
    selones = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a selone selenium (P-64.6.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting ketone C=O or hydroxyl -OH "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for selones"
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
        elif atomic_num == _SELENIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a selenium bonded to more than one heavy atom (e.g. a "
                    "selenoether/selenide) is out of scope; only an isolated "
                    "selone is supported (P-64.6.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a selone selenium must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "a selenium that isn't a selone (C=Se) is out of scope for "
                    "this module (e.g. a selenol, -SeH)"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a selone on an aromatic ring is out of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a selone carbon with fewer than two carbon neighbors "
                    "(a selenoaldehyde) is a different suffix, which "
                    "this module does not attempt to disambiguate"
                )
            selones.add(atom.GetIdx())
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
    if not selones:
        raise UnsupportedStructure("no selone (C=Se) group found; this module only handles selones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return selones


def _suffix_body(ene_locants, yne_locants, selone_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(selone_locants), multiplied_word(len(selone_locants), "selone")))

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


def _name_from_substituents(chain_length, selone_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, selone_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, selone_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    selone_locant_set = lowest_locant_set(selone_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, selone_locants, ene_locants, yne_locants, grouped)
    return (
        (
            selone_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _selone_locants(position_of, selones, graph):
    locants = []
    for s in selones:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, selones):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in selones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_selone(mol, selones, bonds, stereo=None):
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
        if _selone_locants(position_of, selones, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _selone_locants({a: i + 1 for i, a in enumerate(c)}, selones, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every selone-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more C=Se groups is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, selones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, selones)
            key, name = _candidate_key(chain_length, selone_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _substituents_for_ring(graph, ring_order, halogens, selones):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in selones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, selone_locants, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    parent = "cyclo" + alkane_name(ring_size)
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if not has_unsaturation:
        selone_word = multiplied_word(len(selone_locants), "selone")
        elide = selone_word[0] in "aeiouy"
        stem = parent[:-1] if elide else parent
        if total_subs == 0 and len(selone_locants) == 1:
            # P-14.3.3: the sole substituent on an otherwise unsubstituted
            # ring has no locant to distinguish, e.g. 'cyclohexaneselone'.
            return stem + selone_word
        prefix = format_substituent_prefixes(grouped)
        loc_str = ",".join(str(loc) for loc in sorted(selone_locants))
        return f"{prefix}{stem}-{loc_str}-{selone_word}"

    # A competing ring double/triple bond (P-31.1.3) means the selone's
    # locant is never omittable even when it's the sole substituent, e.g.
    # 'cyclohex-2-ene-1-selone' (confirmed via PubChem) -- mirrors
    # `_thione.py`'s identical treatment.
    stem = parent[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    prefix = format_substituent_prefixes(grouped)
    body, elide_stem = _suffix_body(ene_locants, yne_locants, selone_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _ring_candidate_key(ring_size, selone_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    selone_locant_set = lowest_locant_set(selone_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, selone_locants, ene_locants, yne_locants, grouped)
    return selone_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_bond_locant(position_of, bond_atoms, ring_size):
    pa, pb = position_of[bond_atoms[0]], position_of[bond_atoms[1]]
    return ring_size if {pa, pb} == {1, ring_size} else min(pa, pb)


def _ring_bond_locants(position_of, bonds, ring_size):
    """(ene_locants, yne_locants), both sorted, for every ring C=C/C#C bond
    under this ring numbering -- mirrors `_thione.py`'s identical helper."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _ring_bond_locant(position_of, (a, b), ring_size)
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)


def _ring_branch_stereo_display(graph, ring_order, selones, stereo, halogens):
    """Mirrors `_ketone.py`'s identical helper (itself mirroring
    `_aromatic.py`'s `_stereo_display`): if the ring carries exactly one
    specified stereocenter and that stereocenter sits off the ring on the
    ring's own sole substituent branch (P-92), return that branch's
    ring-attachment atom plus its bracketed "[(<locant><R/S>)-<name>]"
    display (P-91.3). Returns None (the caller keeps its existing outright
    rejection) for more than one stereocenter, or the ring having more or
    fewer than one substituent in total."""
    if len(stereo) != 1:
        return None
    stereo_atom, r_or_s = stereo[0]
    ring_set = set(ring_order)
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_order
        for neighbor in graph[ring_atom]
        if neighbor not in ring_set and neighbor not in selones
    ]
    if len(branch_attachments) != 1:
        return None
    ring_atom, branch_root = branch_attachments[0]
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    display = f"[{descriptor}]" if branch_compound else f"({descriptor})"
    return ring_atom, display


def _name_cyclic_selone(mol, selones, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_thione.py`'s `_name_cyclic_thione`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_ketone.py`'s own case):
    exactly one stereocenter on the ring's sole substituent branch instead
    embeds a bracketed descriptor into that substituent's own name, in
    place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, selones, stereo, halogens)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, selones).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "selone is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, selones, graph)
            if selone_locants is None:
                raise UnsupportedStructure(
                    "a selone not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, selones)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = _ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, selone_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_selone(mol, ring_atoms):
    """Name a selone whose C=Se lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 1-phenylpropane-2-selone. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_thione.py`'s `_name_phenyl_chain_thione`. Narrower than the acyclic
    path above: exactly one selone, no chain unsaturation, and no
    specified stereocenter -- each is a separate follow-up rather than
    being combined with the ring case in this first slice. A selone
    carbon directly attached to the ring (an aryl selone) stays out of
    scope, same as the module's existing acyclic-carbonyl check above."""
    selones = _validate_and_collect_selones(mol, aromatic_ring_atoms=ring_atoms)
    if len(selones) != 1:
        raise UnsupportedStructure(
            "more than one selone group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "selone chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in selones and b[1] not in selones and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "selone chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain selone is not supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, selones)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "selone is not supported yet"
        )
    (selone_selenium,) = selones
    (selone_carbon,) = graph[selone_selenium]
    if selone_carbon == chain_root:
        raise UnsupportedStructure(
            "a selone carbon directly attached to the benzene ring (an "
            "aryl selone) is out of scope for this module (see the "
            "separate aromatic-ring module)"
        )
    if selone_carbon not in chain:
        raise UnsupportedStructure(
            "the selone carbon must lie on the chain hanging off the "
            "benzene ring for this benzene-substituent path"
        )

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        selone_locants = _selone_locants(position_of, selones, graph)
        substituents = {
            position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
        }
        key, name = _candidate_key(chain_length, selone_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_selone(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_selone(mol, ring_atoms)
    selones = _validate_and_collect_selones(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in selones and b[1] not in selones]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_selone(mol, selones, bonds, stereo)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic selone "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a selone is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        return _name_cyclic_selone(mol, selones, stereo, bonds)
    raise UnsupportedStructure(
        "polycyclic and spiro selones are not supported yet"
    )
