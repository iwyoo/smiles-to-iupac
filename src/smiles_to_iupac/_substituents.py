"""Naming of compound (branched) substituent groups, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-29.4.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): "A
  compound substituted substituent group is formed by substituting one or
  more simple substituents into another simple substituent that is
  considered as the principal chain." The free valence (the atom attached to
  the parent chain or ring) is always locant 1 of that principal chain;
  remaining branches are cited as nested substituent prefixes, with identical
  branches grouped under a multiplying prefix.
- P-46 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  principal chain of a substituent group is chosen by the same criteria as a
  parent hydride's (P-44.3/P-45.2) — longest chain, then most substituents,
  then lowest locant set, then lowest locants in citation order — except that
  locant 1 is fixed at the free valence, so there is no choice of numbering
  direction the way there is for a parent hydride.
- P-14.2.2 (Chapter P-1): 'bis', 'tris', 'tetrakis', ... multiply identical
  compound substituent prefixes, instead of 'di', 'tri', 'tetra', ..., to
  avoid ambiguity with a substituent's own internal multiplying prefixes.
- P-14.5.2 (Chapter P-1): alphanumerical order is based on a substituent
  prefix's complete name, so a compound substituent like '(1-methylpropyl)'
  alphabetizes under 'm' (from 'methylpropyl'), ignoring the enclosing
  parentheses and locants. Halogeno prefixes (P-35.2.1, see `_common.py`)
  alphabetize the same way, under their own name ('bromo', 'chloro',
  'fluoro', 'iodo') — no special-casing needed.
- P-35.2.1 (Chapter P-3): a halogen atom (F, Cl, Br, I) directly attached to
  a chain/ring atom is itself a simple substituent group ('fluoro', 'chloro',
  'bromo', 'iodo') with no locants or nested prefixes of its own — the base
  case in `name_branch` below.

Cyclic substituent groups (P-29.3.3) are out of scope and raise
UnsupportedStructure, except the minimal case: a plain, unsubstituted
saturated monocyclic ring hanging off the parent chain (e.g. "cyclohexyl" in
cyclohexylmethanol) is recognized by `_simple_ring_substituent` and named
directly ("cyclo" + `alkyl_name`), without walking into
`_longest_chains_from_root`'s cycle-detection rejection.

A ring substituent may also carry one or more named one-atom groups of its
own (-OH, =O as "oxo", -NH2 as "amino", ...), on any ring atom other than
the attachment point itself (e.g. "(4-hydroxycyclohexyl)" in
`1-(4-hydroxycyclohexyl)ethane-1,2-diol`, PubChem CID 21395558) --
`_ring_substituent_with_named_atoms` reuses the same `{atom_idx: "name"}`-
in-`halogens`-dict convention `_alcohol.py`/`_amine.py`/`_ketone.py`
already use for a chain's own named groups, so a caller opts in simply by
including the ring's named atoms in the `halogens` dict passed to
`name_branch` -- every named atom on the ring must share the same name
(a ring mixing two different named groups falls through to the cyclic-
substituent rejection below instead). The attachment point is fixed at
locant 1 (P-29.2's free-valence rule) and the ring-walk direction is
chosen to give the named atoms the lowest locant set (P-14.5.2),
mirroring how `_alcohol.py`'s own plain-ring numbering picks a direction.
Two or more named atoms are cited together with an ordinary "di"/"tri"
multiplying prefix, e.g. "(3,4-dihydroxycyclohexyl)" -- no PubChem-listed
compound was found for this exact multi-hydroxyl shape, so it's a
reviewed (eyeballed), not independently verified, generalization of the
single-hydroxyl mechanism above (for count 1 it produces byte-identical
output). A ring bearing its
own hydroxyl on the attachment atom itself, any other kind of substituent,
an unsaturated ring, or a polycyclic/spiro ring as a substituent all
remain out of scope and still raise `UnsupportedStructure` via the
ordinary cycle-detection path (see `_simple_ring_substituent`'s own
docstring for exactly which zero-substituent shapes it recognizes).
"""

import re

from ._common import UnsupportedStructure, lowest_locant_set, multiplied_word
from ._numerals import alkane_name, alkyl_name, multiplying_prefix

_LEADING_LOCANTS_RE = re.compile(r"^[\d,\-]+")
_ITALIC_PREFIX_RE = re.compile(r"^(tert|sec|iso)-")


def alpha_sort_key(name: str) -> str:
    """P-14.5.2: alphanumerical ordering ignores locants and italicized
    prefixes like 'tert-' -- only the rest of the name counts (so
    'tert-butyl' sorts under 'b', not 't')."""
    stripped = _LEADING_LOCANTS_RE.sub("", name)
    stripped = _ITALIC_PREFIX_RE.sub("", stripped)
    return stripped.lower()


def _locant_sort_key(locant):
    """Numeric locants sort by value; a non-numeric one (e.g. 'N', P-66.4's
    amine-nitrogen locant) always sorts *before* every numeric one within
    a shared substituent's own locant list -- confirmed against PubChem
    (`CNCC(C)C` -> 'N,2-dimethylpropan-1-amine', not
    '2,N-dimethylpropan-1-amine'; corrects a previously unverified
    assumption from `_amine.py`'s N-prefix/halogen interleaving, PR
    #443, which had no coinciding-name test case to catch this)."""
    return (0, str(locant)) if isinstance(locant, str) else (1, locant)


def format_substituent_prefixes(grouped, omit_locants: bool = False) -> str:
    """grouped: {name -> {"locants": [int or 'N', ...], "compound": bool}}.
    Return the assembled, alphanumerically ordered prefix string
    (P-14.5.2), ready to prepend to a parent name; '' if grouped is empty.

    `omit_locants`: for a mononuclear parent hydride (P-14.3.4.2(a)), every
    substituent's locant is always '1' and never cited, no matter how many
    substituents there are — unlike the ordinary case, where a locant is
    droppable only when every substituent's is unambiguous without it. A
    non-numeric locant (e.g. 'N') is never omitted even under
    `omit_locants=True`, since it marks a different atom than the
    mononuclear parent's own carbon (see `_amine.py`'s N-prefix/halogen
    interleaving) -- so the two locant kinds can coexist in one call, with
    a hyphen separating an omitted-locant part from an explicit one but
    never two consecutive omitted-locant parts (still alphabetically
    interleaved either way, since `sorted(grouped, ...)` runs once up
    front over every name regardless of which kind its own locants are).

    A name cited *only* under non-numeric locants (e.g. two identical
    N-substituents, no coinciding numbered one) always multiplies with
    'di'/'tri'/... even when `compound` is True -- confirmed against
    PubChem ('N,N-di(propan-2-yl)methanesulfonamide', CID 284325, not
    '...bis(propan-2-yl)...') -- this project's established convention
    already hardcoded this for pure N,N-/N,N',N''- citation (`_urea.py`,
    `_amide.py`, `_carbamate.py`) before any of them merged into this
    shared helper; a name whose locants mix a numeral with 'N' (a
    genuinely new coinciding-name citation, e.g.
    'N,4-dimethylbenzenesulfonamide') still follows the ordinary
    compound-aware bis/tris rule, unaffected."""
    entries = []
    for name in sorted(grouped, key=alpha_sort_key):
        info = grouped[name]
        locants = sorted(info["locants"], key=_locant_sort_key)
        all_non_numeric = all(isinstance(loc, str) for loc in locants)
        multiplier_compound = info["compound"] and not all_non_numeric
        multiplier = multiplying_prefix(len(locants), compound=multiplier_compound) if len(locants) > 1 else ""
        # P-16.3.3: enclosing marks escalate one level, (), [], {}, ... --
        # a name that already contains its own '(' (e.g. '4-(2-methylpropyl)
        # phenyl') needs the next mark up, or two same-kind marks would abut
        # ambiguously (PubChem '2-[4-(2-methylpropyl)phenyl]propanoic acid').
        if not info["compound"]:
            display_name = name
        elif "(" not in name:
            display_name = f"({name})"
        else:
            display_name = f"[{name}]"
        explicit = (not omit_locants) or any(isinstance(loc, str) for loc in locants)
        if explicit:
            loc_str = ",".join(str(loc) for loc in locants)
            entries.append((f"{loc_str}-{multiplier}{display_name}", True))
        else:
            entries.append((f"{multiplier}{display_name}", False))

    result = ""
    for i, (text, explicit) in enumerate(entries):
        if i == 0:
            result = text
        else:
            sep = "-" if (explicit or entries[i - 1][1]) else ""
            result += sep + text
    return result


def format_mononuclear_prefixes(entries) -> str:
    """Format substituent prefixes for a mononuclear parent hydride whose
    own atom is the sole skeletal atom (e.g. a phosphane/borane central
    atom, P-14.3.4.2(a) locants always omitted): each distinct name gets
    its own ordinary multiplying prefix (di-, tri-) by its own count;
    when two or more distinct names are present, every one is
    parenthesized except the alphabetically first, regardless of that
    name's own count (P-16.5.1.3.1, per the Blue Book's own published
    errata, https://iupac.qmul.ac.uk/bibliog/BBerrors.html) -- but the
    multiplying prefix itself always sits *outside* those parentheses
    ('ethyldi(methyl)phosphane', not 'ethyl(dimethyl)phosphane'), per
    that same rule's own text and confirmed directly by the Blue Book's
    'ethyldi(methyl)phosphane (PIN)' worked example (`tmp/bluebook/
    P1.html`).

    `entries`: a flat list of `(name, is_compound)` tuples (as returned by
    `name_branch`, not a `grouped` dict like `format_substituent_prefixes`
    above takes). A lone compound name (has its own locant, e.g.
    'propan-2-yl') with nothing else to cite is never parenthesized
    ('propan-2-ylphosphane', P-16.5.1.3.2's own "second and subsequent"
    framing implying the sole substituent needs no enclosing marks at
    all) -- but as soon as there's a second, different substituent to
    cite, the literal P-16.5.1.3.1 text ("the first cited substituent
    never has enclosing marks *unless* it is a compound substituent group
    or includes a locant") means a compound first substituent *does* get
    parenthesized too, e.g. '(propan-2-yl)(propyl)phosphane' -- a
    correction from this project's own earlier (undocumented, PubChem-
    trusted) assumption that only non-first compound names needed
    parentheses; PubChem's raw auto-generated names are already
    documented elsewhere in this project as unreliable for omitting
    required parentheses (see e.g. `_azide.py`/`_ether.py`'s benzene-ring
    paths), and this turned out to be another instance of that. A
    *multiplied* compound name needs its own inner parentheses regardless
    of position to avoid ambiguity ('tri(propan-2-yl)phosphane', PubChem
    CID 80969; confirmed independently via the Blue Book's own
    'ethyldi(propan-2-yl)silane (PIN)' worked example, `tmp/bluebook/
    P1.html` P-16.5.1.3.1), mirroring the 'di(...)' rule already
    established for `_carbamate.py`/`_urea.py`. A multiplied compound
    name mixed with a *different* substituent is out of scope
    (`UnsupportedStructure`) -- PubChem's own naming engine is already
    documented elsewhere as unreliable for phosphane/borane cases with
    3+ distinct substituents, and there's no confirmed Blue Book worked
    example settling the resulting punctuation (a hyphen appears to be
    involved, but not reliably enough to encode blind)."""
    counts = {}
    compound_of = {}
    for name, is_compound in entries:
        counts[name] = counts.get(name, 0) + 1
        compound_of[name] = is_compound
    if len(counts) == 1:
        (name, count), = counts.items()
        if count == 1:
            # A compound name that itself begins with a locant digit (e.g.
            # '4-chlorophenyl') needs enclosing marks even alone, unlike a
            # branched name whose own locant isn't leading ('propan-2-yl',
            # left bare above) -- confirmed via PubChem PUG REST:
            # `Clc1ccc(cc1)P` -> '(4-chlorophenyl)phosphane' (CID 17762777).
            if compound_of[name] and name[0].isdigit():
                return f"({name})"
            return name
        # A digit-leading compound name is itself a *substituted*
        # substituent group (e.g. '4-chlorophenyl' = phenyl substituted by
        # chloro), which takes the irregular 'bis'/'tris' series (P-14.2.2)
        # to avoid ambiguity -- unlike a branched name's own locant, which
        # doesn't make it a substituted-substituent-group in this sense
        # ('tri(propan-2-yl)phosphane', PubChem CID 80969, plain 'tri').
        # Confirmed via PubChem PUG REST: three (4-chlorophenyl) groups on
        # one phosphorus -> 'tris(4-chlorophenyl)phosphane' (CID 70874).
        needs_kis = compound_of[name] and name[0].isdigit()
        wrapped = f"({name})" if compound_of[name] else name
        return multiplying_prefix(count, compound=needs_kis) + wrapped

    if any(counts[name] > 1 and compound_of[name] for name in counts):
        raise UnsupportedStructure(
            "a multiplied compound substituent alongside a different "
            "substituent is not supported yet (P-16.5.1.3.1 "
            "parenthesization for this combination is unconfirmed)"
        )

    ordered = sorted(counts, key=alpha_sort_key)
    parts = []
    for i, name in enumerate(ordered):
        count = counts[name]
        if count > 1:
            # P-16.5.1.3.1's own text: "the multiplicative prefixes are
            # not included in the parentheses" -- confirmed via the Blue
            # Book's own 'ethyldi(methyl)phosphane (PIN)' worked example
            # (`tmp/bluebook/P1.html`), so the prefix sits outside the
            # parens at any position, not just the first.
            needs_kis = compound_of[name] and name[0].isdigit()
            parts.append(multiplying_prefix(count, compound=needs_kis) + (name if i == 0 else f"({name})"))
        elif compound_of[name]:
            parts.append(f"({name})")
        else:
            parts.append(name if i == 0 else f"({name})")
    return "".join(parts)


def _group_substituents(entries):
    """entries: iterable of (locant, name, is_compound)."""
    grouped = {}
    for locant, name, is_compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append(locant)
    return grouped


def _longest_chains_from_root(graph, root, coming_from, halogens):
    """All maximum-length simple paths starting at `root`, extending into the
    subtree away from `coming_from`. Used both for the ordinary case where
    `root` is the free valence itself (P-46: fixed at locant 1, so only one
    direction of travel from it) and, via `_branch_point_candidate_chains`,
    for a single side of a chain that runs *through* the free valence
    (P-29.3.2.2). `halogens` atoms are excluded from the walk (P-35.2.1: a
    halogen is a terminal substituent, never a chain-extending atom), the
    same way `carbon_adjacency` excludes them from a parent hydride's chain
    search."""
    best_length = 0
    best_paths = []

    def walk(node, previous, path, path_set):
        nonlocal best_length, best_paths
        neighbors = [n for n in graph[node] if n != previous and n not in halogens]
        if not neighbors:
            if len(path) > best_length:
                best_length, best_paths = len(path), [list(path)]
            elif len(path) == best_length:
                best_paths.append(list(path))
            return
        for neighbor in neighbors:
            if neighbor in path_set:
                raise UnsupportedStructure(
                    "cyclic substituent groups are not supported yet (see "
                    "P-29.3.3, P-46 for cyclic substituent groups)"
                )
            path.append(neighbor)
            path_set.add(neighbor)
            walk(neighbor, node, path, path_set)
            path_set.discard(neighbor)
            path.pop()

    walk(root, coming_from, [root], {root})
    return best_paths


def _candidate_key(grouped):
    """Sort key implementing P-46's analogue of P-45.2.1-P-45.2.3, most
    preferred first (chain length is fixed by the caller, so it is not part
    of this key)."""
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return -total_count, locant_set, citation_locants


def _simple_ring_substituent(graph, root, coming_from, aromatic_atoms=frozenset()):
    """If the branch hanging off `root` (away from `coming_from`) is a
    single, simple, unsubstituted monocyclic ring with `root` as its only
    attachment point, return (ring_size, is_aromatic); else None (a
    non-ring branch, a ring bearing its own substituent, or any
    polycyclic/spiro/fused shape all fall through to the ordinary
    chain-walk in `name_branch`, which raises `UnsupportedStructure` via
    `_longest_chains_from_root`'s cycle-detection check). `is_aromatic` is
    True only for a plain six-membered all-carbon ring whose every atom is
    in `aromatic_atoms` (benzene as a substituent, i.e. 'phenyl' -- see
    `name_branch`); a caller not passing `aromatic_atoms` (the default
    empty set) only ever gets `is_aromatic=False`, matching every existing
    caller's saturated-ring-only scope unchanged. A non-six-membered
    all-aromatic ring is out of scope (`UnsupportedStructure` via the
    ordinary cycle-detection path, same as any other unrecognized ring
    shape) -- no all-carbon monocyclic aromatic exists at another size for
    a neutral hydrocarbon substituent anyway."""
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    if len(ring_neighbors) != 2:
        return None
    visited = {root}
    previous, current = root, ring_neighbors[0]
    while current != root:
        if current in visited:
            return None
        visited.add(current)
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) != 1:
            return None
        previous, current = current, neighbors[0]
    if aromatic_atoms and visited <= aromatic_atoms:
        if len(visited) != 6:
            return None
        return len(visited), True
    return len(visited), False


def _ring_substituent_with_named_atoms(graph, root, coming_from, halogens):
    """Like `_simple_ring_substituent`, but allows one or more non-
    attachment ring atoms to each carry a single one-atom substituent
    found in `halogens` (that dict's `{atom_idx -> prefix name}`
    convention -- see module docstring), so long as every one of them
    shares the same prefix name (e.g. all "hydroxy", or all "oxo", or all
    "amino" -- a ring mixing two different named substituents falls
    through to the ordinary chain-walk's cyclic-substituent rejection,
    same as an unrecognized ring shape). Returns (ring_size, name,
    locants) with the attachment fixed at locant 1 and the ring-walk
    direction chosen to give the named atoms the lowest locant set
    (P-14.5.2); else None (no named atom found in either direction, a
    mix of different names, or any other shape `_simple_ring_substituent`
    itself would already reject).

    Generalized from an earlier version hardcoded to "hydroxy" only
    (`_alcohol.py`'s own P-44.1.1 tie-break was the only caller that
    needed a ring-as-substituent citation) -- once `_amine.py`/
    `_ketone.py` grew the equivalent P-44.1.1 tie-break for "amino"/"oxo",
    their "chain wins" direction hit this same cyclic-substituent
    citation but with a name this function didn't recognize, so it fell
    through to the generic acyclic walker's misleading "cyclic
    substituent groups are not supported" error instead of actually
    working. For the "hydroxy" case this produces byte-identical output
    to the original (same locants, same name), so the already-verified
    case (PubChem CID 21395558, module docstring) is unaffected."""
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    if len(ring_neighbors) != 2:
        return None

    def walk(start):
        visited = {root}
        previous, current = root, start
        position = 1
        name = None
        locants = []
        while current != root:
            if current in visited:
                return None
            visited.add(current)
            position += 1
            neighbors = [n for n in graph[current] if n != previous]
            named = [n for n in neighbors if n in halogens]
            ring_next = [n for n in neighbors if n not in halogens]
            if named:
                if len(named) != 1 or (name is not None and halogens[named[0]] != name):
                    return None
                name = halogens[named[0]]
                locants.append(position)
            if len(ring_next) != 1:
                return None
            previous, current = current, ring_next[0]
        if not locants:
            return None
        return len(visited), name, tuple(sorted(locants))

    results = [r for r in (walk(n) for n in ring_neighbors) if r is not None]
    if not results:
        return None
    ring_sizes = {size for size, _, _ in results}
    names = {name for _, name, _ in results}
    if len(ring_sizes) != 1 or len(names) != 1:
        return None
    best_locants = min(locants for _, _, locants in results)
    return ring_sizes.pop(), names.pop(), best_locants


def halogenated_phenyl_substituent(graph, aromatic_atoms, root, coming_from, halogens):
    """Like `_ring_substituent_with_named_atoms`, but for a benzene ring
    (fixed at locant 1 = `root`, the ring carbon bonded to the parent atom)
    where zero or more of the other five ring atoms each carry a single
    halogen substituent (an exocyclic neighbor found in `halogens`) instead
    of a hydroxyl -- used by `_phosphane.py`/`_borane.py` to extend their
    existing plain-phenyl support to a halogen-substituted phenyl ring.
    Returns the assembled name (e.g. '4-chlorophenyl',
    '2,3,6-trichlorophenyl'), the ring's own atom indices, and the halogen
    atoms' own indices -- as `(name, ring_atoms, halogen_atoms)` -- using
    the ring-walk direction that gives the halogens the lowest locant set
    overall (P-14.5.2); or None if the ring isn't a plain six-membered
    all-aromatic-carbon ring, a non-attachment ring atom carries anything
    other than a single halogen (or nothing), or `root`'s ring has no
    halogen at all (the ordinary unsubstituted-phenyl path already covers
    that case)."""
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    if len(ring_neighbors) != 2 or root not in aromatic_atoms:
        return None

    def walk(start):
        visited = {root}
        previous, current = root, start
        position = 1
        entries = []
        while current != root:
            if current in visited or current not in aromatic_atoms:
                return None
            visited.add(current)
            position += 1
            neighbors = [n for n in graph[current] if n != previous]
            ring_next = [n for n in neighbors if n not in halogens]
            halogen_neighbors = [n for n in neighbors if n in halogens]
            if len(ring_next) != 1 or len(halogen_neighbors) > 1:
                return None
            if halogen_neighbors:
                entries.append((position, halogens[halogen_neighbors[0]], halogen_neighbors[0]))
            previous, current = current, ring_next[0]
        if len(visited) != 6:
            return None
        return visited, entries

    candidates = [r for r in (walk(n) for n in ring_neighbors) if r is not None and r[1]]
    if not candidates:
        return None
    visited, entries = min(candidates, key=lambda vc: [pos for pos, _, _ in vc[1]])
    grouped = {}
    for pos, name, _ in entries:
        # A name starting with its own locant (e.g. '2-methylpropyl') needs
        # enclosing marks here to keep it from reading as a second ring
        # locant butted up against this one (e.g. '4-2-methylpropylphenyl');
        # one that doesn't (e.g. 'propan-2-yl') needs none (PubChem
        # '2-(4-propan-2-ylphenyl)acetic acid' vs
        # '2-[4-(2-methylpropyl)phenyl]propanoic acid', ibuprofen's PIN).
        info = grouped.setdefault(name, {"locants": [], "compound": name[0].isdigit()})
        info["locants"].append(pos)
    full_name = format_substituent_prefixes(grouped) + "phenyl"
    halogen_atoms = {atom_idx for _, _, atom_idx in entries}
    return full_name, frozenset(visited), frozenset(halogen_atoms)


def plain_alkyl_ring_substituents(mol, graph, ring_atoms):
    """{atom_idx -> name} for every ring atom's sole exocyclic substituent
    that is a plain, fully saturated, acyclic alkyl group (any length,
    branched or unbranched) -- a general-algorithm replacement for what
    used to be a narrower terminal-CH3-only helper (`_common.py`'s
    `plain_methyl_ring_substituents`, removed once every caller migrated
    here), reusing `name_branch` itself instead of hand-rolling a second
    walk. Fed into the same {atom_idx -> prefix name} dict
    `halogen_substituents` builds, so
    `ring_chain_attachment_with_halogens`/`halogenated_phenyl_substituent`
    (both halogen-agnostic, just echoing back whatever name a dict value
    gives) recognize e.g. '4-ethylphenyl'/'4-propan-2-ylphenyl' the same
    way they already recognize '4-methylphenyl' -- confirmed via PubChem
    PUG REST IUPACName ('2-(4-ethylphenyl)acetic acid',
    '2-(4-propan-2-ylphenyl)acetic acid': a compound branch name like
    'propan-2-yl' embeds here with no extra inner parens of its own,
    matching this function's plain-string return, since the *whole*
    ring-plus-substituent unit gets its own outer parens from the
    ordinary compound-substituent citation machinery instead).

    `n` itself directly starting another ring (a nested/fused
    substituent, e.g. a cyclohexyl group hanging off the ring) is
    excluded up front via `IsInRing()` -- out of scope for these 19
    phenyl-chain modules, unverified territory this pilot doesn't
    attempt. Anything else `name_branch` itself rejects (unsaturation, a
    *deeper* nested ring several bonds down, hidden heteroatoms) is
    simply skipped via the `UnsupportedStructure` it already raises for
    exactly those shapes -- the atom is left out of the returned dict,
    falling through to each caller's existing "more than one
    non-halogen, non-methyl exocyclic substituent" rejection unchanged,
    so no separate validation is needed here."""
    alkyls = {}
    for atom in ring_atoms:
        for n in graph[atom]:
            if n in ring_atoms:
                continue
            carbon = mol.GetAtomWithIdx(n)
            if carbon.GetAtomicNum() != 6 or carbon.GetIsAromatic() or carbon.IsInRing():
                continue
            if carbon.GetFormalCharge() != 0 or carbon.GetIsotope() != 0:
                continue
            if not _all_carbon_branch(mol, graph, n, atom):
                # `name_branch`'s walk operates on the whole-molecule
                # graph, heteroatoms included -- it has no way to know a
                # branch like -CH2-C(=O)-OH isn't a plain alkyl chain
                # unless this caller filters it out first (unlike every
                # other `name_branch` caller here, which only ever hands
                # it an already-verified all-carbon halogens-dict branch).
                continue
            try:
                name, _ = name_branch(graph, n, atom, {})
            except UnsupportedStructure:
                continue
            alkyls[n] = name
    return alkyls


def _all_carbon_branch(mol, graph, root, coming_from):
    """True if every atom reachable from `root`, away from
    `coming_from`, is a plain (uncharged, non-isotopic) carbon -- the
    pre-check `plain_alkyl_ring_substituents` needs before trusting
    `name_branch` with an arbitrary ring-atom branch."""
    stack = [(root, coming_from)]
    seen = {root}
    while stack:
        node, previous = stack.pop()
        atom = mol.GetAtomWithIdx(node)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return False
        for neighbor in graph[node]:
            if neighbor == previous or neighbor in seen:
                continue
            seen.add(neighbor)
            stack.append((neighbor, node))
    return True


def name_branch(graph, root, coming_from, halogens=None, aromatic_atoms=None):
    """Name the substituent group hanging off `root`, reached from
    `coming_from` (the parent chain/ring atom). Returns (name, is_compound);
    is_compound is True iff the name carries its own locants/nested prefixes
    (P-29.4) and should be parenthesized when cited as a prefix.

    `halogens`: {atom_idx -> prefix name} (see `_common.halogen_substituents`).
    If `root` is itself a halogen atom, it's a simple substituent with no
    locants or nested prefixes of its own (P-35.2.1) — returned directly,
    with no recursion.

    `aromatic_atoms`: the set of every aromatic atom index in the molecule
    (or None/empty, the default, if the caller's scope has no aromatic
    atoms at all) -- opts a plain benzene-ring branch into being named
    'phenyl' instead of falling through to the cyclic-substituent
    rejection (see `_simple_ring_substituent`), and a halogen-substituted
    one into e.g. '4-chlorophenyl' (see `halogenated_phenyl_substituent`)."""
    halogens = halogens or {}
    aromatic_atoms = aromatic_atoms or frozenset()
    if root in halogens:
        return halogens[root], False

    ring_result = _simple_ring_substituent(graph, root, coming_from, aromatic_atoms)
    if ring_result is not None:
        ring_size, is_aromatic = ring_result
        if is_aromatic:
            return "phenyl", False
        return "cyclo" + alkyl_name(ring_size), False

    if aromatic_atoms:
        halophenyl = halogenated_phenyl_substituent(graph, aromatic_atoms, root, coming_from, halogens)
        if halophenyl is not None:
            name, _, _ = halophenyl
            return name, True

    ring_with_named_atoms = _ring_substituent_with_named_atoms(graph, root, coming_from, halogens)
    if ring_with_named_atoms is not None:
        ring_size, name, locants = ring_with_named_atoms
        loc_str = ",".join(str(loc) for loc in locants)
        name_word = multiplied_word(len(locants), name)
        return f"{loc_str}-{name_word}cyclo{alkyl_name(ring_size)}", True

    _, _, name, is_compound = _select_winning_structure(graph, root, coming_from, halogens)
    return name, is_compound


def _substituent_entries_along_chain(graph, chain, first_previous, halogens, extra_exclusions=None):
    """(position, name, is_compound) for every branch hanging off `chain`
    (a chosen principal chain, root/free-valence somewhere on it), position
    1-based along `chain`. `first_previous`: the atom to exclude when
    collecting `chain[0]`'s own branches -- the external parent atom when
    `chain[0]` is the free valence itself, or `None` when `chain[0]` is a
    genuine leaf reached from the far side of a P-29.3.2.2 branch-point
    chain (its only non-chain neighbors, if any, are real substituents,
    since `_longest_chains_from_root` only stops there when nothing longer
    extends past it).

    `extra_exclusions`: {position -> atom}, for a P-29.3.2.2 branch-point
    chain where the free valence sits mid-chain -- its external parent atom
    is a real neighbor there too (not just at `chain[0]`), and, unlike
    `chain[0]`'s neighbor, isn't already screened out by `chain_set` since
    it's outside the chain entirely."""
    extra_exclusions = extra_exclusions or {}
    chain_set = set(chain)
    entries = []
    for position, atom in enumerate(chain, start=1):
        previous = chain[position - 2] if position > 1 else first_previous
        excluded = extra_exclusions.get(position)
        for branch_root in graph[atom]:
            if branch_root == previous or branch_root == excluded or branch_root in chain_set:
                continue
            sub_name, sub_compound = name_branch(graph, branch_root, atom, halogens)
            entries.append((position, sub_name, sub_compound))
    return entries


def _select_winning_chain(graph, root, coming_from, halogens):
    """The ordinary (non-branch-point) P-46 tie-break: `root` is fixed at
    locant 1 (the free valence is always a chain terminus here). Returns
    (chain, root_position, name, is_compound); `root_position` is always 1,
    kept in the return shape so callers can treat this and
    `_branch_point_candidate_chains` interchangeably via
    `_select_winning_structure`."""
    chains = _longest_chains_from_root(graph, root, coming_from, halogens)
    chain_length = len(chains[0])

    best_key = None
    best_chain = None
    best_name = None
    best_compound = None
    for chain in chains:
        entries = _substituent_entries_along_chain(graph, chain, coming_from, halogens)
        grouped = _group_substituents(entries)
        if grouped and chain_length == 1:
            # P-14.3.4.2(a): the branch's own chain is a single (mononuclear)
            # atom, so any substituent on it has no other possible position
            # and its locant is never cited, e.g. '(chloromethyl)', not
            # '(1-chloromethyl)'.
            name, is_compound = format_substituent_prefixes(grouped, omit_locants=True) + alkyl_name(1), True
        elif grouped:
            name, is_compound = format_substituent_prefixes(grouped) + alkyl_name(chain_length), True
        else:
            name, is_compound = alkyl_name(chain_length), False
        key = _candidate_key(grouped) + (name,)
        if best_key is None or key < best_key:
            best_key, best_chain, best_name, best_compound = key, chain, name, is_compound

    return best_chain, 1, best_name, best_compound


def _branch_point_candidate_chains(graph, root, coming_from, halogens):
    """P-29.3.2.2: when the free-valence atom `root` itself forks into two
    or more branches, the principal chain runs *through* it -- `root`
    becomes an internal locant of the chain (e.g. isopropyl's carbon is
    position 2 of propane, giving PIN 'propan-2-yl', not position 1 of the
    pre-PIN CAS-style '1-methylethyl') -- rather than always starting at
    it as `_select_winning_chain` assumes. Mirrors
    `_radical.py::_name_branch_point_radical`'s P-29.3.2.2 derivation, but
    generalized to let each branch be itself further branched (reusing
    `_longest_chains_from_root`'s own recursive walk per branch, instead of
    requiring every branch to be a plain unbranched chain the way the
    radical module's narrower `linear_branch` does).

    Returns (chain, root_position, name, is_compound), choosing among tied
    candidates by the same P-46 criteria `_select_winning_chain` uses
    (`_candidate_key`) -- chain length and `root_position` are fixed before
    that tie-break runs (P-29.2: the free valence must get the lowest
    locant the chain allows), so, like `_candidate_key`'s own docstring
    notes for chain length, they are not themselves part of the key. `None`
    if `root` isn't a branch point (fewer than two non-halogen branches);
    callers fall back to `_select_winning_chain` in that case."""
    branch_roots = [n for n in graph[root] if n != coming_from and n not in halogens]
    if len(branch_roots) < 2:
        return None

    paths = {b: _longest_chains_from_root(graph, b, root, halogens) for b in branch_roots}
    lengths = {b: len(paths[b][0]) for b in branch_roots}

    if len(branch_roots) == 3 and all(length == 1 for length in lengths.values()):
        # P-29.6.1: the retained name 'tert-butyl' is the PIN for the
        # unsubstituted (CH3)3C- group, never the general rule's own
        # '2-methylpropan-2-yl' -- and, being a single retained word with no
        # locant of its own, it is never parenthesized as a compound prefix.
        return [root], 1, "tert-butyl", False

    sorted_lengths = sorted(lengths.values(), reverse=True)
    target_total = sorted_lengths[0] + sorted_lengths[1]
    spine_pairs = [
        (a, b)
        for i, a in enumerate(branch_roots)
        for b in branch_roots[i + 1 :]
        if lengths[a] + lengths[b] == target_total
    ]

    orientations = []
    for a, b in spine_pairs:
        if lengths[a] == lengths[b]:
            orientations.extend([(a, b), (b, a)])
        elif lengths[a] < lengths[b]:
            orientations.append((a, b))
        else:
            orientations.append((b, a))

    chain_length = target_total + 1
    stem = alkane_name(chain_length)[:-1]

    best_key = None
    best_chain = None
    best_position = None
    best_name = None
    for before_branch, after_branch in orientations:
        extra_roots = [r for r in branch_roots if r not in (before_branch, after_branch)]
        for path_before in paths[before_branch]:
            for path_after in paths[after_branch]:
                spine = list(reversed(path_before)) + [root] + path_after
                root_locant = len(path_before) + 1
                entries = _substituent_entries_along_chain(
                    graph, spine, None, halogens, extra_exclusions={root_locant: coming_from}
                )
                for extra in extra_roots:
                    sub_name, sub_compound = name_branch(graph, extra, root, halogens)
                    entries.append((root_locant, sub_name, sub_compound))
                grouped = _group_substituents(entries)
                prefix = format_substituent_prefixes(grouped)
                name = f"{prefix}{stem}-{root_locant}-yl"
                key = _candidate_key(grouped) + (name,)
                if best_key is None or key < best_key:
                    best_key = key
                    best_chain = spine
                    best_position = root_locant
                    best_name = name

    return best_chain, best_position, best_name, True


def _select_winning_structure(graph, root, coming_from, halogens):
    """(chain, root_position, name, is_compound) for the substituent group
    hanging off `root`: the P-29.3.2.2 branch-point chain
    (`_branch_point_candidate_chains`) when `root` forks into two or more
    branches, else the ordinary root-is-locant-1 chain
    (`_select_winning_chain`). Shared by `name_branch` (which only needs
    the name) and `branch_atom_locant` below (which also needs to know
    which chain won and where `root` sits on it), so the two can never
    disagree about which chain was chosen."""
    branch_point = _branch_point_candidate_chains(graph, root, coming_from, halogens)
    if branch_point is not None:
        return branch_point
    return _select_winning_chain(graph, root, coming_from, halogens)


def branch_atom_locant(graph, root, coming_from, atom_idx, halogens=None):
    """The position (1-based) of `atom_idx` on the winning principal chain
    of the substituent group named by
    `name_branch(graph, root, coming_from, halogens)` -- from
    `_select_winning_structure`, so the two always agree on which chain
    that is (and where `root` itself sits on it, per P-29.3.2.2). Used to
    cite a stereodescriptor at the front of a compound substituent prefix
    (P-91.3), e.g. the '1' in '[(1S)-1-chloropropyl]benzene'.

    Raises `UnsupportedStructure` if `atom_idx` isn't on that winning
    chain at all (e.g. it sits on a branch off the chain instead) --
    citing a stereodescriptor in that case would need a locant this
    module doesn't assign to anything, so it's out of scope rather than
    silently wrong."""
    halogens = halogens or {}
    chain, _, _, _ = _select_winning_structure(graph, root, coming_from, halogens)
    if atom_idx not in chain:
        raise UnsupportedStructure(
            "a specified stereocenter that isn't on the substituent's own "
            "principal chain (P-46) is not supported yet"
        )
    return chain.index(atom_idx) + 1
