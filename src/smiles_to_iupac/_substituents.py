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

- P-29.2, P-32.1.1, P-46.1 (Chapters P-2/P-3/P-4): with the molecule's bond
  orders available (`mol`), the principal chain of a substituent group is
  the longest chain through the free-valence atom, then the one with more
  multiple bonds, then lower free-valence, multiple-bond, and substituent
  locants; the free-valence bond itself selects 'yl', 'ylidene', or 'ylidyne'.

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

import contextvars

from rdkit import Chem

from ._multiplicative_text import enclose
from ._free_valence import SUFFIX_OF_ORDER
from ._common import (
    UnsupportedStructure,
    alpha_sort_key,
    group_substituents,
    heteroaromatic_monocycle_name,
    multiplied_word,
    ring_cycle,
    substituent_locant_set_and_citation,
    unsaturation_suffix,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix

# alpha_sort_key lives in _common.py now; re-imported here (not redefined)
# since ~67 modules already import it from this module.


def _locant_sort_key(locant):
    """Numeric locants sort by value; a non-numeric one (e.g. 'N', P-66.4's
    amine-nitrogen locant) always sorts *before* every numeric one within
    a shared substituent's own locant list -- confirmed against PubChem
    (`CNCC(C)C` -> 'N,2-dimethylpropan-1-amine', not
    '2,N-dimethylpropan-1-amine'; corrects a previously unverified
    assumption from `_amine.py`'s N-prefix/halogen interleaving, PR
    #443, which had no coinciding-name test case to catch this)."""
    return (0, str(locant)) if isinstance(locant, str) else (1, locant)


_PLAIN_STEM_PREFIX = None


def is_plain_stem_prefix(name: str) -> bool:
    """True for a prefix that is only a parent stem with its own locants
    ('propan-2-yl', 'prop-2-en-1-yl', 'cyclohex-2-en-1-yl'): such a prefix is
    multiplied with 'di'/'tri' in parentheses, while a substituted prefix
    ('2-chloroethyl', 'bromomethyl') takes 'bis'/'tris' (P-16.3.5, P-16.3.2)."""
    global _PLAIN_STEM_PREFIX
    if _PLAIN_STEM_PREFIX is None:
        import re

        stems = {alkane_name(n)[:-3] for n in range(1, 41) if alkane_name(n).endswith("ane")}
        stems |= {"meth", "eth", "prop", "but", "naphthalen", "anthracen", "phenanthren"}
        escaped = "|".join(sorted(map(re.escape, stems), key=len, reverse=True))
        _PLAIN_STEM_PREFIX = re.compile(
            rf"^(?:\d+H-)?(?:cyclo)?(?:{escaped})(?:a|an)?(?:-[\d,]+-(?:di|tri|tetra)?(?:en|yn))*-?[\d,]*-?(?:(?:di|tri|tetra)?(?:en|yn))?(?:yl|ylidene|ylidyne|diyl)$"
        )
    return bool(_PLAIN_STEM_PREFIX.match(name))


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
        multiplier_compound = (
            info["compound"] and not is_plain_stem_prefix(name)
        ) or (name[:1] in "([{" and not name.startswith("(\u03b7"))
        multiplier = multiplying_prefix(len(locants), compound=multiplier_compound) if len(locants) > 1 else ""
        if "multiplier" in info:
            multiplier = info["multiplier"]
        # P-16.3.3: enclosing marks escalate one level, (), [], {}, ... --
        # a name that already contains its own '(' (e.g. '4-(2-methylpropyl)
        # phenyl') needs the next mark up, or two same-kind marks would abut
        # ambiguously (PubChem '2-[4-(2-methylpropyl)phenyl]propanoic acid').
        display_name = enclose(name) if info["compound"] else name
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


def wrap_marks(name: str) -> str:
    """Enclose `name` with the next mark in the nesting order ( ) [ ] { } (P-16.5.4)."""
    if name[:1] in "([{" and name[-1:] in ")]}":
        return name
    if name.startswith("\x01"):
        return f"({name})"
    if "{" in name:
        return f"({name})"
    if "[" in name:
        return "{" + name + "}"
    if "(" in name:
        return f"[{name}]"
    return f"({name})"


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
            if compound_of[name] and (name[0].isdigit() or "-" not in name):
                return wrap_marks(name)
            return name
        # A digit-leading compound name is itself a *substituted*
        # substituent group (e.g. '4-chlorophenyl' = phenyl substituted by
        # chloro), which takes the irregular 'bis'/'tris' series (P-14.2.2)
        # to avoid ambiguity -- unlike a branched name's own locant, which
        # doesn't make it a substituted-substituent-group in this sense
        # ('tri(propan-2-yl)phosphane', PubChem CID 80969, plain 'tri').
        # Confirmed via PubChem PUG REST: three (4-chlorophenyl) groups on
        # one phosphorus -> 'tris(4-chlorophenyl)phosphane' (CID 70874).
        needs_kis = compound_of[name] and not is_plain_stem_prefix(name)
        wrapped = wrap_marks(name) if compound_of[name] else name
        return multiplying_prefix(count, compound=needs_kis) + wrapped

    ordered = sorted(counts, key=alpha_sort_key)
    if counts[ordered[0]] > 1 and compound_of[ordered[0]]:
        # Unlike a non-first multiplied compound name (confirmed below via
        # 'bromodi(ethenyl)stibane (PIN)', `tmp/bluebook/P6a.txt`
        # ~8608-8609 -- the multiplying prefix and enclosing marks both
        # land correctly on a non-first name in the loop below), no
        # worked example confirms whether a multiplied compound name
        # that sorts *first* still needs its own enclosing marks (the
        # loop below's `name if i == 0 else f"({name})"` bare-first
        # shortcut, written for a plain first name, would silently drop
        # them) -- left unsupported rather than guessed.
        raise UnsupportedStructure(
            "a multiplied compound substituent sorting alphabetically "
            "first, alongside a different substituent, is not supported "
            "yet (P-16.5.1.3.1 parenthesization for this combination is "
            "unconfirmed)"
        )

    parts = []
    for i, name in enumerate(ordered):
        count = counts[name]
        if count > 1:
            # P-16.5.1.3.1's own text: "the multiplicative prefixes are
            # not included in the parentheses" -- confirmed via the Blue
            # Book's own 'ethyldi(methyl)phosphane (PIN)' worked example
            # (`tmp/bluebook/P1.html`), so the prefix sits outside the
            # parens at any position, not just the first.
            needs_kis = compound_of[name] and not is_plain_stem_prefix(name)
            parts.append(multiplying_prefix(count, compound=needs_kis) + (name if i == 0 else wrap_marks(name)))
        elif compound_of[name]:
            parts.append(wrap_marks(name))
        else:
            parts.append(name if i == 0 else wrap_marks(name))
    return "".join(parts)


def _group_substituents(entries):
    """entries: iterable of (locant, name, is_compound)."""
    grouped = {}
    for locant, name, is_compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append(locant)
    return grouped


def _longest_chains_from_root(graph, root, coming_from, halogens, mol=None, aromatic_atoms=frozenset()):
    """All maximum-length simple paths starting at `root`, extending into the
    subtree away from `coming_from`. Used both for the ordinary case where
    `root` is the free valence itself (P-46: fixed at locant 1, so only one
    direction of travel from it) and, via `_branch_point_candidate_chains`,
    for a single side of a chain that runs *through* the free valence
    (P-29.3.2.2). `halogens` atoms are excluded from the walk (P-35.2.1: a
    halogen is a terminal substituent, never a chain-extending atom), the
    same way `carbon_adjacency` excludes them from a parent hydride's chain
    search.

    `mol`: when given, defends against a caller that failed to validate its
    whole branch as carbon-plus-`halogens` before calling in -- every walked
    atom (including `root` itself) must be carbon or already excluded via
    `halogens`, and every bond *within* the branch (not `root`'s own
    attachment bond to `coming_from`, which some callers -- e.g.
    `_hydrazone.py`'s '-ylidene' construction -- legitimately make a double
    bond and account for separately) must be single, else this raises
    rather than silently treating an unvalidated heteroatom (an amine
    nitrogen, an ether oxygen, ...) or an unsaturated bond (name_branch has
    no ene/yne machinery of its own) as if it were an ordinary saturated
    chain-extending carbon (found via several independent
    `smiles-to-iupac-realdata-test` pubchem diffs, PR #487/#488 for the
    atom-type gap and a further diff -- a phenol's ring bearing a plain
    'ethyl' *and* a vinyl substituent came out as '2,3-diethylphenol', the
    double bond silently vanishing -- for this bond-order gap). `None`
    (the default) keeps every not-yet-migrated caller's original behavior
    unchanged.

    `aromatic_atoms`: when `mol` is also given, a neighbor that is itself
    the sole attachment point of a separate, plain, unsubstituted ring
    (checked via `_simple_ring_substituent`) ends the walk there instead
    of descending into it, so a compound substituent like
    '(3-cyclohexylpropyl)' or '(3-phenylpropyl)' is reachable at all --
    without this, the walk would keep going around that ring's own bonds
    and hit its cycle-detection rejection just below."""
    best_length = 0
    best_paths = []

    def walk(node, previous, path, path_set):
        nonlocal best_length, best_paths
        if mol is not None and mol.GetAtomWithIdx(node).GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "a heteroatom in a compound substituent branch, other than "
                "a recognized halogen/named group, is not supported yet"
            )
        if (
            mol is not None
            and previous != coming_from
            and mol.GetBondBetweenAtoms(node, previous).GetBondTypeAsDouble() != 1.0
        ):
            raise UnsupportedStructure(
                "unsaturation in a compound substituent branch is not "
                "supported yet"
            )
        neighbors = [n for n in graph[node] if n != previous and n not in halogens]
        if mol is not None:
            neighbors = [
                n
                for n in neighbors
                if not (
                    mol.GetAtomWithIdx(n).IsInRing()
                    and mol.GetBondBetweenAtoms(node, n).GetBondTypeAsDouble() == 1.0
                    and _is_ring_branch_root(graph, n, node, aromatic_atoms, mol)
                )
            ]
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
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    return -total_count, locant_set, citation_locants


def _is_ring_branch_root(graph, root, coming_from, aromatic_atoms, mol):
    """True when `root` starts a ring substituent that is cited as a unit
    rather than walked as chain: a plain ring, or a monocyclic carbocycle
    (benzene included) carrying substituents of its own."""
    if _simple_ring_substituent(graph, root, coming_from, aromatic_atoms, mol=mol) is not None:
        return True
    ring_info = mol.GetRingInfo()
    ring = next((r for r in ring_info.AtomRings() if root in r), None)
    if ring is None:
        return False
    if any(ring_info.NumAtomRings(a) != 1 for a in ring) or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring):
        return True
    atoms = [mol.GetAtomWithIdx(a) for a in ring]
    aromatic = all(a.GetIsAromatic() for a in atoms)
    return (aromatic and len(ring) == 6) or not any(a.GetIsAromatic() for a in atoms)


def _simple_ring_substituent(graph, root, coming_from, aromatic_atoms=frozenset(), mol=None):
    """(ring_size, is_aromatic) when the branch at `root` is a single unsubstituted monocycle attached only at `root`
    (benzene or a one-heteroatom aromatic monocycle count as aromatic), else None."""
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    if len(ring_neighbors) != 2:
        return None
    order = [root]
    visited = {root}
    previous, current = root, ring_neighbors[0]
    while current != root:
        if current in visited:
            return None
        visited.add(current)
        order.append(current)
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) != 1:
            return None
        previous, current = current, neighbors[0]
    if aromatic_atoms and visited <= aromatic_atoms:
        if mol is not None and heteroaromatic_monocycle_name(mol, order) is not None:
            return len(visited), True
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

    def tiebreak_key(entries):
        # P-14.5.2: lowest locant set first (entries are already in
        # ascending position order from `walk`); when two ring-walk
        # directions give the same locant set (a symmetric halogen
        # pattern), the direction that gives the alphabetically-first
        # substituent name the lower locant wins, mirroring
        # `_candidate_key`'s identical `citation_locants` tiebreak
        # elsewhere in this module.
        locant_set = tuple(pos for pos, _, _ in entries)
        grouped = {}
        for pos, name, _ in entries:
            grouped.setdefault(name, []).append(pos)
        citation_locants = tuple(
            loc for name in sorted(grouped, key=alpha_sort_key) for loc in sorted(grouped[name])
        )
        return locant_set, citation_locants

    candidates = [r for r in (walk(n) for n in ring_neighbors) if r is not None and r[1]]
    if not candidates:
        return None
    visited, entries = min(candidates, key=lambda vc: tiebreak_key(vc[1]))
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


def _ring_of_root_is_all_carbon(mol, root):
    ring = next((r for r in mol.GetRingInfo().AtomRings() if root in r), None)
    return ring is None or all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring)


def _ring_has_other_substituents(graph, mol, root, coming_from):
    ring = next(r for r in mol.GetRingInfo().AtomRings() if root in r)
    return any(n not in ring and not (a == root and n == coming_from) for a in ring for n in graph[a])


_ADAMANTANE_SKELETON = Chem.MolFromSmarts("[#6]12[#6][#6]3[#6][#6]([#6][#6]([#6]3)[#6]1)[#6]2")


def _is_adamantane(mol, root):
    """True when the ring system holding `root` is exactly the adamantane skeleton."""
    from ._diester_ring_diyl import _system_of

    _, atoms = _system_of(mol, root)
    return len(atoms) == 10 and any(set(match) == set(atoms) for match in mol.GetSubstructMatches(_ADAMANTANE_SKELETON))


def _ring_system_branch(graph, mol, root, coming_from):
    """Fused rings and non-aromatic heterocycles are named as substituent
    groups by the ring-system machinery (pyrrolidin-1-yl, naphthalen-2-yl, ...);
    None for the carbocycles and aromatic monocycles handled elsewhere."""
    ring_info = mol.GetRingInfo()
    ring = next(r for r in ring_info.AtomRings() if root in r)
    fused = any(ring_info.NumAtomRings(a) != 1 for a in ring)
    hetero = any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring)
    aromatic = all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
    if not (fused or (hetero and not aromatic)):
        return None
    from ._diester_ring_diyl import ring_substituent_name

    if fused:
        assembly = _fused_assembly_branch(mol, graph, root, coming_from)
        if assembly is not None:
            return assembly
    return ring_substituent_name(mol, graph, root, coming_from)


def _fused_assembly_branch(mol, graph, root, coming_from):
    from ._common import halogen_substituents
    from ._polyfunctional import _arm_atoms, assembly_substituent
    from ._system_assembly import _skeleton_key, _systems

    arm = _arm_atoms(graph, root, coming_from)
    systems = [atoms for _, atoms in _systems(mol) if set(atoms) <= arm]
    own = next((set(a) for a in systems if root in a), None)
    if own is None or len(systems) != 2:
        return None
    other = next(set(a) for a in systems if root not in a)
    if _skeleton_key(mol, own) != _skeleton_key(mol, other) or not any(
        mol.GetBondBetweenAtoms(a, b) is not None for a in own for b in other
    ):
        return None
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    return assembly_substituent(mol, graph, root, coming_from, halogen_substituents(mol), aromatic)


def _hetero_ring_branch(mol, root, coming_from):
    """A ring substituent containing a heteroatom: named through the monocycle
    machinery (pyridinyl, furanyl, ...) or rejected rather than misnamed."""
    from ._multiplicative_groups import classify
    from ._multiplicative_ring import ring_substituent_name

    ring_atoms = next(set(r) for r in mol.GetRingInfo().AtomRings() if root in r)
    if mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("a heterocyclic substituent attached by a multiple bond is not supported yet")
    from ._multiplicative import _bare_key

    own_key = _bare_key(mol, ring_atoms)
    for ring in mol.GetRingInfo().AtomRings():
        joined = set(ring) != ring_atoms and any(
            mol.GetBondBetweenAtoms(a, b) is not None for a in ring_atoms for b in ring
        )
        if joined and not set(ring) & ring_atoms and _bare_key(mol, set(ring)) == own_key:
            from ._common import adjacency, halogen_substituents
            from ._polyfunctional import assembly_substituent

            graph = adjacency(mol)
            halogens = halogen_substituents(mol)
            aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
            assembly = assembly_substituent(mol, graph, root, coming_from, halogens, aromatic)
            if assembly is None:
                raise UnsupportedStructure("this heteroaromatic ring assembly as a substituent is not supported yet")
            return assembly
    try:
        result = ring_substituent_name(mol, ring_atoms, root, coming_from, classify(mol) or [], None)
    except UnsupportedStructure:
        result = None
    if result is None:
        from ._common import adjacency
        from ._diester_ring_diyl import ring_substituent_name as general_ring_substituent_name

        return general_ring_substituent_name(mol, adjacency(mol), root, coming_from)
    return result


def name_branch(graph, root, coming_from, halogens=None, aromatic_atoms=None, mol=None, unsaturated=None):
    try:
        return _name_branch(graph, root, coming_from, halogens, aromatic_atoms, mol, unsaturated)
    except UnsupportedStructure:
        if mol is not None and mol.GetRingInfo().NumRings() >= 3:
            from ._phane_general import phane_substituent

            order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
            fused = phane_substituent(mol, graph, root, coming_from, preferred=False) if order in (1.0, 2.0) else None
            if fused is not None:
                return fused
        raise


def _name_branch(graph, root, coming_from, halogens=None, aromatic_atoms=None, mol=None, unsaturated=None):
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
    one into e.g. '4-chlorophenyl' (see `halogenated_phenyl_substituent`).

    `mol`: see `_longest_chains_from_root` -- pass this through whenever
    the branch hasn't already been fully validated as carbon-plus-
    `halogens` by the caller, so an unrecognized heteroatom in it raises
    instead of being silently treated as carbon.

    `unsaturated` (default: whenever `mol` is given): bond orders are read from `mol`, so the
    branch's own C=C/C#C bonds become 'ene'/'yne' endings and the bond
    joining it to `coming_from` becomes 'yl'/'ylidene'/'ylidyne'
    (P-29.2, P-32.1.1, P-46.1); without `mol` every bond reads as single."""
    halogens = halogens or {}
    aromatic_atoms = aromatic_atoms or frozenset()
    if root in halogens:
        return halogens[root], False
    if unsaturated is None:
        unsaturated = mol is not None
    attach_order = 1.0
    if unsaturated:
        attach_order = _bond_order(mol, root, coming_from)
        aromatic_atoms = aromatic_atoms or frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())

    if mol is not None and mol.GetRingInfo().NumRings() >= 3 and attach_order in (1.0, 2.0):
        from ._phane_general import phane_substituent

        phane = phane_substituent(mol, graph, root, coming_from)
        if phane is not None:
            return phane

    if mol is not None and attach_order == 1.0 and mol.GetAtomWithIdx(root).IsInRing():
        from ._glycosyl import glycosyl_branch

        glycosyl = glycosyl_branch(mol, graph, root, coming_from)
        if glycosyl is not None:
            return glycosyl

    if mol is not None and attach_order == 1.0:
        from ._skeletal_group import skeletal_chain_group, skeletal_ring_group

        skeletal = skeletal_chain_group(mol, graph, root, coming_from) or skeletal_ring_group(mol, graph, root, coming_from)
        if skeletal is not None:
            return skeletal

    if mol is not None and unsaturated:
        from ._hetero_prefixes import hetero_branch_name

        hetero = hetero_branch_name(graph, root, coming_from, halogens, aromatic_atoms, mol)
        if hetero is not None:
            return hetero

    if mol is not None and mol.GetAtomWithIdx(root).IsInRing():
        system_name = _ring_system_branch(graph, mol, root, coming_from)
        if system_name is not None:
            return system_name
    if (
        mol is not None
        and mol.GetAtomWithIdx(root).IsInRing()
        and not _ring_of_root_is_all_carbon(mol, root)
        and _ring_has_other_substituents(graph, mol, root, coming_from)
    ):
        return _hetero_ring_branch(mol, root, coming_from)

    ring_result = _simple_ring_substituent(graph, root, coming_from, aromatic_atoms, mol=mol)
    if ring_result is not None:
        ring_size, is_aromatic = ring_result
        if is_aromatic:
            if mol is not None and not _ring_of_root_is_all_carbon(mol, root):
                from ._diester_ring_diyl import ring_substituent_name

                return ring_substituent_name(mol, graph, root, coming_from)
            return "phenyl", False
        if unsaturated:
            return _unsaturated_ring_branch(graph, root, coming_from, ring_size, attach_order, mol)
        return "cyclo" + alkyl_name(ring_size), False

    if aromatic_atoms and (mol is None or _ring_of_root_is_all_carbon(mol, root)):
        halophenyl = halogenated_phenyl_substituent(graph, aromatic_atoms, root, coming_from, halogens)
        if halophenyl is not None:
            name, _, _ = halophenyl
            return name, True

    ring_with_named_atoms = _ring_substituent_with_named_atoms(graph, root, coming_from, halogens)
    if ring_with_named_atoms is not None:
        ring_size, name, locants = ring_with_named_atoms
        loc_str = ",".join(str(loc) for loc in locants)
        name_word = multiplied_word(len(locants), name)
        ring_name = f"{loc_str}-{name_word}cyclo{alkyl_name(ring_size)}"
        return _free_valence_suffix(ring_name, attach_order), True

    if unsaturated and mol.GetAtomWithIdx(root).IsInRing():
        return _substituted_ring_branch(graph, root, coming_from, halogens, aromatic_atoms, mol, attach_order)

    _, _, name, is_compound = _select_winning_structure(
        graph, root, coming_from, halogens, mol, aromatic_atoms, unsaturated
    )
    return name, is_compound


def substituents_for_chain(
    graph, chain, halogens, excluded=frozenset(), mol=None, aromatic_atoms=frozenset(), unsaturated=None
):
    """{position (1-based) -> [(name, is_compound), ...]} for every branch
    hanging off a candidate principal chain -- shared by every module whose
    parent hydride is an acyclic chain (P-29.2). `excluded`: atom indices to
    skip besides the chain itself, e.g. a functional-group atom the caller
    already accounts for separately (an amine nitrogen, a carbonyl oxygen,
    ...) so it's never mistaken for a substituent branch. `aromatic_atoms`:
    passed through to `name_branch` for a branch rooted on a fused aromatic
    ring."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [
                name_branch(graph, root, atom, halogens, aromatic_atoms, mol=mol, unsaturated=unsaturated)
                for root in branch_roots
            ]
    return substituents


def substituents_for_chain_forced_compound_terminals(graph, chain, terminals, mol=None):
    """Like `substituents_for_chain`, but a `terminals` leaf is always cited
    as a compound (parenthesized) substituent instead of leaving that
    judgment to `name_branch` -- shared by the dichalcogenide modules
    (`_disulfide.py`/`_diselenide.py`/`_ditelluride.py`) whose
    disulfanyl/diselanyl/ditellanyl prefix needs this, unlike the plain
    `terminals` shortcut `_sulfide.py`/`_selenide.py`/`_peroxide.py` share
    via `_acyclic.name_from_carbon_graph`."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        entries = []
        for root in branch_roots:
            if root in terminals:
                entries.append((terminals[root], True))
            else:
                entries.append(name_branch(graph, root, atom, {}, mol=mol))
        substituents[position] = entries
    return substituents


def substituents_for_ring(
    graph, ring_order, halogens, excluded=frozenset(), mol=None, aromatic_atoms=frozenset(), unsaturated=None
):
    """{position (1-based) -> [(name, is_compound), ...]} for every branch
    hanging off a candidate ring numbering -- the ring-parent analogue of
    `substituents_for_chain` (P-29.2). `excluded`: atom indices to skip
    besides the ring itself, e.g. a functional-group atom the caller
    already accounts for separately. `aromatic_atoms`: passed through to
    `name_branch` for a branch rooted on a fused aromatic ring."""
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in excluded]
        if branch_roots:
            substituents[position] = [
                name_branch(graph, root, atom, halogens, aromatic_atoms, mol=mol, unsaturated=unsaturated)
                for root in branch_roots
            ]
    return substituents


def _substituent_entries_along_chain(
    graph, chain, first_previous, halogens, extra_exclusions=None, mol=None, aromatic_atoms=frozenset()
):
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
            sub_name, sub_compound = name_branch(graph, branch_root, atom, halogens, aromatic_atoms, mol=mol)
            entries.append((position, sub_name, sub_compound))
    return entries


def _select_winning_chain(graph, root, coming_from, halogens, mol=None, aromatic_atoms=frozenset()):
    """The ordinary (non-branch-point) P-46 tie-break: `root` is fixed at
    locant 1 (the free valence is always a chain terminus here). Returns
    (chain, root_position, name, is_compound); `root_position` is always 1,
    kept in the return shape so callers can treat this and
    `_branch_point_candidate_chains` interchangeably via
    `_select_winning_structure`."""
    chains = _longest_chains_from_root(graph, root, coming_from, halogens, mol, aromatic_atoms)
    chain_length = len(chains[0])

    best_key = None
    best_chain = None
    best_name = None
    best_compound = None
    for chain in chains:
        entries = _substituent_entries_along_chain(graph, chain, coming_from, halogens, mol=mol, aromatic_atoms=aromatic_atoms)
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


def _branch_point_candidate_chains(graph, root, coming_from, halogens, mol=None, aromatic_atoms=frozenset()):
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

    paths = {b: _longest_chains_from_root(graph, b, root, halogens, mol, aromatic_atoms) for b in branch_roots}
    lengths = {b: len(paths[b][0]) for b in branch_roots}

    if (
        len(branch_roots) == 3
        and all(length == 1 for length in lengths.values())
        and all(len(graph[b]) == 1 for b in branch_roots)
    ):
        # P-29.6.1: the retained name 'tert-butyl' is the PIN for the
        # unsubstituted (CH3)3C- group, never the general rule's own
        # '2-methylpropan-2-yl' -- and, being a single retained word with no
        # locant of its own, it is never parenthesized as a compound prefix.
        # The `len(graph[b]) == 1` check is required alongside the chain-
        # length-1 check above: a branch root with a halogen substituent of
        # its own (e.g. -CH2Br) also has a length-1 chain (halogens don't
        # extend the chain search either), but isn't a bare methyl -- taking
        # this shortcut for it would silently drop the halogen (found via
        # real-data testing: 'CC(C)c1cccc(C(C)(C)CBr)c1' was misnamed
        # '1-tert-butyl-3-(propan-2-yl)benzene', losing the bromine
        # entirely; the correct name keeps it as a substituted
        # '1-bromo-2-methylpropan-2-yl' branch instead).
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
        for path_before in paths[before_branch]:
            for path_after in paths[after_branch]:
                spine = list(reversed(path_before)) + [root] + path_after
                root_locant = len(path_before) + 1
                # Any branch off `root` beyond `before_branch`/`after_branch`
                # (e.g. a third, shorter branch when `root` forks three
                # ways) is already picked up below: `root` always sits at
                # an internal spine locant (`root_locant` >= 2, since
                # `path_before` is never empty), so
                # `_substituent_entries_along_chain`'s own walk over that
                # position's neighbors finds it via the ordinary
                # not-in-chain-set/not-`previous`/not-`excluded` check --
                # a separate pass over it here would double-count it.
                entries = _substituent_entries_along_chain(
                    graph, spine, None, halogens, extra_exclusions={root_locant: coming_from}, mol=mol,
                    aromatic_atoms=aromatic_atoms,
                )
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


_FREE_VALENCE_SUFFIX = {float(order): suffix for order, suffix in SUFFIX_OF_ORDER.items()}


def _bond_order(mol, a, b):
    order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
    if order not in _FREE_VALENCE_SUFFIX and order != 1.5:
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not supported in a substituent group"
        )
    return order


def _free_valence_suffix(alkyl_style_name, order):
    """Swap a trailing 'yl' for 'ylidene'/'ylidyne' per the bond order joining
    the group to its parent (P-29.2)."""
    suffix = _FREE_VALENCE_SUFFIX[order]
    return alkyl_style_name[: -len("yl")] + suffix


def _ring_walk_from(graph, root, coming_from):
    ring_neighbors = [n for n in graph[root] if n != coming_from]
    order = [root]
    previous, current = root, ring_neighbors[0]
    while current != root:
        order.append(current)
        (following,) = [n for n in graph[current] if n != previous]
        previous, current = current, following
    return order


def _ring_base_name(ring_size, ene, yne, attach_order):
    if not ene and not yne:
        return _free_valence_suffix("cyclo" + alkyl_name(ring_size), attach_order)
    body, needs_a = unsaturation_suffix(ene, yne)
    stem = "cyclo" + alkane_name(ring_size)[:-3] + ("a" if needs_a else "")
    return f"{stem}-{body[:-1]}-1-{_FREE_VALENCE_SUFFIX[attach_order]}"


def _ring_multiple_bond_locants(mol, direction):
    ene, yne = [], []
    for i, atom in enumerate(direction):
        order = _bond_order(mol, atom, direction[(i + 1) % len(direction)])
        if order == 1.0:
            continue
        if order not in (2.0, 3.0):
            raise UnsupportedStructure("an unsupported bond order inside a ring substituent")
        (ene if order == 2.0 else yne).append(i + 1)
    return ene, yne


def _unsaturated_ring_branch(graph, root, coming_from, ring_size, attach_order, mol):
    """Name a plain monocyclic ring substituent attached at `root` with the
    free valence fixed at locant 1 (P-32.1.2): 'cyclohexyl',
    'cyclohexylidene', 'cyclohex-2-en-1-yl', ..."""
    ring = _ring_walk_from(graph, root, coming_from)
    best = None
    for direction in (ring, [ring[0]] + ring[:0:-1]):
        ene, yne = _ring_multiple_bond_locants(mol, direction)
        candidate = (sorted(ene + yne), sorted(ene), ene, yne)
        if best is None or candidate < best:
            best = candidate
    _, _, ene, yne = best
    return _ring_base_name(ring_size, ene, yne, attach_order), bool(ene or yne)


BRANCH_STEREO = contextvars.ContextVar("branch_stereo", default=None)


def _branch_stereo_entries(positions, ring=False, record=False):
    """[(locant, code)] for the stereo elements of the enclosing molecule that
    lie on a substituent chain or ring numbered by `positions`
    ({atom: locant}); set by the polyfunctional engine through
    `BRANCH_STEREO` (atoms: {idx: code}, bonds: {(a, b): code}, used: set)."""
    context = BRANCH_STEREO.get()
    if not context:
        return []
    entries = []
    for atom, code in context["atoms"].items():
        if atom in positions:
            entries.append((positions[atom], code))
            if record:
                context["used"].add(("atom", atom))
    if not ring:
        for (a, b), code in context["bonds"].items():
            if a in positions and b in positions and abs(positions[a] - positions[b]) == 1:
                entries.append((min(positions[a], positions[b]), code))
                if record:
                    context["used"].add(("bond", (a, b)))
    return sorted(entries)


def _branch_stereo_rank(entries):
    return tuple(0 if code in "RZr" else 1 for _, code in entries)


def _branch_stereo_prefix(entries):
    return "(" + ",".join(f"{locant}{code}" for locant, code in entries) + ")-" if entries else ""


def _substituted_ring_branch(graph, root, coming_from, halogens, aromatic_atoms, mol, attach_order):
    """A monocyclic carbocyclic (or benzene) substituent group that carries
    substituents of its own (P-29.3.3, P-32.1.2): the free valence is locant
    1, then the ring's multiple bonds and the prefixes take lowest locants,
    e.g. '2-methylcyclohexyl', '2-methylidenecyclohexyl', '4-methylphenyl'."""
    ring_info = mol.GetRingInfo()
    ring_atoms = next((r for r in ring_info.AtomRings() if root in r), None)
    cyclic_error = UnsupportedStructure(
        "cyclic substituent groups are not supported yet (see P-29.3.3, P-46 for cyclic substituent groups)"
    )
    if ring_atoms is None or any(ring_info.NumAtomRings(a) != 1 for a in ring_atoms):
        raise cyclic_error
    from ._multiplicative import _bare_key

    own_key = _bare_key(mol, set(ring_atoms))
    for atom in ring_atoms:
        for neighbor in graph[atom]:
            other = next((r for r in ring_info.AtomRings() if neighbor in r and atom not in r), None)
            if (
                other is not None
                and neighbor not in ring_atoms
                and all(ring_info.NumAtomRings(a) == 1 for a in other)
                and _bare_key(mol, set(other)) == own_key
            ):
                from ._polyfunctional import assembly_substituent

                assembly = assembly_substituent(mol, graph, root, coming_from, halogens, aromatic_atoms)
                if assembly is None:
                    raise UnsupportedStructure(
                        "this ring assembly as a substituent group is not supported yet (P-28)"
                    )
                return assembly
    atoms = [mol.GetAtomWithIdx(a) for a in ring_atoms]
    aromatic = all(a.GetIsAromatic() for a in atoms)
    if any(a.GetAtomicNum() != 6 for a in atoms):
        from ._diester_ring_diyl import ring_substituent_name

        return ring_substituent_name(mol, graph, root, coming_from)
    if not aromatic and any(a.GetIsAromatic() for a in atoms):
        raise cyclic_error
    if aromatic and len(ring_atoms) != 6:
        raise cyclic_error
    order = ring_cycle(graph, list(ring_atoms))
    order = order[order.index(root):] + order[: order.index(root)]
    ring_set = set(ring_atoms)

    best = None
    for direction in (order, [order[0]] + order[:0:-1]):
        ene, yne = ([], []) if aromatic else _ring_multiple_bond_locants(mol, direction)
        entries = []
        for position, atom in enumerate(direction, start=1):
            for neighbor in graph[atom]:
                if neighbor in ring_set or (atom == root and neighbor == coming_from):
                    continue
                sub_name, sub_compound = name_branch(
                    graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True
                )
                entries.append((position, sub_name, sub_compound))
        grouped = _group_substituents(entries)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        prefix = format_substituent_prefixes(grouped)
        base = "phenyl" if aromatic else _ring_base_name(len(order), ene, yne, attach_order)
        name = prefix + base
        stereo_entries = _branch_stereo_entries({atom: i for i, atom in enumerate(direction, start=1)}, ring=True)
        key = (sorted(ene + yne), sorted(ene), locant_set, citation, _branch_stereo_rank(stereo_entries), name)
        if best is None or key < best[0]:
            best = (key, name, bool(prefix) or bool(ene or yne), direction)
    positions = {atom: i for i, atom in enumerate(best[3], start=1)}
    stereo_prefix = _branch_stereo_prefix(_branch_stereo_entries(positions, ring=True, record=True))
    return stereo_prefix + best[1], best[2] or bool(stereo_prefix)


def _chain_children(graph, node, parent, halogens, mol, aromatic_atoms):
    """Neighbors of `node` the substituent's own principal chain may extend
    into: every carbon except atoms named as one-atom prefixes (`halogens`)
    and the root of a plain ring substituent, which are cited as
    substituents of the chain instead."""
    if mol.GetAtomWithIdx(node).GetAtomicNum() != 6:
        raise UnsupportedStructure(
            "a heteroatom in a compound substituent branch, other than "
            "a recognized halogen/named group, is not supported yet"
        )
    from ._hetero_prefixes import is_functional_carbon

    children = []
    for n in graph[node]:
        if n == parent or n in halogens:
            continue
        if mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or is_functional_carbon(mol, n):
            continue
        if (
            mol.GetAtomWithIdx(n).IsInRing()
            and _is_ring_branch_root(graph, n, node, aromatic_atoms, mol)
        ):
            continue
        children.append(n)
    return children


def _longest_arms(graph, node, parent, halogens, mol, aromatic_atoms, visited):
    """Every longest simple path starting at `node` and descending away from
    `parent`."""
    if node in visited:
        raise UnsupportedStructure(
            "cyclic substituent groups are not supported yet (see P-29.3.3, P-46 for cyclic substituent groups)"
        )
    visited = visited | {node}
    best = [[node]]
    for child in _chain_children(graph, node, parent, halogens, mol, aromatic_atoms):
        for arm in _longest_arms(graph, child, node, halogens, mol, aromatic_atoms, visited):
            candidate = [node] + arm
            if len(candidate) > len(best[0]):
                best = [candidate]
            elif len(candidate) == len(best[0]) and len(candidate) > 1:
                best.append(candidate)
    return best


def _unsaturated_candidate_chains(graph, root, coming_from, halogens, mol, aromatic_atoms):
    """Every longest chain through `root` (P-46.1(b)): root-terminated when it
    has a single way onward, running through it along its two longest arms
    otherwise (P-29.3.2.2)."""
    children = _chain_children(graph, root, coming_from, halogens, mol, aromatic_atoms)
    visited = frozenset({root})
    if len(children) <= 1:
        return [list(arm) for arm in _longest_arms(graph, root, coming_from, halogens, mol, aromatic_atoms, frozenset())]
    arms = {c: _longest_arms(graph, c, root, halogens, mol, aromatic_atoms, visited) for c in children}
    total = sorted((len(arms[c][0]) for c in children), reverse=True)
    target = total[0] + total[1]
    chains = []
    for a in children:
        for b in children:
            if a != b and len(arms[a][0]) + len(arms[b][0]) == target:
                for before in arms[a]:
                    for after in arms[b]:
                        chains.append(list(reversed(before)) + [root] + after)
    return chains


def _select_unsaturated_structure(graph, root, coming_from, halogens, mol, aromatic_atoms):
    """`_select_winning_structure`'s bond-order-aware counterpart: the
    principal chain of a substituent group whose chain or free valence
    carries multiple bonds (P-46.1: longest chain, then most multiple
    bonds, then most double bonds, then lowest free-valence, multiple-bond,
    and substituent locants). Returns (chain, root_position, name,
    is_compound)."""
    attach_order = _bond_order(mol, root, coming_from)
    if attach_order not in _FREE_VALENCE_SUFFIX:
        raise UnsupportedStructure("an aromatic bond joining a chain substituent to its parent is not supported")
    suffix = _FREE_VALENCE_SUFFIX[attach_order]

    best = None
    for chain in _unsaturated_candidate_chains(graph, root, coming_from, halogens, mol, aromatic_atoms):
        chain_set = set(chain)
        root_position = chain.index(root) + 1
        ene, yne = [], []
        for i in range(len(chain) - 1):
            order = _bond_order(mol, chain[i], chain[i + 1])
            if order == 2.0:
                ene.append(i + 1)
            elif order == 3.0:
                yne.append(i + 1)
        entries = []
        for position, atom in enumerate(chain, start=1):
            for neighbor in graph[atom]:
                if neighbor in chain_set or (atom == root and neighbor == coming_from):
                    continue
                sub_name, sub_compound = name_branch(
                    graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True
                )
                entries.append((position, sub_name, sub_compound))
        grouped = _group_substituents(entries)
        locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
        name, is_compound = _unsaturated_chain_name(
            len(chain), root_position, suffix, ene, yne, grouped, tert_butyl=_is_tert_butyl(graph, root, coming_from, halogens)
        )
        multiple = sorted(ene + yne)
        stereo_rank = _branch_stereo_rank(_branch_stereo_entries({a: i for i, a in enumerate(chain, start=1)}))
        key = (
            -len(multiple), -len(ene), root_position, multiple, sorted(ene), -total_count, locant_set, citation,
            stereo_rank, name,
        )
        if best is None or key < best[0]:
            best = (key, chain, root_position, name, is_compound)
    _, chain, root_position, name, is_compound = best
    stereo_prefix = _branch_stereo_prefix(
        _branch_stereo_entries({a: i for i, a in enumerate(chain, start=1)}, record=True)
    )
    return chain, root_position, stereo_prefix + name, is_compound or bool(stereo_prefix)


def _is_tert_butyl(graph, root, coming_from, halogens):
    # P-29.6.1: the unsubstituted (CH3)3C- group keeps its retained name.
    others = [n for n in graph[root] if n != coming_from]
    return len(others) == 3 and all(len(graph[n]) == 1 and n not in halogens for n in others)


def _unsaturated_chain_name(length, root_position, suffix, ene, yne, grouped, tert_butyl):
    if tert_butyl and suffix == "yl" and not ene and not yne:
        return "tert-butyl", False
    if length == 1 and list(grouped) == ["phenyl"] and len(grouped["phenyl"]["locants"]) == 1:
        # P-29.6.2.1: 'benzyl'/'benzylidene'/'benzylidyne' are the preferred prefixes when unsubstituted.
        return "benz" + suffix, False
    prefix = ""
    if grouped and length == 1 and len(grouped) >= 2:
        flat = [(name, info["compound"]) for name, info in grouped.items() for _ in info["locants"]]
        try:
            prefix = format_mononuclear_prefixes(flat)
        except UnsupportedStructure:
            prefix = ""
    if not prefix:
        prefix = format_substituent_prefixes(grouped, omit_locants=(length == 1)) if grouped else ""
    multiple = len(ene) + len(yne)
    if not multiple:
        if root_position == 1:
            base = alkyl_name(length)[: -len("yl")] + suffix
        else:
            base = f"{alkane_name(length)[:-1]}-{root_position}-{suffix}"
    else:
        stem = alkane_name(length)[:-3]
        if length == 2 and multiple == 1:
            bond = stem + ("en" if ene else "yn")
            base = f"{bond}-{root_position}-{suffix}" if grouped else bond + suffix
        else:
            body, needs_a = unsaturation_suffix(ene, yne)
            base = f"{stem}{'a' if needs_a else ''}-{body[:-1]}-{root_position}-{suffix}"
    # P-14.3.4.4: only the unsubstituted two-carbon group is unambiguous without locants.
    return prefix + base, bool(prefix) or "-" in base


def _select_winning_structure(
    graph, root, coming_from, halogens, mol=None, aromatic_atoms=frozenset(), unsaturated=None
):
    """(chain, root_position, name, is_compound) for the substituent group
    hanging off `root`: the P-29.3.2.2 branch-point chain
    (`_branch_point_candidate_chains`) when `root` forks into two or more
    branches, else the ordinary root-is-locant-1 chain
    (`_select_winning_chain`). Shared by `name_branch` (which only needs
    the name) and `branch_atom_locant` below (which also needs to know
    which chain won and where `root` sits on it), so the two can never
    disagree about which chain was chosen."""
    if unsaturated:
        return _select_unsaturated_structure(graph, root, coming_from, halogens, mol, aromatic_atoms)
    branch_point = _branch_point_candidate_chains(graph, root, coming_from, halogens, mol, aromatic_atoms)
    if branch_point is not None:
        return branch_point
    return _select_winning_chain(graph, root, coming_from, halogens, mol, aromatic_atoms)


def branch_atom_locant(graph, root, coming_from, atom_idx, halogens=None, mol=None):
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
    chain, _, _, _ = _select_winning_structure(
        graph, root, coming_from, halogens, mol, unsaturated=mol is not None
    )
    if atom_idx not in chain:
        raise UnsupportedStructure(
            "a specified stereocenter that isn't on the substituent's own "
            "principal chain (P-46) is not supported yet"
        )
    return chain.index(atom_idx) + 1


def _ring_branch_stereo_core(graph, ring_order, group_locants, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    """Shared lookup behind `ring_branch_stereo_display` and
    `ring_and_branch_stereo_display` (#776) -- see those two for the
    exact shapes each one exposes. Returns `(ring_atom, display,
    ring_r_or_s)`: `ring_atom`/`display` are the branch's ring-attachment
    atom and its bracketed "[(<locant><R/S>)-<name>]" descriptor (P-91.3)
    exactly as before; `ring_r_or_s` is `None` when the ring atom itself
    isn't a stereocenter, or its own "R"/"S" code when it is (the ring
    atom being additionally a stereocenter is only possible for a ring
    whose own symmetry is already broken -- e.g. by a ring double bond --
    since a plain saturated ring's ring-walk symmetry from a
    singly-substituted atom keeps it from ever being one, see
    `_cyclic_unsaturated.py`). Returns `None` for zero, one-not-matching-
    the-branch, or more than two stereocenters, more than one off-ring
    stereocenter, an off-ring stereocenter not on the ring's sole
    substituent branch, or the ring having other than exactly one
    substituent in total.

    `group_locants`: the calling module's own characteristic-group atom
    indices on the ring (e.g. hydroxyls, amines, ketones) -- excluded
    from the branch-attachment search the same way `ring_order` itself
    is, so the group's own atom is never mistaken for a substituent
    branch."""
    if len(stereo) not in (1, 2):
        return None
    ring_set = set(ring_order)
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_order
        for neighbor in graph[ring_atom]
        if neighbor not in ring_set and neighbor not in group_locants
    ]
    if len(branch_attachments) != 1:
        return None
    ring_atom, branch_root = branch_attachments[0]
    remaining = dict(stereo)
    ring_r_or_s = remaining.pop(ring_atom, None)
    if len(remaining) != 1:
        return None
    stereo_atom, r_or_s = next(iter(remaining.items()))
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens, aromatic_atoms, mol=mol)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens, mol=mol)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    display = f"[{descriptor}]" if branch_compound else f"({descriptor})"
    return ring_atom, display, ring_r_or_s


def ring_branch_stereo_display(graph, ring_order, group_locants, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    """If the ring carries exactly one specified stereocenter and that
    stereocenter sits off the ring on the ring's own sole substituent
    branch (P-92), return that branch's ring-attachment atom plus its
    bracketed "[(<locant><R/S>)-<name>]" display (P-91.3) -- e.g. the '1'
    atom and '[(2S)-butan-2-yl]' in '1-[(2S)-butan-2-yl]cyclohexan-1-ol'
    (PubChem CID confirmed). Returns None (the caller keeps its existing
    outright rejection) for more than one stereocenter, or the ring
    having more or fewer than one substituent in total -- both a
    genuinely more general case this narrow slice doesn't attempt.

    `group_locants`: the calling module's own characteristic-group atom
    indices on the ring (e.g. hydroxyls, amines, ketones) -- excluded
    from the branch-attachment search the same way `ring_order` itself
    is, so the group's own atom is never mistaken for a substituent
    branch."""
    if len(stereo) != 1:
        return None
    result = _ring_branch_stereo_core(graph, ring_order, group_locants, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms)
    if result is None:
        return None
    ring_atom, display, _ring_r_or_s = result
    return ring_atom, display


def ring_and_branch_stereo_display(graph, ring_order, group_locants, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    """Like `ring_branch_stereo_display`, but also accepts the ring's
    branch-attachment atom itself being a second, independent specified
    stereocenter alongside the one on its substituent branch (#776) --
    confirmed real via PubChem (2026-09-21):
    'CC[C@H](C)[C@H]1CCC=CC1' -> '(4S)-4-[(2S)-butan-2-yl]cyclohexene'
    (CID 175915978). The front "(4S)-" is the ring atom's own ordinary
    on-ring locant+R/S (P-91.3, the same mechanism used for an
    all-on-ring stereocenter set); the bracketed "[(2S)-butan-2-yl]" is
    the branch's own descriptor, built exactly as
    `ring_branch_stereo_display` already does.

    Returns `(ring_atom, display, ring_r_or_s)` -- see
    `_ring_branch_stereo_core` for the exact shape and rejection
    conditions."""
    return _ring_branch_stereo_core(graph, ring_order, group_locants, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms)
