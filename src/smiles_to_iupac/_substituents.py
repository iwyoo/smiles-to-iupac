"""Compound substituent groups (P-29.4, P-46, P-14.5.2, P-35.2.1): the free valence is locant 1 of the principal
chain (longest, then most multiple bonds, then lowest locants); the free-valence bond selects yl, ylidene or ylidyne
(P-29.2). Ring and ring-system roots are named by the one general ring-group namer
(`_diester_ring_diyl.ring_substituent_name`, P-29.3.3, P-29.3.4)."""

import contextlib
import contextvars
import re

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._multiplicative_text import enclose, nesting_text
from ._free_valence import SUFFIX_OF_ORDER
from ._common import (
    UnsupportedStructure,
    alpha_sort_key,
    citation_order_key,
    substituent_locant_set_and_citation,
    unsaturation_suffix,
)
from ._locant_omission import omits_all_locants
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
    if not isinstance(locant, str):
        return (1, 0, locant, 0, "")
    digits = re.match(r"\d+", locant)
    if digits is None:
        lettered = re.fullmatch(r"([A-Z][a-z]?)(['′]*)(\d*)", locant)
        if lettered:
            return (0, int(lettered.group(3) or 0), lettered.group(1), len(lettered.group(2)), locant)
        return (0, 0, "", 0, locant)
    return (1, locant.count("\u2032"), int(digits.group()), 0, locant[digits.end():].replace("\u2032", ""))


_PLAIN_STEM_PREFIX = None
# an unsubstituted ring group carries no prefix before its stem, only its own attachment locant (P-16.3.4(a))
_RING_GROUP_PREFIX = re.compile(r"^(?:\d+H-)?[a-z]+-\d+[a-z]?-(?:yl|ylidene|ylidyne)$")


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
    return bool(_PLAIN_STEM_PREFIX.match(name) or _RING_GROUP_PREFIX.match(name))


_CHALCOGEN_HYDRIDE_GROUP = re.compile(r"^(?:di|tri|tetra|penta|hexa)?(?:sulfanyl|selanyl|tellanyl)(?:idene)?$")
_LEADING_NUMERAL = re.compile(r"^(?:di|do|tri|tetra|penta|hexa|hepta|octa|nona|dec)[a-z]*(?:yl|oyl)$")
_PLAIN_CARBONYL = re.compile(r"^(?:cyclo)?[a-z]+(?:ane|an|e)?(?:-[\d,]+)?-?carbon(?:yl|othioyl)$")
_HYDRIDE_ACYL = re.compile(r"^[a-z]{3,}(?:ane|ene)(?:sulfonyl|sulfinyl)$")


_POLYCYCLE_GROUP = re.compile(r"^(?:bi|tri|tetra)?cyclo\[[\d.,^]+\][a-z]+(?:-[\d,]+-(?:en|yn))?-[\d]+[a-z]?-(?:yl|ylidene|ylidyne)$")


def _fully_enclosed(name: str) -> bool:
    """Whether the first enclosing mark of `name` closes at its last character: '(dimethylamino)'."""
    depth = 0
    for index, ch in enumerate(name):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
            if depth == 0:
                return index == len(name) - 1
    return False


# P-16.3.6: 'bis(diazenyl)', not 'di(diazenyl)': the group name itself begins with a multiplying term
_MULTIPLIED_HETERO_STEM = re.compile(r"(?:di|tri|tetra|penta|hexa)(?:az|sil|sulf|selan|tellan|phosph|ars|oxid|germ|stann|plumb|bor)[a-z]*yl")


def prefix_multiplier(count: int, name: str, compound: bool):
    """(multiplier, enclosed) for `name` cited `count` > 1 times as a detachable prefix (P-16.3.3 to P-16.3.6): 'di'
    for simple names; 'di' with enclosing marks for simple names with locants, brackets or a leading numerical term;
    'bis' with marks for substituted names and for the mononuclear groups of a polynuclear chain; 'di-' for tert-butyl."""
    if _isotope_only_prefix(name):
        return multiplying_prefix(count), True
    if name[:1] in "([{" and not name.startswith("(\u03b7") and _fully_enclosed(name):
        return multiplying_prefix(count, compound=True), False
    if (name[:1] in "([{" and not name.startswith("(\u03b7")) or _CHALCOGEN_HYDRIDE_GROUP.match(name):
        return multiplying_prefix(count, compound=True), True
    if name == "tert-butyl":
        return f"{multiplying_prefix(count)}-", False
    if _MULTIPLIED_HETERO_STEM.match(name):
        return multiplying_prefix(count, compound=True), True
    simple = not compound or is_plain_stem_prefix(name) or bool(
        _HYDRIDE_ACYL.match(name) or _POLYCYCLE_GROUP.match(name) or _PLAIN_CARBONYL.match(name)
    )
    if not simple:
        return multiplying_prefix(count, compound=True), True
    enclosed = compound or "[" in name or bool(_LEADING_NUMERAL.match(name) or _HYDRIDE_ACYL.match(name))
    return multiplying_prefix(count), enclosed


def _isotope_only_prefix(name):
    """A plain alkyl or alkoxy group with only its isotopic descriptor in front: cited alone it needs no enclosing
    marks (P-82.2.1)."""
    match = re.fullmatch(r"\(\d[^()]*\)([a-z]+)", name)
    if not match:
        return False
    stem = match.group(1)
    return is_plain_stem_prefix(stem[:-3] + "yl" if stem.endswith("oxy") else stem)


def multiplied_prefix(count: int, name: str, compound: bool) -> str:
    """`name` cited `count` times as a detachable prefix, without locants."""
    if count == 1:
        return enclose(name) if compound else name
    multiplier, enclosed = prefix_multiplier(count, name, compound)
    return multiplier + (enclose(name) if enclosed else name)


def format_substituent_prefixes(grouped, omit_locants: bool = False, omit_all: bool = False) -> str:
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
    if omit_all:
        from ._retained_acids import single_site_prefixes

        return single_site_prefixes(grouped)
    if (
        omit_locants
        and len(grouped) > 1
        and not any("multiplier" in info or any(isinstance(loc, str) for loc in info["locants"]) for info in grouped.values())
    ):
        # P-16.5.1.3.1: with every locant omitted, the second and later prefixes are enclosed
        return format_mononuclear_prefixes(
            [(name, info["compound"]) for name, info in grouped.items() for _ in info["locants"]]
        )
    entries = []
    for name in sorted(grouped, key=citation_order_key):
        info = grouped[name]
        locants = sorted(info["locants"], key=_locant_sort_key)
        if "multiplier" in info:
            display_name = enclose(name) if info["compound"] else name
            body = f"{info['multiplier']}{display_name}"
        else:
            lone = omit_locants and len(grouped) == 1 and len(locants) == 1 and _isotope_only_prefix(name)
            body = multiplied_prefix(len(locants), name, info["compound"] and not lone)
        explicit = (not omit_locants) or any(isinstance(loc, str) for loc in locants)
        if explicit:
            loc_str = ",".join(str(loc) for loc in locants)
            entries.append((f"{loc_str}-{body}", True))
        else:
            entries.append((body, False))

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
    nested = nesting_text(name)
    if "{" in nested:
        return f"({name})"
    if "[" in nested:
        return "{" + name + "}"
    if "(" in nested:
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
    'ethyldi(methyl)phosphane (PIN)' worked example (the Blue Book).

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
    'ethyldi(propan-2-yl)silane (PIN)' worked example, the Blue Book P-16.5.1.3.1), mirroring the 'di(...)' rule already
    established for `_carbamate.py`/`_urea.py`. A multiplied compound
    name mixed with a different substituent keeps its multiplying prefix
    ('bis') outside its own enclosing marks, first position included."""
    counts = {}
    compound_of = {}
    # P-68.3.1.1.1.5: 'aminooxy' is a preselected simple prefix, cited first without enclosing marks
    for name, is_compound in [(n, c and n != "aminooxy") for n, c in entries]:
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
            if compound_of[name] and (name[0] in "0123456789([{" or "-" not in name):
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
        return multiplied_prefix(count, name, compound_of[name])

    ordered = sorted(counts, key=alpha_sort_key)
    parts = []
    for i, name in enumerate(ordered):
        count = counts[name]
        if count > 1:
            # P-16.5.1.3.1's own text: "the multiplicative prefixes are
            # not included in the parentheses" -- confirmed via the Blue
            # Book's own 'ethyldi(methyl)phosphane (PIN)' worked example
            # (the Blue Book), so the prefix sits outside the
            # parens at any position, not just the first.
            needs_kis = compound_of[name] and not is_plain_stem_prefix(name)
            parts.append(multiplying_prefix(count, compound=needs_kis) + (name if i == 0 and not compound_of[name] else wrap_marks(name)))
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


class CompoundPrefix(str):
    """A terminal prefix name (e.g. 'methylsulfanyl') that is itself a compound prefix (P-16.5.1.1), so it is
    cited in enclosing marks when passed through `halogens`."""


ISOTOPE_LABELS = contextvars.ContextVar("isotope_labels", default=None)
FORCED_BRANCH_NAMES = contextvars.ContextVar("forced_branch_names", default=None)


def _label_branch(result, graph, root, coming_from, halogens=None, mol=None, aromatic_atoms=frozenset(), unsaturated=None):
    """`result` with the isotopic descriptor of the labelled atoms of the substituent: inserted before the stem of its
    principal chain with that chain's locants (P-82.2.1, P-82.6.1), before 'phenyl' for a phenyl group, or before
    the 'oxy' of a methoxy group."""
    from ._isotope_labels import descriptor

    context = ISOTOPE_LABELS.get()
    atoms = {root}
    stack = [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != coming_from and n not in atoms:
                atoms.add(n)
                stack.append(n)
    labelled = {a: context["labels"][a] for a in atoms if a in context["labels"] and a not in context["consumed"]}
    if not labelled:
        return result
    name, compound = result
    capacity = _hydrogen_capacity(context, labelled)
    if name == "methoxy":
        text = descriptor(labelled, {a: 1 for a in labelled}, True, capacity=capacity)
        context["consumed"].update(labelled)
        return text + name, False
    if name == "carboxy" and mol is not None and mol.GetAtomWithIdx(root).GetAtomicNum() == 6:
        from ._polyfunctional import _modifications
        from ._isotope_labels import _nuclide_sort_key

        nuclides = sorted((n for e in labelled.values() for n in _modifications(e)), key=_nuclide_sort_key)
        context["consumed"].update(labelled)
        return "(" + ",".join(nuclides) + ")carboxy", False
    if name == "benzyl" and set(labelled) == {root}:
        context["consumed"].add(root)
        return "phenyl" + descriptor(labelled, {root: 1}, True, capacity=capacity) + "methyl", True
    run = _chalcogen_run(mol, graph, root, coming_from) if mol is not None else None
    if run is not None and all(a in run for a in labelled):
        stem = _CHALCOGEN_RUN_STEM.search(name)
        if stem is not None:
            context["consumed"].update(labelled)
            text = descriptor(labelled, {a: i + 1 for i, a in enumerate(run)}, False, capacity=capacity)
            return name[: stem.start()] + text + name[stem.start() :], compound
    positions = _plain_chain_positions(graph, root, coming_from, atoms, mol)
    if positions is not None and name.isalpha():
        context["consumed"].update(labelled)
        omit = len(atoms) == 1 or _fully_modified(labelled, positions, capacity, atoms)
        return descriptor(labelled, positions, omit, capacity=capacity) + name, False
    if mol is None:
        raise UnsupportedStructure("an isotopically modified substituent other than a plain alkyl group is not supported yet")
    if mol.GetAtomWithIdx(root).IsInRing():
        positions = _phenyl_label_positions(mol, root, labelled, name, coming_from)
        if positions is None:
            raise UnsupportedStructure("an isotopically modified ring substituent other than phenyl or cycloalkyl is not supported yet")
        stem_index = name.rfind(_ring_group_stem(mol, root)[0])
        omit = _fully_modified(labelled, positions, capacity, {a for a in positions if a != root})
    else:
        named = {**{a: "x" for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6}, **(halogens or {})}
        chain, _, _, _ = _select_winning_structure(graph, root, coming_from, named, mol, aromatic_atoms, unsaturated)
        if any(a not in chain for a in labelled):
            raise UnsupportedStructure("an isotopically modified atom off the principal chain of a substituent is not supported yet")
        positions = _lowest_chain_positions(mol, chain, root, labelled)
        stem = alkyl_name(len(chain))[:-2]
        stem_index = name.rfind(stem)
        if stem_index < 0:
            raise UnsupportedStructure("the principal chain of this substituent name is not delimited")
        omit = len(chain) == 1 or _fully_modified(labelled, positions, capacity, set(chain))
    context["consumed"].update(labelled)
    text = descriptor(labelled, positions, omit, capacity=capacity)
    return name[:stem_index] + text + name[stem_index:], compound


_CHALCOGEN_RUN_STEM = re.compile(r"(?:di|tri|tetra|penta|hexa)?(?:sulfanyl|selanyl|tellanyl)$")


def _chalcogen_run(mol, graph, root, coming_from):
    """The atoms of the chain of one chalcogen element that starts at `root` (a disulfanyl or trisulfanyl group), else None."""
    element = mol.GetAtomWithIdx(root).GetAtomicNum()
    if element not in (16, 34, 52):
        return None
    run, previous = [root], coming_from
    while True:
        onward = [n for n in graph[run[-1]] if n != previous and mol.GetAtomWithIdx(n).GetAtomicNum() == element]
        if len(onward) != 1:
            break
        previous = run[-1]
        run.append(onward[0])
    return run if len(run) >= 2 else None


def _lowest_chain_positions(mol, chain, root, labelled):
    """{chain atom: locant}; a chain that reads equally from either end (isopropyl) is numbered to give the modified
    atoms the lower locants (P-82.5.2)."""
    forward = {a: i + 1 for i, a in enumerate(chain)}
    from rdkit import Chem

    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    mirrored = all(ranks[a] == ranks[b] for a, b in zip(chain, chain[::-1])) and chain.index(root) == len(chain) - 1 - chain.index(root)
    if not mirrored:
        return forward
    backward = {a: len(chain) - i for i, a in enumerate(chain)}
    return min((forward, backward), key=lambda option: sorted(option[a] for a in labelled))


def _fully_modified(labelled, positions, capacity, expected):
    """P-82.6.1.3: every position of the group modified in the same way, none keeping a hydrogen."""
    if set(labelled) != expected or capacity is None:
        return False
    entries = {repr(sorted(e["H"])) + str(e["skeleton"]) for e in labelled.values()}
    return len(entries) == 1 and all(sum(e["H"].values()) == capacity[a] for a, e in labelled.items())


def _top_level_locants(name):
    depth, flat = 0, []
    for ch in name:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        flat.append(ch if depth == 0 and ch not in "([{" else " ")
    return sorted(int(x) for group in re.findall(r"(\d+(?:,\d+)*)-", "".join(flat)) for x in group.split(","))


def _ring_group_stem(mol, root):
    """('phenyl' or 'cyclo<alkyl>', ring) for an isolated benzene ring or saturated carbocycle attached at `root`."""
    info = mol.GetRingInfo()
    ring = next((r for r in info.AtomRings() if root in r), None)
    if ring is None or any(info.NumAtomRings(a) != 1 or mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in ring):
        return None
    aromatic = [mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring]
    if len(ring) == 6 and all(aromatic):
        return "phenyl", ring
    saturated = not any(aromatic) and all(
        mol.GetBondBetweenAtoms(ring[i], ring[(i + 1) % len(ring)]).GetBondTypeAsDouble() == 1.0 for i in range(len(ring))
    )
    return ("cyclo" + alkyl_name(len(ring)), ring) if saturated else None


def _phenyl_label_positions(mol, root, labelled, name="phenyl", coming_from=None):
    """{ring atom: locant} of a phenyl or cycloalkyl group, numbered from the attachment atom in the direction that
    gives the ring substituents of the name, then the modified atoms, the lower locants (P-82.5.2); None unless the
    ring is an isolated benzene ring or saturated carbocycle."""
    found = _ring_group_stem(mol, root)
    if found is None or not name.endswith(found[0]):
        return None
    ring = found[1]
    if any(a not in ring for a in labelled):
        return None
    start = ring.index(root)
    options = [{ring[(start + step * k) % len(ring)]: k + 1 for k in range(len(ring))} for step in (1, -1)]
    if name != found[0]:
        cited = _top_level_locants(name)
        substituted = [
            a for a in ring
            if any(
                n.GetIdx() not in ring and n.GetIdx() != coming_from and n.GetIdx() not in labelled and n.GetAtomicNum() != 1
                for n in mol.GetAtomWithIdx(a).GetNeighbors()
            )
        ]
        options = [o for o in options if sorted(o[a] for a in substituted) == cited]
        if not options:
            return None
    return min(options, key=lambda option: sorted(option[a] for a in labelled))


def _hydrogen_capacity(context, atoms):
    mol = context.get("mol")
    return {a: mol.GetAtomWithIdx(a).GetTotalNumHs() for a in atoms} if mol is not None else None


def _plain_chain_positions(graph, root, coming_from, atoms, mol=None):
    """{atom: locant} when the branch is an unbranched chain of carbon atoms read from its attachment atom."""
    positions = {}
    previous, current = coming_from, root
    while current is not None:
        if mol is not None and len(atoms) > 1 and mol.GetAtomWithIdx(current).GetAtomicNum() != 6:
            return None
        if len([n for n in graph[current] if n != previous]) > 1:
            return None
        positions[current] = len(positions) + 1
        onward = [n for n in graph[current] if n != previous]
        previous, current = current, (onward[0] if onward else None)
    return positions if set(positions) == atoms else None


def name_branch(graph, root, coming_from, halogens=None, aromatic_atoms=None, mol=None, unsaturated=None):
    token = BRANCH_STEREO_MOL.set(mol)
    try:
        named = _name_branch_isotopes(graph, root, coming_from, halogens, aromatic_atoms, mol, unsaturated)
    finally:
        BRANCH_STEREO_MOL.reset(token)
    context = BRANCH_STEREO.get()
    if not context or mol is None or not context["bonds"]:
        return named
    return cite_heteroatom_double_bonds(named, graph, root, coming_from, mol, context)


_LEADING_DESCRIPTOR = re.compile(r"^[\[{(]*\([^()]*\)-")


def _chain_ylidene_carbon(mol, idx):
    """A chain carbon doubly bonded to a nitrogen: the name of its chain cites the descriptor (P-93.5.1.4.2.1)."""
    from ._hetero_prefixes import is_functional_carbon

    atom = mol.GetAtomWithIdx(idx)
    return atom.GetAtomicNum() == 6 and not atom.IsInRing() and not is_functional_carbon(mol, idx)


def cite_heteroatom_double_bonds(named, graph, root, coming_from, mol, context):
    """P-93.4.2.1.3: the E/Z descriptor of a C=N or N=N bond is cited in front of the smallest group that holds both of
    its atoms, unlocanted because no skeletal locant is present in the name of that group."""
    name, compound = named
    if _LEADING_DESCRIPTOR.match(name):
        return named
    pieces = []
    for child in graph[root]:
        if child != coming_from:
            seen, stack = {child}, [child]
            while stack:
                for n in graph[stack.pop()]:
                    if n not in seen and n != root:
                        seen.add(n)
                        stack.append(n)
            pieces.append(seen)
    inside = {root}.union(*pieces)
    codes = []
    for (a, b), code in sorted(context["bonds"].items()):
        if a not in inside or b not in inside or any(a in p and b in p for p in pieces):
            continue
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond.IsInRing() or 7 not in (mol.GetAtomWithIdx(a).GetAtomicNum(), mol.GetAtomWithIdx(b).GetAtomicNum()):
            continue
        if any(_chain_ylidene_carbon(mol, i) for i in (a, b)):
            continue
        codes.append(code)
        context["used"].add(("bond", (a, b)))
    return (f"({','.join(codes)})-{name}", True) if codes else named


def _name_branch_isotopes(graph, root, coming_from, halogens=None, aromatic_atoms=None, mol=None, unsaturated=None):
    forced = FORCED_BRANCH_NAMES.get()
    if forced and mol is not None and mol.GetNumAtoms() == forced[0] and root in forced[1]:
        return forced[1][root]
    context = ISOTOPE_LABELS.get()
    if not context:
        return _name_branch_with_phane(graph, root, coming_from, halogens, aromatic_atoms, mol, unsaturated)
    cache = context.setdefault("cache", {})
    key = (root, coming_from)
    if key not in cache and mol is not None:
        carboxy = _labelled_carboxy(graph, root, coming_from, mol, context)
        if carboxy is not None:
            cache[key] = (carboxy[0], frozenset(carboxy[1]))
    if key not in cache:
        before = set(context["consumed"])
        result = _name_branch_with_phane(graph, root, coming_from, halogens, aromatic_atoms, mol, unsaturated)
        labelled = _label_branch(
            result, graph, root, coming_from, halogens, mol, aromatic_atoms or frozenset(), unsaturated
        )
        cache[key] = (labelled, frozenset(context["consumed"] - before))
    result, used = cache[key]
    context["consumed"].update(used)
    return result


def _labelled_carboxy(graph, root, coming_from, mol, context):
    """((name, compound), consumed atoms) for an isotopically modified -COOH group, cited as one prefix (P-82.2.2.2)."""
    atom = mol.GetAtomWithIdx(root)
    others = [n for n in graph[root] if n != coming_from]
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or len(others) != 2 or atom.GetFormalCharge():
        return None
    oxo = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0]
    hydroxy = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1 and mol.GetAtomWithIdx(n).GetTotalNumHs() == 1 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 1.0]
    if len(oxo) != 1 or len(hydroxy) != 1:
        return None
    group = [root, oxo[0], hydroxy[0]]
    labelled = {a: context["labels"][a] for a in group if a in context["labels"]}
    if not labelled:
        return None
    from ._isotope_labels import _nuclide_sort_key
    from ._polyfunctional import _modifications

    nuclides = sorted((n for e in labelled.values() for n in _modifications(e)), key=_nuclide_sort_key)
    context["consumed"].update(labelled)
    return ("(" + ",".join(nuclides) + ")carboxy", False), set(labelled)


_CHALCOGEN_OYL = {16: "thioyl", 34: "selenoyl", 52: "telluroyl"}


def _chalcogen_acyl(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """'ethanethioyl' for R-C(=S)- on a carbon R, read off the name of the oxygen acyl group (P-65.1.7.2.2)."""
    if mol is None:
        return None
    atom = mol.GetAtomWithIdx(root)
    link = mol.GetBondBetweenAtoms(root, coming_from)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetTotalNumHs() or atom.GetFormalCharge() or link is None or link.GetBondTypeAsDouble() != 1.0:
        return None
    others = [n for n in graph[root] if n != coming_from]
    chalcogen = [
        n for n in others
        if mol.GetAtomWithIdx(n).GetAtomicNum() in _CHALCOGEN_OYL
        and mol.GetAtomWithIdx(n).GetDegree() == 1
        and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0
    ]
    if len(others) != 2 or len(chalcogen) != 1:
        return None
    (carbon,) = [n for n in others if n != chalcogen[0]]
    if mol.GetAtomWithIdx(carbon).GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(root, carbon).GetBondTypeAsDouble() != 1.0:
        return None
    healed = Chem.RWMol(mol)
    healed.GetAtomWithIdx(chalcogen[0]).SetAtomicNum(8)
    try:
        name, compound = _name_branch_with_phane(graph, root, coming_from, halogens, aromatic_atoms, healed.GetMol(), True)
    except UnsupportedStructure:
        return None
    word = _CHALCOGEN_OYL[mol.GetAtomWithIdx(chalcogen[0]).GetAtomicNum()]
    if name.endswith("acetyl"):
        return name[: -len("acetyl")] + "ethane" + word, True
    if name.endswith("benzoyl"):
        return name[: -len("benzoyl")] + "benzenecarbo" + word, True
    if name.endswith("anoyl"):
        return name[: -len("oyl")] + "e" + word, True
    return None


def _name_branch_with_phane(graph, root, coming_from, halogens, aromatic_atoms, mol, unsaturated):
    acyl = _chalcogen_acyl(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if acyl is not None:
        return acyl
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

    `aromatic_atoms`: the aromatic atom indices of the molecule, for the chain walk (ring roots go to the
    ring-group namer, which reads `mol`).

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
        return halogens[root], isinstance(halogens[root], CompoundPrefix)
    if unsaturated is None:
        unsaturated = mol is not None
    attach_order = 1.0
    if unsaturated:
        attach_order = _bond_order(mol, root, coming_from)
        aromatic_atoms = aromatic_atoms or frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())

    if mol is not None and attach_order == 1.0 and mol.GetRingInfo().NumRings() >= 1:
        from ._appendix3_skeletons import appendix3_group

        natural_product = appendix3_group(mol, graph, root, coming_from)
        if natural_product is not None:
            return natural_product

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

    if mol is not None and attach_order == 1.0 and mol.GetRingInfo().NumRings() >= 1:
        from ._sugar_substituted import sugar_substituent_group

        sugar = sugar_substituent_group(mol, graph, root, coming_from)
        if sugar is not None:
            return sugar

    if mol is not None and attach_order == 1.0 and mol.GetAtomWithIdx(root).GetAtomicNum() in (6, 7) and not mol.GetAtomWithIdx(root).IsInRing():
        from ._formazan import formazan_substituent

        formazan = formazan_substituent(mol, graph, root, coming_from, halogens, aromatic_atoms)
        if formazan is not None:
            return formazan

    if mol is not None and attach_order == 1.0:
        from ._heteroacyclic import skeletal_substituent
        from ._skeletal_group import skeletal_chain_group, skeletal_ring_group

        skeletal = (
            skeletal_substituent(mol, graph, root, coming_from)
            or skeletal_chain_group(mol, graph, root, coming_from)
            or skeletal_ring_group(mol, graph, root, coming_from)
        )
        if skeletal is not None:
            return skeletal

    if mol is not None and unsaturated:
        from ._hetero_prefixes import hetero_branch_name

        hetero = hetero_branch_name(graph, root, coming_from, halogens, aromatic_atoms, mol)
        if hetero is not None:
            return hetero

    if mol is not None and mol.GetAtomWithIdx(root).IsInRing():
        from ._diester_ring_diyl import ring_substituent_name

        return ring_substituent_name(mol, graph, root, coming_from)

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
    """Like `substituents_for_chain`, but a `terminals` leaf is cited as a compound (enclosed) substituent exactly
    when its name is a `CompoundPrefix`; every other branch is named by `name_branch`."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        entries = []
        for root in branch_roots:
            if root in terminals:
                entries.append((terminals[root], isinstance(terminals[root], CompoundPrefix)))
            else:
                entries.append(name_branch(graph, root, atom, terminals, mol=mol))
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


def _chain_isotope_key(chain):
    context = ISOTOPE_LABELS.get()
    if not context:
        return ()
    from ._isotope_labels import modification_key

    return modification_key(context["labels"], {a: i for i, a in enumerate(chain, start=1)})


def _labels_consumed():
    context = ISOTOPE_LABELS.get()
    return context["consumed"] if context else set()


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
    consumed = _labels_consumed()
    base = set(consumed)
    best_gained = base
    for chain in chains:
        consumed.clear()
        consumed.update(base)
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
        key = _candidate_key(grouped) + (_chain_isotope_key(chain), name)
        if best_key is None or key < best_key:
            best_key, best_chain, best_name, best_compound = key, chain, name, is_compound
            best_gained = set(consumed)

    consumed.clear()
    consumed.update(best_gained)
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
        and not _carries_label(root, branch_roots)
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
    consumed = _labels_consumed()
    base = set(consumed)
    best_gained = base
    for before_branch, after_branch in orientations:
        for path_before in paths[before_branch]:
            for path_after in paths[after_branch]:
                consumed.clear()
                consumed.update(base)
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
                key = _candidate_key(grouped) + (_chain_isotope_key(spine), name)
                if best_key is None or key < best_key:
                    best_key = key
                    best_chain = spine
                    best_position = root_locant
                    best_name = name
                    best_gained = set(consumed)

    consumed.clear()
    consumed.update(best_gained)
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


BRANCH_STEREO = contextvars.ContextVar("branch_stereo", default=None)
BRANCH_STEREO_MOL = contextvars.ContextVar("branch_stereo_mol", default=None)


@contextlib.contextmanager
def cited_branch_stereo(mol, graph, blocked, roots):
    """Substituents named inside cite the CIP descriptors of their own stereo elements; raises when one is left uncited."""
    inside, stack = set(), [r for r in roots if r not in blocked]
    while stack:
        a = stack.pop()
        if a not in inside:
            inside.add(a)
            stack.extend(n for n in graph[a] if n not in blocked and n not in inside)
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    context = {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in inside and a.HasProp("_CIPCode")},
        "bonds": {
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
            for b in probe.GetBonds()
            if b.GetBeginAtomIdx() in inside and b.GetEndAtomIdx() in inside and b.HasProp("_CIPCode")
        },
        "used": set(),
    }
    token = BRANCH_STEREO.set(context)
    try:
        yield
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element inside a substituent is not cited by any supported name")


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
            elif (a in positions) != (b in positions) and _ylidene_to_heteroatom(a, b, positions):
                # P-93.5.1.4.2.1: the double bond to an ylidene group takes the locant of the parent atom
                entries.append((positions[a] if a in positions else positions[b], code))
                if record:
                    context["used"].add(("bond", (a, b)))
    return sorted(entries)


def _ylidene_to_heteroatom(a, b, positions):
    mol = BRANCH_STEREO_MOL.get()
    if mol is None:
        return False
    outside = mol.GetAtomWithIdx(b if a in positions else a)
    return outside.GetAtomicNum() == 7 and not mol.GetBondBetweenAtoms(a, b).IsInRing()


def _branch_stereo_rank(entries):
    return tuple(0 if code in "RZr" else 1 for _, code in entries)


def _branch_stereo_prefix(entries, single_atom=False):
    if single_atom:
        entries = [("", code) for _, code in entries]
    return "(" + ",".join(f"{locant}{code}" for locant, code in entries) + ")-" if entries else ""


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
        if mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or (is_functional_carbon(mol, n) and not _terminal_amide(graph, mol, n, node)):
            continue
        if mol.GetAtomWithIdx(n).IsInRing():
            continue
        children.append(n)
    return children


def _terminal_amide(graph, mol, carbon, parent):
    """A -C(=O)-NR2 carbon ending a carbon chain of more than one carbon, cited as 'amino' and 'oxo' on the chain
    (P-66.1.1.4.1.1)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.IsInRing() or atom.GetFormalCharge() or mol.GetAtomWithIdx(parent).GetAtomicNum() != 6 or len(graph[carbon]) != 3:
        return False
    others = [n for n in graph[carbon] if n != parent]
    oxo = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() == 2.0]
    amino = [n for n in others if (mol.GetAtomWithIdx(n).GetAtomicNum() == 7 or (mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 2)) and not mol.GetAtomWithIdx(n).GetFormalCharge() and not mol.GetAtomWithIdx(n).IsInRing() and mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() == 1.0]
    return len(oxo) == 1 and len(amino) == 1 and mol.GetBondBetweenAtoms(carbon, parent).GetBondTypeAsDouble() == 1.0


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
    consumed = _labels_consumed()
    base = set(consumed)
    for chain in _unsaturated_candidate_chains(graph, root, coming_from, halogens, mol, aromatic_atoms):
        consumed.clear()
        consumed.update(base)
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
            len(chain),
            root_position,
            suffix,
            ene,
            yne,
            grouped,
            tert_butyl=_is_tert_butyl(graph, root, coming_from, halogens),
            omit_locants=omits_all_locants(mol, chain, grouped, single_kind=False, free_atoms={coming_from}),
        )
        multiple = sorted(ene + yne)
        stereo_rank = _branch_stereo_rank(_branch_stereo_entries({a: i for i, a in enumerate(chain, start=1)}))
        key = (
            -len(multiple), -len(ene), root_position, multiple, sorted(ene), -total_count, locant_set, citation,
            _chain_isotope_key(chain), stereo_rank, name,
        )
        if best is None or key < best[0]:
            best = (key, chain, root_position, name, is_compound, set(consumed))
    _, chain, root_position, name, is_compound, gained = best
    consumed.clear()
    consumed.update(gained)
    stereo_prefix = _branch_stereo_prefix(
        _branch_stereo_entries({a: i for i, a in enumerate(chain, start=1)}, record=True), single_atom=len(chain) == 1
    )
    return chain, root_position, stereo_prefix + name, is_compound or bool(stereo_prefix)


def _is_tert_butyl(graph, root, coming_from, halogens):
    # P-29.6.1: the unsubstituted (CH3)3C- group keeps its retained name.
    others = [n for n in graph[root] if n != coming_from]
    return (
        len(others) == 3
        and all(len(graph[n]) == 1 and n not in halogens for n in others)
        and not _carries_label(root, others)
    )


def _carries_label(root, atoms):
    """A modified group cannot keep a retained name (P-82.2)."""
    context = ISOTOPE_LABELS.get()
    return bool(context) and any(a in context["labels"] for a in [root, *atoms])


def _unsaturated_chain_name(length, root_position, suffix, ene, yne, grouped, tert_butyl, omit_locants=False):
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
        prefix = format_substituent_prefixes(
            grouped, omit_locants=(length == 1), omit_all=omit_locants and root_position == 1 and not (ene or yne)
        ) if grouped else ""
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
    chain, _, _, _ = _select_winning_structure(graph, branch_root, ring_atom, halogens or {}, mol, unsaturated=mol is not None)
    descriptor = f"({'' if len(chain) == 1 else site_locant}{r_or_s})-{branch_name}"
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
