"""Naming of saturated and mancude (maximally unsaturated) monocyclic
rings containing exactly one heteroatom (O, S, Se, Te, or N) and no
substituents, using Hantzsch-Widman-system and retained names, per the
IUPAC 2013 Recommendations ("the Blue Book"):

Saturated rings (3- to 7-membered):
- P-22.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf,
  Table 22.1): the Hantzsch-Widman stem for a saturated ring is
  '-irane'/'-irene' (3), '-etane' (4), '-olane' (5), '-ane' (6), '-epane'
  (7) for O and S ('oxa'/'thia' + stem), giving oxirane, oxetane,
  oxolane, oxane, oxepane and thiirane, thietane, thiolane, thiane,
  thiepane.
- For N, the saturated-ring stem is '-iridine' (3), '-etidine' (4),
  '-olidine' (5), '-inane' (6), '-epane' (7) ('aza' + stem) -- but four of
  these five (all but 7-membered azepane) are retained names that are
  themselves the preferred IUPAC name in place of the literal 'aza' +
  stem form: aziridine, azetidine, pyrrolidine, and piperidine (not
  azirane/azetane/azolidine/azinane) -- confirmed via IUPAC's P-22.2.1
  Table 2.3 listing piperidine/pyrrolidine (among others) as retained
  names that are PINs. Aziridine/azetidine happen to coincide with the
  literal stem-based construction; pyrrolidine/piperidine don't.

Mancude (aromatic) rings, single heteroatom, 5- and 6-membered only
(P-22.2.1 Table 2.2): furan/thiophene/selenophene/tellurophene (5-membered
O/S/Se/Te) and pyridine (6-membered N) are retained names that are PINs
outright; pyrrole (5-membered N) needs an indicated-hydrogen prefix --
its PIN is '1H-pyrrole', not bare 'pyrrole' (confirmed via Table 2.2 and
PubChem).

Mancude rings, two heteroatoms (also P-22.2.1 Table 2.2), unsubstituted
only: the 5-membered N+N/N+O/N+S rings need an indicated-hydrogen prefix
on the NH ring member -- 1H-imidazole (1,3-diazole) and 1H-pyrazole
(1,2-diazole); 1,3-oxazole and 1,2-oxazole (isoxazole); 1,3-thiazole and
1,2-thiazole (isothiazole) don't carry indicated hydrogen since neither
heteroatom bears an H. The 6-membered N+N rings (pyridazine, pyrimidine,
pyrazine) are, like pyridine, fully mancude without indicated hydrogen.
All nine confirmed as PubChem's IUPACName for the exact SMILES: CID
795/1048/9255/9254/9256/67515/9259/9260/9261.

Se/Te analogues of the N+S mancude pair (Table 2.2 lists them alongside
thiazole/isothiazole as a symmetric O->S->Se->Te chalcogen series):
1,3-selenazole and 1,2-selenazole, 1,3-tellurazole and 1,2-tellurazole --
none carry indicated hydrogen, same as their S analogues. Three of the
four confirmed via PubChem's IUPACName for the exact SMILES: 1,3-selenazole
(CID 11686913), 1,2-selenazole (CID 13224788), 1,2-tellurazole (CID
102212476); 1,3-tellurazole has no PubChem record (CID 0 for the exact
SMILES) so it's confirmed from Table 2.2's text alone, matching the same
symmetric pattern as the other three. The N+O pair (oxazole/isoxazole)
has no listed Se/Te analogue in Table 2.2, so that combination stays out
of scope.

Three-or-more heteroatom rings (triazole, tetrazole, etc.), 6-membered
O/S/Se/Te rings (pyran/thiopyran/selenopyran/telluropyran, which need an
indicated-hydrogen prefix themselves since they aren't fully mancude with
a single chalcogen), and 2-or-more substituents are out of scope --
separate future tasks.

Since this module's core job (`has_hetero_monocyclic_name`/
`name_hetero_monocyclic`) is recognizing the exact unsubstituted parent
for a fixed, small (element(s), ring size, saturation) table -- no
locants to assign, no substituent numbering -- an exact whole-molecule
canonical-SMILES match against each name's structure is both sufficient
and simplest, mirroring `_peri_fused_aromatic.py`'s approach.

Formulas cross-checked (all well-known compounds): oxirane C2H4O,
piperidine C5H11N, thiane C5H10S, furan C4H4O, pyridine C5H5N, imidazole
C3H4N2, pyrimidine C4H4N2, etc. -- see `_RETAINED_NAME_SMILES` and
`_MANCUDE_NAME_SMILES`.

Saturated rings, two heteroatoms, unsubstituted only: the 6-membered
1,4-related pairs (N+O morpholine, N+N piperazine, N+S thiomorpholine,
O+O 1,4-dioxane, O+S 1,4-oxathiane, S+S 1,4-dithiane), the 5-membered
1,3-related pairs (N+N imidazolidine, N+O 1,3-oxazolidine, N+S
1,3-thiazolidine, O+O 1,3-dioxolane, O+S 1,3-oxathiolane, S+S
1,3-dithiolane), and the 5-membered 1,2-related pairs (N+N pyrazolidine,
N+O 1,2-oxazolidine, N+S 1,2-thiazolidine, O+O 1,2-dioxolane, O+S
1,2-oxathiolane, S+S 1,2-dithiolane) all have their own retained/
systematic name -- per Table 2.3 and P-22.2.2.1.2/.1.3, a pair with no
retained name (all but the N+N ones, which are simply irregular retained
names carrying no locants at all, 1,2- or 1,3- alike) always cites its
heteroatom locants; PubChem's own computed names are unreliable here
(they drop the locants for several of the O/S-only 1,2-pairs, contradicted
by Table 2.3's explicit '1,2-oxazolidine (PIN)'/'1,2-thiazolidine (PIN)'
entries and P-22.2.2.1.3's own '1,2-oxathiolane (PIN)' worked example),
so the 1,2-pairs are sourced from the primary text directly rather than
from PubChem.

Explicitly out of scope for the unsubstituted-only functions above: any
substituent, partially-saturated indicated-hydrogen forms other than the
ones listed above, three or more heteroatoms, heteroatoms other than
O/S/Se/Te/N, saturated two-heteroatom ring sizes/relationships other than
the three listed above, and ring sizes outside the single-heteroatom
tables above. `has_hetero_monocyclic_name`
returns False for all of these, so `core.py`'s existing dispatch (which
already rejects heteroatoms outside a few specific recognized shapes)
continues to raise `UnsupportedStructure` for them, unchanged.

Substituent naming (`has_hetero_monocyclic_substituent_name`/
`name_hetero_monocyclic_substituent`), one axis further than the
unsubstituted-only functions above:

- Scope: one or more substituents (a halogen, or any group `name_branch`
  can name, repeats and mixed kinds both allowed) on one of the 19
  mancude parents listed in `_MANCUDE_NAME_SMILES`/
  `_TWO_HETEROATOM_MANCUDE_NAME_SMILES` above (single-heteroatom 5/6-
  membered plus two-heteroatom 5/6-membered, Se/Te analogues included)
  -- everything else about the parent (ring size, unsaturation,
  heteroatom set) stays exactly as recognized by the unsubstituted
  table. Saturated/partially saturated rings and polycyclic systems
  remain out of scope (a second ring anywhere -- including one folded
  into a substituent itself, like a cyclopropyl group -- is rejected via
  `mol.GetRingInfo().NumRings() == 1`).
- Multiple substituents (extending the original single-substituent axis):
  each substituted ring
  atom must still carry exactly one exocyclic branch (no gem-disubstitution
  on a ring atom, which mancude/aromatic ring carbons can't have anyway).
  Locants are chosen the same way `_alcohol.py`'s `_name_cyclic_alcohol`
  numbers a substituted ring: try every rotation/direction of the ring
  that satisfies the parent's fixed heteroatom role sequence, and among
  the valid ones pick the lowest locant *set* for the substituted atoms
  as a whole (P-14.4), breaking ties by giving the alphabetically-first
  substituent name the lower individual locant (P-14.5.2) -- mirrors
  `_alcohol.py`'s `_ring_candidate_key` exactly (locant-set, then
  citation-locants-in-alphabetical-name-order, then the assembled name
  itself as a final tiebreak). Substituents are cited alphabetically with
  ordinary di-/tri- multiplying prefixes for repeats
  (`format_substituent_prefixes`, same helper `_alcohol.py`/
  `_carboxylic_acid.py` already use). PubChem verification: `2,5-
  dimethylfuran` (CID 12266), `2-chloro-5-methylthiophene` (CID 140208,
  confirming alphabetical mixed-kind ordering), `3,4-dichloropyridine`
  (CID 2736081), `2,4-dimethyl-1,3-oxazole` (CID 138961, heteroatom-
  asymmetric parent so the locants are structurally fixed, not a lowest-
  set choice), `2,3-dichloropyrazine` (CID 78575).
- A chalcogen ring atom (O/S/Se/Te) or a "pyridine-type" ring nitrogen
  (no N-H in the unsubstituted parent: pyridine's N, and the second
  nitrogen of imidazole/pyrazole/oxazole/isoxazole/thiazole/isothiazole/
  selenazole/isoselenazole/tellurazole/isotellurazole, plus both
  nitrogens of pyridazine/pyrimidine/pyrazine) can never itself bear the
  substituent -- there's no spare valence to replace, unlike a ring
  carbon's H or, uniquely among the nitrogens here, pyrrole/imidazole/
  pyrazole's N-H (position 1 in each), which is a real substitutable
  position (P-22.2.1) -- confirmed by PubChem's IUPACName for `Cn1cccc1`
  ("1-methylpyrrole", CID 7304): the indicated-hydrogen prefix
  disappears once the position it marked is substituted directly,
  exactly as the analogous case already works for `1H-pyrrole` itself
  losing its `1H-` prefix nowhere else in this project's ring-naming
  modules, since the locant `1-` alone already pins the position.
- `_ROLE_SEQUENCES` fixes each parent's own heteroatom locants (already
  implied by the unsubstituted dictionaries above -- e.g. `1,3-thiazole`
  means S=1/N=3) as an explicit (element, substitutable) tuple per
  position 1..ring_size, in one direction. Matching a candidate molecule
  tries every rotation and both directions of its own ring traversal
  (mirroring `_alcohol.py`'s `_name_cyclic_alcohol` candidate-rotation
  approach) against every parent's role sequence; a parent with a real
  reflection symmetry (furan/thiophene/selenophene/tellurophene/pyridine/
  pyrrole's own carbons, and the three 6-membered two-nitrogen rings
  pyridazine/pyrimidine/pyrazine) yields two or more valid alignments
  with different substituent locants, and the lowest one wins (P-14.4) --
  this falls out of the brute-force search automatically, with no need to
  hand-classify which parents are symmetric. The other, heteroatom-
  asymmetric parents (thiazole/oxazole/imidazole/pyrazole and their Se/Te
  analogues) yield exactly one valid alignment, so the substituent locant
  is structurally fixed.
- Broad PubChem verification of the exact SMILES (representative single-
  substituent case per parent, methyl or chloro): 3-methylfuran (CID
  13587), 3-methylthiophene (CID 12024), 3-methylselenophene (CID
  13022371), 3-methyltellurophene (CID 13022372), 2-methylpyridine (CID
  7975) / 4-chloropyridine (CID 12288), 1-methylpyrrole (CID 7304),
  5-methyl-1H-imidazole (CID 13195), 5-methyl-1H-pyrazole (CID 15073),
  3-methyl-1,2-oxazole (CID 96098), 4-methyl-1,3-thiazole (CID 12748),
  3-methyl-1,2-thiazole (CID 12747), 4-methylpyridazine (CID 136882),
  2-methylpyrimidine (CID 78748), 2-methylpyrazine (CID 7976). The Se/Te
  two-heteroatom analogues (selenazole/tellurazole family) aren't
  separately re-verified substituted here -- their role sequences mirror
  thiazole/isothiazole's exactly (same asymmetric heteroatom pattern,
  just O->S->Se->Te), so the same locant logic applies unchanged.
"""

from rdkit import Chem

from ._common import adjacency, halogen_substituents, lowest_locant_set, ring_cycle
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_RETAINED_NAME_SMILES = {
    ("O", 3): ("oxirane", "C1CO1"),
    ("O", 4): ("oxetane", "C1CCO1"),
    ("O", 5): ("oxolane", "C1CCCO1"),
    ("O", 6): ("oxane", "C1CCCCO1"),
    ("O", 7): ("oxepane", "C1CCCCCO1"),
    ("S", 3): ("thiirane", "C1CS1"),
    ("S", 4): ("thietane", "C1CCS1"),
    ("S", 5): ("thiolane", "C1CCCS1"),
    ("S", 6): ("thiane", "C1CCCCS1"),
    ("S", 7): ("thiepane", "C1CCCCCS1"),
    ("N", 3): ("aziridine", "C1CN1"),
    ("N", 4): ("azetidine", "C1CCN1"),
    ("N", 5): ("pyrrolidine", "C1CCCN1"),
    ("N", 6): ("piperidine", "C1CCCCN1"),
    ("N", 7): ("azepane", "C1CCCCCN1"),
}
_MANCUDE_NAME_SMILES = {
    ("O", 5): ("furan", "c1ccoc1"),
    ("S", 5): ("thiophene", "c1ccsc1"),
    ("Se", 5): ("selenophene", "c1cc[se]c1"),
    ("Te", 5): ("tellurophene", "c1cc[te]c1"),
    ("N", 5): ("1H-pyrrole", "c1cc[nH]c1"),
    ("N", 6): ("pyridine", "c1ccncc1"),
}
_TWO_HETEROATOM_SATURATED_NAME_SMILES = {
    frozenset(("N", "O")): ("morpholine", "C1COCCN1"),
    frozenset(("N", "N")): ("piperazine", "C1CNCCN1"),
    frozenset(("N", "S")): ("thiomorpholine", "C1CSCCN1"),
    frozenset(("O", "O")): ("1,4-dioxane", "C1COCCO1"),
    frozenset(("O", "S")): ("1,4-oxathiane", "C1COCCS1"),
    frozenset(("S", "S")): ("1,4-dithiane", "C1CSCCS1"),
}
_FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("imidazolidine", "C1CNCN1"),
    frozenset(("N", "O")): ("1,3-oxazolidine", "C1CNCO1"),
    frozenset(("N", "S")): ("1,3-thiazolidine", "C1CSCN1"),
    frozenset(("O", "O")): ("1,3-dioxolane", "C1COCO1"),
    frozenset(("O", "S")): ("1,3-oxathiolane", "C1CSCO1"),
    frozenset(("S", "S")): ("1,3-dithiolane", "C1CSCS1"),
}
_FIVE_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("pyrazolidine", "C1CCNN1"),
    frozenset(("N", "O")): ("1,2-oxazolidine", "C1CCON1"),
    frozenset(("N", "S")): ("1,2-thiazolidine", "C1CCSN1"),
    frozenset(("O", "O")): ("1,2-dioxolane", "C1CCOO1"),
    frozenset(("O", "S")): ("1,2-oxathiolane", "C1CCOS1"),
    frozenset(("S", "S")): ("1,2-dithiolane", "C1CCSS1"),
}
_TWO_HETEROATOM_MANCUDE_NAME_SMILES = {
    ("1,3-diazole", 5): ("1H-imidazole", "c1cnc[nH]1"),
    ("1,2-diazole", 5): ("1H-pyrazole", "c1cc[nH]n1"),
    ("1,3-oxazole", 5): ("1,3-oxazole", "c1cocn1"),
    ("1,2-oxazole", 5): ("1,2-oxazole", "c1ccon1"),
    ("1,3-thiazole", 5): ("1,3-thiazole", "c1cscn1"),
    ("1,2-thiazole", 5): ("1,2-thiazole", "c1ccsn1"),
    ("1,2-diazine", 6): ("pyridazine", "c1ccnnc1"),
    ("1,3-diazine", 6): ("pyrimidine", "c1ccncn1"),
    ("1,4-diazine", 6): ("pyrazine", "c1cnccn1"),
    ("1,3-selenazole", 5): ("1,3-selenazole", "c1cnc[se]1"),
    ("1,2-selenazole", 5): ("1,2-selenazole", "c1ccn[se]1"),
    ("1,3-tellurazole", 5): ("1,3-tellurazole", "c1cnc[te]1"),
    ("1,2-tellurazole", 5): ("1,2-tellurazole", "c1ccn[te]1"),
}
_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name
    for name, smiles in (
        *_RETAINED_NAME_SMILES.values(),
        *_MANCUDE_NAME_SMILES.values(),
        *_TWO_HETEROATOM_SATURATED_NAME_SMILES.values(),
        *_FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.values(),
        *_FIVE_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES.values(),
        *_TWO_HETEROATOM_MANCUDE_NAME_SMILES.values(),
    )
}


def has_hetero_monocyclic_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_hetero_monocyclic(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]


def saturated_ring_name(element: str, size: int):
    """The retained/Hantzsch-Widman name for the unsubstituted saturated
    monocyclic ring with a single `element` heteroatom (O/S/N) and `size`
    ring atoms (e.g. ('N', 6) -> 'piperidine'), or None if that combination
    isn't in P-22.2.1 Table 2.3's scope (3-7 membered O/S/N only). Exposed
    for `_ketone.py`'s hetero-ring ketone naming, which needs the bare stem
    name rather than a full unsubstituted-molecule match."""
    entry = _RETAINED_NAME_SMILES.get((element, size))
    return entry[0] if entry else None


def saturated_two_heteroatom_1_4_ring_name(elements):
    """The retained/systematic name for the unsubstituted, 6-membered,
    1,4-related two-heteroatom saturated ring whose heteroatom elements
    are `elements` (an (element, element) pair or frozenset, e.g.
    ('N', 'O') -> 'morpholine'), or None if that element pair isn't one
    of the six in P-22.2.1's scope (morpholine/piperazine/
    thiomorpholine/1,4-dioxane/1,4-oxathiane/1,4-dithiane -- other
    element pairs (Se/Te included), and other ring sizes/relationships,
    have no retained name and are out of scope). Exposed for
    `_ketone.py`'s hetero-ring ketone naming, which needs the bare stem
    name rather than a full unsubstituted-molecule match."""
    entry = _TWO_HETEROATOM_SATURATED_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_five_membered_1_3_two_heteroatom_ring_name(elements):
    """The retained/systematic name for the unsubstituted, 5-membered,
    1,3-related two-heteroatom saturated ring whose heteroatom elements
    are `elements` (an (element, element) pair or frozenset, e.g.
    ('N', 'N') -> 'imidazolidine'), or None if that element pair isn't
    one of the six in P-22.2.1's scope (imidazolidine/1,3-oxazolidine/
    1,3-thiazolidine/1,3-dioxolane/1,3-oxathiolane/1,3-dithiolane --
    other element pairs, the 1,2-relationship, and other ring sizes have
    no retained name here and are out of scope). Exposed for
    `_ketone.py`'s hetero-ring ketone naming, which needs the bare stem
    name rather than a full unsubstituted-molecule match."""
    entry = _FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def _role(element, has_h=False):
    """A ring position's fixed (element, has_h) role in its unsubstituted
    parent -- `has_h` is whether that position carries a replaceable H
    (True for every ring carbon, and for the single N-H nitrogen of
    pyrrole/imidazole/pyrazole; False for every other heteroatom, which
    has no spare valence to give up). Matching checks `has_h` two ways:
    at the substituted position, it gates whether a substituent could
    have landed there at all; at every other position, the *actual*
    H-count must still equal it -- otherwise element symbol alone can't
    tell imidazole's two chemically different nitrogens apart (both are
    plain 'N'), and two rotations would wrongly look equally valid."""
    return (element, has_h)


_C = _role("C", has_h=True)
_ROLE_SEQUENCES = {
    "furan": (_role("O"), _C, _C, _C, _C),
    "thiophene": (_role("S"), _C, _C, _C, _C),
    "selenophene": (_role("Se"), _C, _C, _C, _C),
    "tellurophene": (_role("Te"), _C, _C, _C, _C),
    "1H-pyrrole": (_role("N", has_h=True), _C, _C, _C, _C),
    "pyridine": (_role("N"), _C, _C, _C, _C, _C),
    "1H-imidazole": (_role("N", has_h=True), _C, _role("N"), _C, _C),
    "1H-pyrazole": (_role("N", has_h=True), _role("N"), _C, _C, _C),
    "1,3-oxazole": (_role("O"), _C, _role("N"), _C, _C),
    "1,2-oxazole": (_role("O"), _role("N"), _C, _C, _C),
    "1,3-thiazole": (_role("S"), _C, _role("N"), _C, _C),
    "1,2-thiazole": (_role("S"), _role("N"), _C, _C, _C),
    "1,3-selenazole": (_role("Se"), _C, _role("N"), _C, _C),
    "1,2-selenazole": (_role("Se"), _role("N"), _C, _C, _C),
    "1,3-tellurazole": (_role("Te"), _C, _role("N"), _C, _C),
    "1,2-tellurazole": (_role("Te"), _role("N"), _C, _C, _C),
    "pyridazine": (_role("N"), _role("N"), _C, _C, _C, _C),
    "pyrimidine": (_role("N"), _C, _role("N"), _C, _C, _C),
    "pyrazine": (_role("N"), _C, _C, _role("N"), _C, _C),
}

# imidazole/pyrazole's N-H can migrate to the *other* ring nitrogen
# (a real prototropic tautomer, not a naming choice) -- this only
# reshuffles which ring carbon counts as adjacent to N1, so a substituent
# set that includes N1 (locant 1, replacing the H directly) is unaffected
# and always safe; one that doesn't touch N1 at all is excluded below (see
# `_match_hetero_monocyclic_substituents`).
_TAUTOMER_AMBIGUOUS_UNLESS_N1 = {"1H-imidazole", "1H-pyrazole"}


def _ring_alignments(ring_order):
    """Every (rotation, direction) reading of a ring traversal, as a
    generic stand-in for "try both numbering directions from every
    starting atom" -- mirrors `_alcohol.py`'s `_name_cyclic_alcohol`
    candidate search, but tries every rotation (not just the one aligning
    a fixed atom to position 1) since a substituent's own position isn't
    known in advance here."""
    n = len(ring_order)
    for base in (ring_order, list(reversed(ring_order))):
        for start in range(n):
            yield base[start:] + base[:start]


def _group(substituents):
    """{ring_atom: (name, is_compound)} -> {name: {"locants": [...],
    "compound": bool}}, mirroring `_alcohol.py`'s own `_group` helper."""
    grouped = {}
    for position, (name, is_compound) in substituents.items():
        info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
        info["locants"].append(position)
    return grouped


def _candidate_key(grouped):
    """(locant_set, citation_locants, name) sort key -- lowest locant set
    for the substituted atoms as a whole wins first (P-14.4), ties broken
    by which alphabetically-ordered substituent gets the lower individual
    locant (P-14.5.2), mirroring `_alcohol.py`'s `_ring_candidate_key`."""
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return locant_set, citation_locants


def _match_hetero_monocyclic_substituents(mol):
    """(parent_name, grouped, graph) for a mancude ring carrying one or more
    substituents that matches one of `_ROLE_SEQUENCES`'s 19 parents; None if
    the molecule doesn't fit that shape at all (a second ring anywhere, a
    ring atom bearing more than one exocyclic branch, a ring size/
    heteroatom pattern outside the table, a substituent sitting on a
    non-substitutable heteroatom, or -- for imidazole/pyrazole only -- no
    substituent at the tautomer-fixing N-H position)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    if len(ring_atoms) not in (5, 6):
        return None
    if not all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms):
        return None

    graph = adjacency(mol)
    ring_set = set(ring_atoms)
    exo_by_atom = {}
    for atom in ring_atoms:
        exo = [n for n in graph[atom] if n not in ring_set]
        if not exo:
            continue
        if len(exo) != 1:
            return None
        exo_by_atom[atom] = exo[0]
    if not exo_by_atom:
        return None
    substituted_atoms = set(exo_by_atom)

    ring_order = ring_cycle(graph, ring_atoms)
    elements = {atom: mol.GetAtomWithIdx(atom).GetSymbol() for atom in ring_atoms}
    h_counts = {atom: mol.GetAtomWithIdx(atom).GetTotalNumHs() for atom in ring_atoms}
    halogens = halogen_substituents(mol)
    names_by_atom = {
        atom: name_branch(graph, root, atom, halogens) for atom, root in exo_by_atom.items()
    }

    for parent_name, role_sequence in _ROLE_SEQUENCES.items():
        if len(role_sequence) != len(ring_atoms):
            continue
        best_key = None
        best_grouped = None
        for candidate in _ring_alignments(ring_order):
            position_of = {atom: position for position, atom in enumerate(candidate, start=1)}
            valid = True
            for atom, position in position_of.items():
                role_element, role_has_h = role_sequence[position - 1]
                if elements[atom] != role_element:
                    valid = False
                    break
                if atom in substituted_atoms:
                    if not role_has_h:
                        valid = False
                        break
                elif h_counts[atom] != (1 if role_has_h else 0):
                    valid = False
                    break
            if not valid:
                continue
            if parent_name in _TAUTOMER_AMBIGUOUS_UNLESS_N1 and all(
                position_of[atom] != 1 for atom in substituted_atoms
            ):
                # None of the substituents sits at the N-H-derived locant
                # 1 -- a real prototropic-tautomer ambiguity (see
                # `_TAUTOMER_AMBIGUOUS_UNLESS_N1`'s note), so this
                # alignment (and, since the ring is otherwise rigid, every
                # alignment for this parent) can't be trusted.
                continue
            grouped = _group(
                {position_of[atom]: names_by_atom[atom] for atom in substituted_atoms}
            )
            key = _candidate_key(grouped)
            if best_key is None or key < best_key:
                best_key, best_grouped = key, grouped
        if best_grouped is not None:
            return parent_name, best_grouped, graph
    return None


def has_hetero_monocyclic_substituent_name(mol) -> bool:
    return _match_hetero_monocyclic_substituents(mol) is not None


def name_hetero_monocyclic_substituent(mol) -> str:
    parent_name, grouped, _graph = _match_hetero_monocyclic_substituents(mol)
    all_locants = {loc for info in grouped.values() for loc in info["locants"]}
    # The parent's own indicated-hydrogen prefix (pyrrole/imidazole/
    # pyrazole's '1H-') marks position 1 as the ring's one substitutable
    # N-H; once a substituent sits there instead, the locant '1-' alone
    # already pins the position, so the '1H-' becomes redundant and is
    # dropped (confirmed via PubChem: `Cn1cccc1` -> "1-methylpyrrole", not
    # "1-methyl-1H-pyrrole").
    if 1 in all_locants and parent_name.startswith("1H-"):
        parent_name = parent_name[len("1H-") :]
    prefix = format_substituent_prefixes(grouped)
    separator = "-" if parent_name[0].isdigit() else ""
    return f"{prefix}{separator}{parent_name}"
