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

Three-or-more heteroatom rings, unsubstituted only (P-22.2.2.1.7):
1,2,3-triazole and tetrazole -- both real, structurally distinct
N-H tautomers of each (confirmed via RDKit: `c1c[nH]nn1`/`c1cn[nH]n1`
and `c1[nH]nnn1`/`c1n[nH]nn1` are four separate canonical structures,
though PubChem itself collapses each ring's pair into a single CID --
67516 and 67519 respectively). P-22.2.2.1.2's "lowest possible locants
are assigned to heteroatoms, locant 1 being assigned to one of the
heteroatoms" gives the ring position count with no choice left (3 or 4
N's fill every position but one), so the indicated-hydrogen locant alone
distinguishes the tautomers: the "outer" (carbon-adjacent) N-H tautomer
is locant 1 (`1H-1,2,3-triazole`, `1H-tetrazole`, the latter without a
heteroatom-locant set at all per the worked example "1H-tetrazole (not
1H-1,2,3,4-tetrazole)"), the "inner" (nitrogen-flanked) one locant 2
(`2H-1,2,3-triazole`, `2H-tetrazole`). 1,2,4-triazole and every other
three-or-more heteroatom ring, plus 2-or-more substituents on any ring
here, are out of scope -- separate future tasks.

Pyran (6-membered, one O, five C -- P-25.7.1.3.1's indicated-hydrogen
case, `has_pyran_indicated_hydrogen_name`/`name_pyran_indicated_hydrogen`):
unlike furan/thiophene, RDKit does not treat this ring as aromatic at all
-- it Kekulizes to one sp3 ring carbon (2 H, no double bond) and four
other ring carbons each carrying one double bond, O itself always
single-bonded on both sides. O is fixed at locant 1 (P-2.txt ~522-529's
'2H-pyran'/'4H-pyran' worked examples, same numbering convention as
furan/thiophene); which locant needs the indicated 'H' depends on where
the sp3 carbon actually sits, so both ring directions from O are tried
and the lower locant wins (P-14.4), same tie-breaking principle as
`_ring_alignments` above. Confirmed against PubChem: 2H-pyran (CID
186148) and 4H-pyran (CID 136135). Unsubstituted only, and the S/Se/Te
analogues (thiopyran/selenopyran/telluropyran) are out of scope for this
function -- separate future task.

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
O+O 1,4-dioxane, O+S 1,4-oxathiane, S+S 1,4-dithiane, N+Se
selenomorpholine, N+Te telluromorpholine, O+Se 1,4-oxaselenane, O+Te
1,4-oxatellurane, S+Se 1,4-thiaselenane, S+Te 1,4-thiatellurane, Se+Se
1,4-diselenane, Te+Te 1,4-ditellurane -- Se+Te has no name registered
in PubChem and stays out of scope), the 5-membered
1,3-related pairs (N+N imidazolidine, N+O 1,3-oxazolidine, N+S
1,3-thiazolidine, O+O 1,3-dioxolane, O+S 1,3-oxathiolane, S+S
1,3-dithiolane, N+Se 1,3-selenazolidine, N+Te 1,3-tellurazolidine, O+Se
1,3-oxaselenolane, O+Te 1,3-oxatellurolane, S+Se 1,3-thiaselenolane, S+Te
1,3-thiatellurolane, Se+Se 1,3-diselenolane, Te+Te 1,3-ditellurolane --
Se+Te has no name registered in PubChem here either and stays out of
scope, same as the 1,4-ring), and the 5-membered 1,2-related pairs (N+N pyrazolidine,
N+O 1,2-oxazolidine, N+S 1,2-thiazolidine, O+O 1,2-dioxolane, O+S
1,2-oxathiolane, S+S 1,2-dithiolane, N+Se 1,2-selenazolidine, N+Te
1,2-tellurazolidine, O+Se 1,2-oxaselenolane, O+Te 1,2-oxatellurolane, S+Se
1,2-thiaselenolane, S+Te 1,2-thiatellurolane, Se+Se 1,2-diselenolane,
Te+Te 1,2-ditellurolane -- Se+Te not attempted, same reasoning as the
1,3-/1,4-rings' Se+Te gap) all have their own retained/
systematic name -- per Table 2.3 and P-22.2.2.1.2/.1.3, a pair with no
retained name (all but the N+N ones, which are simply irregular retained
names carrying no locants at all, 1,2- or 1,3- alike) always cites its
heteroatom locants; PubChem's own computed names are unreliable here
(they drop the locants for several of the O/S/Se/Te-only 1,2-pairs,
contradicted by Table 2.3's explicit '1,2-oxazolidine (PIN)'/
'1,2-thiazolidine (PIN)'/'1,2-selenazolidine (PIN)'/'1,2-tellurazolidine
(PIN)' entries and P-22.2.2.1.3's own '1,2-oxathiolane (PIN)' worked
example -- the N-containing pairs are the exception where PubChem's
computed name does happen to agree with the primary text), so the
1,2-pairs are sourced from the primary text directly rather than
from PubChem.

The 7-membered 1,4-related pairs (N+N 1,4-diazepane, N+O 1,4-oxazepane,
N+S 1,4-thiazepane, O+O 1,4-dioxepane, O+S 1,4-oxathiepane, S+S
1,4-dithiepane -- Se/Te not attempted, same reasoning as the 6-/5-membered
axes' Se+Te gap) are confirmed via PubChem's IUPACName for the exact
SMILES, all six matching the systematic Hantzsch-Widman-style pattern
directly (no retained irregular name here, unlike piperazine et al.).

The 7-membered 1,3-related pairs (N+N 1,3-diazepane, N+O 1,3-oxazepane,
N+S 1,3-thiazepane, O+O 1,3-dioxepane, O+S 1,3-oxathiepane, S+S
1,3-dithiepane -- Se/Te not attempted, same reasoning as the other axes'
Se+Te gap) are confirmed via PubChem's IUPACName the same way.

The 7-membered 1,2-related pairs (N+N 1,2-diazepane, N+O 1,2-oxazepane,
N+S 1,2-thiazepane, O+O 1,2-dioxepane, O+S 1,2-oxathiepane, S+S
1,2-dithiepane -- Se/Te not attempted, same reasoning as the other axes'
Se+Te gap): PubChem's raw IUPACName drops the locant for all six
('diazepane', 'oxazepane', etc.), but P-22.2.2.1.7 only omits
Hantzsch-Widman locants "if there is no ambiguity if locants are
omitted" -- and there is ambiguity here, since the 1,3-related pair above
already uses the identical bare stem (e.g. both 1,2-oxazepane and
1,3-oxazepane would collapse to plain 'oxazepane'). So, mirroring the
5-membered 1,2-axis's own PubChem-vs-primary-text correction above, this
axis keeps its locants explicit rather than trusting PubChem's bare
name.

Explicitly out of scope for the unsubstituted-only functions above: any
substituent, partially-saturated indicated-hydrogen forms other than the
ones listed above, three or more heteroatoms, heteroatoms other than
O/S/Se/Te/N, saturated two-heteroatom ring sizes/relationships other than
the four listed above, and ring sizes outside the single-heteroatom
tables above. `has_hetero_monocyclic_name`
returns False for all of these, so `core.py`'s existing dispatch (which
already rejects heteroatoms outside a few specific recognized shapes)
continues to raise `UnsupportedStructure` for them, unchanged.

Substituent naming (`has_hetero_monocyclic_substituent_name`/
`name_hetero_monocyclic_substituent`), one axis further than the
unsubstituted-only functions above:

- Scope: one or more substituents (a halogen, or a plain hydrocarbon
  group `name_branch` can name -- repeats and mixed kinds both allowed)
  on one of the 19 mancude parents listed in `_MANCUDE_NAME_SMILES`/
  `_TWO_HETEROATOM_MANCUDE_NAME_SMILES` above (single-heteroatom 5/6-
  membered plus two-heteroatom 5/6-membered, Se/Te analogues included)
  -- everything else about the parent (ring size, unsaturation,
  heteroatom set) stays exactly as recognized by the unsubstituted
  table. Saturated/partially saturated rings and polycyclic systems
  remain out of scope (a second ring anywhere -- including one folded
  into a substituent itself, like a cyclopropyl group -- is rejected via
  `mol.GetRingInfo().NumRings() == 1`). Every atom outside the ring must
  be carbon or a halogen (`_SUBSTITUENT_ATOMIC_NUMS`) -- `name_branch`'s
  chain-walking fallback doesn't actually check element types as it
  counts chain length, so an unvalidated branch containing e.g. an amine
  nitrogen or a carboxylic-acid oxygen would silently get counted as
  carbon and produce a wrong plain-alkyl name instead of being rejected
  (found via real-data testing: 'NC(CCc1ncccc1Cl)C(=O)O', an amino acid
  with a pyridine side chain, was misnamed
  '3-chloro-2-(3,4-dimethylpentyl)pyridine', silently dropping the amino
  and carboxylic acid groups). Any other heteroatom anywhere in the
  molecule makes this function return `None` so `core.py` falls through
  to a module that actually understands it.
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

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    non_single_bonds,
    ring_cycle,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, name_branch

_SUBSTITUENT_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}

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
    frozenset(("N", "Se")): ("selenomorpholine", "C1CNCC[Se]1"),
    frozenset(("N", "Te")): ("telluromorpholine", "C1CNCC[Te]1"),
    frozenset(("O", "Se")): ("1,4-oxaselenane", "C1COCC[Se]1"),
    frozenset(("O", "Te")): ("1,4-oxatellurane", "C1COCC[Te]1"),
    frozenset(("S", "Se")): ("1,4-thiaselenane", "C1CSCC[Se]1"),
    frozenset(("S", "Te")): ("1,4-thiatellurane", "C1CSCC[Te]1"),
    frozenset(("Se", "Se")): ("1,4-diselenane", "C1C[Se]CC[Se]1"),
    frozenset(("Te", "Te")): ("1,4-ditellurane", "C1C[Te]CC[Te]1"),
}
_FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("imidazolidine", "C1CNCN1"),
    frozenset(("N", "O")): ("1,3-oxazolidine", "C1CNCO1"),
    frozenset(("N", "S")): ("1,3-thiazolidine", "C1CSCN1"),
    frozenset(("O", "O")): ("1,3-dioxolane", "C1COCO1"),
    frozenset(("O", "S")): ("1,3-oxathiolane", "C1CSCO1"),
    frozenset(("S", "S")): ("1,3-dithiolane", "C1CSCS1"),
    frozenset(("N", "Se")): ("1,3-selenazolidine", "C1CNC[Se]1"),
    frozenset(("N", "Te")): ("1,3-tellurazolidine", "C1CNC[Te]1"),
    frozenset(("O", "Se")): ("1,3-oxaselenolane", "C1C[Se]CO1"),
    frozenset(("O", "Te")): ("1,3-oxatellurolane", "C1C[Te]CO1"),
    frozenset(("S", "Se")): ("1,3-thiaselenolane", "C1C[Se]CS1"),
    frozenset(("S", "Te")): ("1,3-thiatellurolane", "C1C[Te]CS1"),
    frozenset(("Se", "Se")): ("1,3-diselenolane", "C1C[Se]C[Se]1"),
    frozenset(("Te", "Te")): ("1,3-ditellurolane", "C1C[Te]C[Te]1"),
}
_FIVE_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("pyrazolidine", "C1CCNN1"),
    frozenset(("N", "O")): ("1,2-oxazolidine", "C1CCON1"),
    frozenset(("N", "S")): ("1,2-thiazolidine", "C1CCSN1"),
    frozenset(("O", "O")): ("1,2-dioxolane", "C1CCOO1"),
    frozenset(("O", "S")): ("1,2-oxathiolane", "C1CCOS1"),
    frozenset(("S", "S")): ("1,2-dithiolane", "C1CCSS1"),
    frozenset(("N", "Se")): ("1,2-selenazolidine", "C1CC[Se]N1"),
    frozenset(("N", "Te")): ("1,2-tellurazolidine", "C1CC[Te]N1"),
    frozenset(("O", "Se")): ("1,2-oxaselenolane", "C1CCO[Se]1"),
    frozenset(("O", "Te")): ("1,2-oxatellurolane", "C1CCO[Te]1"),
    frozenset(("S", "Se")): ("1,2-thiaselenolane", "C1CCS[Se]1"),
    frozenset(("S", "Te")): ("1,2-thiatellurolane", "C1CCS[Te]1"),
    frozenset(("Se", "Se")): ("1,2-diselenolane", "C1CC[Se][Se]1"),
    frozenset(("Te", "Te")): ("1,2-ditellurolane", "C1CC[Te][Te]1"),
}
_SEVEN_MEMBERED_1_4_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("1,4-diazepane", "C1CNCCNC1"),
    frozenset(("N", "O")): ("1,4-oxazepane", "C1CNCCOC1"),
    frozenset(("N", "S")): ("1,4-thiazepane", "C1CNCCSC1"),
    frozenset(("O", "O")): ("1,4-dioxepane", "C1COCCOC1"),
    frozenset(("O", "S")): ("1,4-oxathiepane", "C1COCCSC1"),
    frozenset(("S", "S")): ("1,4-dithiepane", "C1CSCCSC1"),
}
_SEVEN_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("1,3-diazepane", "N1CNCCCC1"),
    frozenset(("N", "O")): ("1,3-oxazepane", "O1CNCCCC1"),
    frozenset(("N", "S")): ("1,3-thiazepane", "S1CNCCCC1"),
    frozenset(("O", "O")): ("1,3-dioxepane", "O1COCCCC1"),
    frozenset(("O", "S")): ("1,3-oxathiepane", "C1CCCOCS1"),
    frozenset(("S", "S")): ("1,3-dithiepane", "S1CSCCCC1"),
}
_SEVEN_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES = {
    frozenset(("N", "N")): ("1,2-diazepane", "C1CCCCNN1"),
    frozenset(("N", "O")): ("1,2-oxazepane", "C1CCCCON1"),
    frozenset(("N", "S")): ("1,2-thiazepane", "C1CCCCSN1"),
    frozenset(("O", "O")): ("1,2-dioxepane", "C1CCCCOO1"),
    frozenset(("O", "S")): ("1,2-oxathiepane", "C1CCCCOS1"),
    frozenset(("S", "S")): ("1,2-dithiepane", "C1CCCCSS1"),
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
_THREE_OR_FOUR_HETEROATOM_MANCUDE_NAME_SMILES = {
    "1H-1,2,3-triazole": ("1H-1,2,3-triazole", "c1c[nH]nn1"),
    "2H-1,2,3-triazole": ("2H-1,2,3-triazole", "c1cn[nH]n1"),
    "1H-tetrazole": ("1H-tetrazole", "c1[nH]nnn1"),
    "2H-tetrazole": ("2H-tetrazole", "c1n[nH]nn1"),
}
_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name
    for name, smiles in (
        *_RETAINED_NAME_SMILES.values(),
        *_MANCUDE_NAME_SMILES.values(),
        *_TWO_HETEROATOM_SATURATED_NAME_SMILES.values(),
        *_FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.values(),
        *_FIVE_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES.values(),
        *_SEVEN_MEMBERED_1_4_TWO_HETEROATOM_NAME_SMILES.values(),
        *_SEVEN_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.values(),
        *_SEVEN_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES.values(),
        *_TWO_HETEROATOM_MANCUDE_NAME_SMILES.values(),
        *_THREE_OR_FOUR_HETEROATOM_MANCUDE_NAME_SMILES.values(),
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
    of the fourteen in scope (morpholine/piperazine/thiomorpholine/
    1,4-dioxane/1,4-oxathiane/1,4-dithiane, plus their Se/Te analogues
    selenomorpholine/telluromorpholine/1,4-oxaselenane/1,4-oxatellurane/
    1,4-thiaselenane/1,4-thiatellurane/1,4-diselenane/1,4-ditellurane --
    Se+Te has no PubChem-registered name and stays out of scope, and
    other ring sizes/relationships are out of scope too). Exposed for
    `_ketone.py`'s hetero-ring ketone naming, which needs the bare stem
    name rather than a full unsubstituted-molecule match."""
    entry = _TWO_HETEROATOM_SATURATED_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_five_membered_1_3_two_heteroatom_ring_name(elements):
    """The retained/systematic name for the unsubstituted, 5-membered,
    1,3-related two-heteroatom saturated ring whose heteroatom elements
    are `elements` (an (element, element) pair or frozenset, e.g.
    ('N', 'N') -> 'imidazolidine'), or None if that element pair isn't
    one of the fourteen in scope (imidazolidine/1,3-oxazolidine/
    1,3-thiazolidine/1,3-dioxolane/1,3-oxathiolane/1,3-dithiolane, plus
    their Se/Te analogues 1,3-selenazolidine/1,3-tellurazolidine/
    1,3-oxaselenolane/1,3-oxatellurolane/1,3-thiaselenolane/
    1,3-thiatellurolane/1,3-diselenolane/1,3-ditellurolane -- Se+Te itself
    has no name registered in PubChem and stays out of scope, same as the
    1,4-ring; other element pairs and other ring sizes have no retained
    name here and are out of scope; the 1,2-relationship has its own six, see
    `saturated_five_membered_1_2_two_heteroatom_ring_name`). Exposed for
    `_ketone.py`'s hetero-ring ketone naming, which needs the bare stem
    name rather than a full unsubstituted-molecule match."""
    entry = _FIVE_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_five_membered_1_2_two_heteroatom_ring_name(elements):
    """The retained/systematic name for the unsubstituted, 5-membered,
    1,2-related two-heteroatom saturated ring whose heteroatom elements
    are `elements` (an (element, element) pair or frozenset, e.g.
    ('N', 'N') -> 'pyrazolidine'), or None if that element pair isn't one
    of the fourteen in scope (pyrazolidine/1,2-oxazolidine/
    1,2-thiazolidine/1,2-dioxolane/1,2-oxathiolane/1,2-dithiolane, plus
    their Se/Te analogues 1,2-selenazolidine/1,2-tellurazolidine/
    1,2-oxaselenolane/1,2-oxatellurolane/1,2-thiaselenolane/
    1,2-thiatellurolane/1,2-diselenolane/1,2-ditellurolane -- Se+Te itself
    was not attempted, same reasoning as the 1,3-/1,4-rings' Se+Te gap;
    other element pairs and other ring sizes have no retained name here
    and are out of scope). Exposed for `_ketone.py`'s hetero-ring ketone
    naming, which needs the bare stem name rather than a full
    unsubstituted-molecule match."""
    entry = _FIVE_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_seven_membered_1_4_two_heteroatom_ring_name(elements):
    """The systematic name for the unsubstituted, 7-membered, 1,4-related
    two-heteroatom saturated ring whose heteroatom elements are `elements`
    (an (element, element) pair or frozenset, e.g. ('N', 'N') ->
    '1,4-diazepane'), or None if that element pair isn't one of the six in
    scope (1,4-diazepane/1,4-oxazepane/1,4-thiazepane/1,4-dioxepane/
    1,4-oxathiepane/1,4-dithiepane -- Se/Te analogues weren't attempted,
    same as the 6-membered/5-membered two-heteroatom axes; other element
    pairs, other ring sizes, and the 1,2-/1,3-relationships have no
    retained name here and are out of scope). Exposed for the same reason
    as `saturated_two_heteroatom_1_4_ring_name`, should a future ketone
    suffix or N-alkyl substituent PR need the bare stem name."""
    entry = _SEVEN_MEMBERED_1_4_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_seven_membered_1_3_two_heteroatom_ring_name(elements):
    """The systematic name for the unsubstituted, 7-membered, 1,3-related
    two-heteroatom saturated ring whose heteroatom elements are `elements`
    (an (element, element) pair or frozenset, e.g. ('N', 'N') ->
    '1,3-diazepane'), or None if that element pair isn't one of the six in
    scope (1,3-diazepane/1,3-oxazepane/1,3-thiazepane/1,3-dioxepane/
    1,3-oxathiepane/1,3-dithiepane -- Se/Te analogues weren't attempted,
    same as the other two-heteroatom axes; other element pairs, other ring
    sizes, and the 1,2-/1,4-relationships have no retained name here and
    are out of scope). Exposed for the same reason as
    `saturated_seven_membered_1_4_two_heteroatom_ring_name`, should a
    future ketone suffix or N-alkyl substituent PR need the bare stem
    name."""
    entry = _SEVEN_MEMBERED_1_3_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
    return entry[0] if entry else None


def saturated_seven_membered_1_2_two_heteroatom_ring_name(elements):
    """The systematic name for the unsubstituted, 7-membered, 1,2-related
    two-heteroatom saturated ring whose heteroatom elements are `elements`
    (an (element, element) pair or frozenset, e.g. ('N', 'N') ->
    '1,2-diazepane'), or None if that element pair isn't one of the six in
    scope (1,2-diazepane/1,2-oxazepane/1,2-thiazepane/1,2-dioxepane/
    1,2-oxathiepane/1,2-dithiepane -- Se/Te analogues weren't attempted,
    same as the other two-heteroatom axes; other element pairs, other ring
    sizes, and the 1,3-/1,4-relationships have no retained name here and
    are out of scope). Unlike PubChem's own computed name for these six
    (which drops the locant), this module keeps '1,2-' explicit -- see the
    module docstring's note on P-22.2.2.1.7 and the clash with the
    1,3-axis's identical bare stem. Exposed for the same reason as
    `saturated_seven_membered_1_4_two_heteroatom_ring_name`, should a
    future ketone suffix or N-alkyl substituent PR need the bare stem
    name."""
    entry = _SEVEN_MEMBERED_1_2_TWO_HETEROATOM_NAME_SMILES.get(frozenset(elements))
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

# imidazole/pyrazole's own N-H is a real prototropic tautomer -- for a
# *single* ring this doesn't create any actual naming ambiguity (see
# `_match_hetero_monocyclic_substituents`'s own role-sequence search,
# which derives the one structurally-valid alignment directly from the
# input molecule's real N-H position), but `_ring_assembly_chain.py`'s
# P-28.3 multi-ring assembly numbering has its own separate, still-open
# question of how per-ring tautomer choice interacts with an assembly-
# wide locant search, so that module excludes these two names from its
# own ring-parent table pending its own dedicated research.
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
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    return locant_set, citation_locants


def _match_hetero_monocyclic_substituents(mol):
    """(parent_name, grouped, graph) for a mancude ring carrying one or more
    substituents that matches one of `_ROLE_SEQUENCES`'s 19 parents; None if
    the molecule doesn't fit that shape at all (a second ring anywhere, a
    ring atom bearing more than one exocyclic branch, a ring size/
    heteroatom pattern outside the table, or a substituent sitting on a
    non-substitutable heteroatom)."""
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
    if any(atom.GetIdx() not in ring_set and atom.GetAtomicNum() not in _SUBSTITUENT_ATOMIC_NUMS for atom in mol.GetAtoms()):
        # `name_branch`'s chain-walking fallback (`_longest_chains_from_root`)
        # treats every non-halogen neighbor as chain-extending regardless of
        # element, so an exocyclic branch containing e.g. an amine nitrogen
        # or a carboxylic-acid oxygen would silently get counted as if it
        # were carbon and named as a plain alkyl substituent -- a wrong
        # name, not a raised rejection (confirmed via real-data testing:
        # 'NC(CCc1ncccc1Cl)C(=O)O' was misnamed
        # '3-chloro-2-(3,4-dimethylpentyl)pyridine', dropping the amino and
        # carboxylic acid groups entirely). Bail out to `None` so `core.py`
        # falls through to a module that actually understands the other
        # heteroatoms, rather than claim a shape this module can't safely
        # name.
        return None
    if any(a not in ring_set or b not in ring_set for a, b, _ in non_single_bonds(mol)):
        # Same blind spot as the atom-type check above, but for bond order:
        # `name_branch`'s chain-walking fallback doesn't check bond order
        # either, so a C=C in an exocyclic branch would silently be
        # counted as if it were a saturated chain (confirmed via real-data
        # testing: 'C=CCc1ccnc(Cl)c1Cl' was misnamed
        # '2,3-dichloro-4-propylpyridine', dropping the branch's own
        # double bond -- the correct 'prop-2-enyl' substituent never even
        # gets considered). Ring-internal bonds (aromatic, bond order 1.5)
        # are fine and excluded by the `a not in ring_set or b not in
        # ring_set` check.
        return None
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
        atom: name_branch(graph, root, atom, halogens, mol=mol) for atom, root in exo_by_atom.items()
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
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
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


def _match_pyran_indicated_hydrogen(mol):
    """The indicated-hydrogen locant for an unsubstituted 6-membered,
    one-oxygen, five-carbon mancude ring (2H-pyran/4H-pyran), or None if
    `mol` doesn't fit that shape at all."""
    if mol.GetNumAtoms() != 6:
        return None
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    if len(ring_atoms) != 6:
        return None
    if any(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms):
        return None
    oxygens = [atom for atom in ring_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() == 8]
    if len(oxygens) != 1 or mol.GetAtomWithIdx(oxygens[0]).GetTotalNumHs() != 0:
        return None
    o_atom = oxygens[0]
    carbons = [atom for atom in ring_atoms if atom != o_atom]
    if any(mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 for atom in carbons):
        return None
    sp3_carbons = [atom for atom in carbons if mol.GetAtomWithIdx(atom).GetTotalNumHs() == 2]
    if len(sp3_carbons) != 1:
        return None
    sp3_atom = sp3_carbons[0]
    if any(mol.GetAtomWithIdx(atom).GetTotalNumHs() != 1 for atom in carbons if atom != sp3_atom):
        return None

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, ring_atoms)
    start = ring_order.index(o_atom)
    ring_order = ring_order[start:] + ring_order[:start]
    forward_locant = ring_order.index(sp3_atom) + 1
    backward_locant = len(ring_order) - forward_locant + 2
    return min(forward_locant, backward_locant)


def has_pyran_indicated_hydrogen_name(mol) -> bool:
    return _match_pyran_indicated_hydrogen(mol) is not None


def name_pyran_indicated_hydrogen(mol) -> str:
    locant = _match_pyran_indicated_hydrogen(mol)
    return f"{locant}H-pyran"
