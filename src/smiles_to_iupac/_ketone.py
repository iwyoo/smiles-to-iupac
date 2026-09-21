"""Naming of ketones (the '-one' suffix, C=O with two carbon substituents) on
acyclic saturated or unsaturated carbon chains and on simple monocyclic
saturated rings, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-33.4, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'one' is the preselected suffix for a ketone carbonyl, cited in a combined
  chain/suffix name the same way 'ol'/'amine' are for -OH/-NH2
  (`_alcohol.py`/`_amine.py`), e.g. 'propan-2-one', 'cyclohexanone'. Table
  3.3 ranks 'one' senior to 'ol' and 'amine' but junior to the carboxylic
  acid/ester/amide/nitrile/aldehyde suffixes; a carbonyl carbon shaped like
  an aldehyde (only one carbon neighbor) or a carboxylic acid/ester/amide
  (a second oxygen on the same carbon) is rejected outright rather than
  silently named as if it were a ketone. A coexisting hydroxyl (-OH), being
  junior to 'one', is *not* rejected: it is cited as the 'hydroxy'
  substituent prefix instead (P-41), e.g. 'CC(=O)CCO' ->
  '4-hydroxybutan-2-one'.
- P-44.4.1.8 / P-45.2: suffix locants are minimized before 'ene'/'yne'
  locants, which are minimized before substituent-prefix locants — same
  ordering as `_alcohol.py`/`_amine.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  ketone suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- P-91.3/P-92: a molecule with one or
  more *specified* tetrahedral stereocenters -- every one on the principal
  chain/ring itself, no unspecified one alongside them, and no C=C/C#N
  double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix, ascending
  locant order, e.g. '(3R)-3-chloropentan-2-one',
  '(2R)-2-chlorocyclohexan-1-one', same pattern as `_carboxylic_acid.py`/
  `_aldehyde.py` (chain) and `_alcohol.py`'s `_name_cyclic_alcohol` (ring).
- P-31.1.3: a monocyclic ring bearing a ketone and exactly one C=C ring
  double bond -- e.g. 'cyclohex-2-en-1-one', 'cyclohex-3-en-1-one', both
  confirmed via PubChem PUG REST. The ketone's own locant is never
  omittable here even as the ring's sole substituent (P-44.4.1.8: suffix
  locant chosen first, but a competing ring double bond means '1' must
  still be cited, unlike the fully saturated 'cyclohexanone' case above);
  the numbering direction is then chosen to minimize the double bond's own
  locant, mirroring `_cyclic_unsaturated.py`'s `name_cyclic_unsaturated_yl`
  (which fixes its free valence at '1' the same way this fixes the
  ketone). This is deliberately the narrowest slice of the axis: a ring
  triple bond, and any other substituent alongside the ring double bond,
  are both still explicitly rejected pending further verification.

Unlike -OH/-NH2, a ketone carbonyl carbon always has exactly two carbon
neighbors and no hydrogens, so it can never be the sole substituent on a
mononuclear (P-14.3.4.2(a)) or homogeneous two-carbon (P-14.3.4.2(b)) chain
— those locant-omission special cases from `_alcohol.py` don't apply here
and are simply absent below; even the seemingly unambiguous 'propan-2-one'
(acetone) still cites its locant, confirmed against PubChem.

This module otherwise mirrors `_alcohol.py`'s scope restrictions:

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated carbonyl oxygen or a
  singly-bonded hydroxyl (ethers, and any other oxygen shape).
- A carbonyl carbon with fewer than two carbon neighbors (aldehyde) or an
  aromatic carbonyl carbon (aryl ketone) — a separate module's territory.
- Any other heteroatom (N, S, ...) — this module only resolves the 'one'
  vs. 'ol' seniority competition (Table 3.3) between a ketone carbonyl and a
  coexisting hydroxyl; only C, halogen, and ketone/hydroxyl oxygen atoms are
  accepted at all.
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- -one on a von Baeyer polycyclic or spiro skeleton — deferred, same as
  `_alcohol.py`.

A narrow extra path handles one or more ketone carbonyls on an otherwise
unsubstituted, saturated, single- or two-heteroatom (N/O/S) monocyclic
ring of size 5-7 (single-heteroatom) or the 1,4-related 6-membered
two-heteroatom rings with their own retained/systematic name (morpholine
N+O, piperazine N+N, thiomorpholine N+S, 1,4-dioxane O+O, 1,4-oxathiane
O+S, 1,4-dithiane S+S -- `_hetero_monocyclic.py`'s
`saturated_ring_name`/`saturated_two_heteroatom_1_4_ring_name`) -- the
ring heteroatom(s) always fixed at locant 1 (only their own H, if any,
may be present -- no other exocyclic substituent), and the ketone locant
set minimized over the possible numbering directions. For two different
heteroatoms, the higher-priority one (P-22.2.1 element seniority O > S >
N, `_TWO_HETERO_PRIORITY`) always takes locant 1; for two identical
heteroatoms, either may be locant 1 so both are tried. A ketone carbonyl
directly bonded to the ring heteroatom (a lactam/lactone, e.g.
piperidin-2-one, 1,4-oxathian-3-one) fits this same shape and is named
correctly here; `core.py` routes `has_hetero_ring_ketone_shape` ahead of
both `has_amide_shape` and `has_ester_shape` so this module claims it
before either of those modules' ring-always-out-of-scope guards would. A
symmetric, unsubstituted cyclic imide/dione (e.g. succinimide,
`O=C1CCC(=O)N1`) fits the same shape too and is named here as a plain
ring dione rather than via `_imide.py`'s acyclic-only "N-acyl amide"
construction (P-66.6.3), which was never meant to cover the cyclic case.
This path is deliberately narrow: 3+ heteroatom rings, Se/Te heteroatoms,
mancude (aromatic) rings, 3-/4-membered rings, non-1,4- two-heteroatom
relationships, and any substituent other than the ring heteroatom's own
indicated hydrogen and the ketone carbonyl(s) themselves (alkyl, halogen,
hydroxyl, ...) are all out of scope -- these fall through to this
module's existing carbocyclic-only validation, which raises its own
(more general) error.

PubChem-confirmed: `O=C1CCNCC1` -> 'piperidin-4-one' (CID 33721),
`O=C1CCCNC1` -> 'piperidin-3-one' (CID 33722), `O=C1CCOCC1` ->
'oxan-4-one' (CID 121599), `O=C1CCSCC1` -> 'thian-4-one' (CID 66173),
`O=C1CCCCN1` -> 'piperidin-2-one' (CID 12665), `O=C1CCCN1` ->
'pyrrolidin-2-one' (CID 12025), `O=C1CCCCCN1` -> 'azepan-2-one' (CID
7768), `O=C1CCC(=O)N1` -> 'pyrrolidine-2,5-dione' (CID 11439),
`O=C1COCCN1` -> 'morpholin-3-one' (CID 66953), `O=C1CNCCN1` ->
'piperazin-2-one' (CID 231360), `O=C1CNC(=O)CN1` -> 'piperazine-2,5-dione'
(CID 7817), `O=C1CSCCN1` -> 'thiomorpholin-3-one' (CID 88402),
`O=C1COCCO1` -> '1,4-dioxan-2-one' (CID 18233), `O=C1COC(=O)CO1` ->
'1,4-dioxane-2,5-dione' (CID 65432), `O=C1COCCS1` -> '1,4-oxathian-3-one'
(CID 15238324), `O=C1CSCCS1` -> '1,4-dithian-2-one' (CID 542724),
`O=C1CSC(=O)CS1` -> '1,4-dithiane-2,5-dione' (CID 319007). The element
seniority itself (O > S > N) is confirmed via PubChem's own
'4-methylmorpholine'/'4-methylthiomorpholine' (N always at locant 4) and
'1,4-oxathian-3-one' (O at locant 1, S at locant 4).

The 7-membered analogue of this same 1,4-two-heteroatom shape
(1,4-diazepane/1,4-oxazepane/1,4-thiazepane/1,4-dioxepane/
1,4-oxathiepane/1,4-dithiepane -- `_hetero_ring_seven_membered_1_4`) is
handled separately, since its two heteroatom-to-heteroatom arcs differ in
length (2 carbons vs. 3) rather than being equal as in the 6-membered
case -- see `_seven_membered_1_4_numbering`. PubChem-confirmed:
`C1CNCCNC1=O` -> '1,4-diazepan-5-one' (CID 2737264), `C1CNCC(=O)NC1` ->
'1,4-diazepan-2-one' (CID 13428532), `C1COCCNC1=O` ->
'1,4-oxazepan-5-one' (CID 7023020), `C1CSCCNC1=O` ->
'1,4-thiazepan-5-one' (CID 286337), `C1COCC(=O)OC1` ->
'1,4-dioxepan-2-one' (CID 10219408), `C1COC(=O)CSC1` ->
'1,4-oxathiepan-2-one' (CID 12391259). The dione form and N-/ring-carbon
alkyl substituents are not attempted here, unlike the 6-membered case --
no PubChem-registered example was found to confirm either.

A separate, narrower path handles one or two ketone carbonyls on a
five-membered 1,3-related saturated ring (N+N/N+O/N+S/O+O/O+S/S+S,
`_hetero_ring_five_membered_1_3`) -- one carbonyl always sits between the
two heteroatoms (locant 2, the sole ring atom on the short arc, fixed
regardless of numbering direction); a second, optional carbonyl sits on
the long arc (locant 4 or 5). When the two heteroatoms differ, P-22.2.1's
element seniority (O > S > N, `_TWO_HETERO_PRIORITY`) fixes which one is
locant 1 -- this outranks locant minimization entirely, so the second
carbonyl's locant is then whatever that forced numbering gives it, not
whichever is lower (confirmed by contrast: `O=C1OC(=O)CN1` names as
'1,3-oxazolidine-2,5-dione', not the lower-locant '...-2,4-dione' that
naive minimization would produce, because O must stay locant 1). When the
two heteroatoms are identical, either may be locant 1, so both are tried
and the lower ketone locant set wins (e.g. imidazolidine-2,4-dione).
PubChem-confirmed unsubstituted/single-ketone forms: `O=C1NCCN1` ->
'imidazolidin-2-one' (CID 8453), `O=C1OCCN1` -> '1,3-oxazolidin-2-one'
(CID 73949), `O=C1SCCN1` -> '1,3-thiazolidin-2-one' (CID 97431),
`O=C1OCCO1` -> '1,3-dioxolan-2-one' (CID 7303), `O=C1SCCO1` ->
'1,3-oxathiolan-2-one' (CID 72822), `O=C1SCCS1` -> '1,3-dithiolan-2-one'
(CID 123140). PubChem-confirmed dione forms: `O=C1NC(=O)CN1` ->
'imidazolidine-2,4-dione' (CID 10006), `O=C1OC(=O)CN1` ->
'1,3-oxazolidine-2,5-dione' (CID 75136), `O=C1SC(=O)CN1` ->
'1,3-thiazolidine-2,5-dione' (CID 542718), `O=C1OC(=O)CO1` ->
'1,3-dioxolane-2,4-dione' (CID 12793796), `O=C1SC(=O)CO1` ->
'1,3-oxathiolane-2,4-dione' (CID 67415624), `O=C1SC(=O)CS1` ->
'1,3-dithiolane-2,4-dione' (CID 637793).

Either ring nitrogen of this five-membered shape (single-ketone or dione
alike) may also carry a single plain, unbranched, unsubstituted alkyl
substituent -- its own locant is whatever the already-fixed ketone-locant
numbering above gives it, never separately minimized (P-44.4.1.8: suffix
locants are decided first and win outright). PubChem-confirmed:
`O=C1N(C)C(=O)CN1` -> '3-methylimidazolidine-2,4-dione' (CID 138851),
`O=C1NC(=O)CN1C` -> '1-methylimidazolidine-2,4-dione' (CID 69217),
`O=C1N(C)C(=O)CN1C` -> '1,3-dimethylimidazolidine-2,4-dione' (CID
123410).

For the dione (hydantoin) form specifically, the sole remaining plain
ring carbon may likewise carry a single plain, unbranched, unsubstituted
alkyl substituent -- its locant is, again, whatever the already-fixed
ketone numbering gives it. The single-ketone form is excluded: two plain
ring carbons remain there, and which one carries the substituent would be
a new locant tie-break input this module doesn't implement.
PubChem-confirmed: `O=C1NC(=O)C(C)N1` -> '5-methylimidazolidine-2,4-dione'
(CID 69216), `O=C1OC(=O)C(C)N1` -> '4-methyl-1,3-oxazolidine-2,5-dione'
(CID 70938), `O=C1OC(=O)C(C)O1` -> '5-methyl-1,3-dioxolane-2,4-dione'
(CID 22227598).

The 7-membered analogue of this same 1,3-two-heteroatom shape
(1,3-diazepane/1,3-oxazepane/1,3-thiazepane/1,3-dioxepane/
1,3-oxathiepane/1,3-dithiepane -- `_hetero_ring_seven_membered_1_3`) is
handled separately and much more narrowly: only a single ketone directly
on the bridging carbon (locant 2, same as the 5-membered case, since the
bridge is still the sole ring atom between the two heteroatoms on the
short arc) is supported -- no PubChem-registered example was found for a
second ketone on the 7-membered ring's longer (4-carbon) far arc, or for
a heteroatom/ring-carbon alkyl substituent, so both stay out of scope
(unlike the 5-membered case above). PubChem-confirmed: `C1CCNC(=O)NC1`
-> '1,3-diazepan-2-one' (CID 29396), `C1CCOC(=O)NC1` ->
'1,3-oxazepan-2-one' (CID 12218833), `C1CCSC(=O)NC1` ->
'1,3-thiazepan-2-one' (CID 20253446), `C1CCOC(=O)OC1` ->
'1,3-dioxepan-2-one' (CID 10197665), `C1CCSC(=O)OC1` ->
'1,3-oxathiepan-2-one' (CID 59128151), `C1CCSC(=O)SC1` ->
'1,3-dithiepan-2-one' (CID 59128158).

A parallel path handles the five-membered 1,2-related two-heteroatom
shape (the two heteroatoms directly bonded, locants 1/2 fixed rather than
separated by a bridging carbon; the ketone lands on one of the three
remaining ring carbons, locants 3/4/5) -- otherwise the same mechanics as
the 1,3-case above (element seniority forces locant 1 when the pair
differs, minimization decides it when identical; N-substituents allowed).
PubChem-confirmed for the N-containing pairs: `O=C1CCNN1` ->
'pyrazolidin-3-one' (CID 151497), `O=C1CCON1` -> '1,2-oxazolidin-3-one'
(CID 192737), `O=C1CCSN1` -> '1,2-thiazolidin-3-one' (CID 21878697). For
the O/S-only pairs, PubChem's own computed name drops the '1,2-' locant
(e.g. 'dioxolan-3-one'), contradicted by the same Table 2.3/
P-22.2.2.1.3 primary-text evidence already used for the unsubstituted
parent names (`saturated_five_membered_1_2_two_heteroatom_ring_name`) --
so `O=C1CCOO1`/`O=C1CCOS1`/`O=C1CCSS1` are named
'1,2-dioxolan-3-one'/'1,2-oxathiolan-3-one'/'1,2-dithiolan-3-one' here.

The dione form of this 1,2-shape likewise allows a single plain,
unbranched, unsubstituted alkyl substituent on its sole remaining plain
ring carbon, mirroring the 1,3-case's identical hydantoin rule above --
the single-ketone form stays excluded for the same reason (two plain
ring carbons remain, and picking which one carries the substituent would
need a new locant tie-break rule this module doesn't implement).
PubChem-confirmed: `O=C1C(C)C(=O)NN1` -> '4-methylpyrazolidine-3,5-dione'
(CID 12391681).

The 7-membered analogue of this same 1,2-two-heteroatom shape
(1,2-diazepane/1,2-oxazepane/1,2-thiazepane/1,2-dioxepane/
1,2-oxathiepane/1,2-dithiepane -- `_hetero_ring_seven_membered_1_2`) is
handled separately and, like the 7-membered 1,3-case above, much more
narrowly: only the single-ketone form (locant 3, the ring carbon
immediately after the second heteroatom) is supported -- no
PubChem-registered example was found for a second ketone on the longer
(4-carbon) far arc or for a ring-carbon alkyl substituent. Unlike the
5-membered case, PubChem's own computed name drops the '1,2-' locant for
every one of the six element pairs here, including the N-containing ones
(not just O/S-only) -- '1,2-' is kept anyway, for the same reason the
unsubstituted parent name already does
(`saturated_seven_membered_1_2_two_heteroatom_ring_name`: P-22.2.2.1.7 and
a stem clash with the 1,3-axis). PubChem-confirmed by connectivity:
`N1NC(=O)CCCC1` -> '1,2-diazepan-3-one' (CID 21962094), `O1NC(=O)CCCC1`
-> '1,2-oxazepan-3-one' (CID 12116845), `S1NC(=O)CCCC1` ->
'1,2-thiazepan-3-one' (CID 22346532), `O1OC(=O)CCCC1` ->
'1,2-dioxepan-3-one' (CID 17988591), `O1SC(=O)CCCC1` ->
'1,2-oxathiepan-3-one' (CID 123366628), `S1SC(=O)CCCC1` ->
'1,2-dithiepan-3-one' (CID 129723503).

Unlike -OH/-NH2, a ketone carbon can never itself also be a C=C/C#C alkene
carbon (its two remaining bonds, after the C=O double bond, are already
committed to its two required carbon substituents — a ketone carbon with a
third, double-bonded C=C neighbor would have five bonds), so there is no
enol-analogous "enone" case to scope out here: an alpha,beta-unsaturated
ketone like 'CC(=O)C=CC' (the double bond adjacent to, but not on, the
carbonyl carbon) is a perfectly nameable 'pent-3-en-2-one' and is supported
via the same suffix-locant-priority mechanism as `_alcohol.py`'s
'pent-4-en-1-ol'.

A further narrow path (`_name_phenyl_chain_ketone`) allows the sole ketone
to sit on an otherwise plain, unbranched, saturated chain hanging off one
atom of an unsubstituted benzene ring -- the ring is then cited as a
'phenyl' substituent prefix on the chain (P-44/P-52.2.8), mirroring
`_carboxylic_acid.py`'s identical first slice
(`_name_phenyl_chain_carboxylic_acid`). A carbonyl carbon directly
attached to the ring (an aryl ketone) still stays out of scope.
PubChem-confirmed: `CC(=O)Cc1ccccc1` -> '1-phenylpropan-2-one' (CID 7678),
`CC(=O)CCc1ccccc1` -> '4-phenylbutan-2-one' (CID 17355).

The same chain-substituent shape also extends to a simple heteroaromatic
monocycle (pyridine/pyrrole/furan/thiophene, P-29.3.4.1) in place of
benzene (P-616 M2 step 7) -- e.g. `CC(=O)Cc1cccnc1` ->
'1-(pyridin-3-yl)propan-2-one' (PubChem CID 238414), `CC(=O)Cc1ccco1` ->
'1-(furan-2-yl)propan-2-one' (CID 228583), `CC(=O)Cc1cccs1` ->
'1-(thiophen-2-yl)propan-2-one' (CID 529394), `CC(=O)Cc1cc[nH]c1` ->
'1-(1H-pyrrol-3-yl)propan-2-one' (CID 18383973), mirroring `_thiol.py`'s/
`_selenol.py`'s/`_tellurol.py`'s identical generalization. A ketone
carbon directly attached to a heteroaromatic ring stays out of scope,
same as the benzene case.
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
    heteroaromatic_monocycle_name,
    is_plain_benzene_ring,
    linear_branch,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ordered_chain,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereo_elements,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    two_separate_rings_with_plain_aromatic_substituent,
)
from ._hetero_monocyclic import (
    saturated_five_membered_1_2_two_heteroatom_ring_name,
    saturated_five_membered_1_3_two_heteroatom_ring_name,
    saturated_ring_name,
    saturated_seven_membered_1_2_two_heteroatom_ring_name,
    saturated_seven_membered_1_3_two_heteroatom_ring_name,
    saturated_seven_membered_1_4_two_heteroatom_ring_name,
    saturated_two_heteroatom_1_4_ring_name,
)
from ._numerals import alkyl_name
from ._substituents import (
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
    ring_branch_stereo_display,
    substituents_for_chain,
    substituents_for_ring,
)

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}
_HETERO_RING_ELEMENTS = {7: "N", 8: "O", 16: "S"}
_HETERO_RING_SIZES = (5, 6, 7)
_TWO_HETERO_RING_ELEMENT_PAIRS = {
    frozenset(("N", "O")),
    frozenset(("N", "N")),
    frozenset(("N", "S")),
    frozenset(("O", "O")),
    frozenset(("O", "S")),
    frozenset(("S", "S")),
}
_TWO_HETERO_RING_SIZE = 6
# P-22.2.1 element seniority for locant 1 (Table 2.8's replacement-nomenclature
# order O > S > N, confirmed via PubChem's own 'oxathian-3-one' citing O at
# locant 1 and S at locant 4): the strictly lower value always wins, so O
# outranks S which outranks N. Se/Te are deliberately absent here -- the
# unsubstituted-ring names for Se/Te-containing 1,4-pairs are supported
# (`_hetero_monocyclic.py`), but PubChem has no registered name for any
# ketone on those rings to confirm the locant this priority table would
# pick, so the ketone axis for Se/Te pairs stays out of scope.
_TWO_HETERO_PRIORITY = {"O": 0, "S": 1, "N": 2}
_FIVE_MEMBERED_1_3_RING_ELEMENT_PAIRS = {
    frozenset(("N", "N")),
    frozenset(("N", "O")),
    frozenset(("N", "S")),
    frozenset(("O", "O")),
    frozenset(("O", "S")),
    frozenset(("S", "S")),
}
_FIVE_MEMBERED_1_2_RING_ELEMENT_PAIRS = {
    frozenset(("N", "N")),
    frozenset(("N", "O")),
    frozenset(("N", "S")),
    frozenset(("O", "O")),
    frozenset(("O", "S")),
    frozenset(("S", "S")),
}
# Same six N/O/S pairs as the 6-membered 1,4-axis above -- Se/Te stay out of
# scope for the same reason (no PubChem-registered ketone name to confirm
# the locant either priority table or numbering direction would pick).
_SEVEN_MEMBERED_1_4_RING_ELEMENT_PAIRS = _TWO_HETERO_RING_ELEMENT_PAIRS
_SEVEN_MEMBERED_1_4_RING_SIZE = 7
# Same six N/O/S pairs as the 5-membered 1,3-axis above -- Se/Te stay out of
# scope for the same reason (no PubChem-registered ketone name found).
_SEVEN_MEMBERED_1_3_RING_ELEMENT_PAIRS = _FIVE_MEMBERED_1_3_RING_ELEMENT_PAIRS
_SEVEN_MEMBERED_1_3_RING_SIZE = 7
# Same six N/O/S pairs as the 5-membered 1,2-axis below -- Se/Te stay out of
# scope for the same reason (no PubChem-registered ketone name found).
_SEVEN_MEMBERED_1_2_RING_ELEMENT_PAIRS = _FIVE_MEMBERED_1_2_RING_ELEMENT_PAIRS
_SEVEN_MEMBERED_1_2_RING_SIZE = 7


def _validate_and_collect_ketones(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (ketones, hydroxyls): the set of carbonyl-oxygen atom indices,
    and the set of any coexisting hydroxyl-oxygen atom indices. A hydroxyl is
    junior to 'one' in Table 3.3's suffix seniority order, so it is cited as
    the 'hydroxy' substituent prefix instead of competing for the suffix
    (P-41).

    `aromatic_ring_atoms`: atom indices already independently verified
    elsewhere as a single plain benzene ring or heteroaromatic monocycle
    (pyridine/furan/thiophene/pyrrole) cited as a substituent prefix (see
    `_name_phenyl_chain_ketone`, and `name_ketone`'s `two_separate_rings_
    with_plain_aromatic_substituent` shape, #622) -- exempted here
    wholesale from the per-atomic-number checks below (already
    independently verified by that shape check itself) so those callers
    can still share this one validator for everything outside the ring."""
    ketones = set()
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in aromatic_ring_atoms:
            if atomic_num == 6:
                has_carbon = True
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a ketone carbonyl oxygen (P-33.4) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; only an isolated ketone "
                    "carbonyl or hydroxyl is supported (Table 3.3, P-33.4)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a ketone/hydroxyl oxygen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() == 1.0:
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        "an oxygen that isn't a carbonyl (=O) or hydroxyl "
                        "(-OH) is out of scope for this module"
                    )
                hydroxyls.add(atom.GetIdx())
                continue
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring (an aryl ketone) is out "
                    "of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a carbonyl carbon with fewer than two carbon neighbors "
                    "(an aldehyde or terminal carbonyl) is a more senior "
                    "characteristic group than a plain ketone (Table 3.3), "
                    "which this module does not attempt to disambiguate"
                )
            ketones.add(atom.GetIdx())
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
    if not ketones:
        raise UnsupportedStructure("no ketone (C=O) group found; this module only handles ketones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return ketones, hydroxyls


def _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(len(one_locants), "one"), one_locants
    )


def _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    one_locant_set = lowest_locant_set(one_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, one_locants, ene_locants, yne_locants, grouped)
    return (
        (
            one_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _one_locants(position_of, ketones, graph):
    locants = []
    for o in ketones:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants

def _best_acyclic_ketone_candidate(
    mol, ketones, hydroxyls, bonds, stereo=None, extra_names=None, required_atoms=frozenset(), carbon_graph=None
):
    """(best_name, best_position_of) -- the winning chain numbering and its
    fully-formatted name, factored out of `_name_acyclic_ketone` (which
    just adds the stereo-descriptor prefix on top) so `_isotope_ketone.py`
    can reuse the identical numbering decision: P-82.5.1 states the
    presence of isotopic nuclides is "considered last" among numbering
    criteria, so an isotope label never changes which candidate wins here,
    and that module only needs `best_position_of` to independently compute
    its own isotope descriptor's locant.

    `extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group demoted to a substituent prefix by
    `_seniority.senior_class` (e.g. a demoted amine's 'amino'), reused by
    `_coexisting_groups.py` so a pairwise module doesn't have to
    reimplement this function's chain search/candidate selection.
    `required_atoms`: additional carbon atoms (e.g. every demoted amine's
    own carbon neighbor) that a candidate chain must also carry -- both
    empty/None by default so existing callers are unaffected.
    `carbon_graph`: the carbon-only graph to search for the principal
    chain -- defaults to `carbon_adjacency(mol)` (unchanged behavior);
    `_ether_ketone.py` passes one with a coexisting ether's alkoxy-branch
    component already removed, mirroring `_thiol.py`'s identical
    `carbon_graph` parameter (PR #428) and for the same reason: an ether
    oxygen isn't itself a carbon, so its alkoxy branch would otherwise
    form a separate component that could wrongly outrank the real
    ketone-bearing chain in `longest_chains`' global-diameter search."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}, **(extra_names or {})}
    chains = longest_chains(carbon_graph if carbon_graph is not None else carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [idx for kind, idx, _ in stereo if kind == "atom"] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _one_locants(position_of, ketones, graph) is None:
            continue
        if not required_atoms <= set(chain):
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _one_locants({a: i + 1 for i, a in enumerate(c)}, ketones, graph) is not None
            and required_atoms <= set(c)
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every ketone-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more C=O groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = substituents_for_chain(graph, candidate, halogens, ketones, mol=mol)
            key, name = _candidate_key(chain_length, one_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    return best_name, best_position_of


def _name_acyclic_ketone(
    mol, ketones, hydroxyls, bonds, stereo=None, extra_names=None, required_atoms=frozenset(), carbon_graph=None
):
    """`stereo`: None, or a list of ("atom"/"bond", idx, "R"/"S"/"E"/"Z")
    from `specified_stereo_elements`/`specified_stereocenters` -- if given,
    only chain candidates that include every tetrahedral stereocenter are
    eligible (P-92: a stereocenter on a substituent branch rather than the
    principal chain is out of scope, mirroring
    `_carboxylic_acid.py`/`_aldehyde.py`'s identical treatment; a
    double-bond E/Z element's atoms are already required to lie on the
    chain via `bonds`, so no separate check is needed for those), and the
    winning candidate's own locants are used to format a
    "(<locant><R/S/E/Z>,...)-" prefix onto the name, ascending locant order
    (P-91.3, including when both kinds coexist -- only reachable from the
    acyclic caller, see `name_ketone`). See `_best_acyclic_ketone_candidate`
    for `extra_names`/`required_atoms`/`carbon_graph`."""
    best_name, best_position_of = _best_acyclic_ketone_candidate(
        mol, ketones, hydroxyls, bonds, stereo, extra_names, required_atoms, carbon_graph
    )
    if stereo is not None:
        labels = []
        for kind, idx, code in stereo:
            if kind == "atom":
                locant = best_position_of[idx]
            else:
                bond = mol.GetBondWithIdx(idx)
                locant = min(best_position_of[bond.GetBeginAtomIdx()], best_position_of[bond.GetEndAtomIdx()])
            labels.append((locant, code))
        labels.sort()
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_ketone(mol, ring_atoms):
    """Name a ketone whose C=O lies entirely on a single unbranched chain
    hanging off one atom of an otherwise-plain, unsubstituted benzene ring
    -- e.g. 1-phenylpropan-2-one. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`. Narrower
    than the acyclic path above: exactly one ketone, no coexisting
    standalone hydroxyl, no chain unsaturation, and no specified
    stereocenter -- each is a separate follow-up rather than being combined
    with the ring case in this first slice. A ketone carbon directly
    attached to the ring (an aryl ketone) stays out of scope, same as the
    module's existing acyclic-carbonyl check above."""
    ketones, hydroxyls = _validate_and_collect_ketones(mol, aromatic_ring_atoms=ring_atoms)
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a benzene-ring-substituent "
            "ketone chain is not supported yet"
        )
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "more than one ketone group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "ketone chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in ketones
        and b[1] not in ketones
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "ketone chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain ketone is not "
            "supported yet"
        )
    ring_atom, chain_root = attachment
    (ketone_oxygen,) = ketones
    (ketone_carbon,) = graph[ketone_oxygen]
    if ketone_carbon == chain_root:
        raise UnsupportedStructure(
            "a carbonyl carbon directly attached to the benzene ring (an "
            "aryl ketone) is out of scope for this module (see the "
            "separate aromatic-ring module)"
        )

    chain, branches = longest_branched_chain_through(graph, ketone_carbon, ring_atoms, ketones, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        one_locants = _one_locants(position_of, ketones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, one_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_ketone(mol, ketones):
    """Name a ketone whose C=O lies entirely on a single branched chain
    hanging off one atom of an otherwise-plain saturated monocyclic ring
    (the ring itself bears no ketone) -- e.g. 1-cyclohexylethanone. The
    ring is cited as a "cyclo..." substituent prefix (P-29.3.3) on the
    chain, which is the parent hydride, mirroring `_name_phenyl_chain_
    ketone` above and `_alcohol.py`'s `_name_ring_substituent_chain_
    alcohol`. Narrower than the benzene-ring case: exactly one ketone, no
    coexisting standalone hydroxyl -- each a separate follow-up."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, ketones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    (ketone_oxygen,) = ketones
    (ketone_carbon,) = graph[ketone_oxygen]

    chain, branches = longest_branched_chain_through(graph, ketone_carbon, ring_atoms, ketones, halogens=halogen_substituents(mol))
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
        one_locants = _one_locants(position_of, ketones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, one_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_ketone_chain_ketone(mol, ketones):
    """Name a ketone compound where the ring itself bears at least as many
    ketones as a single unbranched chain hanging off exactly one ring atom
    does (P-44.1.1: the candidate with the greater count of the principal
    characteristic group 'one' is senior; P-44.1.2.2 resolves an exact tie
    in the ring's favor) -- e.g. 2-(2-oxopropyl)cyclohexan-1-one. The ring
    is the parent; the chain is cited as an '(oxo...alkyl)' substituent
    prefix, mirroring `_alcohol.py`'s `_name_ring_with_hydroxy_chain_
    alcohol` and `_amine.py`'s `_name_ring_with_amine_chain_amine`. A
    ring's own carbonyl carbon can never be the chain-attachment atom
    (it's sp2 with two ring bonds plus the C=O, leaving no room for a
    fifth, exocyclic bond), so no aryl-ketone-style exception is needed
    here."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, ketones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, ketones)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_ketones = {o for o in ketones if next(iter(graph[o])) in chain_set}
    ring_ketones = ketones - chain_ketones
    if len(ring_ketones) < len(chain_ketones):
        # P-44.1.1: the chain captures strictly more ketones, so it's the
        # senior parent and the ring (with its own one or more ketones)
        # is cited as a substituent instead -- mirrors
        # `_name_ring_substituent_chain_ketone` exactly, substituting
        # the ring's own name_branch-computed name for the plain
        # "cyclo..." one that function uses.
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{o: "oxo" for o in ring_ketones}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, chain_ketones, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, one_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{o: "oxo" for o in chain_ketones}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ring_ketones, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, one_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def _ring_name_from_substituents(ring_size, one_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size, ene_locants, yne_locants, prefix, total_subs, multiplied_word(len(one_locants), "one"), one_locants
    )


def _ring_candidate_key(ring_size, one_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    one_locant_set = lowest_locant_set(one_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, one_locants, ene_locants, yne_locants, grouped)
    return one_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, ketones, stereo, halogens, mol=None, aromatic_atoms=frozenset()):
    return ring_branch_stereo_display(graph, ring_order, ketones, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms)


def _name_cyclic_ketone(mol, ketones, hydroxyls, stereo=None, bonds=(), ring_atoms=None, aromatic_atoms=frozenset()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_alcohol.py`'s `_name_cyclic_alcohol`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_alcohol.py`'s own case):
    exactly one stereocenter on the ring's sole substituent branch instead
    embeds a bracketed descriptor into that substituent's own name, in
    place of the usual ring-locant prefix.

    `ring_atoms`/`aromatic_atoms`: when the ketone-bearing ring reaches
    this function as one half of a `two_separate_rings_with_plain_
    aromatic_substituent` shape (P-25 M1 step 2, #622), the caller passes
    the ketone-bearing ring's own atoms explicitly (RDKit's SSSR would
    otherwise list either of the two disjoint rings first) along with the
    other, aromatic ring's atoms, cited as a plain substituent (phenyl or
    a heteroaromatic monocycle) via `name_branch` the same way an ordinary
    alkyl ring substituent already is. `ring_atoms=None` (default)
    preserves the original single-ring dispatch unchanged."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    if ring_atoms is None:
        ring_atoms = list(mol.GetRingInfo().AtomRings()[0])
    else:
        ring_atoms = list(ring_atoms)
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(
            graph, ring_order, ketones, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms
        )
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(
        substituents_for_ring(graph, ring_order, halogens, ketones, mol=mol, aromatic_atoms=aromatic_atoms).values()
    ):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "ketone is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            if one_locants is None:
                raise UnsupportedStructure(
                    "a ketone not on the ring itself (e.g. on a substituent "
                    "branch) is not supported yet"
                )
            substituents = substituents_for_ring(
                graph, candidate, halogens, ketones, mol=mol, aromatic_atoms=aromatic_atoms
            )
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, one_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _hetero_ring_heteroatom(mol):
    """The ring-atom index of the sole non-carbon heteroatom in this
    molecule's one ring, if it looks like a saturated single-heteroatom
    monocyclic ketone shape (N/O/S in a 5-, 6-, or 7-membered ring) --
    None if it doesn't match that shape at all, in which case the caller
    falls through to the existing carbocyclic-only path (which raises its
    own, more general error for whatever doesn't fit)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) not in _HETERO_RING_SIZES:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    (heteroatom,) = heteroatoms
    if mol.GetAtomWithIdx(heteroatom).GetAtomicNum() not in _HETERO_RING_ELEMENTS:
        return None
    return heteroatom


def _linear_alkyl_substituent(mol, full_carbon_graph, graph, root, attachment_point):
    """alkyl_name(length) for the plain, unbranched, unsubstituted alkyl
    chain hanging off `root` (reached from `attachment_point`, a ring
    atom), or raise `UnsupportedStructure` if it's branched, unsaturated,
    or itself carries any further substituent -- shared by the ring
    heteroatom (N-alkyl) and ring carbon substituent paths in
    `_validate_and_collect_hetero_ring_ketone` below, which only differ in
    which atom `attachment_point` is."""
    length = linear_branch(full_carbon_graph, root, attachment_point)
    if length is None:
        raise UnsupportedStructure("a branched substituent is not supported yet")
    atoms = set()
    previous, current = attachment_point, root
    while current is not None:
        atoms.add(current)
        neighbors = [n for n in full_carbon_graph[current] if n != previous]
        previous, current = current, (neighbors[0] if neighbors else None)
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated substituent is not supported yet")
    if any(set(graph[atom]) - atoms - {attachment_point} for atom in atoms):
        raise UnsupportedStructure(
            "a substituted substituent is not supported yet; only a "
            "plain, unsubstituted alkyl substituent is in scope"
        )
    return alkyl_name(length)


def _validate_and_collect_hetero_ring_ketone(
    mol, heteroatoms, allow_n_substituent=False, allow_ring_carbon_substituent=False
):
    """(ring_order, ketones, elements_by_atom, substituents) for a
    hetero-ring ketone -- validates that every ring heteroatom
    (`heteroatoms`, 1 or 2 atom indices) carries no substituent beyond its
    own indicated hydrogen, every other ring atom is a plain CH2 or an
    unsubstituted ketone carbonyl carbon, and the ring itself is fully
    saturated (P-22.2.1's plain retained-name ring shape). `substituents`
    is {} unless `allow_n_substituent`/`allow_ring_carbon_substituent` is
    True, in which case a ring nitrogen (mirroring `_amide.py`/
    `_hydrazide.py`'s identical N-substituent restriction) and/or a plain
    ring carbon may instead carry a single plain, unbranched,
    unsubstituted alkyl substituent -- collected as
    {atom_idx: (alkyl_name, is_compound=False)}; O/S heteroatoms never get
    this option (P-22.2.1: no spare valence to give up). A single-ketone
    ring (as opposed to a dione/hydantoin shape) is equally in scope here
    -- the caller's own locant computation is responsible for the
    P-14.5.2 substituent-locant tie-break this can introduce (see
    `_five_membered_1_3_numbering`)."""
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    graph = adjacency(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_set = set(ring_atoms)
    full_carbon_graph = (
        carbon_adjacency(mol) if allow_n_substituent or allow_ring_carbon_substituent else None
    )

    elements_by_atom = {}
    substituents = {}
    for heteroatom in heteroatoms:
        hetero_atom_obj = mol.GetAtomWithIdx(heteroatom)
        element = _HETERO_RING_ELEMENTS[hetero_atom_obj.GetAtomicNum()]
        elements_by_atom[heteroatom] = element
        if hetero_atom_obj.GetFormalCharge() != 0 or hetero_atom_obj.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        expected_h = 1 if element == "N" else 0
        if hetero_atom_obj.GetDegree() == 2 and hetero_atom_obj.GetTotalNumHs() == expected_h:
            continue
        if (
            allow_n_substituent
            and element == "N"
            and hetero_atom_obj.GetDegree() == 3
            and hetero_atom_obj.GetTotalNumHs() == 0
        ):
            (root,) = [n for n in graph[heteroatom] if n not in ring_set]
            if mol.GetAtomWithIdx(root).GetAtomicNum() == 6:
                substituents[heteroatom] = (
                    _linear_alkyl_substituent(mol, full_carbon_graph, graph, root, heteroatom),
                    False,
                )
                continue
        raise UnsupportedStructure(
            "a ring heteroatom bearing a substituent is out of scope for "
            "this module's hetero-ring ketone path (see "
            "P-22.2.1)"
        )

    ketones = set()
    carbon_substituents = {}
    for atom_idx in ring_atoms:
        if atom_idx in elements_by_atom:
            continue
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        exo = [n for n in graph[atom_idx] if n not in ring_set]
        if not exo:
            if atom.GetTotalNumHs() != 2:
                raise UnsupportedStructure(
                    "a substituent other than a ketone carbonyl is out of "
                    "scope for this module's hetero-ring ketone path"
                )
            continue
        if (
            allow_ring_carbon_substituent
            and len(exo) == 1
            and atom.GetTotalNumHs() == 1
            and mol.GetAtomWithIdx(exo[0]).GetAtomicNum() == 6
        ):
            carbon_substituents[atom_idx] = (
                _linear_alkyl_substituent(mol, full_carbon_graph, graph, exo[0], atom_idx),
                False,
            )
            continue
        if len(exo) != 1:
            raise UnsupportedStructure(
                "a substituent other than a ketone carbonyl is out of "
                "scope for this module's hetero-ring ketone path"
            )
        (other,) = exo
        other_atom = mol.GetAtomWithIdx(other)
        bond = mol.GetBondBetweenAtoms(atom_idx, other)
        if other_atom.GetAtomicNum() != 8 or other_atom.GetDegree() != 1 or bond.GetBondTypeAsDouble() != 2.0:
            raise UnsupportedStructure(
                "a substituent other than a ketone carbonyl is out of "
                "scope for this module's hetero-ring ketone path"
            )
        ketones.add(other)
    if not ketones:
        raise UnsupportedStructure("no ketone (C=O) group found; this module only handles ketones")
    substituents.update(carbon_substituents)

    for i in range(len(ring_order)):
        a, b = ring_order[i], ring_order[(i + 1) % len(ring_order)]
        if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "an unsaturated hetero ring is out of scope for this "
                "module's hetero-ring ketone path"
            )

    return ring_order, ketones, elements_by_atom, substituents


def _best_one_locants(graph, ring_order, ketones, starts):
    """Lowest ketone locant set over every (start, direction) candidate in
    `starts` (each start rotated to position 1, tried both directions) --
    shared by the single- and two-heteroatom hetero-ring ketone paths."""
    best_locants = None
    for start in starts:
        rotated_start = ring_order.index(start)
        rotated = ring_order[rotated_start:] + ring_order[:rotated_start]
        for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            one_locants = _one_locants(position_of, ketones, graph)
            if one_locants is None:
                raise UnsupportedStructure(
                    "a ketone not on the ring itself (e.g. on a substituent "
                    "branch) is not supported yet"
                )
            locant_set = lowest_locant_set(one_locants)
            if best_locants is None or locant_set < best_locants:
                best_locants = locant_set
    return best_locants


def _hetero_ring_ketone_name(stem, best_locants):
    one_word = multiplied_word(len(best_locants), "one")
    elide = one_word[0] in "aeiouy"
    base = stem[:-1] if elide else stem
    loc_str = ",".join(str(loc) for loc in best_locants)
    return f"{base}-{loc_str}-{one_word}"


def _name_hetero_cyclic_ketone(mol, heteroatom):
    graph = adjacency(mol)
    ring_order, ketones, elements_by_atom, substituents = _validate_and_collect_hetero_ring_ketone(
        mol, {heteroatom}, allow_n_substituent=True
    )
    ring_size = len(ring_order)
    element = elements_by_atom[heteroatom]
    stem = saturated_ring_name(element, ring_size)
    if stem is None:
        raise UnsupportedStructure(
            f"no retained/Hantzsch-Widman name for a {ring_size}-membered "
            f"{element}-heteroatom saturated ring (P-22.2.1)"
        )
    best_locants = _best_one_locants(graph, ring_order, ketones, [heteroatom])
    name = _hetero_ring_ketone_name(stem, best_locants)
    if not substituents:
        return name
    ((sub_name, is_compound),) = substituents.values()
    prefix = format_substituent_prefixes({sub_name: {"locants": [1], "compound": is_compound}})
    separator = "-" if name[0].isdigit() else ""
    return f"{prefix}{separator}{name}"


def _hetero_ring_two_heteroatoms(mol):
    """(het1, het2) ring-atom indices for a saturated, 6-membered,
    1,4-related two-heteroatom ketone shape (morpholine/piperazine/
    thiomorpholine's element pairs only -- see
    `_TWO_HETERO_RING_ELEMENT_PAIRS`), or None if it doesn't match that
    shape at all."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _TWO_HETERO_RING_SIZE:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _TWO_HETERO_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    if abs(ring_order.index(het1) - ring_order.index(het2)) != _TWO_HETERO_RING_SIZE // 2:
        # Not the 1,4 (antipodal) relationship morpholine/piperazine/
        # thiomorpholine need -- e.g. a 1,2- or 1,3-diheteroatom ring,
        # which has no retained name and is out of scope.
        return None
    return het1, het2


def _name_two_hetero_cyclic_ketone(mol, het1, het2):
    graph = adjacency(mol)
    ring_order, ketones, elements_by_atom, _ = _validate_and_collect_hetero_ring_ketone(mol, {het1, het2})
    stem = saturated_two_heteroatom_1_4_ring_name((elements_by_atom[het1], elements_by_atom[het2]))
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this two-heteroatom saturated ring "
            "(P-22.2.1)"
        )
    if elements_by_atom[het1] == elements_by_atom[het2]:
        starts = [het1, het2]
    else:
        starts = [het1 if _TWO_HETERO_PRIORITY[elements_by_atom[het1]] < _TWO_HETERO_PRIORITY[elements_by_atom[het2]] else het2]
    best_locants = _best_one_locants(graph, ring_order, ketones, starts)
    return _hetero_ring_ketone_name(stem, best_locants)


def _hetero_ring_seven_membered_1_4(mol):
    """(het1, het2) ring-atom indices for a saturated, 7-membered,
    1,4-related two-heteroatom ketone shape (the 1,4-diazepan-2-one/
    1,4-diazepan-5-one family -- see
    `_SEVEN_MEMBERED_1_4_RING_ELEMENT_PAIRS`), or None if it doesn't match
    that shape at all. Unlike the 6-membered 1,4-ring (whose two arcs
    between the heteroatoms are both length 3, so either numbering
    direction reaches locant 4), a 7-membered ring's two arcs are length 2
    and 3 -- only the short (length-2) arc gives a valid '1,4-' numbering,
    so a heteroatom pair must sit exactly 3 ring bonds apart on one side
    (equivalently 4 the other way) to qualify at all; see
    `_seven_membered_1_4_numbering` for how the single valid direction per
    start is then picked out."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _SEVEN_MEMBERED_1_4_RING_SIZE:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _SEVEN_MEMBERED_1_4_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    diff = abs(ring_order.index(het1) - ring_order.index(het2))
    if min(diff, _SEVEN_MEMBERED_1_4_RING_SIZE - diff) != _SEVEN_MEMBERED_1_4_RING_SIZE // 2:
        return None
    return het1, het2


def _seven_membered_1_4_numbering(ring_order, het1, het2, elements_by_atom, ketones, graph):
    """Lowest ketone locant set for the 7-membered 1,4-two-heteroatom
    shape. When the two heteroatoms differ, P-22.2.1 element seniority
    (`_TWO_HETERO_PRIORITY`) fixes the higher-priority one at locant 1;
    when identical, either may be locant 1, so both are tried. Whichever
    start is used, only one of the two rotation directions actually lands
    the other heteroatom on locant 4 (the ring's short arc) -- the other
    direction would misnumber it to locant 5, contradicting the ring's own
    '1,4-' name -- so (unlike the symmetric 6-membered case in
    `_best_one_locants`) there is no free direction choice to minimize the
    ketone locant against; each start contributes exactly one candidate."""
    starts = [het1, het2] if elements_by_atom[het1] == elements_by_atom[het2] else [
        het1 if _TWO_HETERO_PRIORITY[elements_by_atom[het1]] < _TWO_HETERO_PRIORITY[elements_by_atom[het2]] else het2
    ]
    best_locants = None
    for start in starts:
        other = het2 if start == het1 else het1
        rotated_start = ring_order.index(start)
        rotated = ring_order[rotated_start:] + ring_order[:rotated_start]
        if rotated.index(other) != _SEVEN_MEMBERED_1_4_RING_SIZE // 2:
            rotated = [rotated[0]] + list(reversed(rotated[1:]))
        position_of = {atom: i + 1 for i, atom in enumerate(rotated)}
        one_locants = _one_locants(position_of, ketones, graph)
        if one_locants is None:
            raise UnsupportedStructure(
                "a ketone not on the ring itself (e.g. on a substituent "
                "branch) is not supported yet"
            )
        locant_set = lowest_locant_set(one_locants)
        if best_locants is None or locant_set < best_locants:
            best_locants = locant_set
    return best_locants


def _name_seven_membered_1_4_ring_ketone(mol, het1, het2):
    graph = adjacency(mol)
    ring_order, ketones, elements_by_atom, _ = _validate_and_collect_hetero_ring_ketone(mol, {het1, het2})
    stem = saturated_seven_membered_1_4_two_heteroatom_ring_name(
        (elements_by_atom[het1], elements_by_atom[het2])
    )
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this seven-membered two-heteroatom "
            "saturated ring (P-22.2.1)"
        )
    best_locants = _seven_membered_1_4_numbering(ring_order, het1, het2, elements_by_atom, ketones, graph)
    return _hetero_ring_ketone_name(stem, best_locants)


def _hetero_ring_five_membered_1_3(mol):
    """(het1, het2) ring-atom indices for a saturated, 5-membered,
    1,3-related two-heteroatom ketone shape (the imidazolidin-2-one/
    imidazolidine-2,4-dione family -- see
    `_FIVE_MEMBERED_1_3_RING_ELEMENT_PAIRS`), or None if it doesn't match
    that shape at all (heteroatoms adjacent to each other rather than
    1,3-related, an unsupported element pair, or not a 5-membered
    saturated ring with exactly two heteroatoms)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 5:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    ring_set = set(ring_atoms)
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _FIVE_MEMBERED_1_3_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    if het2 in graph[het1]:
        return None
    bridging = (set(graph[het1]) & set(graph[het2]) & ring_set) - {het1, het2}
    if len(bridging) != 1:
        return None
    return het1, het2


def _five_membered_1_3_numbering(mol, het1, het2, elements_by_atom, ketones, substituent_atoms=frozenset()):
    """(best_locants, position_of) for the five-membered 1,3-two-heteroatom
    ring shape -- `position_of` is the full {atom_idx: locant} mapping for
    whichever numbering direction produced `best_locants`, exposed so a
    caller can also look up a substituent's own locant. Locant 1 goes to
    the higher-priority heteroatom (P-22.2.1's O > S > N element
    seniority, `_TWO_HETERO_PRIORITY`) when the two elements differ -- the
    heteroatom-seniority rule outranks locant minimization entirely, so
    there is no freedom to choose the lower-locant direction instead
    (confirmed against PubChem: 'oxazolidine-2,4-dione' would minimize the
    ketone locant to 4, but the real compound is
    '1,3-oxazolidine-2,5-dione' because O must stay locant 1). When the
    two heteroatoms are identical, either may be locant 1, so both are
    tried and the one minimizing the ketone locant set wins (e.g.
    imidazolidine-2,4-dione, not -2,5-dione). With a dione this always
    settles the tie outright (the two candidate ketone-locant sets always
    differ, e.g. {2,4} vs {2,5}); with a single ketone the lone ketone
    always sits on the bridging carbon (locant 2) regardless of
    direction, so both candidates tie on the ketone locant and the choice
    falls through to P-14.5.2's next criterion -- lowest locants to the
    full set of substituents cited by prefix (`substituent_atoms`, ring
    carbon and/or N alike), e.g. '4-methylimidazolidin-2-one' (PubChem
    CID 97832), not '5-methyl...'. Only the direction that reaches the
    other heteroatom via the bridging carbon (locant 2) is considered in
    either case -- the opposite direction would misnumber the other
    heteroatom to locant 4 instead of 3, contradicting the ring's own
    '1,3-' name (unlike the symmetric six-membered 1,4-case in
    `_best_one_locants`, the two directions are not equivalent here)."""
    graph = adjacency(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])
    (bridge,) = (set(graph[het1]) & set(graph[het2]) & ring_atoms) - {het1, het2}
    long_near = {}
    for het, other in ((het1, het2), (het2, het1)):
        (long_near[het],) = (set(graph[het]) & ring_atoms) - {bridge, other}

    def candidate(start, other):
        position_of = {start: 1, bridge: 2, other: 3, long_near[other]: 4, long_near[start]: 5}
        locants = sorted(position_of[next(iter(graph[oxygen]))] for oxygen in ketones)
        return locants, position_of

    def sort_key(result):
        locants, position_of = result
        return locants, sorted(position_of[atom] for atom in substituent_atoms)

    if elements_by_atom[het1] == elements_by_atom[het2]:
        a, b = candidate(het1, het2), candidate(het2, het1)
        return min(a, b, key=sort_key)
    start, other = (
        (het1, het2)
        if _TWO_HETERO_PRIORITY[elements_by_atom[het1]] < _TWO_HETERO_PRIORITY[elements_by_atom[het2]]
        else (het2, het1)
    )
    return candidate(start, other)


def _name_five_membered_1_3_ring_ketone(mol, het1, het2):
    ring_order, ketones, elements_by_atom, n_substituents = _validate_and_collect_hetero_ring_ketone(
        mol, {het1, het2}, allow_n_substituent=True, allow_ring_carbon_substituent=True
    )
    stem = saturated_five_membered_1_3_two_heteroatom_ring_name(
        (elements_by_atom[het1], elements_by_atom[het2])
    )
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this five-membered two-heteroatom "
            "saturated ring (P-22.2.1)"
        )
    best_locants, position_of = _five_membered_1_3_numbering(
        mol, het1, het2, elements_by_atom, ketones, set(n_substituents)
    )
    name = _hetero_ring_ketone_name(stem, best_locants)
    if not n_substituents:
        return name
    grouped = {}
    for atom_idx, (sub_name, is_compound) in n_substituents.items():
        info = grouped.setdefault(sub_name, {"locants": [], "compound": is_compound})
        info["locants"].append(position_of[atom_idx])
    prefix = format_substituent_prefixes(grouped)
    separator = "-" if name[0].isdigit() else ""
    return f"{prefix}{separator}{name}"


def _hetero_ring_seven_membered_1_3(mol):
    """(het1, het2) ring-atom indices for a saturated, 7-membered,
    1,3-related two-heteroatom ketone shape (the 1,3-diazepan-2-one
    family -- see `_SEVEN_MEMBERED_1_3_RING_ELEMENT_PAIRS`), or None if it
    doesn't match that shape at all. Structurally identical to the
    5-membered 1,3-case above (a single bridging carbon between the two
    heteroatoms, the rest of the ring plain CH2) except for the longer
    4-carbon far arc; unlike that case, only the narrow single-ketone
    shape (ketone on the bridging carbon, always locant 2) is attempted
    here -- no PubChem-registered example was found for a second ketone
    on the far arc or for a heteroatom/ring-carbon alkyl substituent, so
    those stay unsupported (see `_name_seven_membered_1_3_ring_ketone`)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _SEVEN_MEMBERED_1_3_RING_SIZE:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    ring_set = set(ring_atoms)
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _SEVEN_MEMBERED_1_3_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    if het2 in graph[het1]:
        return None
    bridging = (set(graph[het1]) & set(graph[het2]) & ring_set) - {het1, het2}
    if len(bridging) != 1:
        return None
    return het1, het2


def _name_seven_membered_1_3_ring_ketone(mol, het1, het2):
    graph = adjacency(mol)
    ring_set = set(mol.GetRingInfo().AtomRings()[0])
    (bridge,) = (set(graph[het1]) & set(graph[het2]) & ring_set) - {het1, het2}
    _, ketones, elements_by_atom, _ = _validate_and_collect_hetero_ring_ketone(mol, {het1, het2})
    stem = saturated_seven_membered_1_3_two_heteroatom_ring_name(
        (elements_by_atom[het1], elements_by_atom[het2])
    )
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this seven-membered two-heteroatom "
            "saturated ring (P-22.2.1)"
        )
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "a second ketone on the far arc of this ring shape is not "
            "supported yet; only a single ketone on the bridging carbon "
            "between the two heteroatoms is in scope"
        )
    (ketone_oxygen,) = ketones
    (carbon,) = graph[ketone_oxygen]
    if carbon != bridge:
        raise UnsupportedStructure(
            "a ketone anywhere but the bridging carbon between the two "
            "heteroatoms is not supported yet for this ring shape"
        )
    return _hetero_ring_ketone_name(stem, [2])


def _hetero_ring_five_membered_1_2(mol):
    """(het1, het2) ring-atom indices for a saturated, 5-membered,
    1,2-related (directly bonded) two-heteroatom ketone shape (the
    pyrazolidin-3-one family -- see
    `_FIVE_MEMBERED_1_2_RING_ELEMENT_PAIRS`), or None if it doesn't match
    that shape at all (heteroatoms not directly bonded, an unsupported
    element pair, or not a 5-membered saturated ring with exactly two
    heteroatoms)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 5:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _FIVE_MEMBERED_1_2_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    if het2 not in graph[het1]:
        return None
    return het1, het2


def _five_membered_1_2_numbering(mol, het1, het2, elements_by_atom, ketones):
    """(best_locants, position_of) for the five-membered 1,2-two-heteroatom
    ring shape -- the two heteroatoms are always locants 1 and 2 (directly
    bonded), and the remaining three ring carbons follow as locants 3, 4,
    5 in the one direction that reaches the other heteroatom immediately
    (the opposite direction from a given start reaches a ring carbon
    first instead, misnumbering the heteroatom pair as e.g. {1, 5} rather
    than the required lowest set {1, 2}, so it's never a valid candidate
    -- mirrors `_five_membered_1_3_numbering`'s identical "only one
    direction is geometrically valid" reasoning). Locant 1 goes to the
    higher-priority heteroatom (P-22.2.1's O > S > N, `_TWO_HETERO_PRIORITY`)
    when the two elements differ -- no minimization freedom, same
    seniority-over-locants rule as the 1,3-case. When identical, both
    starts are tried and the one minimizing the ketone locant set wins."""
    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(mol.GetRingInfo().AtomRings()[0]))

    def candidates(start, other):
        rotated_start = ring_order.index(start)
        rotated = ring_order[rotated_start:] + ring_order[:rotated_start]
        for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
            if candidate[1] != other:
                continue
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locants = sorted(position_of[next(iter(graph[oxygen]))] for oxygen in ketones)
            yield locants, position_of

    if elements_by_atom[het1] == elements_by_atom[het2]:
        options = list(candidates(het1, het2)) + list(candidates(het2, het1))
        return min(options, key=lambda option: option[0])
    start, other = (
        (het1, het2)
        if _TWO_HETERO_PRIORITY[elements_by_atom[het1]] < _TWO_HETERO_PRIORITY[elements_by_atom[het2]]
        else (het2, het1)
    )
    (result,) = candidates(start, other)
    return result


def _name_five_membered_1_2_ring_ketone(mol, het1, het2):
    ring_order, ketones, elements_by_atom, n_substituents = _validate_and_collect_hetero_ring_ketone(
        mol, {het1, het2}, allow_n_substituent=True, allow_ring_carbon_substituent=True
    )
    stem = saturated_five_membered_1_2_two_heteroatom_ring_name(
        (elements_by_atom[het1], elements_by_atom[het2])
    )
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this five-membered two-heteroatom "
            "saturated ring (P-22.2.1)"
        )
    best_locants, position_of = _five_membered_1_2_numbering(mol, het1, het2, elements_by_atom, ketones)
    name = _hetero_ring_ketone_name(stem, best_locants)
    if not n_substituents:
        return name
    grouped = {}
    for atom_idx, (sub_name, is_compound) in n_substituents.items():
        info = grouped.setdefault(sub_name, {"locants": [], "compound": is_compound})
        info["locants"].append(position_of[atom_idx])
    prefix = format_substituent_prefixes(grouped)
    separator = "-" if name[0].isdigit() else ""
    return f"{prefix}{separator}{name}"


def _hetero_ring_seven_membered_1_2(mol):
    """(het1, het2) ring-atom indices for a saturated, 7-membered,
    1,2-related (directly bonded) two-heteroatom ketone shape (the
    1,2-diazepan-3-one family -- see
    `_SEVEN_MEMBERED_1_2_RING_ELEMENT_PAIRS`), or None if it doesn't match
    that shape at all. Structurally identical to the 5-membered 1,2-case
    above (heteroatoms at locants 1/2, the ketone on the ring carbon
    immediately after locant 2) except for the longer far arc (4 carbons
    instead of 2, at locants 4-7); unlike that case, only the narrow
    single-ketone shape (locant 3) is attempted here -- no
    PubChem-registered example was found for a second ketone on the far
    arc or for a ring-carbon alkyl substituent, so both stay unsupported
    (see `_name_seven_membered_1_2_ring_ketone`)."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != _SEVEN_MEMBERED_1_2_RING_SIZE:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 2:
        return None
    het1, het2 = heteroatoms
    atomic_num1 = mol.GetAtomWithIdx(het1).GetAtomicNum()
    atomic_num2 = mol.GetAtomWithIdx(het2).GetAtomicNum()
    if atomic_num1 not in _HETERO_RING_ELEMENTS or atomic_num2 not in _HETERO_RING_ELEMENTS:
        return None
    elements = frozenset((_HETERO_RING_ELEMENTS[atomic_num1], _HETERO_RING_ELEMENTS[atomic_num2]))
    if elements not in _SEVEN_MEMBERED_1_2_RING_ELEMENT_PAIRS:
        return None
    graph = adjacency(mol)
    if het2 not in graph[het1]:
        return None
    return het1, het2


def _name_seven_membered_1_2_ring_ketone(mol, het1, het2):
    _, ketones, elements_by_atom, _ = _validate_and_collect_hetero_ring_ketone(mol, {het1, het2})
    stem = saturated_seven_membered_1_2_two_heteroatom_ring_name(
        (elements_by_atom[het1], elements_by_atom[het2])
    )
    if stem is None:
        raise UnsupportedStructure(
            "no retained name for this seven-membered two-heteroatom "
            "saturated ring (P-22.2.1)"
        )
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "a second ketone on the far arc of this ring shape is not "
            "supported yet; only a single ketone adjacent to the second "
            "heteroatom is in scope"
        )
    # `_five_membered_1_2_numbering` makes no assumption about ring size --
    # it derives locants entirely from the molecule's own ring traversal --
    # so it's reused as-is here rather than duplicated.
    best_locants, _ = _five_membered_1_2_numbering(mol, het1, het2, elements_by_atom, ketones)
    return _hetero_ring_ketone_name(stem, best_locants)


def has_hetero_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the narrow hetero-ring-ketone shape
    (single heteroatom, see `_hetero_ring_heteroatom`, or the 1,4
    two-heteroatom morpholine/piperazine/thiomorpholine/1,4-dioxane
    shape, see `_hetero_ring_two_heteroatoms`) -- used by `core.py` to
    route ahead of both `has_ester_shape` and `has_amide_shape`, since a
    ketone directly bonded to a ring oxygen/sulfur (a lactone, e.g.
    oxan-2-one/1,4-dioxan-2-one) or ring nitrogen (a lactam, e.g.
    piperidin-2-one/morpholin-3-one) would otherwise look ester- or
    amide-shaped to those checks and get rejected by `_ester.py`'s
    acyclic-only construction or `_amide.py`'s ring-always-out-of-scope
    guard before ever reaching this module. Deliberately excludes the
    five-membered 1,3-two-heteroatom shape (see
    `has_five_membered_1_3_ring_ketone_shape` below) -- that shape needs
    to be routed much earlier (ahead of ether/acetal/carbamate too), while
    this single-/1,4-two-heteroatom shape must stay routed after
    `has_anhydride_shape`: a single-heteroatom ring flanked by a ketone on
    both neighboring carbons (e.g. 'O=C1CCC(=O)O1') is itself a cyclic
    anhydride shape, and moving this check any earlier would silently
    reclassify it before `_anhydride.py`'s own explicit ring rejection."""
    return (
        _hetero_ring_heteroatom(mol) is not None
        or _hetero_ring_two_heteroatoms(mol) is not None
        or _hetero_ring_seven_membered_1_4(mol) is not None
    )


def has_five_membered_1_3_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the five-membered 1,3-two-heteroatom
    ring-ketone shape (see `_hetero_ring_five_membered_1_3`, one or two
    ketones -- e.g. 1,3-dioxolan-2-one, or the hydantoin family like
    imidazolidine-2,4-dione) -- kept as its own separate check from
    `has_hetero_ring_ketone_shape` above because it must be routed much
    earlier in `core.py` (ahead of ether/acetal/carbamate/alkoxide/
    hydroxylamine/cyanate, which would otherwise misname its various O/N/S
    element-pair combinations), and unlike the single-/1,4-two-heteroatom
    shapes, it can never collide with the cyclic anhydride shape (that
    needs exactly one ring heteroatom flanked by two ketones; this shape
    always has exactly two ring heteroatoms), so the earlier routing is
    safe."""
    return _hetero_ring_five_membered_1_3(mol) is not None


def has_seven_membered_1_3_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the 7-membered 1,3-two-heteroatom
    ring-ketone shape (see `_hetero_ring_seven_membered_1_3`, e.g.
    1,3-diazepan-2-one) -- kept alongside
    `has_five_membered_1_3_ring_ketone_shape` for the identical reason: an
    O+N pair here looks carbamate-shaped (N-C(=O)-O) and must be routed
    ahead of `has_carbamate_shape` (and the other checks that shape's
    docstring lists), and this shape can never collide with the cyclic
    anhydride shape either (always exactly two ring heteroatoms)."""
    return _hetero_ring_seven_membered_1_3(mol) is not None


def has_five_membered_1_2_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the five-membered 1,2-two-heteroatom
    ring-ketone shape (see `_hetero_ring_five_membered_1_2`, e.g.
    pyrazolidin-3-one) -- kept separate from
    `has_five_membered_1_3_ring_ketone_shape` for the same reasons: it
    must be routed just as early in `core.py` (this shape, too, is
    otherwise misnamed by the aldehyde/ether checks depending on the
    element pair), and can never collide with the cyclic anhydride shape
    (which needs exactly one ring heteroatom)."""
    return _hetero_ring_five_membered_1_2(mol) is not None


def has_seven_membered_1_2_ring_ketone_shape(mol) -> bool:
    """True if this molecule fits the 7-membered 1,2-two-heteroatom
    ring-ketone shape (see `_hetero_ring_seven_membered_1_2`, e.g.
    1,2-diazepan-3-one) -- kept alongside
    `has_five_membered_1_2_ring_ketone_shape` for the identical reason
    (misnamed by the aldehyde/ether checks depending on the element pair
    if not claimed here first), and likewise can never collide with the
    cyclic anhydride shape."""
    return _hetero_ring_seven_membered_1_2(mol) is not None


def name_ketone(mol) -> str:
    hetero_atom = _hetero_ring_heteroatom(mol)
    if hetero_atom is not None:
        return _name_hetero_cyclic_ketone(mol, hetero_atom)
    five_membered_1_3 = _hetero_ring_five_membered_1_3(mol)
    if five_membered_1_3 is not None:
        return _name_five_membered_1_3_ring_ketone(mol, *five_membered_1_3)
    seven_membered_1_3 = _hetero_ring_seven_membered_1_3(mol)
    if seven_membered_1_3 is not None:
        return _name_seven_membered_1_3_ring_ketone(mol, *seven_membered_1_3)
    five_membered_1_2 = _hetero_ring_five_membered_1_2(mol)
    if five_membered_1_2 is not None:
        return _name_five_membered_1_2_ring_ketone(mol, *five_membered_1_2)
    seven_membered_1_2 = _hetero_ring_seven_membered_1_2(mol)
    if seven_membered_1_2 is not None:
        return _name_seven_membered_1_2_ring_ketone(mol, *seven_membered_1_2)
    two_heteroatoms = _hetero_ring_two_heteroatoms(mol)
    if two_heteroatoms is not None:
        return _name_two_hetero_cyclic_ketone(mol, *two_heteroatoms)
    seven_membered_1_4 = _hetero_ring_seven_membered_1_4(mol)
    if seven_membered_1_4 is not None:
        return _name_seven_membered_1_4_ring_ketone(mol, *seven_membered_1_4)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        is_benzene = is_plain_benzene_ring(mol, ring_atoms)
        is_heteroaromatic = not is_benzene and (
            heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
        )
        if is_benzene or is_heteroaromatic:
            return _name_phenyl_chain_ketone(mol, ring_atoms)

    # Two separate simple monocycles joined by one direct bond, one a plain
    # benzo/heteroaromatic ring with no substituent of its own (P-25 M1 step
    # 2, #622) -- e.g. 2-phenylcyclohexan-1-one -- reuses the existing
    # single-ring `_name_cyclic_ketone` dispatch below with the aromatic
    # ring's atoms passed through as an exemption/substituent, rather than
    # falling into the generic "polycyclic" rejection.
    aromatic_shape = None
    if ring_info.NumRings() == 2:
        aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
    aromatic_atoms = aromatic_shape[1] if aromatic_shape is not None else frozenset()

    ketones, hydroxyls = _validate_and_collect_ketones(mol, aromatic_ring_atoms=aromatic_atoms)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the ketone oxygen, never part of any carbon chain). Also exclude a
    # bond entirely inside the aromatic-substituent ring above -- its
    # aromatic bond order (1.5) isn't a real chain/ring 'ene'/'yne' bond
    # either, and is already accounted for by naming that ring via
    # `name_branch` instead.
    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in ketones
        and b[1] not in ketones
        and not (b[0] in aromatic_atoms and b[1] in aromatic_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior ketone/aldehyde "
                "form and is out of scope for this module (P-31.1.4.2.4)"
            )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        # No ring is ever involved past this point, so a specified C=C
        # double-bond E/Z element may coexist with a specified tetrahedral
        # stereocenter (P-91.3), mirroring `_alcohol.py`'s identical
        # acyclic-only combined-stereo split.
        return _name_acyclic_ketone(mol, ketones, hydroxyls, bonds, specified_stereo_elements(mol))
    stereo = specified_stereocenters(mol)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic ketone is "
                "not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a ketone is not "
                "supported yet -- only a ring double bond is in scope for "
                "this first pass (see P-31.1.3)"
            )
        ring_ketones = {o for o in ketones if next(iter(graph[o])) in ring_atoms}
        if not bonds and not hydroxyls and not ring_ketones and len(ketones) == 1:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_ketone(mol, ketones)
        if not bonds and not hydroxyls and ring_ketones and ring_ketones != ketones:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain ketone "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_ketone_chain_ketone(mol, ketones)
        return _name_cyclic_ketone(mol, ketones, hydroxyls, stereo, bonds)
    if num_rings == 2 and aromatic_shape is not None:
        ring_atoms, _, _, _ = aromatic_shape
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic ketone is "
                "not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a ketone is not "
                "supported yet -- only a ring double bond is in scope for "
                "this first pass (see P-31.1.3)"
            )
        return _name_cyclic_ketone(
            mol, ketones, hydroxyls, stereo, bonds, ring_atoms=ring_atoms, aromatic_atoms=aromatic_atoms
        )
    raise UnsupportedStructure(
        "polycyclic and spiro ketones are not supported yet (P-23/P-24/P-25 "
        "numbering integration with a suffix group is future work)"
    )
