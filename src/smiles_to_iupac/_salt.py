"""Naming of the simplest salts (P-77, Blue Book P-7,
https://iupac.qmul.ac.uk/BlueBook/P7.html): a name formed by citing the
cation's name followed by the anion's name, e.g. 'sodium ethanoate'
(PubChem's own computed name for the same structure is 'sodium acetate';
this project's `_carboxylate.py` already deliberately uses the systematic
'-oate' stem over the retained one for consistency with
`_carboxylic_acid.py`, so this module inherits that same divergence -- see
`_carboxylate.py`'s own docstring).

This first pass is deliberately the narrowest possible slice, chosen
because the general case needs infrastructure this project doesn't have at
all yet (multi-fragment SMILES are not otherwise recognized anywhere else
in this codebase -- every other module assumes a single connected
molecule and would reject a foreign atom like sodium outright):

- The cation must be a single bare monoatomic metal with a fixed formal
  charge (Group 1 alkali metals at +1, Group 2 alkaline earth metals at
  +2, aluminium at +3) -- its own name is a fixed lookup, no substitutive
  nomenclature involved.
- The anion(s) must all be the same carboxylate, recognized by
  `_carboxylate.py` (`has_carboxylate_shape`/`name_carboxylate`), reused
  unchanged. Exactly as many anion fragments as the cation's charge are
  required (one cation, charge-balanced), each cited once with a
  multiplying prefix (P-16.3.3) when there is more than one, e.g.
  'calcium diethanoate' (P-65.6.2.1's own worked example, 'calcium
  diacetate' PIN, confirms this multiplier-on-the-anion pattern).

The cation may also be a single ammonium-shaped fragment recognized by
`_ammonium.py` (`has_ammonium_shape`/`name_ammonium`, reused unchanged) --
NH4+ or an N-alkyl-substituted ammonium, always singly charged, so exactly
one anion fragment is required and no multiplying prefix ever applies.
PubChem-verified: `[NH4+].CC(=O)[O-]` -> 'azanium acetate',
`C[NH3+].CC(=O)[O-]` -> 'methylazanium acetate' -- this project's own
`_ammonium.py` PIN convention names the latter cation 'methanaminium'
rather than PubChem's 'methylazanium' (a divergence already established
and accepted, see `_ammonium.py`'s own docstring), so this module's output
for that example is 'methanaminium acetate'.

The anion may also be a single alkoxide, thioate, or selenoate, recognized
by `_alkoxide.py`/`_thioate.py`/`_selenoate.py` (each module's own
`has_*_shape`/`name_*`, reused unchanged), but -- unlike carboxylate --
only paired with a singly-charged cation (alkali metal or ammonium):
PubChem-verified `[Na+].C[O-]` -> 'sodium methanolate',
`[NH4+].C[O-]` -> 'azanium methanolate',
`[Na+].CC(=O)[S-]` -> 'sodium ethanethioate' (an exact match with this
project's own thioate naming, unlike the alkoxide/carboxylate retained-
name divergences),
`[NH4+].CC(=O)[S-]` -> 'azanium ethanethioate',
`[Na+].CC(=O)[Se-]` -> 'sodium ethaneselenoate' -- all exactly one cation
to one anion. A 2+/3+ metal cation with two or three such anions is *not*
attempted here -- PubChem's own generator inconsistently omits the
multiplying prefix there for all three (`[Ca+2].C[O-].C[O-]` -> 'calcium
methanolate', `[Ca+2].CC(=O)[S-].CC(=O)[S-]` -> 'calcium ethanethioate',
neither 'di...', unlike the carboxylate case's confirmed 'calcium
diacetate'), so there is no reliable worked example to verify the
multivalent-cation form of any of these salts against yet.

The anion may also be a single monoatomic halide (F-/Cl-/Br-/I-, a fixed
name lookup off `_common.py`'s own `HALOGEN_PREFIXES`), with no
singly-charged restriction -- unlike alkoxide/thioate/selenoate above,
PubChem's own generator *does* apply the multiplying prefix consistently
here: `[Na+].[Cl-]` -> 'sodium chloride', `[NH4+].[Cl-]` -> 'azanium
chloride', `[Ca+2].[Cl-].[Cl-]` -> 'calcium dichloride'. This is a
distinct shape from `_hydrohalide_salt.py`'s bare *neutral* HX fragment
(P-77.1.3(3) general nomenclature for an organic base) -- no overlap.

The anion may also be a single fixed-formula polyatomic inorganic ion --
sulfate (SO4^2-), carbonate (CO3^2-), nitrate (NO3^-), or phosphate
(PO4^3-) (`_is_bare_polyatomic_ion`, accepting any resonance depiction of
the charge split across the oxygens) -- each with its own plain retained
name (no organic-acid-derived namer to reuse, unlike carboxylate/
alkoxide/thioate/selenoate above). This anion shape is the *opposite* of
every other one above: exactly *one* anion fragment, balanced by
`magnitude / cation_charge` copies of *one* cation type (monoatomic or
ammonium) when that division is exact, with the multiplying prefix now on
the *cation* side instead -- e.g. 2 K+ for 1 SO4^2- -> 'dipotassium
sulfate', 2 Na+ for 1 CO3^2- -> 'disodium carbonate' (P-65.6.2's own
worked example), 1 K+ for 1 NO3- -> 'potassium nitrate' (no prefix). A
charge ratio that doesn't divide evenly (e.g. Al3+ with SO4^2-, needing
mixed 2:3 cation:anion counts) is out of scope, not attempted here.

A polyatomic-inorganic anion (only) may also be balanced by two or more
*different* cation types together instead of several copies of one type
-- e.g. 1 NH4+ + 1 K+ for 1 SO4^2- -> 'azanium potassium sulfate', or a
mix of charges, e.g. 1 Ca2+ + 1 Na+ for 1 PO4^3- -> 'calcium sodium
phosphate' -- each cation cited once (no multiplying prefix, since the
cations differ) in alphabetical order, as long as their charges sum
exactly to the anion's magnitude (P-65.6.2's own worked examples:
'potassium sodium butanedioate', 'ammonium potassium hexanedioate' --
those two are themselves dicarboxylate examples this module still can't
reach, since a 'dioate' multi-carboxylate anion doesn't exist anywhere in
this project yet, but the same charge-balancing mechanism they illustrate
is what this covers for the polyatomic-anion side). A repeated cation
type is always resolved by the multiplying-prefix shape above instead,
even alongside a different second cation type (e.g. two Na+ and one K+
for a -3 anion is not attempted -- only "all one type" or "one copy each
of all-different types" are).

Explicitly out of scope (raise `UnsupportedStructure`, or -- for
`has_salt_shape` -- simply return False so the shape falls through to
every other branch's own, usually less helpful, rejection): any metal
cation other than the fixed-valence ones above (transition metals need
Stock/oxidation-number disambiguation, not attempted here), any anion
other than a plain carboxylate, (singly-charged-cation-only) alkoxide/
thioate/selenoate, halide, or the four polyatomic inorganic anions above,
mixed/different anions on the same cation, more than one cation *type*
balancing anything other than one of the four polyatomic inorganic
anions above (a multi-cation carboxylate/alkoxide/thioate/selenoate/
halide salt is still unattempted here), a 2+/3+ metal cation paired with
one of the singly-charged-cation-only organic anions (see above), a
repeated cation type alongside a different second type (see above), and
a polyatomic-anion charge magnitude that no combination of the
fixed-valence cations above sums to exactly."""

from rdkit import Chem

from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._ammonium import has_ammonium_shape, name_ammonium
from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._numerals import multiplying_prefix
from ._selenoate import has_selenoate_shape, name_selenoate
from ._thioate import has_thioate_shape, name_thioate


def _has_halide_anion_shape(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetFormalCharge() == -1


def _name_halide_anion(frag) -> str:
    prefix = HALOGEN_PREFIXES[frag.GetAtomWithIdx(0).GetAtomicNum()]
    return prefix[:-1] + "ide"


_ANION_KINDS = [
    (has_carboxylate_shape, name_carboxylate),
    (has_alkoxide_shape, name_alkoxide),
    (has_thioate_shape, name_thioate),
    (has_selenoate_shape, name_selenoate),
    (_has_halide_anion_shape, _name_halide_anion),
]


def _is_bare_polyatomic_ion(frag, center_atomic_num, num_oxygens, total_charge):
    """True iff `frag` is a single central atom (`center_atomic_num`, formal
    charge 0) bonded to exactly `num_oxygens` monovalent oxygens and
    nothing else, with the fragment's combined formal charge exactly
    `total_charge` -- accepts any resonance depiction of the charge split
    across the oxygens (e.g. sulfate's two anionic + two neutral
    doubly-bonded oxygens), since only the stoichiometry and net charge
    matter for naming these fixed-formula ions (P-65.6.2)."""
    if frag.GetNumAtoms() != 1 + num_oxygens:
        return False
    centers = [a for a in frag.GetAtoms() if a.GetAtomicNum() == center_atomic_num]
    if len(centers) != 1:
        return False
    (center,) = centers
    # The central atom's own formal charge varies by resonance depiction
    # (e.g. nitrate's N is routinely written charge +1, '[O-][N+](=O)[O-]',
    # to satisfy valence -- unlike carbonate's neutral-C depiction) --
    # only the fragment's *combined* charge below is checked, not any
    # individual atom's.
    if center.GetIsotope() != 0 or center.GetDegree() != num_oxygens:
        return False
    others = [a for a in frag.GetAtoms() if a.GetIdx() != center.GetIdx()]
    if any(a.GetAtomicNum() != 8 or a.GetDegree() != 1 or a.GetIsotope() != 0 for a in others):
        return False
    return sum(a.GetFormalCharge() for a in frag.GetAtoms()) == total_charge


def has_sulfate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 16, 4, -2)


def name_sulfate(frag) -> str:
    return "sulfate"


def has_carbonate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 6, 3, -2)


def name_carbonate(frag) -> str:
    return "carbonate"


def has_nitrate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 7, 3, -1)


def name_nitrate(frag) -> str:
    return "nitrate"


def has_phosphate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 15, 4, -3)


def name_phosphate(frag) -> str:
    return "phosphate"


_POLYATOMIC_ANION_KINDS = [
    (has_sulfate_shape, name_sulfate),
    (has_carbonate_shape, name_carbonate),
    (has_nitrate_shape, name_nitrate),
    (has_phosphate_shape, name_phosphate),
]


def _polyatomic_anion(frag):
    """(namer, charge_magnitude) for `frag` if it matches one of the fixed-
    formula polyatomic anions above, else None."""
    for has_shape, namer in _POLYATOMIC_ANION_KINDS:
        if has_shape(frag):
            magnitude = -sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
            return namer, magnitude
    return None

_SINGLY_CHARGED_CATION_ONLY_ANIONS = {name_alkoxide, name_thioate, name_selenoate}

_MONOATOMIC_CATION_NAMES = {
    ("Li", 1): "lithium",
    ("Na", 1): "sodium",
    ("K", 1): "potassium",
    ("Rb", 1): "rubidium",
    ("Cs", 1): "cesium",
    ("Be", 2): "beryllium",
    ("Mg", 2): "magnesium",
    ("Ca", 2): "calcium",
    ("Sr", 2): "strontium",
    ("Ba", 2): "barium",
    ("Al", 3): "aluminium",
}


def _monoatomic_cation(frag):
    if frag.GetNumAtoms() != 1:
        return None
    atom = frag.GetAtomWithIdx(0)
    charge = atom.GetFormalCharge()
    name = _MONOATOMIC_CATION_NAMES.get((atom.GetSymbol(), charge))
    if name is None:
        return None
    return name, charge


def _cation(frag):
    cation = _monoatomic_cation(frag)
    if cation is not None:
        return cation
    if has_ammonium_shape(frag):
        # `has_ammonium_shape` only looks at the charged nitrogen itself,
        # ignoring the rest of the fragment -- when called (as below) on a
        # would-be *single*-fragment "cation" that is actually a whole
        # zwitterion, `name_ammonium` then tries to name the entire
        # fragment as if it were just the ammonium compound and fails
        # deeper in (e.g. on a coexisting carboxylate oxygen). That failure
        # means this fragment was never a real bare cation, not that
        # `_split_cation_anions` itself should crash -- `has_zwitterion_shape`
        # (see `_zwitterion.py`) is routed ahead of this module for the
        # zwitterion case specifically, but this catch stays regardless as
        # a defensive fallback for any other shape that superficially
        # matches on the nitrogen alone.
        try:
            return name_ammonium(frag), 1
        except UnsupportedStructure:
            return None
    return None


def _split_cation_anions(mol):
    """(cation_parts, namer, anion_frag, anion_count), or None -- where
    `cation_parts` is the already-ordered list of cation display strings
    (each with its own multiplying prefix already baked in when it's
    multiple copies of one cation type) to join with spaces before the
    anion name.

    Three distinct charge-balancing shapes are tried:
    - One cation fragment (charge `c`) balanced by `c` singly-charged (-1)
      anion fragments of the same organic-anion kind (the original shape
      this module supported, e.g. 'calcium diethanoate').
    - One polyatomic-inorganic anion fragment (charge magnitude `k`,
      P-65.6.2's sulfate/carbonate/nitrate/phosphate) balanced by `k / c`
      copies of one cation type of charge `c`, when that division is exact
      (e.g. 2 K+ for 1 SO4^2-, 'dipotassium sulfate') -- a mismatched
      charge ratio (e.g. Al3+ with SO4^2-) would need mixed cation/anion
      counts and stays out of scope, not attempted here.
    - One polyatomic-inorganic anion fragment balanced by a combination of
      *different* cation types (each appearing once, no multiplying prefix
      needed since the cations differ), cited in alphabetical order by
      cation name (P-65.6.2's own worked examples: 'potassium sodium
      butanedioate', 'ammonium potassium hexanedioate') -- tried only
      after the single-cation-type shape above fails, so a molecule with
      several copies of the *same* cation type still prefers the
      multiplying-prefix form."""
    frags = Chem.GetMolFrags(mol, asMols=True)

    for i, cation_frag in enumerate(frags):
        cation = _cation(cation_frag)
        if cation is None:
            continue
        cation_name, charge = cation
        anion_frags = frags[:i] + frags[i + 1 :]
        if len(anion_frags) != charge:
            continue
        anion_charges = [
            sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
            for frag in anion_frags
        ]
        if any(c != -1 for c in anion_charges):
            continue
        for has_shape, namer in _ANION_KINDS:
            if namer in _SINGLY_CHARGED_CATION_ONLY_ANIONS and charge != 1:
                continue
            if not all(has_shape(frag) for frag in anion_frags):
                continue
            smiles = {Chem.MolToSmiles(frag) for frag in anion_frags}
            if len(smiles) != 1:
                continue
            return [cation_name], namer, anion_frags[0], len(anion_frags)

    for i, anion_frag in enumerate(frags):
        polyatomic = _polyatomic_anion(anion_frag)
        if polyatomic is None:
            continue
        namer, magnitude = polyatomic
        cation_frags = frags[:i] + frags[i + 1 :]
        if not cation_frags:
            continue
        cation = _cation(cation_frags[0])
        if cation is None:
            continue
        cation_name, charge = cation
        if magnitude % charge != 0:
            continue
        required_count = magnitude // charge
        if len(cation_frags) != required_count:
            continue
        if any(_cation(frag) != cation for frag in cation_frags):
            continue
        smiles = {Chem.MolToSmiles(frag) for frag in cation_frags}
        if len(smiles) != 1:
            continue
        cation_part = f"{multiplying_prefix(required_count)}{cation_name}" if required_count > 1 else cation_name
        return [cation_part], namer, anion_frag, 1

    for i, anion_frag in enumerate(frags):
        polyatomic = _polyatomic_anion(anion_frag)
        if polyatomic is None:
            continue
        namer, magnitude = polyatomic
        cation_frags = frags[:i] + frags[i + 1 :]
        cations = [_cation(frag) for frag in cation_frags]
        if any(c is None for c in cations):
            continue
        names = [name for name, _ in cations]
        if len(set(names)) != len(names):
            # A repeated cation type here is the previous shape's
            # territory (a multiplying prefix), not this one's.
            continue
        if sum(charge for _, charge in cations) != magnitude:
            continue
        return sorted(names), namer, anion_frag, 1
    return None


def has_salt_shape(mol) -> bool:
    return _split_cation_anions(mol) is not None


def name_salt(mol) -> str:
    cation_parts, namer, anion_frag, anion_count = _split_cation_anions(mol)
    anion_name = namer(anion_frag)

    # A cation name is never itself a compound/substituted prefix (it's
    # always a fixed element/retained-cation-name lookup), so a repeated
    # cation type always multiplies with the plain 'di'/'tri' word, never
    # 'bis'/'tris' -- confirmed against PubChem's own 'disodium carbonate'
    # (P-65.6.2). Several *different* cation types are cited side by side
    # instead, each already alphabetically ordered by `_split_cation_
    # anions` with no multiplying prefix of its own (P-65.6.2's own
    # 'potassium sodium butanedioate'-style worked examples).
    cation_part = " ".join(cation_parts)

    if anion_count == 1:
        return f"{cation_part} {anion_name}"
    # P-16.3.4(a)/P-16.3.5(a): an anion name carrying its own substituent
    # prefix (always locant-leading in this project's convention, since
    # the anion's own suffix carbon is always C1 and never cited) is a
    # compound prefix and multiplies with 'bis'/'tris'/... in parentheses,
    # not the plain 'di'/'tri' used for an unsubstituted anion name (e.g.
    # 'calcium diethanoate', P-65.6.2.1's own worked example) -- confirmed
    # against PubChem's 'calcium bis(2-methyl-2-phenylhexanoate)'.
    is_compound = anion_name[0].isdigit()
    prefix = multiplying_prefix(anion_count, compound=is_compound)
    if is_compound:
        return f"{cation_part} {prefix}({anion_name})"
    return f"{cation_part} {prefix}{anion_name}"
