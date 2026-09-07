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

The anion may also be a single alkoxide, recognized by `_alkoxide.py`
(`has_alkoxide_shape`/`name_alkoxide`, reused unchanged), but -- unlike
carboxylate -- only paired with a singly-charged cation (alkali metal or
ammonium): PubChem-verified `[Na+].C[O-]` -> 'sodium methanolate' and
`[NH4+].C[O-]` -> 'azanium methanolate', both exactly one cation to one
anion. A 2+/3+ metal cation with two or three alkoxide anions is *not*
attempted here -- PubChem's own generator inconsistently omits the
multiplying prefix there (`[Ca+2].C[O-].C[O-]` -> 'calcium methanolate',
not 'calcium dimethanolate', unlike the carboxylate case's confirmed
'calcium diacetate'), so there is no reliable worked example to verify
the multivalent-cation form of an alkoxide salt against yet.

Explicitly out of scope (raise `UnsupportedStructure`, or -- for
`has_salt_shape` -- simply return False so the shape falls through to
every other branch's own, usually less helpful, rejection): any metal
cation other than the fixed-valence ones above (transition metals need
Stock/oxidation-number disambiguation, not attempted here), any anion
other than a plain carboxylate or (singly-charged-cation-only) alkoxide,
mixed/different anions on the same cation, a 2+/3+ metal cation paired
with an alkoxide anion (see above), and more than one cation fragment
(multi-cation salts, e.g. 'potassium sodium butanedioate', are
unattempted here)."""

from rdkit import Chem

from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._ammonium import has_ammonium_shape, name_ammonium
from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._numerals import multiplying_prefix

_ANION_KINDS = [
    (has_carboxylate_shape, name_carboxylate),
    (has_alkoxide_shape, name_alkoxide),
]

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
        return name_ammonium(frag), 1
    return None


def _split_cation_anions(mol):
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
            if namer is name_alkoxide and charge != 1:
                continue
            if not all(has_shape(frag) for frag in anion_frags):
                continue
            smiles = {Chem.MolToSmiles(frag) for frag in anion_frags}
            if len(smiles) != 1:
                continue
            return cation_name, namer, anion_frags[0], len(anion_frags)
    return None


def has_salt_shape(mol) -> bool:
    return _split_cation_anions(mol) is not None


def name_salt(mol) -> str:
    cation_name, namer, anion_frag, anion_count = _split_cation_anions(mol)
    anion_name = namer(anion_frag)
    if anion_count == 1:
        return f"{cation_name} {anion_name}"
    return f"{cation_name} {multiplying_prefix(anion_count)}{anion_name}"
