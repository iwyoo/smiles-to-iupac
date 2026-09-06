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

Explicitly out of scope (raise `UnsupportedStructure`, or -- for
`has_salt_shape` -- simply return False so the shape falls through to
every other branch's own, usually less helpful, rejection): any metal
cation other than the fixed-valence ones above (transition metals need
Stock/oxidation-number disambiguation, not attempted here), a polyatomic
cation (ammonium etc. -- those name fine standalone already, via
`_ammonium.py`, but combining them into a salt name is unattempted here),
any anion other than a plain carboxylate, mixed/different anions on the
same cation, and more than one cation fragment (multi-cation salts, e.g.
'potassium sodium butanedioate', are unattempted here)."""

from rdkit import Chem

from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._numerals import multiplying_prefix

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


def _split_cation_anions(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    for i, cation_frag in enumerate(frags):
        cation = _monoatomic_cation(cation_frag)
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
        if not all(has_carboxylate_shape(frag) for frag in anion_frags):
            continue
        smiles = {Chem.MolToSmiles(frag) for frag in anion_frags}
        if len(smiles) != 1:
            continue
        return cation_name, anion_frags[0], len(anion_frags)
    return None


def has_salt_shape(mol) -> bool:
    return _split_cation_anions(mol) is not None


def name_salt(mol) -> str:
    cation_name, anion_frag, anion_count = _split_cation_anions(mol)
    anion_name = name_carboxylate(anion_frag)
    if anion_count == 1:
        return f"{cation_name} {anion_name}"
    return f"{cation_name} {multiplying_prefix(anion_count)}{anion_name}"
