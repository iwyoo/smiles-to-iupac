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

- Exactly two disconnected fragments (P-77.1.1's general 'cation anion'
  citation order collapses to a single instance of each here).
- The cation must be a single bare Group 1 alkali metal atom with formal
  charge +1 (Li/Na/K/Rb/Cs) -- monoatomic, so its own name is a fixed
  lookup, no substitutive nomenclature involved.
- The anion must be a carboxylate recognized by `_carboxylate.py`
  (`has_carboxylate_shape`/`name_carboxylate`), reused unchanged.

Explicitly out of scope (raise `UnsupportedStructure`, or -- for
`has_salt_shape` -- simply return False so the shape falls through to
every other branch's own, usually less helpful, rejection): more than two
fragments, any 2+/3+ metal cation (needs charge-balancing stoichiometry,
e.g. calcium diacetate), a polyatomic cation (ammonium etc. -- those name
fine standalone already, via `_ammonium.py`, but combining them into a
salt name is unattempted here), any anion other than a plain carboxylate,
and any stoichiometry other than exactly one cation to one anion (a
2:1/1:2 salt needs a multiplying prefix on the appropriate ion's name,
P-16.3.3, not attempted here)."""

from rdkit import Chem

from ._carboxylate import has_carboxylate_shape, name_carboxylate

_MONOATOMIC_CATION_NAMES = {
    "Li": "lithium",
    "Na": "sodium",
    "K": "potassium",
    "Rb": "rubidium",
    "Cs": "cesium",
}


def _monoatomic_cation_name(frag):
    if frag.GetNumAtoms() != 1:
        return None
    atom = frag.GetAtomWithIdx(0)
    if atom.GetFormalCharge() != 1:
        return None
    return _MONOATOMIC_CATION_NAMES.get(atom.GetSymbol())


def _split_cation_anion(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 2:
        return None
    a, b = frags
    for cation_frag, anion_frag in ((a, b), (b, a)):
        cation_name = _monoatomic_cation_name(cation_frag)
        if cation_name is None:
            continue
        anion_charge = sum(atom.GetFormalCharge() for atom in anion_frag.GetAtoms())
        if anion_charge != -1:
            continue
        if not has_carboxylate_shape(anion_frag):
            continue
        return cation_name, anion_frag
    return None


def has_salt_shape(mol) -> bool:
    return _split_cation_anion(mol) is not None


def name_salt(mol) -> str:
    cation_name, anion_frag = _split_cation_anion(mol)
    return f"{cation_name} {name_carboxylate(anion_frag)}"
