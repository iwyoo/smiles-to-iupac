"""Naming of di-/polynuclear noncarbon oxoacids with their own P-67.2.1
preselected names, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.2.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf,
  ~4784-4930): 'diphosphoric acid' (pyrophosphoric acid) is the
  preselected name for (HO)2P(O)-O-P(O)(OH)2 -- like the mononuclear
  noncarbon oxoacids `_phosphate.py`/etc. already cover, di-/polynuclear
  ones are retained/preselected names, not built compositionally
  ("Substitutive or additive names are not recommended"), so this is an
  exact whole-molecule-shape match (mirroring `_hetero_monocyclic.py`'s
  own retained-name-table approach) rather than a generalization of
  `_phosphate.py`'s existing single-center substitutive-ester logic.
  Confirmed against a real PubChem structure, CID 1023 (`OP(=O)(O)OP(=O)
  (O)O`) -- PubChem's own computed name for it ('phosphono dihydrogen
  phosphate') uses a different, non-retained ester-style construction,
  not trusted here (see this project's prior PubChem-vs-Blue-Book
  divergence findings elsewhere).

Scope: exact-match unsubstituted diphosphoric acid only -- the first of
~15 similar P-67.2.1 preselected dinuclear-acid names (hypodiphosphoric/
diphosphonic/hypodiphosphonic/diphosphorous/hypodiphosphorous acids for
phosphorus alone, plus the parallel As/Sb/S series and the boron/silicon
pair); the rest are each their own follow-up step, batch-extending this
same table once the recognition mechanism is proven here.
"""

from rdkit import Chem

_CANONICAL_TO_NAME = {
    Chem.CanonSmiles("OP(=O)(O)OP(=O)(O)O"): "diphosphoric acid",
}


def has_dinuclear_oxoacid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_dinuclear_oxoacid(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
