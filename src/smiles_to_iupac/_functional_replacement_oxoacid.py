"""Naming of mononuclear noncarbon oxoacids modified by functional-
replacement (infix) nomenclature, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-67.1.2.3.1/P-67.1.4.1.1.4 (Chapter P-6a, https://iupac.qmul.ac.uk/
  BlueBook/PDF/P6a.pdf, ~3520-3531, ~4245-4256): thiophosphoric acid
  (`P(=S)(OH)3`, chalcogen substitution only at the original =O position,
  every original -OH kept) is named with the 'thio' infix inserted into
  the acid stem, plus an italic tautomer-locant prefix ('O,O,O-') citing
  that all three non-substituted positions are the ordinary -OH tautomer
  -- confirmed directly, `phosphorothioic O,O,O-acid (preselected name)`,
  the Blue Book (not just by analogy to arsenic's own
  `arsorothioic O,O,O-acid`, ~3681, which uses the identical suffix
  pattern too). Like the mononuclear noncarbon oxoacids `_phosphate.py`/
  etc. and the dinuclear ones in `_dinuclear_oxoacid.py` already cover,
  this is an exact whole-molecule-shape match (mirroring
  `_hetero_monocyclic.py`'s own retained-name-table approach), not a
  generalization of any substitutive-naming module. Confirmed against a
  real PubChem structure, CID 167254 (`OP(=S)(O)O`) -- PubChem's own
  computed name for it (`trihydroxy(sulfanylidene)-lambda5-phosphane`)
  uses the λ-convention instead of this retained-plus-infix construction,
  not trusted here (see this project's prior PubChem-vs-Blue-Book
  divergence findings elsewhere).

Scope: exact-match unsubstituted thiophosphoric acid only -- the first of
many P-67.1.2.1 infix-modifiable mononuclear noncarbon oxoacids (boric/
azonic/phosphoric/phosphorous/arsoric/stiboric/sulfuric/selenic/telluric
acid families, each combinable with 7 chalcogen-replacement infixes); the
rest are each their own follow-up step, batch-extending this same table
once the recognition mechanism is proven here.
"""

from rdkit import Chem

_CANONICAL_TO_NAME = {
    Chem.CanonSmiles("OP(=S)(O)O"): "phosphorothioic O,O,O-acid",
}


def has_functional_replacement_oxoacid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_functional_replacement_oxoacid(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
