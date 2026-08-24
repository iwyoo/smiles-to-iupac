"""Naming of the two simplest retained-name heteroaromatic fused ring
systems -- quinoline and indole -- by hardcoded recognition of the
unsubstituted parent hydride only, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-25.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  'quinoline' and 'indole' are retained names for these two benzo-fused
  heteroaromatic ring systems (quinoline = benzo + pyridine, indole =
  benzo + pyrrole).
- P-25.7.1.3 / indicated hydrogen: indole's pyrrole-type nitrogen carries
  an aromatic N-H that isn't structurally forced to one position by the
  ring skeleton alone -- a second, non-aromatic tautomer (3H-indole,
  "indolenine", sp3 at C3) is a distinct real compound -- so the PIN cites
  it explicitly as '1H-indole', not bare 'indole' (confirmed via
  Wikipedia's Indole article, which gives '1H-indole' as the IUPAC name).
  Quinoline needs no such prefix: its nitrogen is pyridine-type (no H,
  degree 2, fully specified by the mancude ring skeleton alone), so
  there's no tautomeric ambiguity to disambiguate.
- Deliberately not built on `_aromatic.py`'s all-carbon ortho-fused-chain
  machinery: `find_aromatic_fused_core` requires every ring atom to be
  carbon, which is false for both compounds by definition, and indole's
  five-membered pyrrole ring also fails that function's "every SSSR ring
  is 6-membered" precondition -- so neither would ever reach that
  module's dispatch at all. As with `_peri_fused_aromatic.py` and
  `_branched_fused_aromatic.py`, an exact whole-molecule canonical-SMILES
  match against each retained name's unsubstituted structure is both
  sufficient and simplest for this narrow scope (no locants to assign).

Formulas cross-checked: quinoline C9H7N, indole (1H-indole) C8H7N; both
reference SMILES independently verified against PubChem's canonical
SMILES for CID 7047 (quinoline) and CID 798 (indole).

Explicitly out of scope: any substituted derivative, any other retained-
name heteroaromatic fused system (purine, carbazole, etc. -- P-25.2.1
lists many more), a second ring heteroatom, and non-aromatic tautomers
(e.g. 3H-indole) or partially saturated forms (e.g. indoline).
`has_retained_heteroaromatic_fused_name` returns False for all of these,
so `core.py`'s existing dispatch continues to raise `UnsupportedStructure`
for them, unchanged.
"""

from rdkit import Chem

_RETAINED_NAME_SMILES = {
    "quinoline": "c1ccc2ncccc2c1",
    "1H-indole": "c1ccc2[nH]ccc2c1",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for name, smiles in _RETAINED_NAME_SMILES.items()}


def has_retained_heteroaromatic_fused_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_retained_heteroaromatic_fused(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
