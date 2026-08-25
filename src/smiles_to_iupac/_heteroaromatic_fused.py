"""Naming of retained-name heteroaromatic fused ring systems -- quinoline,
indole, and the benzofuran/benzothiophene isomer pairs -- by hardcoded
recognition of the unsubstituted parent hydride only, per the IUPAC 2013
Recommendations ("the Blue Book"):

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
- benzofuran/benzothiophene (P-25.2.1's benzo + furan/thiophene fusion):
  originally investigated as a general "compute the fusion locant letter"
  case (P-25.3.1.3/FR-4.1's italic-letter mechanism, e.g. 'benzo[b]furan'
  vs the isomeric 'benzo[c]furan'), but the actual PIN for this specific
  retained-name pair turns out to use a different, simpler convention: a
  leading peripheral *numeral* (not a bracketed letter) disambiguates the
  two isomers -- '1-benzofuran' (the common isomer, O adjacent to the
  ring-fusion carbon) vs '2-benzofuran' ("isobenzofuran", O at the "meso"
  position between the two fusion carbons), and analogously
  '1-benzothiophene'/'2-benzothiophene' for the sulfur analogues
  (confirmed against PubChem CID 9223's canonical SMILES for
  1-benzofuran, and cross-checked against Wikipedia's Benzofuran/
  Isobenzofuran/Benzothiophene articles, all of which give these as PINs
  -- 'benzo[b]furan'/'benzo[c]furan' are valid general-nomenclature names
  but not the PINs). Since there's no general letter-locant computation
  involved after all, and only these four whole molecules are verified,
  they're recognized the same hardcoded way as quinoline/indole rather
  than via a general fusion-locant algorithm.
- Deliberately not built on `_aromatic.py`'s all-carbon ortho-fused-chain
  machinery: `find_aromatic_fused_core` requires every ring atom to be
  carbon, which is false for all four compounds by definition, and the
  five-membered second ring in each also fails that function's "every
  SSSR ring is 6-membered" precondition -- so none would ever reach that
  module's dispatch at all. As with `_peri_fused_aromatic.py` and
  `_branched_fused_aromatic.py`, an exact whole-molecule canonical-SMILES
  match against each retained name's unsubstituted structure is both
  sufficient and simplest for this narrow scope (no locants to assign).

Formulas cross-checked: quinoline C9H7N, indole (1H-indole) C8H7N,
benzofuran/isobenzofuran C8H6O, benzothiophene/isobenzothiophene C8H6S;
quinoline/indole reference SMILES independently verified against
PubChem's canonical SMILES for CID 7047 (quinoline) and CID 798 (indole);
1-benzofuran's reference SMILES independently verified against PubChem's
canonical SMILES for CID 9223.

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
    "1-benzofuran": "c1ccc2occc2c1",
    "2-benzofuran": "c1ccc2cocc2c1",
    "1-benzothiophene": "c1ccc2sccc2c1",
    "2-benzothiophene": "c1ccc2cscc2c1",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for name, smiles in _RETAINED_NAME_SMILES.items()}


def has_retained_heteroaromatic_fused_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_retained_heteroaromatic_fused(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
