# smiles-to-iupac

A Python library that converts SMILES strings to IUPAC names using only the
nomenclature rules from the *Nomenclature of Organic Chemistry: Recommendations
and Preferred Names 2013* (the "Blue Book") — no machine learning, no lookup
tables of known names.

## Status

Early-stage / pre-alpha. Currently supported: acyclic and simple monocyclic
saturated hydrocarbons (alkanes and cycloalkanes), monospiro saturated
hydrocarbons (two carbocyclic rings sharing exactly one atom), plus acyclic
hydrocarbons with one or more carbon-carbon double and/or triple bonds on the
principal chain (alkenes, alkynes, dienes/trienes, diynes/triynes, and mixed
enynes), including branched ("compound") substituent groups (P-29.4), e.g.
`CC(C)C(CC)CCC` → `3-ethyl-2-methylhexane`, `CC1CCCCC1` →
`methylcyclohexane`, `C1CCCC12CCCCC2` → `spiro[4.5]decane`,
`CCCCC(C(C)CC)CCCCC` → `5-(1-methylpropyl)decane`,
`CCCC(C(C)C)CC=C` → `4-(1-methylethyl)hept-1-ene`, `C=CC=C` →
`buta-1,3-diene`, and `C=CC#C` → `but-1-en-3-yne`. Fused/bridged/polyspiro
ring systems, unsaturated rings, a multiple bond that isn't on any candidate
principal chain, heteroatoms, and cyclic substituent groups raise
`NotImplementedError` and are future work. See
[REFERENCES.md](REFERENCES.md) for the source of the rules applied.

## Installation

```bash
pip install -e ".[dev]"
```

## Usage

```python
from smiles_to_iupac import smiles_to_iupac

name = smiles_to_iupac("CC(C)C(CC)CCC")  # "3-ethyl-2-methylhexane"
```

## CLI

```bash
smiles-to-iupac CC(C)C(CC)CCC  # "3-ethyl-2-methylhexane"
```

Multiple SMILES print one `<smiles>\t<name>` line each. Also runnable as
`python3 -m smiles_to_iupac <smiles>`.

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License

MIT — see [LICENSE](LICENSE).
