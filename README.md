# chemonym

A Python library that converts SMILES strings to IUPAC names using only the
nomenclature rules from the *Nomenclature of Organic Chemistry: Recommendations
and Preferred Names 2013* (the "Blue Book") — no machine learning, no lookup
tables of known names.

## Status

Early-stage / pre-alpha. Currently supported: acyclic and simple monocyclic
saturated hydrocarbons (alkanes and cycloalkanes) whose substituents are
themselves unbranched, e.g. `CC(C)C(CC)CCC` → `3-ethyl-2-methylhexane` and
`CC1CCCCC1` → `methylcyclohexane`. Fused/bridged/spiro ring systems,
unsaturation, heteroatoms, and branched ("compound") substituent groups raise
`NotImplementedError` and are future work. See [REFERENCES.md](REFERENCES.md)
for the source of the rules applied.

## Installation

```bash
pip install -e ".[dev]"
```

## Usage

```python
from chemonym import smiles_to_iupac

name = smiles_to_iupac("CC(C)C(CC)CCC")  # "3-ethyl-2-methylhexane"
```

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License

MIT — see [LICENSE](LICENSE).
