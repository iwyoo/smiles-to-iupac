# chemonym

A Python library that converts SMILES strings to IUPAC names using only the
nomenclature rules from the *Nomenclature of Organic Chemistry: Recommendations
and Preferred Names 2013* (the "Bluebook") — no machine learning, no lookup
tables of known names.

## Status

Early-stage / pre-alpha. The naming engine is not implemented yet.

## Installation

```bash
pip install -e ".[dev]"
```

## Usage

```python
from chemonym import smiles_to_iupac

name = smiles_to_iupac("CCO")
```

## Development

```bash
pip install -e ".[dev]"
pytest
```

## License

MIT — see [LICENSE](LICENSE).
