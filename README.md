# smiles-to-iupac

A Python library that converts SMILES strings to IUPAC names using only the
nomenclature rules from the *Nomenclature of Organic Chemistry: Recommendations
and Preferred Names 2013* (the "Blue Book") — no machine learning, no lookup
tables of known names.

## Status

Early-stage / pre-alpha. Currently supported: acyclic and simple monocyclic
saturated hydrocarbons (alkanes and cycloalkanes), plus acyclic hydrocarbons
with a single carbon-carbon double or triple bond (alkenes and alkynes),
including branched ("compound") substituent groups (P-29.4), e.g.
`CC(C)C(CC)CCC` → `3-ethyl-2-methylhexane`, `CC1CCCCC1` →
`methylcyclohexane`, `CCCCC(C(C)CC)CCCCC` → `5-(1-methylpropyl)decane`, and
`CCCC(C(C)C)CC=C` → `4-(1-methylethyl)hept-1-ene`. Fused/bridged/spiro ring
systems, unsaturated rings, more than one multiple bond, heteroatoms, and
cyclic substituent groups raise `NotImplementedError` and are future work.
See [REFERENCES.md](REFERENCES.md) for the source of the rules applied.

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
