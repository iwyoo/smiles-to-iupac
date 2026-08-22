# smiles-to-iupac

A Python library that converts SMILES strings to IUPAC names using only the
nomenclature rules from the *Nomenclature of Organic Chemistry: Recommendations
and Preferred Names 2013* (the "Blue Book") — no machine learning, no lookup
tables of known names.

## Status

Early-stage / pre-alpha. Currently supported: acyclic and simple monocyclic
saturated hydrocarbons (alkanes and cycloalkanes), monospiro saturated
hydrocarbons (two carbocyclic rings sharing exactly one atom), saturated
bicyclic hydrocarbons — both fused and bridged (von Baeyer nomenclature,
two carbocyclic rings sharing two or more atoms) — a subset of saturated
tricyclic hydrocarbons (von Baeyer systems with exactly four skeletal atoms
of degree 3, forming a main ring plus a main bridge plus one independent
secondary bridge, P-23.2.5 — including both the "K4" case (adamantane,
twistane) and the "doubled main bridgeheads/secondary bridgeheads" case
(the Blue Book's own tricyclo[4.2.2.2²,⁵]dodecane worked example, and
ortho-fused ring chains like perhydroanthracene)), plus acyclic hydrocarbons
with one or more carbon-carbon double and/or triple bonds on the principal
chain (alkenes, alkynes, dienes/trienes, diynes/triynes, and mixed enynes),
including branched ("compound") substituent groups (P-29.4), e.g.
`CC(C)C(CC)CCC` → `3-ethyl-2-methylhexane`, `CC1CCCCC1` →
`methylcyclohexane`, `C1CCCC12CCCCC2` → `spiro[4.5]decane`,
`C1CC2CCC1C2` → `bicyclo[2.2.1]heptane`, `C1CCC2CCCCC2C1` →
`bicyclo[4.4.0]decane`, `C1C2CC3CC1CC(C2)C3` (adamantane) →
`tricyclo[3.3.1.1^3,7]decane`, `CCCCC(C(C)CC)CCCCC` →
`5-(1-methylpropyl)decane`, `CCCC(C(C)C)CC=C` →
`4-(1-methylethyl)hept-1-ene`, `C=CC=C` → `buta-1,3-diene`, and `C=CC#C` →
`but-1-en-3-yne`. Fluoro, chloro, bromo, and iodo substituents (P-35.2.1)
are supported on any of the above parent hydrides, e.g. `CCCF` →
`1-fluoropropane`, `C(Cl)(Cl)(Cl)Cl` → `tetrachloromethane`, `BrCC(Cl)C(F)CI`
→ `1-bromo-2-chloro-3-fluoro-4-iodobutane`, `ClC1CCCCC1` →
`chlorocyclohexane`, and `C=C(Cl)CC` → `2-chlorobut-1-ene`,
`C1CC2CCC1C1CCC2CC1` → `tricyclo[4.2.2.2^2,5]dodecane`,
`C1CCC2CC3CCCCC3CC2C1` (perhydroanthracene) → `tricyclo[8.4.0.0^3,8]tetradecane`,
plus propellane-type tricyclics — two skeletal branch atoms of degree 4,
directly bonded to each other, joined by three further bridges, with that
direct bond treated as a zero-length independent secondary bridge — e.g.
`C1C23CC12C3` ([1.1.1]propellane) → `tricyclo[1.1.1.0^1,3]pentane`.
Other tricyclic topologies, tetracyclic-and-higher/polyspiro ring systems,
unsaturated rings, a multiple bond that isn't on any candidate principal
chain, heteroatoms other than the four halogens above, and cyclic
substituent groups raise `NotImplementedError` and are future work. See
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
