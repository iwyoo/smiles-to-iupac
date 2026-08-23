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
ortho-fused ring chains like perhydroanthracene)), a subset of saturated
tetracyclic hydrocarbons (von Baeyer systems with exactly six skeletal
atoms of degree 3, forming a main bicyclic system plus two independent
secondary bridges, P-23.2.6, e.g. quadricyclane), plus acyclic hydrocarbons
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
`C1C23CC12C3` ([1.1.1]propellane) → `tricyclo[1.1.1.0^1,3]pentane`, and
`C1C2C3C2C4C1C34` (quadricyclane) → `tetracyclo[3.2.0.0^2,7.0^4,6]heptane`,
and linear (unbranched) polyspiro saturated hydrocarbons — three or more
carbocyclic rings connected in a chain by spiro atoms, each internal ring
sharing exactly one spiro atom with each of its two neighbors (P-24.2.2) —
e.g. `C1CCC12CCC3(CC2)CCC3` → `dispiro[3.2.3^7.2^4]dodecane`.
Also supported: ortho-fused aromatic (mancude) six-membered all-carbon ring
chains carrying one of six retained names — `benzene` (one ring),
`naphthalene` (two rings), `anthracene`/`phenanthrene` (three rings,
straight/angular, with their traditional fixed numbering per P-25.3.3), and
`tetracene`/`pentacene` (four/five rings, straight) — with halogen and
alkyl substituents, e.g. `c1ccccc1` → `benzene`,
`c1ccc2ccccc2c1` → `naphthalene`, `C1=CC=C2C=C(C=CC2=C1)Cl` →
`2-chloronaphthalene`.

Alcohols (P-33.2, the '-ol' suffix) are also supported on acyclic and
simple monocyclic saturated skeletons, and on acyclic skeletons with
existing double/triple-bond support, e.g. `CCO` → `ethanol`, `CC(O)C` →
`propan-2-ol`, `OCCO` → `ethane-1,2-diol`, `OC1CCCCC1` → `cyclohexanol`,
and `OCCCC=C` → `pent-4-en-1-ol`. Primary amines (P-33.1, the '-amine'
suffix) are supported on the same scope of skeletons, e.g. `CCN` →
`ethanamine`, `CC(N)C` → `propan-2-amine`, `NCCN` → `ethane-1,2-diamine`,
`NC1CCCCC1` → `cyclohexanamine`, and `NCCCC=C` → `pent-4-en-1-amine`.
Ketones (P-33.4, the '-one' suffix) are supported on the same scope of
skeletons, e.g. `CC(=O)C` → `propan-2-one`, `CCCC(=O)CC` → `hexan-3-one`,
`CC(=O)CC(=O)C` → `pentane-2,4-dione`, `O=C1CCCCC1` → `cyclohexanone`, and
`CC(=O)C=CC` → `pent-3-en-2-one`. Aldehydes (P-33.3, the '-al' suffix) are
supported on acyclic skeletons only (a ring-bound -CHO uses the different
'carbaldehyde' suffix pattern, out of scope) — since the -CHO carbon is
always the chain terminus, its own locant is never cited, e.g. `CCC=O` →
`propanal`, `CC(C)C=O` → `2-methylpropanal`, `O=CCCCC=O` → `pentanedial`,
and `C=CCCC=O` → `pent-4-enal`. Carboxylic acids (P-65.1.1, the '-oic
acid' suffix) are supported on acyclic saturated or unsaturated carbon
chains, e.g. `CC(=O)O` → `ethanoic acid`, `CC(C)CC(=O)O` →
`3-methylbutanoic acid`, `OC(=O)CCCCC(=O)O` → `hexanedioic acid`, and
`CC=CC(=O)O` → `but-2-enoic acid`.

Other tricyclic topologies, other tetracyclic topologies (any branch atom
of degree 4, or fewer/more than six skeletal atoms of degree 3),
pentacyclic-and-higher ring systems, branched polyspiro ring systems
(P-24.2.3, a spiro atom shared by three or more rings) and heterocyclic
spiro ring systems (P-24.2.4), unsaturated non-aromatic rings, a multiple
bond that isn't on any candidate principal chain, heteroatoms other than
the four halogens, hydroxyl oxygen, primary amine nitrogen, ketone carbonyl
oxygen, aldehyde carbonyl oxygen, and carboxylic-acid oxygens above,
secondary/tertiary amines, characteristic groups more senior than a plain
alcohol/amine/ketone/aldehyde/carboxylic acid (e.g. esters, amides,
nitriles), a ring-bound aldehyde (the 'carbaldehyde' suffix), a carboxylic
acid on/in a ring, aryl ketones/aldehydes, cyclic substituent groups,
peri-fused or branched aromatic ring systems
(e.g. pyrene, triphenylene), heteroaromatic rings, and any other
ortho-fused aromatic ring chain (angular chains of four or more rings,
chains of six or more rings, or anything needing genuine
`benzo[x,y-z]fusion[...]` name construction) raise `NotImplementedError`
and are future work. See [REFERENCES.md](REFERENCES.md) for the source of
the rules applied.

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
