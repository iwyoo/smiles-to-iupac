# smiles-to-iupac

Rule-based conversion of SMILES strings to IUPAC names, implemented directly
from the *Nomenclature of Organic Chemistry: Recommendations and Preferred
Names 2013* (the "Blue Book"). No machine learning, no lookup table of known
names: every name is constructed by applying the written rules.

```python
>>> from smiles_to_iupac import smiles_to_iupac
>>> smiles_to_iupac("CC(C)C(CC)CCC")
'3-ethyl-2-methylhexane'
>>> smiles_to_iupac("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
'2-[4-(2-methylpropyl)phenyl]propanoic acid'
>>> smiles_to_iupac("CN1CCC[C@H]1c1cccnc1")
'3-[(2S)-1-methylpyrrolidin-2-yl]pyridine'
```

## Guiding principles

Two principles define the scope of this project and every change made to it.

### 1. Anything a SMILES string can express is in scope

The input is a SMILES string, so the goal is to name every molecule that
SMILES can represent: acyclic and cyclic skeletons, fused, bridged and spiro
ring systems, heterocycles, all common characteristic groups, charged and
radical species, isotopes, stereochemistry, and organometallic and
coordination compounds. What SMILES cannot express (for example, a
free-standing locant set, a mixture description, or a polymer repeat unit
with no defined structure) is out of scope by definition.

### 2. Names follow the rules written in the Blue Book

Behavior is derived from the text of the Blue Book rules, not from example
outputs, third-party name generators, or database lookups. Worked examples
and PubChem entries are used only to *verify* a result, never to define it.

- If the Blue Book defines a **preferred IUPAC name (PIN)** for the structure,
  that PIN is returned.
- If the Blue Book defines no PIN (for example, P-69.0 for most
  organometallics), a single rule-valid name is returned and a
  `NonPreferredNameWarning` is emitted, so the caller always knows when the
  result is valid but not preferred.
- Where the Blue Book is itself contradictory, or a structure cannot be
  expressed unambiguously in SMILES, the structure is reported as
  unsupported rather than guessed.

## What it covers

Rule coverage is organized by Blue Book chapter:

| Blue Book chapter | Covered |
|---|---|
| P-1 | General principles: numerals, multiplying prefixes, locants, alphanumerical ordering, enclosing marks |
| P-2 | Parent hydrides: acyclic and cyclic hydrocarbons, heteroatom skeletons (replacement "a" nomenclature), von Baeyer, spiro, fused (ortho-/peri-fused, heterocyclic), bridged, ring assemblies, phanes, fullerenes |
| P-3 | Characteristic groups, substituent prefixes and suffixes, seniority of classes |
| P-4 | Selection of the principal chain/ring, numbering and lowest locants, name assembly |
| P-5 | Constructing names: functional replacement, added/indicated hydrogen, hydro/dehydro prefixes |
| P-6 | Characteristic groups in names: alcohols, amines, carbonyls, acids and their derivatives (esters, amides, anhydrides, halides, nitriles), sulfur/selenium/tellurium and phosphorus acids, organometallic compounds |
| P-7 | Radicals, ions and zwitterions, salts |
| P-8 | Isotopically modified compounds |
| P-9 | Stereodescriptors (`R`/`S`, `E`/`Z`, and related) and stereo-aware retained names |
| P-10 | Natural products: carbohydrates, amino acids, nucleosides and nucleotides, lipids, and the Appendix 3 parent structures (alkaloids, steroids, terpenoids, carotenoids, tetrapyrroles, flavans) with their substituents, unsaturation and configuration |

Beyond organic skeletons it also handles coordination and organometallic
compounds (with `NonPreferredNameWarning` where the Blue Book defines no
PIN), organoelement hydrides (B, Si, P, ...) and onium, carbocation and
carbanion species.

Coverage grows rule by rule and is audited per rule. Structures that fall
outside what is currently implemented raise `NotImplementedError` (see
below); they are never silently mis-named.

## Installation

Requires Python 3.9 or newer. Dependencies (RDKit and NetworkX) are installed
automatically.

```bash
pip install smiles-to-iupac
```

From a clone of the repository:

```bash
pip install -e ".[dev]"
```

## Usage

### Python

```python
import warnings
from smiles_to_iupac import smiles_to_iupac, NonPreferredNameWarning

smiles_to_iupac("c1ccc2[nH]ccc2c1")        # '1H-indole'
smiles_to_iupac("C[C@H](N)C(=O)O")         # 'L-alanine'
smiles_to_iupac("C1CC2CCC1C2")             # 'bicyclo[2.2.1]heptane'
smiles_to_iupac("CS(=O)(=O)O")             # 'methanesulfonic acid'
```

A name returned with `NonPreferredNameWarning` is valid under the Blue Book
rules but is not a preferred IUPAC name:

```python
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    name = smiles_to_iupac("[Fe]")          # 'iron'
    is_pin = not any(
        issubclass(w.category, NonPreferredNameWarning) for w in caught
    )
```

### Errors

| Exception | Meaning |
|---|---|
| `ValueError` | The input is not a valid SMILES string. |
| `NotImplementedError` | The SMILES is valid but uses a rule that is not implemented yet. |

### Command line

```bash
smiles-to-iupac "CC(C)C(CC)CCC"            # 3-ethyl-2-methylhexane
smiles-to-iupac CCO c1ccncc1               # one "<smiles>\t<name>" line each
python3 -m smiles_to_iupac CCO             # equivalent module form
```

With a single SMILES only the name is printed. With several, each line is
`<smiles>\t<name>`. Errors are written to stderr and the exit status is 1 if
any input failed.

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

When adding or changing behavior, derive it from the Blue Book rule text and
cite the section number; use worked examples and PubChem only to check the
result.

## References

See [REFERENCES.md](REFERENCES.md) for the Blue Book source and chapter PDFs.

## License

MIT — see [LICENSE](LICENSE).
