# Changelog

## Unreleased

- `smiles_to_alternative_names()` and `smiles-to-iupac --alternatives` return the
  general-nomenclature alternatives to the preferred name, starting with
  functional-modifier names (oximes, hydrazones, azines, semicarbazones,
  acetals, ketals and hemiacetals) and ionic skeletal replacement prefixes
  (azonia, thionia, azanylia, boranuida, phosphanida, ...).

## 0.1.0 (alpha)

First public release.

- Rule-based SMILES to IUPAC name conversion implemented from the 2013 IUPAC
  Recommendations (Blue Book), chapters P-1 to P-10.
- `smiles_to_iupac()` Python API, `NonPreferredNameWarning` for names that are
  valid but not preferred IUPAC names, and the `smiles-to-iupac` command line.
- Structures outside the implemented rules raise `NotImplementedError`.
