# Changelog

## Unreleased

- The library's own errors are exported and share the base class `SmilesToIupacError`:
  `UnsupportedStructure` (also a `NotImplementedError`) is a rule that is not implemented yet;
  `RuleUndefinedStructure` (no longer an `UnsupportedStructure`) is a structure for which the
  Blue Book defines no single name (undefined or self-contradictory); `InternalInvariantError`
  (also a `RuntimeError`) is a bug in the naming pipeline; `InvalidSmilesError`
  (also a `ValueError`) is an input that is not a valid SMILES string.
- `smiles_to_alternative_names()` and `smiles-to-iupac --alternatives` return the
  general-nomenclature alternatives to the preferred name, starting with
  functional-modifier names (oximes, hydrazones, azines, semicarbazones,
  acetals, ketals and hemiacetals) and ionic skeletal replacement prefixes
  (azonia, thionia, azanylia, boranuida, phosphanida, ...) and conjunctive
  names (benzeneacetic acid, naphthalene-2-propanol, ...).

## 0.1.0 (alpha)

First public release.

- Rule-based SMILES to IUPAC name conversion implemented from the 2013 IUPAC
  Recommendations (Blue Book), chapters P-1 to P-10.
- `smiles_to_iupac()` Python API, `NonPreferredNameWarning` for names that are
  valid but not preferred IUPAC names, and the `smiles-to-iupac` command line.
- Structures outside the implemented rules raise `NotImplementedError`.
