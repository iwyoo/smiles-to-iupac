"""Seniority-of-suffix ranking (P-41 Table 4.1 / P-43 Table 4.4) for the
characteristic-group classes already implemented as their own suffix
module in this codebase.

- P-41 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): "The
  selection of a preferred parent structure is based on the seniority of
  classes ... which gives priority first to characteristic groups
  expressed as suffixes." Table 4.1 ranks the general compound classes;
  Table 4.4 (P-43) gives the complete numbered seniority order of the
  individual suffixes within/across those classes. When two different
  characteristic groups that could each be cited as a suffix are both
  present on one molecule, only the more senior one actually is -- the
  rest are cited as substituent prefixes instead (P-41.1).
- This module is a single, reusable rank table for that suffix-vs-suffix
  comparison, extracted from the numbered Table 4.4 list (see
  `tmp/bluebook/P41.txt` for the cached excerpt this was built from).
  Lower rank number = more senior. Only classes this codebase already has
  a dedicated suffix module for are listed -- a class with no module here
  yet has no suffix to demote *to* a prefix in the first place, so it
  would be meaningless to rank it.
- This does NOT yet replace the five existing pairwise "X+Y coexist"
  modules (`_aldehyde_ketone.py`, `_ketone_ester.py`, `_ketone_amide.py`,
  `_aldehyde_carboxylic_acid.py`, `_carboxylic_acid_amine.py`), which each
  still hardcode their own seniority conclusion directly in prose (their
  module docstrings each cite "Table 3.3 ranks 'X' senior to 'Y'"
  verbatim) -- refactoring them to consult this table instead is a
  separate follow-up, not done here to avoid regression risk on modules
  that already work. `_sulfonic_acid_thiol.py` is the first new pairwise
  module built directly on top of this table instead.
- Chalcogen analogues of a listed class (e.g. thioamide/selenoamide for
  amide, thiol/selenol/tellurol for alcohol, thione/selone/tellone for
  ketone) share their parent class's rank -- Table 4.4 groups them as a
  single numbered entry (P-43.0's functional-replacement note), and this
  codebase already treats them as structurally interchangeable elsewhere
  (e.g. `_thiol.py` mirrors `_alcohol.py`).
"""

# Table 4.4's own numbering (not renumbered/compacted here, so a rank gap
# between two entries is not itself meaningful -- only relative order is).
SUFFIX_CLASS_RANK = {
    "carboxylic_acid": 1,
    "sulfonic_acid": 4,
    "sulfinic_acid": 9,
    "selenonic_acid": 12,
    "seleninic_acid": 13,
    "telluronic_acid": 14,
    "tellurinic_acid": 15,
    "amide": 16,
    "sulfonamide": 19,
    "sulfinamide": 24,
    "hydrazide": 31,
    "nitrile": 46,
    "aldehyde": 47,
    "ketone": 48,  # thione/selone/tellone share this rank (chalcogen analogues)
    "alcohol": 49,  # thiol/selenol/tellurol share this rank (chalcogen analogues)
    "hydroperoxide": 50,
    "amine": 51,
    "imine": 52,
}


def senior_class(a: str, b: str) -> str:
    """Return whichever of `a`/`b` is senior per Table 4.1/4.4 (P-41/P-43).
    Both must be keys of `SUFFIX_CLASS_RANK` -- raises `KeyError`
    otherwise, since an unranked class has no suffix here to compare."""
    return a if SUFFIX_CLASS_RANK[a] < SUFFIX_CLASS_RANK[b] else b
