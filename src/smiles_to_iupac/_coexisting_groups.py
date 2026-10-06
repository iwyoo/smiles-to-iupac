"""Shared acyclic-chain naming for two coexisting characteristic-group
classes where one is demoted to a substituent prefix per P-41/P-43
(`_seniority.py`), per the IUPAC 2013 Recommendations ("the Blue Book").

- Each existing pairwise module (`_sulfonic_acid_thiol.py`,
  `_carboxylic_acid_amine.py`, ...) already hardcodes its own seniority
  conclusion and reimplements the senior group's whole chain-search/
  numbering/candidate-selection logic from scratch, copy-pasted from the
  senior group's own single-group module (`_sulfonic_acid.py`/
  `_carboxylic_acid.py`). That logic doesn't depend on which junior group
  is being demoted -- only the {atom_idx -> prefix name} map does. This
  module lets a pairwise module instead call straight into the senior
  module's own acyclic-naming function (`_name_acyclic_sulfonic_acid`,
  `_name_acyclic_carboxylic_acid`), passing the junior group's atoms as
  `extra_names`, so the numbering logic is never duplicated.
- Deliberately narrow first pilot (2026-09-05): only the acyclic (plain
  chain, no benzene-ring substituent) path of `_sulfonic_acid_thiol.py`/
  `_carboxylic_acid_amine.py` is rebuilt on top of this module, to prove
  the pattern without touching the six other existing pairwise modules or
  either pilot's separate benzene-ring-substituent function (regression
  risk). All six pairwise modules' acyclic branches were later migrated
  the same way (PR #507/#508/#509/#510/#511/#512).
- The three pairwise modules that also have a benzene-ring-substituent
  ("phenyl-chain") path (`_carboxylic_acid_sulfinic_acid.py`/
  `_carboxylic_acid_sulfonamide.py`/`_carboxylic_acid_sulfonic_acid.py`)
  still reimplement that path's own chain-search/numbering logic too --
  `name_via_senior_phenyl_chain` below extends the same pattern to it,
  reusing the senior module's own phenyl-chain namer (e.g.
  `_name_phenyl_chain_carboxylic_acid`). A generic N-way dispatcher across
  all ranked classes remains a follow-up.
"""

from ._seniority import senior_class


def name_via_senior_acyclic(senior_acyclic_namer, senior_key, junior_key, senior_args, extra_names, **kwargs):
    """Assert `senior_key` actually outranks `junior_key` per
    `_seniority.senior_class` (loud failure if a future rank-table edit
    ever reversed the pair's assumed direction, mirroring the standalone
    `assert` each existing pairwise module already makes at import time),
    then call `senior_acyclic_namer(*senior_args, extra_names=extra_names,
    **kwargs)` -- the senior single-group module's own acyclic-chain
    naming function, unmodified except for accepting this `extra_names`
    injection point."""
    if senior_class(senior_key, junior_key) != senior_key:
        raise AssertionError(f"expected {senior_key!r} senior to {junior_key!r} per _seniority.py")
    return senior_acyclic_namer(*senior_args, extra_names=extra_names, **kwargs)


def name_via_senior_phenyl_chain(senior_phenyl_chain_namer, senior_key, junior_key, mol, ring_atoms, extra_names, **kwargs):
    """Same idea as `name_via_senior_acyclic`, for the benzene-ring-
    substituent ("phenyl-chain") path instead -- calls
    `senior_phenyl_chain_namer(mol, ring_atoms, extra_names=extra_names,
    **kwargs)` (e.g. `_name_phenyl_chain_carboxylic_acid`), the senior
    single-group module's own phenyl-chain naming function."""
    if senior_class(senior_key, junior_key) != senior_key:
        raise AssertionError(f"expected {senior_key!r} senior to {junior_key!r} per _seniority.py")
    return senior_phenyl_chain_namer(mol, ring_atoms, extra_names=extra_names, **kwargs)
