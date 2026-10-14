"""Shared chain naming for two coexisting characteristic-group classes where one is demoted to a prefix (P-41, P-43; `_seniority.py`).

The senior group's own namer is called with the junior group's atoms injected as `extra_names`, so its numbering logic is never
duplicated per pair; `name_via_senior_phenyl_chain` does the same for a chain hanging off a benzene ring.
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
