"""Shared mechanism for citing a single principal-characteristic-group
suffix locant on an already-numbered von Baeyer bicyclic/polycyclic
(P-23.2.1) or monospiro (P-24.2.1) parent hydride.

Extracted from `_alcohol.py`'s original `_name_von_baeyer_alcohol`/
`_name_monospiro_alcohol` (#782, M1), which already implemented this end
to end for -OH: given the parent's already-established skeleton
numbering, the only per-suffix pieces are (1) which atom carries the
suffix, (2) the excluded-heteroatom set passed to `substituents_for_ring`
so the suffix atom itself is never also treated as a substituent branch,
(3) the suffix word appended to the parent stem, and (4) the noun used in
the "not on the ring system itself" error message. Everything else --
candidate iteration, substituent enumeration, ranking via each shape's own
`_candidate_key`'s `suffix_locant` parameter -- is suffix-agnostic (#822,
M2 step 1).

Lives in its own module rather than `_common.py` because `_bicyclic.py`/
`_polycyclic.py`/`_spiro.py` already import from `_common.py` -- importing
them back from there would be circular.
"""

from ._bicyclic import (
    _candidate_key as _bicyclic_candidate_key,
    bicyclic_parent_name,
    iter_bicyclic_numberings,
)
from ._common import (
    ENE_BOND_ORDER,
    YNE_BOND_ORDER,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    stereo_locants_prefix,
    von_baeyer_unsaturation_citations,
)
from ._polycyclic import (
    _candidate_key as _polycyclic_candidate_key,
    iter_polycyclic_candidates,
)
from ._spiro import _candidate_key as _spiro_candidate_key, iter_monospiro_numberings
from ._substituents import substituents_for_ring
from ._unsaturated import _unsaturation_suffix_from_citations


def von_baeyer_core_atoms(bicyclic_core, polycyclic_core):
    if bicyclic_core is not None:
        bh1, bh2, bridges = bicyclic_core
        return {bh1, bh2} | {atom for bridge in bridges for atom in bridge}
    branch_atoms, bridges = polycyclic_core
    return set(branch_atoms) | {atom for _, _, path in bridges for atom in path}


def _suffixed_parent(base_parent, locant, suffix_word, elide_e):
    """P-16.3.3: the parent stem's final 'e' is elided before a
    vowel-initial suffix ('heptane' + '-2-ol' -> 'heptan-2-ol') but kept
    before a consonant-initial one ('heptane' + '-2-thiol' ->
    'heptane-2-thiol', confirmed by `_thiol.py`'s own module docstring
    and existing acyclic/cyclic naming, PubChem CID 13487780). `elide_e`
    defaults to True (`-ol`/`-amine`/`-one`, every suffix piloted so far)
    so existing callers are unaffected."""
    stem = base_parent[:-1] if elide_e else base_parent
    return stem + f"-{locant}-{suffix_word}"


def _validate_bond_orders(bonds):
    invalid_orders = [order for _, _, order in bonds if order not in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if invalid_orders:
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )


def _suffixed_parent_and_key(base_parent, locant, suffix_word, elide_e, bonds, full_order):
    """Builds the suffixed parent name and its `_candidate_key` tie-break
    slot, composing ring unsaturation (P-31.1.4.1/.2) with the suffix
    locant when `bonds` is given (#1081, M6 step 1) -- mirrors
    `_bicyclic._name_bicyclic_unsaturated`/`_polycyclic._name_polycyclic_
    unsaturated`'s own unsaturated-stem construction and tie-break exactly,
    except the suffix atom's own locant is compared first and the
    unsaturation-citation quality only breaks ties among numberings that
    already agree on it (P-31.1.4.3's "lowest locant to the principal
    characteristic group suffix first" precedent, `tests/test_ketone_
    hydroxyl_monospiro.py`). `base_parent`: the plain saturated parent
    name ending in '...ane' (`bicyclic_parent_name`'s or `iter_polycyclic_
    candidates`'s own output)."""
    if not bonds:
        return _suffixed_parent(base_parent, locant, suffix_word, elide_e), locant

    position = {atom: i + 1 for i, atom in enumerate(full_order)}
    ene_citations, yne_citations, compound_count, primary_locants, full_locants = (
        von_baeyer_unsaturation_citations(position, bonds)
    )
    primary_locant_set = lowest_locant_set(primary_locants)
    full_locant_set = lowest_locant_set(full_locants)
    ene_locant_set = lowest_locant_set([bond_locant for bond_locant, _ in ene_citations])
    body, needs_stem_a = _unsaturation_suffix_from_citations(ene_citations, yne_citations)
    unsaturated_stem = base_parent[:-3] + ("a" if needs_stem_a else "") + "-" + body
    parent = _suffixed_parent(unsaturated_stem, locant, suffix_word, elide_e)
    key = (locant, compound_count, primary_locant_set, full_locant_set, ene_locant_set)
    return parent, key


def _with_stereo_prefix(stereo, full_order, name):
    """Prefix `name` with its P-92 stereodescriptor(s), once the winning
    candidate's own `full_order` numbering is known -- mirrors `_spiro.py`'s
    `name_monospiro` (the un-suffixed monospiro parent naming), which
    already proved this exact pattern for a bare ring system. `stereo`:
    the caller's own `specified_stereocenters(mol)` result, or None."""
    if stereo is None:
        return name
    position_of = {atom: i + 1 for i, atom in enumerate(full_order)}
    return stereo_locants_prefix(stereo, position_of) + name


def name_von_baeyer_suffix(
    mol,
    suffix_carbon,
    excluded,
    suffix_word,
    noun,
    bicyclic_core,
    polycyclic_core,
    ring_count,
    elide_e=True,
    extra_substituents=None,
    stereo=None,
    bonds=None,
):
    """`suffix_carbon`: the single ring atom the suffix is attached to.
    `excluded`: heteroatom indices to keep out of `substituents_for_ring`'s
    own substituent enumeration (the suffix group's own atom(s)).
    `suffix_word`: appended after the parent stem (see `_suffixed_parent`
    for the `elide_e` rule, e.g. 'ol' -> 'bicyclo[2.2.1]heptan-2-ol',
    'thiol' -> 'bicyclo[2.2.1]heptane-2-thiol'). `noun`: used only in the
    not-on-the-ring-system error message (e.g. 'hydroxyl', 'amine').
    `extra_substituents`: an optional `{atom_idx: name}` dict merged
    alongside the usual halogen substituents (e.g. a coexisting hydroxyl
    cited as a prefix, mirroring `_name_cyclic_ketone`'s identical
    monocyclic pattern, #1029) -- unlike `excluded`, these atoms are NOT
    excluded from `substituents_for_ring`'s own branch-root search; they
    must still be found as a branch root there so this dict's name gets
    attached to them. `stereo`: optional `specified_stereocenters(mol)`
    result, cited via `_with_stereo_prefix` against the winning candidate's
    own numbering (#1078, M5 step 1). `bonds`: optional `non_single_bonds
    (mol)` result, composed into the parent stem and numbering tie-break
    via `_suffixed_parent_and_key` (#1081, M6 step 1)."""
    core_atoms = von_baeyer_core_atoms(bicyclic_core, polycyclic_core)
    if suffix_carbon not in core_atoms:
        raise UnsupportedStructure(
            f"a {noun} not on the bicyclic/polycyclic ring system itself "
            "(e.g. on a substituent branch) is not supported yet"
        )
    if bonds:
        _validate_bond_orders(bonds)

    halogens = {**halogen_substituents(mol), **(extra_substituents or {})}
    graph = adjacency(mol)
    best_key = None
    best_order = None
    if bicyclic_core is not None:
        base_parent = bicyclic_parent_name(bicyclic_core)
        for full_order in iter_bicyclic_numberings(bicyclic_core):
            locant = full_order.index(suffix_carbon) + 1
            substituents = substituents_for_ring(graph, full_order, halogens, excluded, mol=mol)
            parent, suffix_key = _suffixed_parent_and_key(base_parent, locant, suffix_word, elide_e, bonds, full_order)
            key = _bicyclic_candidate_key(parent, substituents, suffix_locant=suffix_key)
            if best_key is None or key < best_key:
                best_key, best_order = key, full_order
        return _with_stereo_prefix(stereo, best_order, best_key[-1])

    for full_order, parent, outer_key in iter_polycyclic_candidates(polycyclic_core, ring_count):
        locant = full_order.index(suffix_carbon) + 1
        substituents = substituents_for_ring(graph, full_order, halogens, excluded, mol=mol)
        suffixed_parent, suffix_key = _suffixed_parent_and_key(parent, locant, suffix_word, elide_e, bonds, full_order)
        key = outer_key + _polycyclic_candidate_key(suffixed_parent, substituents, suffix_locant=suffix_key)
        if best_key is None or key < best_key:
            best_key, best_order = key, full_order
    if best_key is None:
        raise UnsupportedStructure(
            "this polycyclic topology is not supported yet (disjoint ring "
            "systems joined only by an acyclic linker are out of scope; "
            "see _polycyclic.py's name_polycycloalkane for the analogous "
            "non-suffix guard)"
        )
    return _with_stereo_prefix(stereo, best_order, best_key[-1])


def name_monospiro_suffix(
    mol,
    suffix_carbon,
    excluded,
    suffix_word,
    noun,
    spiro_atom,
    elide_e=True,
    extra_substituents=None,
    stereo=None,
):
    """Same mechanism as `name_von_baeyer_suffix`, for a monospiro
    skeleton (`_spiro.py`). `extra_substituents`: see that function's own
    docstring (#1029) -- an optional `{atom_idx: name}` dict merged
    alongside the usual halogen substituents, e.g. a coexisting hydroxyl
    cited as a prefix. `stereo`: see `name_von_baeyer_suffix`'s own
    docstring (#1078, M5 step 1)."""
    ring_atoms = {atom for ring in mol.GetRingInfo().AtomRings() for atom in ring}
    if suffix_carbon not in ring_atoms:
        raise UnsupportedStructure(
            f"a {noun} not on the monospiro ring system itself (e.g. on a "
            "substituent branch) is not supported yet"
        )

    halogens = {**halogen_substituents(mol), **(extra_substituents or {})}
    graph = adjacency(mol)
    best_key = None
    best_order = None
    for parent, full_order in iter_monospiro_numberings(mol, spiro_atom):
        locant = full_order.index(suffix_carbon) + 1
        substituents = substituents_for_ring(graph, full_order, halogens, excluded, mol=mol)
        suffixed_parent = _suffixed_parent(parent, locant, suffix_word, elide_e)
        key = _spiro_candidate_key(suffixed_parent, substituents, suffix_locant=locant)
        if best_key is None or key < best_key:
            best_key, best_order = key, full_order
    return _with_stereo_prefix(stereo, best_order, best_key[-1])
