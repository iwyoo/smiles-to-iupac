"""Auto-derive a fused aromatic base's P-25.3.3.1.1 peripheral numbering
and P-25.3.1.3 `letter_by_pair` mapping from its reference molecule plus
one or more real "anchor" structures -- replacing the manual process
every `_X_fusion.py` module's own docstring describes doing by hand so
far (walk the SSSR periphery, then match real monosubstituted PubChem
structures by substructure to fix the numbering direction, e.g.
`_aceanthrylene_fusion.py` used four monomethylaceanthrylene anchors,
`_fluoranthene_fusion.py` used five). Two anchor kinds are supported,
combinable in the same call: a substituent-position anchor (`anchors`,
"this real molecule is `ref` plus one substituent at known locant N")
and a known-compound anchor (`compound_anchors`, "this real, already-
named compound is `ref` plus one benzo ring fused at known letter L").

Reuses `_aromatic.py`'s existing graph primitives unchanged
(`find_aromatic_fused_core`, `_periphery_cycle`, `_atom_ring_membership`,
`_walk_locants`) -- this module only adds the anchor-matching and
letter-assignment layer on top.

**Anchors alone cannot always pin down a *unique* numbering, no matter
how many are given -- this is a real mathematical limit, not an
engineering gap.** If the reference molecule has a nontrivial
automorphism (e.g. phenanthrene's own mirror symmetry, or the D2h
long-axis-mirror-plus-ring-swap symmetry every plain straight acene from
anthracene through heptacene has), then for *every* group element other
than the identity, applying it to any anchor-consistent numbering
produces a second, equally anchor-consistent numbering -- because the
automorphism, by definition, maps the whole molecule onto itself, locant
labels and all. No amount of "atom X is locant N" anchor data can ever
distinguish between these; only an additional rule external to the
periphery-walk-plus-anchors model (e.g. a fixed convention for breaking
ties among otherwise-equivalent numberings) could. This module does not
implement such a rule -- `derive_letter_by_pair` raises
`UnsupportedStructure` naming how many candidates remain whenever more
than one survives, rather than silently guessing, and that count will
never drop below the reference molecule's own automorphism group size
regardless of how many anchors are supplied. In practice this means: a
reference molecule with trivial symmetry (like aceanthrylene) converges
to one answer from just enough anchors to break its accidental periphery-
walk ambiguities; a symmetric one (like phenanthrene, or any of the
straight acenes) never converges past its symmetry-group size, and the
historically-hardcoded numbering (whichever member of that symmetric
family the Blue Book happened to depict) will always be *among* the
survivors, verified below, but is not picked out uniquely by this module
alone.
"""

from ._common import UnsupportedStructure
from ._aromatic import (
    _atom_ring_membership,
    _periphery_cycle,
    _walk_locants,
    find_aromatic_fused_core,
)


def _anchor_ref_locant_candidates(ref, anchor_mol):
    """Ref atom indices that could be "the substituted position" in
    `anchor_mol`, one per valid embedding of `ref` as a substructure of
    `anchor_mol`'s own ring-plus-substituent skeleton. More than one
    entry means `ref`'s own symmetry doesn't yet distinguish which
    physical ref atom this is, from this anchor alone."""
    candidates = set()
    for match in anchor_mol.GetSubstructMatches(ref, useChirality=False, uniquify=False):
        core = set(match)
        attach_positions = set()
        for ref_idx, anchor_idx in enumerate(match):
            has_extra_neighbor = any(
                n.GetIdx() not in core for n in anchor_mol.GetAtomWithIdx(anchor_idx).GetNeighbors()
            )
            if has_extra_neighbor:
                attach_positions.add(ref_idx)
        if len(attach_positions) == 1:
            candidates.add(next(iter(attach_positions)))
    return candidates


def _anchor_ref_bond_candidates(ref, compound_mol):
    """Ref bonds (as `frozenset({a, b})` of ref atom indices) that could
    be the fusion bond where an already-named real compound has an extra
    ring ortho-fused onto `ref`, one per valid embedding of `ref` as a
    substructure of `compound_mol`. This is the "known-compound" anchor
    counterpart to `_anchor_ref_locant_candidates`'s "position-substituent"
    anchor -- same substructure-match skeleton, but it looks for a bond
    whose *both* ref atoms have an extra neighbor outside the match (the
    two fusion-bond atoms shared with the new ring), not a single
    substituted atom."""
    candidates = set()
    for match in compound_mol.GetSubstructMatches(ref, useChirality=False, uniquify=False):
        core = set(match)
        fused_ref_atoms = set()
        for ref_idx, compound_idx in enumerate(match):
            has_extra_neighbor = any(
                n.GetIdx() not in core
                for n in compound_mol.GetAtomWithIdx(compound_idx).GetNeighbors()
            )
            if has_extra_neighbor:
                fused_ref_atoms.add(ref_idx)
        for bond in ref.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in fused_ref_atoms and b in fused_ref_atoms:
                candidates.add(frozenset({a, b}))
    return candidates


def _numbering_candidates(cycle, fusion_atoms):
    """Every (start_index, step, locants) triple `_walk_locants` can
    produce from this periphery cycle, starting only at a nonfusion atom
    (P-25.3.3.1.1 numbering never starts on a fusion atom)."""
    n = len(cycle)
    for start in range(n):
        if cycle[start] in fusion_atoms:
            continue
        for step in (1, -1):
            yield start, step, _walk_locants(cycle, fusion_atoms, start, step)


def _assign_letters(cycle, start, step, fusion_atoms):
    """P-25.3.1.3's continuous peripheral lettering: 'a' for the bond
    leaving position 1, then 'b', 'c', ... all the way around -- but only
    bonds whose two atoms are both nonfusion make it into the returned
    dict (the only bonds a plain benzo ring can ortho-fuse at)."""
    n = len(cycle)
    letter_by_pair = {}
    for k in range(n):
        a = cycle[(start + step * k) % n]
        b = cycle[(start + step * (k + 1)) % n]
        if a not in fusion_atoms and b not in fusion_atoms:
            letter_by_pair[frozenset({a, b})] = chr(ord("a") + k)
    return letter_by_pair


def letter_by_pair_candidates(ref, anchors, compound_anchors=()):
    """Every `letter_by_pair` dict consistent with `ref`'s periphery graph
    and every anchor given. `anchors` is a list of (anchor_mol,
    claimed_locant) pairs, each a real molecule that is `ref` plus exactly
    one plain substituent at a known IUPAC locant. `compound_anchors` is
    the same idea one level up: a list of (compound_mol, claimed_letter)
    pairs, each an already-named real compound that is `ref` plus exactly
    one plain benzo ring ortho-fused at a known P-25.3.1.3 letter --
    useful when the compound's *name* (hence its fusion letter) is known
    but no single substituent-locant anchor is available or needed.
    Ordinarily a list of 1 (a fully anchor-determined numbering) or more
    (a residual symmetry -- see module docstring); never empty unless the
    anchors are mutually inconsistent, in which case this raises
    `UnsupportedStructure` directly rather than returning an empty list,
    since that's always a bad-input error, not a valid "zero candidates"
    answer.

    Raises `UnsupportedStructure` if `ref` isn't a supported aromatic
    core, if an anchor doesn't cleanly match `ref` plus one substituent
    (or one fused ring), or if the anchors are mutually inconsistent (no
    valid numbering satisfies all of them at once)."""
    core = find_aromatic_fused_core(ref)
    if core is None:
        raise UnsupportedStructure(
            "not a plain all-carbon aromatic mancude ring system (see P-25.3.1.3)"
        )
    atom_rings, ring_atom_sets, fusion_bond_idxs = core
    cycle = _periphery_cycle(ref, ring_atom_sets, fusion_bond_idxs)
    fusion_atoms = {atom for atom, rings in _atom_ring_membership(atom_rings).items() if len(rings) >= 2}

    anchor_ref_candidates = []
    for anchor_mol, claimed_locant in anchors:
        ref_candidates = _anchor_ref_locant_candidates(ref, anchor_mol)
        if not ref_candidates:
            raise UnsupportedStructure(
                "an anchor structure does not match the reference molecule plus "
                "exactly one plain substituent (see P-25.3.1.3)"
            )
        anchor_ref_candidates.append((ref_candidates, claimed_locant))

    compound_anchor_bond_candidates = []
    for compound_mol, claimed_letter in compound_anchors:
        ref_bond_candidates = _anchor_ref_bond_candidates(ref, compound_mol)
        if not ref_bond_candidates:
            raise UnsupportedStructure(
                "a compound anchor does not match the reference molecule plus "
                "exactly one plain fused ring (see P-25.3.1.3)"
            )
        compound_anchor_bond_candidates.append((ref_bond_candidates, claimed_letter))

    unique_letterings = []
    for start, step, locants in _numbering_candidates(cycle, fusion_atoms):
        if not all(
            any(locants.get(c) == claimed_locant for c in ref_candidates)
            for ref_candidates, claimed_locant in anchor_ref_candidates
        ):
            continue
        letters = _assign_letters(cycle, start, step, fusion_atoms)
        if not all(
            any(letters.get(c) == claimed_letter for c in ref_bond_candidates)
            for ref_bond_candidates, claimed_letter in compound_anchor_bond_candidates
        ):
            continue
        if letters not in unique_letterings:
            unique_letterings.append(letters)

    if not unique_letterings:
        raise UnsupportedStructure(
            "the given anchors are not consistent with any valid numbering of "
            "this reference molecule (see P-25.3.1.3)"
        )
    return unique_letterings


def derive_letter_by_pair(ref, anchors, compound_anchors=()):
    """Like `letter_by_pair_candidates`, but requires the anchors to
    narrow the result down to exactly one `letter_by_pair` dict, raising
    `UnsupportedStructure` (naming how many candidates remain) otherwise
    -- use this once enough anchors are in hand to fully determine a
    base's numbering; use `letter_by_pair_candidates` directly to inspect
    the candidate set for a symmetric base where that will never happen
    (see module docstring)."""
    candidates = letter_by_pair_candidates(ref, anchors, compound_anchors)
    if len(candidates) > 1:
        raise UnsupportedStructure(
            f"the given anchors leave {len(candidates)} numbering "
            "candidates undecided (the reference molecule's own symmetry isn't "
            "fully broken yet) - supply another, more discriminating anchor"
        )
    return candidates[0]
