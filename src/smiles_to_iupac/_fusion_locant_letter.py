"""Shared engine for P-25.3.1.3's fusion-locant-letter mechanism: a plain,
unsubstituted benzo ring ortho-fused onto an all-carbon retained-name
aromatic base at one of its periphery bonds not touching a ring-fusion
atom (https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf).

`_anthracene_fusion.py`/`_phenanthrene_fusion.py`/`_pyrene_fusion.py` each
reimplemented this exact matching procedure verbatim, differing only in
the base's own reference SMILES, its atom/ring count, its side-letter
table, and which letters are excluded (already produced via a different,
more senior base component, or already a distinct retained name) --
extracted here once a third copy made the duplication impossible to miss
(2026-09-12 refactor, no behavior change; see each caller's own module
docstring for the base-specific chemistry/citations/PubChem CIDs).

`nonaromatic_ref_atoms` (2026-09-12, added for `_aceanthrylene_fusion.py`):
some bases (aceanthrylene, and per `_peri_fused_aromatic.py`'s own
docstring, acenaphthylene/fluoranthene/acephenanthrylene too) have a
genuine localized -CH=CH- vinylene bridge in their five-membered ring
that RDKit's own aromaticity perception marks non-aromatic even in the
plain, unsubstituted parent -- real chemistry, not a construction bug.
The other five bases so far happen to be fully aromatic throughout, so
this defaults to empty (identical behavior to before this parameter
existed).
"""

from rdkit import Chem


def find_fusion_letter(
    mol,
    ref,
    letter_by_pair,
    excluded_letters,
    num_atoms,
    num_rings,
    nonaromatic_ref_atoms=frozenset(),
    num_extra_atoms=4,
    nonaromatic_extra_atom_count=0,
):
    """`ref`: the base's reference molecule (all-carbon, aromatic aside
    from any `nonaromatic_ref_atoms`). `letter_by_pair`: {frozenset of the
    two `ref`-index fusion atoms -> P-25.3.1.3 side letter}, covering only
    the periphery bonds this base supports plain ortho-fusion at.
    `excluded_letters`: letters that match structurally but must not be
    claimed here (already a distinct retained name, or already produced
    via a different, more senior base component) -- returns None for
    those, same as a non-match, so the caller falls through to whichever
    module does claim that shape. `num_atoms`/`num_rings`: the exact
    fused molecule's expected atom and ring count, ruling out every other
    size before the more expensive substructure search below runs.
    `nonaromatic_ref_atoms`: `ref`-index atoms expected to be non-aromatic
    (see module docstring). `num_extra_atoms`/`nonaromatic_extra_atom_count`:
    the attached ring's own atom count/non-aromatic count beyond the two
    shared fusion atoms (default 4/0, a plain aromatic benzo ring; a
    non-aromatic attachment like cyclopenta passes a smaller/nonzero pair).

    Returns the winning letter (alphabetically lowest among every
    automorphism of `ref` that matches, per P-25.3.1.3's own tie-break),
    or None if the molecule isn't this exact shape."""
    if mol.GetNumAtoms() != num_atoms:
        return None
    if mol.GetRingInfo().NumRings() != num_rings:
        return None
    expected_nonaromatic = len(nonaromatic_ref_atoms) + nonaromatic_extra_atom_count
    if sum(1 for atom in mol.GetAtoms() if not atom.GetIsAromatic()) != expected_nonaromatic:
        return None
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if atom.GetAtomicNum() != 6:
            return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    matches = mol.GetSubstructMatches(ref, useChirality=False, uniquify=False)
    best = None
    for match in matches:
        core = set(match)
        extra = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
        if len(extra) != num_extra_atoms:
            continue
        if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in extra):
            continue
        if sum(1 for a in extra if not mol.GetAtomWithIdx(a).GetIsAromatic()) != nonaromatic_extra_atom_count:
            continue
        target_to_ref_idx = {match[i]: i for i in range(len(match))}
        if any(
            mol.GetAtomWithIdx(target).GetIsAromatic() == (ref_idx in nonaromatic_ref_atoms)
            for target, ref_idx in target_to_ref_idx.items()
        ):
            continue
        fusion_ref_idxs = set()
        bad = False
        for a in extra:
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                if n.GetIdx() in target_to_ref_idx:
                    fusion_ref_idxs.add(target_to_ref_idx[n.GetIdx()])
                elif n.GetIdx() not in extra:
                    bad = True
        if bad or len(fusion_ref_idxs) != 2:
            continue
        letter = letter_by_pair.get(frozenset(fusion_ref_idxs))
        if letter is not None and (best is None or letter < best):
            best = letter
    if best in excluded_letters:
        return None
    return best
