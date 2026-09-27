"""Shared strip-then-match idiom for recognizing a retained parent-hydride
skeleton underneath one or more exocyclic substituents: remove a set of
substituent atoms, then let the caller canonicalize the bare ring system
and match it against its own name lookup.

Generalizes `_steroid_parent_hydrides.py`'s original single-substituent
`steroid_suffix_name` (which stripped exactly one suffix heteroatom) to an
arbitrary set of exocyclic substituent atoms removed simultaneously --
needed for Appendix 3 alkaloid retained names (#1065), where a real
skeleton like morphine carries several substituents (two O-substituents,
an N-methyl, a ring-closing ether bridge) on one recognized parent at
once, not just one.
"""

from rdkit import Chem


def strip_substituents(mol, atoms_to_remove):
    """Remove every atom index in `atoms_to_remove` from `mol`, recomputing
    each surviving neighbor's implicit H count from scratch (so a single-
    bonded substituent, e.g. an -OH oxygen or a ring-closing bridge atom,
    leaves correct valence behind on its ring-carbon neighbor -- a bracket
    atom like `[C@H]` normally has its H count fixed as *explicit* by the
    SMILES parser, so sanitizing alone leaves it under-valent unless this
    is reset first; a double-bonded substituent's carbon has no such fixed
    explicit count, so this is a no-op for it, matching prior behavior).

    Returns `(stripped_mol, old_to_new)`, where `old_to_new` maps every
    surviving atom's original index in `mol` to its index in
    `stripped_mol` -- or `(None, None)` if the result doesn't sanitize
    (e.g. stripping breaks aromaticity)."""
    to_remove = set(atoms_to_remove)
    neighbors_to_fix = set()
    for idx in to_remove:
        for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
            if neighbor.GetIdx() not in to_remove:
                neighbors_to_fix.add(neighbor.GetIdx())

    remaining_sorted = sorted(set(range(mol.GetNumAtoms())) - to_remove)
    old_to_new = {old: new for new, old in enumerate(remaining_sorted)}

    rw = Chem.RWMol(mol)
    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)
    stripped = rw.GetMol()

    for old_neighbor in neighbors_to_fix:
        atom = stripped.GetAtomWithIdx(old_to_new[old_neighbor])
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)

    try:
        Chem.SanitizeMol(stripped)
    except Chem.rdchem.KekulizeException:
        return None, None
    return stripped, old_to_new
