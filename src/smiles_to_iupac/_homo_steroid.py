"""Naming of the 'homo' skeletal-modification prefix (P-101.3.2,
https://iupac.qmul.ac.uk/BlueBook/PDF/P10.pdf) on the bare steroid parent
hydrides already recognized by `_steroid_parent_hydrides.py`: one extra
methylene (-CH2-) inserted into an otherwise-unmodified parent skeleton.

Two locant forms, per P-101.3.2.2.1/.2.2.2:

- Insertion into a ring bond between two adjacent numbered ring atoms
  (a "bond connector") is cited with both termini, the higher one in
  parentheses, plus a letter: `13(17)a-homo-<parent>` (P-101.3.2.2.2's own
  worked example, expanding ring D from five to six atoms).
- Insertion between a ring atom and one of the two angular methyls
  (C10-C19 or C13-C18), extending the methyl to an ethyl bridge (a
  "terminal segment"), is cited with just the methyl's own locant plus a
  letter: `19a-homo-<parent>` (P-101.3.2.2.1's own worked example).

Only single-methylene insertion is handled (the `dihomo` multiplicative
case and the indicated-hydrogen edge case for mancude ring systems, which
doesn't apply to these fully-saturated steroid skeletons, are out of
scope). Reuses `_nor_steroid.py`'s per-skeleton steroid-numbering
infrastructure (`_locant_map`, `_RAW_PARENTS`) rather than duplicating it.

As with `_nor_steroid.py`, an insertion that happens to reproduce an
already-listed retained-name parent is excluded from the lookup below
(none found empirically across all 20 ring bonds and both angular methyls
on all 7 raw parents) -- `core.py` dispatches to the exact-parent and
`nor` checks before this module regardless, so this is defensive rather
than load-bearing.

Several distinct ring bonds can produce the exact same product structure
-- confirmed on gonane, whose plain unsubstituted rings have no angular
methyls to break the symmetry, so inserting a methylene anywhere around
an unsubstituted ring segment (e.g. any of ring A's five non-fusion-free
bonds) gives the same molecule regardless of which bond was split. When
this happens, standard lowest-locant tie-breaking picks the citation: of
all bond pairs producing one product, the pair compared lexicographically
smallest wins. This is also why the Blue Book's own `13(17)a-homo-`
worked example on ring D is the *lowest* of that ring's four
symmetry-equivalent bonds ((13,17), (14,15), (15,16), (16,17)) in a
methyl-free ring D, not an arbitrary pick.
"""

from rdkit import Chem
from rdkit.Chem import Atom, BondType, RWMol

from ._nor_steroid import _RAW_PARENTS, _locant_map
from ._steroid_parent_hydrides import _CANONICAL_TO_NAME

_RING_BONDS = (
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
    (10, 1), (5, 10), (9, 11), (11, 12), (12, 13), (13, 14), (14, 8),
    (13, 17), (17, 16), (16, 15), (15, 14),
)
_ANGULAR_METHYLS = ((18, 13), (19, 10))


def _insert_methylene(mol, a_idx, b_idx):
    rw = RWMol(mol)
    new_idx = rw.AddAtom(Atom(6))
    rw.RemoveBond(a_idx, b_idx)
    rw.AddBond(a_idx, new_idx, BondType.SINGLE)
    rw.AddBond(new_idx, b_idx, BondType.SINGLE)
    candidate = rw.GetMol()
    try:
        Chem.SanitizeMol(candidate)
    except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    return Chem.MolToSmiles(candidate)


def _add(lookup, candidate, sort_key, label, name):
    if candidate is None or candidate in _CANONICAL_TO_NAME:
        return
    existing = lookup.get(candidate)
    if existing is None:
        lookup[candidate] = (sort_key, label, name)
        return
    existing_key, _, existing_name = existing
    assert existing_name == name, f"{label}-homo-{name} collides with an unrelated homo entry"
    if sort_key < existing_key:
        lookup[candidate] = (sort_key, label, name)


def _build_homo_lookup():
    raw = {}
    for name, smiles in _RAW_PARENTS.items():
        mol = Chem.MolFromSmiles(smiles)
        locants = _locant_map(name, mol)
        for a, b in _RING_BONDS:
            lo, hi = sorted((a, b))
            candidate = _insert_methylene(mol, locants[a], locants[b])
            _add(raw, candidate, (lo, hi), f"{lo}({hi})a", name)
        for methyl_locant, ring_locant in _ANGULAR_METHYLS:
            if methyl_locant not in locants:
                continue
            candidate = _insert_methylene(mol, locants[ring_locant], locants[methyl_locant])
            _add(raw, candidate, (methyl_locant,), f"{methyl_locant}a", name)
    return {smiles: (label, name) for smiles, (_, label, name) in raw.items()}


_HOMO_LOOKUP = _build_homo_lookup()


def has_homo_steroid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _HOMO_LOOKUP


def name_homo_steroid(mol) -> str:
    label, parent = _HOMO_LOOKUP[Chem.MolToSmiles(mol)]
    return f"{label}-homo-{parent}"
