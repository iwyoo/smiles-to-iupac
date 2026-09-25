"""Naming of the 'seco' skeletal-modification prefix (P-101.3.4.1,
https://iupac.qmul.ac.uk/BlueBook/PDF/P10.pdf) on the bare steroid parent
hydrides already recognized by `_steroid_parent_hydrides.py`: one ring
bond cleaved, with one extra hydrogen added at each of the two former
bond termini, reducing the ring count by one while leaving the atom count
and the original numbering unchanged. Named `<locant>,<locant>-seco-
<parent>`, e.g. `2,3-secohopane`/`3,4-secocuran` (the Blue Book's own
worked examples) -- unlike `nor` (P-101.3.1), there is no atomic-
connector/bond-connector distinction here and no renumbering: P-101.3.4.1
states plainly "the original numbering is retained".

No locant-collision or degeneracy issues were found empirically across
all 20 ring bonds on all 7 raw parents (unlike `_homo_steroid.py`, which
needs lowest-locant tie-breaking for bare gonane's symmetric rings --
cleaving a bond leaves the two termini as distinguishable degree-1
branch points elsewhere on the same fused ring system, so different
cleavage points always give different products even without any
methyl substitution to break symmetry). Cleaving a ring-fusion bond
(shared between two rings, e.g. locants 9,10) and cleaving a plain
non-fusion ring bond are handled identically here -- both simply remove
one ring bond and reduce the whole-molecule ring count by one; no special
casing is needed despite the general von Baeyer engine failing
differently on the two shapes when it encounters them unclaimed (a
9,10-cleaved structure hits an outright "polycyclic ring systems are not
supported yet" rejection there, a plain-ring-bond cleavage instead gets
silently misnamed as an unrelated tricyclic system) -- both are simply
absent from this lookup below before this module runs.

Reuses `_nor_steroid.py`'s existing steroid-numbering infrastructure
(`_locant_map`, `_RAW_PARENTS`) rather than duplicating it.
"""

from rdkit import Chem
from rdkit.Chem import RWMol

from ._nor_steroid import _RAW_PARENTS, _locant_map
from ._steroid_parent_hydrides import _CANONICAL_TO_NAME

_RING_BONDS = (
    (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
    (10, 1), (5, 10), (9, 11), (11, 12), (12, 13), (13, 14), (14, 8),
    (13, 17), (17, 16), (16, 15), (15, 14),
)


def _cleave_bond(mol, a_idx, b_idx):
    rw = RWMol(mol)
    rw.RemoveBond(a_idx, b_idx)
    candidate = rw.GetMol()
    try:
        Chem.SanitizeMol(candidate)
    except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    return Chem.MolToSmiles(candidate)


def _build_seco_lookup():
    lookup = {}
    for name, smiles in _RAW_PARENTS.items():
        mol = Chem.MolFromSmiles(smiles)
        locants = _locant_map(name, mol)
        for a, b in _RING_BONDS:
            lo, hi = sorted((a, b))
            candidate = _cleave_bond(mol, locants[a], locants[b])
            if candidate is None or candidate in _CANONICAL_TO_NAME:
                continue
            label = f"{lo},{hi}"
            if candidate in lookup:
                assert lookup[candidate] == (label, name), (
                    f"{label}-seco-{name} collides with an unrelated seco entry"
                )
                continue
            lookup[candidate] = (label, name)
    return lookup


_SECO_LOOKUP = _build_seco_lookup()


def has_seco_steroid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _SECO_LOOKUP


def name_seco_steroid(mol) -> str:
    label, parent = _SECO_LOOKUP[Chem.MolToSmiles(mol)]
    return f"{label}-seco-{parent}"
