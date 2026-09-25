"""Naming of the 'cyclo' skeletal-modification prefix (P-101.3.3,
https://iupac.qmul.ac.uk/BlueBook/PDF/P10.pdf) on the bare steroid parent
hydrides already recognized by `_steroid_parent_hydrides.py`: one extra
ring formed by a direct bond between two non-adjacent ring atoms (locants
1-17 only -- not through the angular methyls 18/19, a different,
unresearched shape) of an otherwise-unmodified parent skeleton. Named
`<locant>,<locant>-cyclo-<parent>` (P-101.3.3's own worked example,
`3alpha,5-cyclo-5alpha-pregnane`).

Unlike `_nor_steroid.py`/`_dinor_steroid.py`, a new bond may legitimately
terminate at a ring-fusion atom (the Blue Book's own worked example bonds
to C5, a fusion carbon) -- no degree-based eligibility filter is applied
here. Instead, every non-adjacent locant pair is tried and simply
rejected if bond formation pushes either atom's valence too high for
carbon (a fusion atom already carrying an angular methyl, e.g. C10 or
C13, is already tetravalent and can't take a new bond) -- caught via
`Chem.SanitizeMol`'s own valence check rather than precomputed degree
rules, since which pairs are actually inaccessible depends on the whole
local structure, not a single atom's degree in isolation.

No locant-collision or degeneracy issues were found empirically across
all viable non-adjacent locant-1-17 pairs on all 7 raw parents (unlike
`_homo_steroid.py`/`_dinor_steroid.py`, no lowest-locant tie-breaking was
needed here -- every successfully-formed bond gives a distinct product).

Reuses `_nor_steroid.py`'s existing steroid-numbering infrastructure
(`_locant_map`, `_RAW_PARENTS`) rather than duplicating it.
"""

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import BondType, RWMol

from ._nor_steroid import _RAW_PARENTS, _locant_map
from ._steroid_parent_hydrides import _CANONICAL_TO_NAME


def _form_bond(mol, a_idx, b_idx):
    if mol.GetBondBetweenAtoms(a_idx, b_idx) is not None:
        return None
    rw = RWMol(mol)
    rw.AddBond(a_idx, b_idx, BondType.SINGLE)
    candidate = rw.GetMol()
    try:
        Chem.SanitizeMol(candidate)
    except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    return Chem.MolToSmiles(candidate)


def _build_cyclo_lookup():
    lookup = {}
    for name, smiles in _RAW_PARENTS.items():
        mol = Chem.MolFromSmiles(smiles)
        locants = _locant_map(name, mol)
        ring_locants = sorted(loc for loc in locants if loc <= 17)
        for a, b in combinations(ring_locants, 2):
            candidate = _form_bond(mol, locants[a], locants[b])
            if candidate is None or candidate in _CANONICAL_TO_NAME:
                continue
            label = f"{a},{b}"
            if candidate in lookup:
                assert lookup[candidate] == (label, name), (
                    f"{label}-cyclo-{name} collides with an unrelated cyclo entry"
                )
                continue
            lookup[candidate] = (label, name)
    return lookup


_CYCLO_LOOKUP = _build_cyclo_lookup()


def has_cyclo_steroid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CYCLO_LOOKUP


def name_cyclo_steroid(mol) -> str:
    label, parent = _CYCLO_LOOKUP[Chem.MolToSmiles(mol)]
    return f"{label}-cyclo-{parent}"
