"""Naming of the 'dinor' skeletal-modification prefix (P-101.3.1.1,
https://iupac.qmul.ac.uk/BlueBook/PDF/P10.pdf) on the bare steroid parent
hydrides already recognized by `_steroid_parent_hydrides.py`: two
non-fusion ring atoms or angular methyls -- each individually a valid
single-atom `nor` target (see `_nor_steroid.py`) -- removed together from
an otherwise-unmodified parent skeleton, with their attached hydrogens.
Named `<locant>,<locant>-dinor-<parent>` (the numerical multiplicative
extension P-101.3.1.1 gives to `nor` for two or more removed atoms).

Removing two individually-eligible atoms is not always well-formed on its
own: if both removed atoms are the only two ring neighbors of some third
atom, that third atom is stranded as its own disconnected fragment (e.g.
removing androstane's locants 1 and 3 together isolates locant 2 as a
free methane molecule) -- rejected here via a post-removal connectivity
check (`Chem.GetMolFrags`), not by reasoning about locant adjacency
directly.

As with `_homo_steroid.py`, several distinct locant pairs can produce the
same product (confirmed empirically: 21 of the 463 candidate pairs across
all 7 raw parents collide with another pair on the same parent) -- lowest
locants win the citation, same tie-break rule. One pair (removing both
angular methyls, 18 and 19) reproduces gonane exactly and is excluded via
the same `_CANONICAL_TO_NAME` collision check `_nor_steroid.py` uses.

Reuses `_nor_steroid.py`'s existing steroid-numbering infrastructure
(`_locant_map`, `_RAW_PARENTS`) rather than duplicating it.
"""

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import RWMol

from ._nor_steroid import _RAW_PARENTS, _locant_map
from ._steroid_parent_hydrides import _CANONICAL_TO_NAME


def _eligible_locants(mol, locants):
    return [loc for loc, idx in locants.items() if mol.GetAtomWithIdx(idx).GetDegree() <= 2]


def _remove_two(mol, a_idx, b_idx):
    rw = RWMol(mol)
    for idx in sorted((a_idx, b_idx), reverse=True):
        rw.RemoveAtom(idx)
    candidate = rw.GetMol()
    try:
        Chem.SanitizeMol(candidate)
    except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    if len(Chem.GetMolFrags(candidate)) != 1:
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
    assert existing_name == name, f"{label}-dinor-{name} collides with an unrelated dinor entry"
    if sort_key < existing_key:
        lookup[candidate] = (sort_key, label, name)


def _build_dinor_lookup():
    raw = {}
    for name, smiles in _RAW_PARENTS.items():
        mol = Chem.MolFromSmiles(smiles)
        locants = _locant_map(name, mol)
        eligible = sorted(_eligible_locants(mol, locants))
        for a, b in combinations(eligible, 2):
            candidate = _remove_two(mol, locants[a], locants[b])
            _add(raw, candidate, (a, b), f"{a},{b}", name)
    return {smiles: (label, name) for smiles, (_, label, name) in raw.items()}


_DINOR_LOOKUP = _build_dinor_lookup()


def has_dinor_steroid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _DINOR_LOOKUP


def name_dinor_steroid(mol) -> str:
    label, parent = _DINOR_LOOKUP[Chem.MolToSmiles(mol)]
    return f"{label}-dinor-{parent}"
