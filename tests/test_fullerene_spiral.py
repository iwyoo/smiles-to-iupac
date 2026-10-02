from rdkit import Chem

from smiles_to_iupac._fullerene import (
    _FULLERENE_C60_SMILES,
    _FULLERENE_C70_SMILES,
    _FULLERENE_C76_SMILES,
)
from smiles_to_iupac._fullerene_spiral import (
    _C84_IPR_ISOMERS,
    _canonical_spiral,
    _pentagon_positions,
    _planar_faces,
)


def _pentagons_of(smiles):
    mol = Chem.MolFromSmiles(smiles)
    faces = _planar_faces(mol)
    return _pentagon_positions(_canonical_spiral(faces))


def test_canonical_spiral_matches_published_c60_code():
    # buckminsterfullerene's own well-known canonical spiral (e.g. Fowler &
    # Manolopoulos's atlas): 12 pentagons among 32 faces.
    assert _pentagons_of(_FULLERENE_C60_SMILES) == (1, 7, 9, 11, 13, 15, 18, 20, 22, 24, 26, 32)


def test_canonical_spiral_matches_published_c70_code():
    # cross-checked against the CTCP "Program Fullerene" IPR database.
    assert _pentagons_of(_FULLERENE_C70_SMILES) == (1, 7, 9, 11, 13, 15, 27, 29, 31, 33, 35, 37)


def test_canonical_spiral_matches_published_c76_code():
    # cross-checked against the CTCP "Program Fullerene" IPR database.
    assert _pentagons_of(_FULLERENE_C76_SMILES) == (1, 7, 9, 11, 13, 18, 26, 31, 33, 35, 37, 39)


def test_c84_isomer_table_has_24_distinct_entries():
    # the full, known set of C84 IPR isomers -- no duplicates, no missing
    # index, every spiral exactly 12 pentagons among 44 faces.
    assert len(_C84_IPR_ISOMERS) == 24
    indices = sorted(index for index, _point_group in _C84_IPR_ISOMERS.values())
    assert indices == list(range(1, 25))
    for pentagons in _C84_IPR_ISOMERS:
        assert len(pentagons) == 12
        assert all(1 <= p <= 44 for p in pentagons)
