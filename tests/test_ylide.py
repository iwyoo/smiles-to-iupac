import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_trimethylammonium_methylide():
    # Real PubChem structure, CID 12436722.
    assert smiles_to_iupac("[CH2-][N+](C)(C)C") == "(N,N-dimethylmethanaminiumyl)methanide"


def test_mixed_substituent_ammonium_ylide():
    # The longest substituent (ethyl) becomes the aminium's own parent
    # chain, the other two methyls cited as an 'N,N-dimethyl' prefix.
    assert smiles_to_iupac("[CH2-][N+](C)(C)CC") == "(N,N-dimethylethanaminiumyl)methanide"


def test_branched_anion_carbon_raises():
    # A carbanion bonded to more than just the ylide nitrogen (a
    # branched/chain carbanion parent, e.g. the worked example's own
    # 'propan-2-ide' case) is a narrower deferred slice -- see module
    # docstring.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[CH-][N+](C)(C)C")


def test_tertiary_ammonium_ylide_raises():
    # An ammonium nitrogen with degree other than 4 (here 3: the ylide
    # bond plus two methyls) isn't yet supported.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[CH2-][N+](C)C")
