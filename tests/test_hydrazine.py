import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_plain_hydrazine_name():
    # H2N-NH2, PubChem CID 713 (auto-generated name itself is "hydrazine",
    # matching the Blue Book PIN exactly -- an unusually complete match
    # for this project).
    assert smiles_to_iupac("NN") == "hydrazine"


def test_methylhydrazine_name():
    # CH3-NH-NH2, PubChem auto-generated name matches the PIN exactly
    # ("methylhydrazine") -- no locant needed, the two nitrogens are
    # interchangeable by symmetry with a single substituent (mirrors
    # `_diazene.py`'s identical no-locant rule for its own single-
    # substituent case).
    assert smiles_to_iupac("CNN") == "methylhydrazine"


def test_ethylhydrazine_name():
    assert smiles_to_iupac("CCNN") == "ethylhydrazine"


def test_1_1_dimethylhydrazine_name():
    # (CH3)2N-NH2, Blue Book P-68.3.1.2.1's own worked example:
    # "1,1-dimethylhydrazine (PIN)". PubChem auto-generated name matches
    # exactly.
    assert smiles_to_iupac("CN(C)N") == "1,1-dimethylhydrazine"


def test_1_2_dimethylhydrazine_name():
    # CH3-NH-NH-CH3, PubChem auto-generated name matches the PIN exactly
    # ("1,2-dimethylhydrazine") -- confirms the locant-assignment logic
    # picks the opposite numbering direction from the 1,1- case above.
    assert smiles_to_iupac("CNNC") == "1,2-dimethylhydrazine"


def test_1_1_2_trimethylhydrazine_name():
    # PubChem auto-generated name matches the PIN exactly
    # ("1,1,2-trimethylhydrazine") -- confirms an uneven 2+1 substituent
    # split is numbered with the doubly-substituted nitrogen as locant 1
    # (lower locant set).
    assert smiles_to_iupac("CN(C)NC") == "1,1,2-trimethylhydrazine"


def test_1_1_2_2_tetramethylhydrazine_name():
    assert smiles_to_iupac("CN(C)N(C)C") == "1,1,2,2-tetramethylhydrazine"


def test_1_1_diethylhydrazine_name():
    assert smiles_to_iupac("CCN(CC)N") == "1,1-diethylhydrazine"


def test_branched_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)NN")


def test_aromatic_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1NN")


def test_halogen_on_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClNN")


def test_chloroethylhydrazine_name():
    # PubChem auto-generated name matches exactly.
    assert smiles_to_iupac("ClCCNN") == "2-chloroethylhydrazine"


def test_bis_chloroethylhydrazine_name():
    # Two identical compound (halogen-bearing) substituents combine with
    # the ordinary 'bis' multiplying prefix, same as elsewhere in this
    # project. PubChem auto-generated name matches exactly.
    assert smiles_to_iupac("ClCCN(N)CCCl") == "1,1-bis(2-chloroethyl)hydrazine"


def test_chloromethylhydrazine_name():
    # The sole-substituent case omits its own hydrazine locant even
    # though the substituent itself is a compound (halogen-bearing) name.
    # PubChem auto-generated name matches exactly.
    assert smiles_to_iupac("ClCNN") == "chloromethylhydrazine"


def test_halogenated_branched_substituent_raises():
    # A real carbon fork (not just a halogen leaf) must still be rejected,
    # even with a halogen elsewhere in the molecule.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC(C)(C)NN")
