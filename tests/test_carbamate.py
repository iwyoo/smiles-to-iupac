import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 'R carbamate' (P-65.6.5) for a plain unbranched R -- a common,
        # unambiguous naming pattern (e.g. 'methyl carbamate' is the
        # well-known name of H2NCOOCH3) with no locant/alphabetization
        # choice, so these are trivially checkable without a per-case
        # PubChem lookup.
        ("COC(N)=O", "methyl carbamate"),
        ("CCOC(N)=O", "ethyl carbamate"),
        ("CCCOC(N)=O", "propyl carbamate"),
        ("CCCCOC(N)=O", "butyl carbamate"),
    ],
)
def test_carbamate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_branched_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)OC(N)=O")


def test_methyl_n_methylcarbamate():
    # PubChem CID 81151.
    assert smiles_to_iupac("CNC(=O)OC") == "methyl N-methylcarbamate"


def test_ethyl_n_ethylcarbamate():
    # PubChem CID 12195.
    assert smiles_to_iupac("CCNC(=O)OCC") == "ethyl N-ethylcarbamate"


def test_ethyl_n_methylcarbamate():
    # PubChem CID 7752.
    assert smiles_to_iupac("CNC(=O)OCC") == "ethyl N-methylcarbamate"


def test_n_n_disubstituted_carbamate_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)N(C)C")


def test_branched_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)NC(C)C")


def test_free_carbamic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(N)=O")


def test_cyclic_r_group_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)OC1CCCC1")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOC(N)=O")
