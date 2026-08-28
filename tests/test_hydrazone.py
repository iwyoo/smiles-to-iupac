import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methylidenehydrazine_name():
    # Formaldehyde hydrazone, C=N-NH2. PubChem PUG REST CID 81125,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("C=NN") == "methylidenehydrazine"


def test_ethylidenehydrazine_name():
    # Acetaldehyde hydrazone. PubChem PUG REST CID 53651639, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC=NN") == "ethylidenehydrazine"


def test_propylidenehydrazine_name():
    # Propionaldehyde hydrazone -- the double-bond carbon lands at
    # position 1 of its 3-carbon chain (aldehyde-shaped attachment), so no
    # locant is cited, matching the free-valence rule `name_branch` already
    # applies for an ordinary substituent. PubChem PUG REST CID 53726800,
    # auto-generated name matches exactly.
    assert smiles_to_iupac("CCC=NN") == "propylidenehydrazine"


def test_propan_2_ylidenehydrazine_name():
    # Acetone hydrazone -- the double-bond carbon has two carbon neighbors
    # (ketone-shaped attachment), so it lands at position 2 of the
    # 3-carbon chain and the systematic locanted form is used instead.
    # PubChem PUG REST CID 78937, auto-generated name matches exactly.
    assert smiles_to_iupac("CC(C)=NN") == "propan-2-ylidenehydrazine"


def test_butan_2_ylidenehydrazine_name():
    # 2-Butanone hydrazone. PubChem PUG REST CID 71356052, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CCC(C)=NN") == "butan-2-ylidenehydrazine"


def test_butan_2_ylidenehydrazine_name_alternate_smiles():
    # Same compound as above, written the other way round -- same CID
    # 71356052 confirms both SMILES describe the same structure/name.
    assert smiles_to_iupac("CC(=NN)CC") == "butan-2-ylidenehydrazine"


def test_2_methylpropylidenehydrazine_name():
    # Isobutyraldehyde hydrazone, a branched aldehyde hydrazone -- the
    # double-bond carbon is still position-1 (aldehyde-shaped), so
    # `name_branch`'s own old-style branched-substituent name
    # ("2-methylpropyl") is reused with "idene" appended. PubChem PUG REST
    # CID 71356051, auto-generated name matches exactly.
    assert smiles_to_iupac("CC(C)C=NN") == "2-methylpropylidenehydrazine"


def test_2_chloroethylidenehydrazine_name():
    # A halogen substituent on the carbon chain coexists freely (P-35.2.1).
    # PubChem PUG REST CID 163551743, auto-generated name matches exactly.
    assert smiles_to_iupac("ClCC=NN") == "2-chloroethylidenehydrazine"


def test_3_chloropropylidenehydrazine_name():
    # PubChem PUG REST CID 174934626, auto-generated name matches exactly.
    assert smiles_to_iupac("ClCCC=NN") == "3-chloropropylidenehydrazine"


def test_n_substituted_hydrazone_raises():
    # =N-NH-R (the terminal nitrogen bearing an alkyl group) follows a
    # completely different naming scheme (PubChem PUG REST CID 86520,
    # "N-(ethylideneamino)methanamine") that this module doesn't attempt.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=NNC")


def test_azine_raises():
    # R2C=N-N=CR2 (both hydrazine nitrogens double-bonded to their own
    # carbon) is a different naming scheme (PubChem PUG REST CID 79085,
    # "N-(propan-2-ylideneamino)propan-2-imine") -- more than one C=N
    # double bond, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)=NN=C(C)C")


def test_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1=NN")


def test_aromatic_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C=NN")
