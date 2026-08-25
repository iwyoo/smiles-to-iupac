import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_thieno_2_3_b_thiophene():
    # C6H4S2, cross-checked against PubChem CID 136062's IUPACName field
    # and ConnectivitySMILES ("C1=CSC2=C1C=CS2").
    assert smiles_to_iupac("C1=CSC2=C1C=CS2") == "thieno[2,3-b]thiophene"


def test_thieno_3_2_b_thiophene():
    # C6H4S2, cross-checked against PubChem CID 136063's IUPACName field
    # and ConnectivitySMILES ("C1=CSC2=C1SC=C2") -- the other attachment
    # direction of the same fusion bond as the [2,3-b] isomer above.
    assert smiles_to_iupac("C1=CSC2=C1SC=C2") == "thieno[3,2-b]thiophene"


def test_furo_2_3_b_furan():
    # C6H4O2, cross-checked against PubChem CID 17945614's IUPACName field
    # and ConnectivitySMILES ("C1=COC2=C1C=CO2").
    assert smiles_to_iupac("C1=COC2=C1C=CO2") == "furo[2,3-b]furan"


def test_furo_3_2_b_furan():
    # C6H4O2, cross-checked against PubChem (SMILES lookup for
    # "C1=COC2=C1OC=C2") CID 22416599's IUPACName field.
    assert smiles_to_iupac("C1=COC2=C1OC=C2") == "furo[3,2-b]furan"


def test_c_lettered_fusion_raises():
    # thieno[2,3-c]thiophene (PubChem CID 520171's own IUPACName for this
    # ConnectivitySMILES): the fusion bond doesn't touch either ring's own
    # heteroatom, a structurally distinct isomer this module's narrow
    # 'b'-lettered scope must not silently misname.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CSC2=CSC=C21")


def test_thieno_2_3_b_furan():
    # C6H4OS, cross-checked against PubChem CID 21868935's IUPACName field
    # and ConnectivitySMILES ("C1=COC2=C1C=CS2"). A mixed thiophene+furan
    # pair: P-25.3.2.4(a)'s heteroatom seniority (O > S, confirmed against
    # the Blue Book's own "furan is senior to dithiepine" worked example)
    # always makes furan the base component.
    assert smiles_to_iupac("C1=CSC2=C1C=CO2") == "thieno[2,3-b]furan"


def test_thieno_3_2_b_furan():
    # C6H4OS, cross-checked against PubChem CID 18436918's IUPACName field
    # and ConnectivitySMILES ("C1=COC2=C1SC=C2").
    assert smiles_to_iupac("C1=COC2=C1SC=C2") == "thieno[3,2-b]furan"


def test_substituted_thieno_thiophene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1csc2ccsc12")


def test_selenophene_thiophene_raises():
    # Se is out of scope -- only the O-vs-S slice of P-25.3.2.4(a)'s
    # seniority order is implemented (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=C[Se]C2=C1C=CS2")
