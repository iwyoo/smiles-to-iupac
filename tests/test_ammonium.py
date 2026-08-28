import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_unsubstituted_ammonium_name():
    # NH4+, PubChem CID 223 -- the Blue Book PIN and PubChem's own
    # auto-generated name agree exactly on "azanium" (P-73.1.1.2).
    assert smiles_to_iupac("[NH4+]") == "azanium"


def test_methanaminium_name():
    # CH3-NH3+, structure-verified via PubChem CID 644041
    # (ConnectivitySMILES "[NH3+]C"; PubChem's own auto-generated name
    # "methylazanium" isn't PIN format). Name derived from the Blue Book's
    # own P-73.1.1.2 rule ("final 'e' -> 'ium'") applied to this project's
    # existing, independently verified "methanamine" (P-33.1).
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_propan_1_aminium_name():
    # CH3-CH2-CH2-NH3+, structure-verified via PubChem CID 3483736
    # ("propylazanium"). Derived from this project's existing
    # "propan-1-amine" the same way as above.
    assert smiles_to_iupac("CCC[NH3+]") == "propan-1-aminium"


def test_propan_2_aminium_name():
    # (CH3)2CH-NH3+, structure-verified via PubChem CID 3364502
    # ("propan-2-ylazanium"). Derived from this project's existing
    # "propan-2-amine" the same way as above.
    assert smiles_to_iupac("CC(C)[NH3+]") == "propan-2-aminium"


def test_chloromethanaminium_name():
    # Cl-CH2-NH3+, structure-verified via PubChem CID 20462735
    # ("chloromethylazanium"). Confirms halogen-prefix coexistence,
    # derived from this project's existing "chloromethanamine" the same
    # way as above.
    assert smiles_to_iupac("C(Cl)[NH3+]") == "chloromethanaminium"


def test_secondary_ammonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+]C")


def test_quaternary_ammonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[N+](C)(C)C")


def test_doubly_charged_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3++]")
