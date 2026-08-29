import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Blue Book P-72.2.2.2.3's own worked example: CH3-NH(-) ->
        # 'methanaminide (PIN)'. This project's `_amine.py` already names
        # CH3-NH2 as 'methanamine', confirming the 'amine'->'aminide'
        # suffix-replacement rule. PubChem CID 21952893 (structure match;
        # PubChem's own generated name, 'methylazanide', uses a different
        # naming system).
        ("C[NH-]", "methanaminide"),
        # PubChem CID 187872 (structure match; 'ethylazanide' there too) --
        # P-14.3.4.2(b) locant omission on an unsubstituted 2-carbon chain,
        # mirroring `_amine.py`'s own 'ethanamine'.
        ("CC[NH-]", "ethanaminide"),
        # PubChem CID 23176953 (structure match) -- matches this project's
        # existing `_amine.py` locant convention exactly
        # (smiles_to_iupac("CCCN") == "propan-1-amine").
        ("CCC[NH-]", "propan-1-aminide"),
        ("CC(C)[NH-]", "propan-2-aminide"),
        # PubChem CID 20624633 (structure match; 'tert-butylazanide' there)
        # -- the branched skeleton is handled the same way `_amine.py`'s
        # own chain search + branch-substituent naming handles
        # 2-methylpropan-2-amine.
        ("CC(C)(C)[NH-]", "2-methylpropan-2-aminide"),
    ],
)
def test_aminide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aminide_with_halogen_substituent():
    # Mirrors `_amine.py`'s own halogen-coexistence support; mononuclear
    # parent locant omitted (P-14.3.4.2(a)).
    assert smiles_to_iupac("FC[NH-]") == "fluoromethanaminide"


def test_aminide_with_halogen_on_two_carbon_chain():
    # Mirrors this project's established '2-fluoroethan-1-ol'/
    # '2-fluoroethan-1-olate' locant-citation convention (a substituent
    # breaks the 2-carbon chain's symmetry, so the locant is cited).
    assert smiles_to_iupac("FCC[NH-]") == "2-fluoroethan-1-aminide"


def test_aromatic_aminide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH-]c1ccccc1")


def test_two_aminide_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH-]CC[NH-]")


def test_enamine_aminide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C[NH-]")


def test_second_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC[NH-]")
