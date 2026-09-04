import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_selenourea():
    # PubChem structure match: NC(=[Se])N is CID 6327594, synonym
    # "Selenourea"; Blue Book P-66.1.6.1.3.1 confirms 'selenourea' as PIN.
    assert smiles_to_iupac("NC(=[Se])N") == "selenourea"


def test_n_methylselenourea():
    assert smiles_to_iupac("CNC(=[Se])N") == "N-methylselenourea"


def test_n_butan_2_yl_selenourea():
    # Blue Book P-66.1.6.1.3.1's own worked example, confirmed verbatim:
    # "N-(butan-2-yl)selenourea (PIN)" (`tmp/bluebook/P6a.txt` line 1143).
    assert smiles_to_iupac("CCC(C)NC(=[Se])N") == "N-(butan-2-yl)selenourea"


def test_n_n_dimethylselenourea_same_nitrogen():
    assert smiles_to_iupac("CN(C)C(=[Se])N") == "N,N-dimethylselenourea"


def test_n_ethyl_n_methylselenourea_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=[Se])N") == "N-ethyl-N-methylselenourea"


def test_n_n_prime_dimethylselenourea_different_nitrogens():
    assert smiles_to_iupac("CNC(=[Se])NC") == "N,N'-dimethylselenourea"


def test_different_substituents_on_different_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCNC(=[Se])NC")


def test_branched_n_substituent():
    assert smiles_to_iupac("CC(C)NC(=[Se])N") == "N-(propan-2-yl)selenourea"


def test_n_tert_butylselenourea():
    assert smiles_to_iupac("CC(C)(C)NC(=[Se])N") == "N-tert-butylselenourea"


def test_two_identical_branched_n_substituents_parenthesized_when_compound():
    assert smiles_to_iupac("CC(C)NC(=[Se])NC(C)C") == "N,N'-di(propan-2-yl)selenourea"


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=[Se])N")


def test_thiourea_not_confused_with_selenourea():
    assert smiles_to_iupac("NC(=S)N") == "thiourea"


def test_urea_not_confused_with_selenourea():
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_phenylselenourea():
    # PubChem structure match: CID 6329315, synonym "Phenylselenourea"
    # (IUPACName unavailable for this selenium compound).
    assert smiles_to_iupac("NC(=[Se])Nc1ccccc1") == "N-phenylselenourea"


def test_n_n_prime_diphenylselenourea():
    # PubChem structure match: CID 6327858, synonym
    # "N,N'-Diphenylselenourea".
    assert smiles_to_iupac("c1ccc(NC(=[Se])Nc2ccccc2)cc1") == "N,N'-diphenylselenourea"


def test_substituted_phenyl_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=[Se])Nc1ccc(C)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=[Se])N")
