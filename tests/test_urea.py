import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_urea():
    # PubChem structure match: "urea".
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_methylurea():
    # PubChem structure match ("methylurea"); this project uses the Blue
    # Book's own letter-locant style confirmed directly from source text
    # (tmp/bluebook/P6.txt lines 630, 1373-1378), not PubChem's numeric
    # locants.
    assert smiles_to_iupac("CNC(=O)N") == "N-methylurea"


def test_n_n_dimethylurea_same_nitrogen():
    # PubChem structure match ("1,1-dimethylurea" -- both methyls on the
    # same nitrogen).
    assert smiles_to_iupac("CN(C)C(=O)N") == "N,N-dimethylurea"


def test_n_ethyl_n_methylurea_same_nitrogen():
    # PubChem structure match ("1-ethyl-1-methylurea" -- both
    # substituents on the same nitrogen).
    assert smiles_to_iupac("CCN(C)C(=O)N") == "N-ethyl-N-methylurea"


def test_n_n_prime_dimethylurea_different_nitrogens():
    # PubChem structure match ("1,3-dimethylurea" -- one methyl on each of
    # the two different nitrogens); Blue Book P6.txt lines 1373-1378
    # confirm the primed N,N' convention for two distinct nitrogens.
    assert smiles_to_iupac("CNC(=O)NC") == "N,N'-dimethylurea"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # One different substituent on each nitrogen -- the alphabetically
        # first substituent name (P-14.5.2, locants ignored) becomes 'N-',
        # the other becomes 'N''-', confirmed via the Blue Book's own PIN
        # worked example 'N-[1-cyano-3-(methylsulfanyl)propyl]-
        # N'-methylurea'. PubChem structure match (its own numeric-locant
        # style): 'ethyl' < 'methyl' -> `CCNC(=O)NC` -> '1-ethyl-3-methylurea'
        # (CID 206567); 'methyl' < 'propyl' -> `CCCNC(=O)NC` ->
        # '1-methyl-3-propylurea' (CID 217014).
        ("CCNC(=O)NC", "N-ethyl-N'-methylurea"),
        ("CCCNC(=O)NC", "N-methyl-N'-propylurea"),
    ],
)
def test_different_substituents_on_different_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_different_substituent_counts_on_different_nitrogens_not_supported():
    # One nitrogen with two substituents, the other with one -- still out
    # of scope: no confirmed worked example settles this locant tie-break.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCN(C)C(=O)NC")


def test_branched_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)NC(=O)N")


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=O)N")


def test_thiourea_not_confused_with_urea():
    # Thiourea is supported by its own module (`_thiourea.py`, see
    # test_thiourea.py) but must not be mistaken for plain urea.
    assert smiles_to_iupac("NC(=S)N") == "thiourea"


def test_semicarbazide():
    # PubChem structure match: "aminourea" (CID 5196).
    assert smiles_to_iupac("NC(=O)NN") == "aminourea"


def test_semicarbazide_with_n_alkyl_substituent_not_supported():
    # PubChem's own name for this combination ("1-amino-3-methylurea",
    # CID 256005) uses a numeric-locant style this module doesn't
    # otherwise follow for urea -- out of scope until that's resolved.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CNC(=O)NN")


def test_double_amino_substituted_urea_not_supported():
    # Carbonohydrazide (H2N-NH-C(=O)-NH-NH2) -- amino on both nitrogens is
    # out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NNC(=O)NN")
