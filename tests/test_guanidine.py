import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_guanidine():
    # PubChem structure match: "guanidine".
    assert smiles_to_iupac("NC(=N)N") == "guanidine"


def test_n_methylguanidine():
    # PubChem structure match ("2-methylguanidine", deprecated numeral
    # locant); this project uses the Blue Book's own letter-locant style
    # established by `_urea.py`/`_thiourea.py`.
    assert smiles_to_iupac("CNC(=N)N") == "N-methylguanidine"


def test_n_ethylguanidine():
    assert smiles_to_iupac("CCNC(=N)N") == "N-ethylguanidine"


def test_n_n_dimethylguanidine_same_nitrogen():
    # PubChem structure match: "1,1-dimethylguanidine".
    assert smiles_to_iupac("CN(C)C(=N)N") == "N,N-dimethylguanidine"


def test_n_ethyl_n_methylguanidine_same_nitrogen():
    assert smiles_to_iupac("CCN(C)C(=N)N") == "N-ethyl-N-methylguanidine"


def test_n_n_prime_dimethylguanidine_different_nitrogens():
    # Blue Book P-66.4.1.2.1.2 worked example: "N,N′-dimethylguanidine
    # (PIN) (not N,N′′-dimethylguanidine)". PubChem structure match:
    # "1,2-dimethylguanidine" (deprecated numeral locants).
    assert smiles_to_iupac("CNC(=N)NC") == "N,N'-dimethylguanidine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # One different substituent on each amino nitrogen -- the
        # alphabetically first substituent name becomes 'N-', the other
        # 'N''-', same rule as `_urea.py`/`_thiourea.py` (PubChem
        # structure match, its own numeric-locant style): `CCNC(=N)NC` ->
        # '1-ethyl-2-methylguanidine' (CID 17814701), `CCCNC(=N)NC` ->
        # '2-methyl-1-propylguanidine' (CID 20383691).
        ("CCNC(=N)NC", "N-ethyl-N'-methylguanidine"),
        ("CCCNC(=N)NC", "N-methyl-N'-propylguanidine"),
    ],
)
def test_different_substituents_on_different_amino_nitrogens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_different_substituent_counts_on_amino_nitrogens_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(C)C(=N)NC")


def test_n_double_prime_methylguanidine():
    # PubChem structure match: "2-methylguanidine" (deprecated numeral
    # locant for the imino nitrogen).
    assert smiles_to_iupac("NC(=NC)N") == "N''-methylguanidine"


def test_n_double_prime_ethylguanidine():
    # PubChem structure match: "2-ethylguanidine".
    assert smiles_to_iupac("NC(=NCC)N") == "N''-ethylguanidine"


def test_imino_plus_one_amino_substituted():
    # PubChem structure match: "1,2-dimethylguanidine" (deprecated numeral
    # locants for one amino + the imino nitrogen).
    assert smiles_to_iupac("CNC(=NC)N") == "N,N''-dimethylguanidine"


def test_imino_plus_both_amino_substituted_same_name():
    assert smiles_to_iupac("CNC(=NC)NC") == "N,N',N''-trimethylguanidine"


def test_imino_plus_amino_different_names_alphabetized():
    # Different substituent names are cited alphabetically (P-14.5.2)
    # regardless of which nitrogen (letter) they sit on: "ethyl" sorts
    # before "methyl", so the N''-ethyl group is cited first even though
    # its letter has more primes than the amino N-methyl group.
    assert smiles_to_iupac("CNC(=NCC)N") == "N''-ethyl-N-methylguanidine"


def test_branched_imino_n_substituent():
    # Same structure as CC(C)NC(=N)N (PubChem CID 11491919), just written
    # with the substituent on the double-bonded nitrogen -- this module's
    # existing convention (see module docstring) always labels the
    # explicit C=N nitrogen 'N''', regardless of which tautomer the input
    # SMILES happens to spell out (matches the already-passing unbranched
    # 'N''-methylguanidine' case above).
    assert smiles_to_iupac("NC(=NC(C)C)N") == "N''-(propan-2-yl)guanidine"


def test_unsaturated_imino_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=NC=C)N")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 11491919, 12830400.
        ("CC(C)NC(=N)N", "N-(propan-2-yl)guanidine"),
        ("CC(C)(C)NC(=N)N", "N-tert-butylguanidine"),
    ],
)
def test_branched_n_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_identical_branched_substituents_different_amino_nitrogens_parenthesized_when_compound():
    # PubChem CID 198192.
    assert smiles_to_iupac("CC(C)NC(=N)NC(C)C") == "N,N'-di(propan-2-yl)guanidine"


def test_two_identical_branched_substituents_different_amino_nitrogens_not_parenthesized_when_retained():
    # PubChem CID 23103888.
    assert smiles_to_iupac("CC(C)(C)NC(=N)NC(C)(C)C") == "N,N'-ditert-butylguanidine"


def test_unsaturated_n_substituent_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CNC(=N)N")


def test_amidine_not_confused_with_guanidine():
    assert smiles_to_iupac("CC(=N)N") == "ethanimidamide"


def test_urea_not_confused_with_guanidine():
    assert smiles_to_iupac("NC(=O)N") == "urea"


def test_n_phenylguanidine():
    # PubChem structure match: "2-phenylguanidine".
    assert smiles_to_iupac("c1ccccc1NC(=N)N") == "N-phenylguanidine"


def test_n_n_prime_diphenylguanidine():
    # PubChem structure match: "1,2-diphenylguanidine".
    assert smiles_to_iupac("c1ccccc1NC(=N)Nc1ccccc1") == "N,N'-diphenylguanidine"


def test_n_methyl_n_prime_phenylguanidine_different_nitrogens():
    # PubChem structure match: "2-methyl-1-phenylguanidine".
    assert smiles_to_iupac("CNC(=N)Nc1ccccc1") == "N-methyl-N'-phenylguanidine"


def test_n_double_prime_phenylguanidine():
    assert smiles_to_iupac("NC(=Nc1ccccc1)N") == "N''-phenylguanidine"


def test_substituted_phenyl_n_substituent_not_supported():
    # Only a plain, unsubstituted benzene ring is supported, mirroring
    # `_urea.py`'s identical restriction.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc(NC(=N)N)cc1")


def test_phenyl_alongside_another_substituent_on_same_nitrogen_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CN(c1ccccc1)C(=N)N")


def test_ring_fused_guanidine_still_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N1CCCCC1=N")
