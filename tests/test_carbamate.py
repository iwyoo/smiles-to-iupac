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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 15628, 77922, 10973, 10984.
        ("CC(C)OC(N)=O", "propan-2-yl carbamate"),
        ("CC(C)(C)OC(N)=O", "tert-butyl carbamate"),
        ("CC(C)COC(N)=O", "2-methylpropyl carbamate"),
        ("CC(C)CCOC(N)=O", "3-methylbutyl carbamate"),
    ],
)
def test_branched_r(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methyl_n_methylcarbamate():
    # PubChem CID 81151.
    assert smiles_to_iupac("CNC(=O)OC") == "methyl N-methylcarbamate"


def test_ethyl_n_ethylcarbamate():
    # PubChem CID 12195.
    assert smiles_to_iupac("CCNC(=O)OCC") == "ethyl N-ethylcarbamate"


def test_ethyl_n_methylcarbamate():
    # PubChem CID 7752.
    assert smiles_to_iupac("CNC(=O)OCC") == "ethyl N-methylcarbamate"


def test_methyl_n_n_dimethylcarbamate():
    # PubChem structure match: "methyl N,N-dimethylcarbamate".
    assert smiles_to_iupac("COC(=O)N(C)C") == "methyl N,N-dimethylcarbamate"


def test_methyl_n_ethyl_n_methylcarbamate():
    # PubChem structure match: "methyl N-ethyl-N-methylcarbamate".
    assert smiles_to_iupac("COC(=O)N(C)CC") == "methyl N-ethyl-N-methylcarbamate"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 568334, 575684.
        ("COC(=O)NC(C)C", "methyl N-propan-2-ylcarbamate"),
        ("COC(=O)NC(C)(C)C", "methyl N-tert-butylcarbamate"),
    ],
)
def test_branched_n_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_two_identical_branched_n_substituents_parenthesized_when_compound():
    # PubChem CID 568417: a compound name (has its own locant) is
    # parenthesized under the 'di' multiplying prefix to avoid ambiguity.
    assert smiles_to_iupac("COC(=O)N(C(C)C)C(C)C") == "methyl N,N-di(propan-2-yl)carbamate"


def test_two_identical_branched_n_substituents_not_parenthesized_when_retained():
    # PubChem CID 12567954: a non-compound (retained) name is not
    # parenthesized under 'di'.
    assert smiles_to_iupac("COC(=O)N(C(C)(C)C)C(C)(C)C") == "methyl N,N-ditert-butylcarbamate"


def test_two_different_n_substituents_alphabetized_ignoring_italic_prefix():
    # PubChem CID 87089877: alphanumerical ordering (P-14.5.2) ignores
    # italicized 'tert-', so 'tert-butyl' sorts under 'b', ahead of 'ethyl'.
    assert smiles_to_iupac("COC(=O)N(CC)C(C)(C)C") == "methyl N-tert-butyl-N-ethylcarbamate"


def test_free_carbamic_acid_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(N)=O")


def test_cyclic_r_group_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)OC1CCCC1")


def test_unsaturated_r_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=CCOC(N)=O")
