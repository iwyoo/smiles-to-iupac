import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # cross-checked against PubChem PUG REST IUPACName:
        # CID 6331 (B), CID 5326078 (CB), CID 15879358 (CBC),
        # CID 68979 (CB(C)C), CID 7357 (CCB(CC)CC), CID 143352201 (CCBC),
        # CID 517955 (CCB(C)C), CID 543198 (CCB(CC)C).
        ("B", "borane"),
        ("CB", "methylborane"),
        ("CBC", "dimethylborane"),
        ("CB(C)C", "trimethylborane"),
        ("CCB(CC)CC", "triethylborane"),
        ("CCBC", "ethyl(methyl)borane"),
        # the multiplying prefix sits *outside* the parentheses (Blue
        # Book "chlorodi(methyl)borane (PIN)" worked example,
        # `tmp/bluebook/P6.txt`), not PubChem's own raw
        # "ethyl(dimethyl)borane" for CID 517955's structure.
        ("CCB(C)C", "ethyldi(methyl)borane"),
        ("CCB(CC)C", "diethyl(methyl)borane"),
        # halogen bonded directly to boron: structure cross-checked
        # against PubChem PUG REST, CID 140714 (ClB), CID 137221
        # (ClB(C)C) -- the second name matches the Blue Book's own
        # "chlorodi(methyl)borane (PIN)" directly, not PubChem's own raw
        # "chloro(dimethyl)borane".
        ("ClB", "chloroborane"),
        ("ClB(C)C", "chlorodi(methyl)borane"),
    ],
)
def test_smiles_to_iupac_simple_borane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogen_substituted_alkyl_chain():
    assert smiles_to_iupac("ClCCB") == "(2-chloroethyl)borane"


def test_borane_chain_raises():
    # diborane (B-B bond): two boron atoms, a separate "chain" nomenclature
    # problem this module doesn't handle.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("BB")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: CID 101871698, 5326178.
        ("CC(C)B", "propan-2-ylborane"),
        ("CC(C)(C)B", "tert-butylborane"),
    ],
)
def test_branched_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_substituent():
    assert smiles_to_iupac("C=CB") == "ethenylborane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: `c1ccccc1B` -> "phenylborane" (CID 101697465),
        # `c1ccccc1B(c1ccccc1)c1ccccc1` -> "triphenylborane" (CID 70400).
        # Unlike `_phosphanone.py`'s identical extension, boron's normal
        # valence is exactly 3, so three phenyls never trigger a
        # lambda-convention label.
        ("c1ccccc1B", "phenylborane"),
        ("c1ccccc1B(c1ccccc1)c1ccccc1", "triphenylborane"),
    ],
)
def test_phenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_mixed_alkyl_and_phenyl_substituents():
    # Not itself a Blue Book/PubChem worked example, but a direct
    # generalization of the two confirmed shapes above (mixed substituent
    # identity is already handled generically by
    # `format_mononuclear_prefixes`).
    assert smiles_to_iupac("CB(c1ccccc1)C") == "dimethyl(phenyl)borane"


def test_substituted_phenyl():
    assert smiles_to_iupac("Cc1ccccc1B") == "(2-methylphenyl)borane"


def test_non_aromatic_ring():
    assert smiles_to_iupac("C1CCCCC1B") == "cyclohexylborane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified: `Clc1ccc(cc1)B` -> "(4-chlorophenyl)borane" (CID
        # 102445216), `Fc1ccc(cc1)B` -> "(4-fluorophenyl)borane" (CID
        # 101045919) -- see `_phosphane.py`'s identical, more extensively
        # verified extension.
        ("Clc1ccc(cc1)B", "(4-chlorophenyl)borane"),
        ("Fc1ccc(cc1)B", "(4-fluorophenyl)borane"),
    ],
)
def test_halogenated_phenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_halogenated_phenyl_mixed_with_alkyl():
    assert smiles_to_iupac("CB(c1ccc(Cl)cc1)") == "(4-chlorophenyl)(methyl)borane"
