import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # -oyl halide is a mechanical stem+halide-word construction (P-65.1.5,
        # see _acyl_halide.py docstring) with a single, non-ambiguous
        # substitution point; no locant tie-break or alphabetization choice
        # is involved, so these are cross-checked directly against common
        # usage (acetyl chloride == 'ethanoyl chloride' etc.) rather than a
        # per-case PubChem lookup.
        ("C(=O)Cl", "methanoyl chloride"),
        ("CC(=O)Cl", "ethanoyl chloride"),
        ("CCC(=O)Cl", "propanoyl chloride"),
        ("CCCC(=O)F", "butanoyl fluoride"),
        ("CCCC(=O)Br", "butanoyl bromide"),
        ("CCCC(=O)I", "butanoyl iodide"),
    ],
)
def test_saturated_acyl_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_acyl_halide():
    assert smiles_to_iupac("C=CC(=O)Cl") == "prop-2-enoyl chloride"


def test_other_halogen_substituent_coexists():
    assert smiles_to_iupac("CC(Cl)C(=O)Cl") == "2-chloropropanoyl chloride"


def test_branched_substituent():
    assert smiles_to_iupac("CC(C)CC(=O)Cl") == "3-methylbutanoyl chloride"


def test_two_acyl_halides_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC(=O)CC(=O)Cl")


def test_cyclic_acyl_halide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClC(=O)C1CCCCC1")


def test_coexisting_ketone_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)CC(=O)Cl")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified (CID 7157125/7006443/7568783): the acyl halide
        # carbon is always C1 (see module docstring), so the stereocenter
        # prefix mirrors _carboxylic_acid.py's P-91.3 mechanism exactly.
        ("C[C@H](Cl)C(=O)Cl", "(2S)-2-chloropropanoyl chloride"),
        ("C[C@@H](Cl)C(=O)Cl", "(2R)-2-chloropropanoyl chloride"),
        ("CC[C@H](Cl)C(=O)Cl", "(2S)-2-chlorobutanoyl chloride"),
    ],
)
def test_stereocenter_on_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CC(C)C(Cl)C(=O)Cl") == "2-chloro-3-methylbutanoyl chloride"


def test_branch_stereocenter_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C[C@H](C)Cl)C(=O)Cl")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/
        # #271/#272/#273/#274/#275's carboxylic-acid/ketone/alcohol/ester/
        # aldehyde/amide/nitrile chains): the ring is cited as a "phenyl"
        # substituent prefix. Cross-checked against PubChem CID 61529
        # ("3-phenylpropanoyl chloride", hydrocinnamoyl chloride).
        ("c1ccccc1CCC(=O)Cl", "3-phenylpropanoyl chloride"),
        # Two-carbon chain: this module always uses the systematic
        # 'ethanoyl' stem (see test_saturated_acyl_halide above) rather
        # than a retained "phenylacetyl" form, so this is an accepted,
        # reviewed result -- same policy as the aldehyde module's
        # '2-phenylethanal' (PR #273).
        ("c1ccccc1CC(=O)Cl", "2-phenylethanoyl chloride"),
        # A different halide word coexists the same way (P-35.2.1's
        # module-level analog for the acyl halogen itself).
        ("c1ccccc1CCC(=O)Br", "3-phenylpropanoyl bromide"),
    ],
)
def test_phenyl_chain_acyl_halide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phenyl_directly_attached_acyl_halide_raises():
    # -C(=O)X directly on the ring uses a separate naming construction
    # (benzoyl-style), out of scope for this acyclic-chain-parent module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)Cl")


def test_phenyl_substituted_benzene_ring_acyl_halide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)Cl")


def test_phenyl_chain_acyl_halide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)Cl")
