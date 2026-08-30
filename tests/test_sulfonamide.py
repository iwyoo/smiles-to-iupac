import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanesulfonamide():
    # PubChem structure match: "methanesulfonamide".
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_ethanesulfonamide():
    # PubChem structure match: "ethanesulfonamide".
    assert smiles_to_iupac("CCS(=O)(=O)N") == "ethanesulfonamide"


def test_propane_1_sulfonamide():
    # PubChem structure match: "propane-1-sulfonamide".
    assert smiles_to_iupac("CCCS(=O)(=O)N") == "propane-1-sulfonamide"


def test_propane_2_sulfonamide():
    assert smiles_to_iupac("CC(S(=O)(=O)N)C") == "propane-2-sulfonamide"


def test_chlorobutanesulfonamide():
    assert smiles_to_iupac("ClCCCCS(=O)(=O)N") == "4-chlorobutane-1-sulfonamide"


def test_2_chloroethane_1_sulfonamide():
    # Locant cited even on a 2-carbon chain once a substituent (here, the
    # halogen) is present -- same project-wide convention as
    # '2-chloroethane-1-selenol'/'2-chloroethane-1-thiol' (PubChem's own
    # generated name omits the locant here; that divergence is accepted
    # project-wide, see test_selenol.py's identical case).
    assert smiles_to_iupac("ClCCS(=O)(=O)N") == "2-chloroethane-1-sulfonamide"


def test_pent_4_ene_1_sulfonamide():
    assert smiles_to_iupac("C=CCCCS(=O)(=O)N") == "pent-4-ene-1-sulfonamide"


def test_disulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CCS(=O)(=O)N")


def test_cyclohexanesulfonamide():
    # PubChem structure match: "cyclohexanesulfonamide".
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1") == "cyclohexanesulfonamide"


def test_2_methylcyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1C") == "2-methylcyclohexane-1-sulfonamide"


def test_cyclopentanesulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCC1") == "cyclopentanesulfonamide"


def test_2_chlorocyclohexane_1_sulfonamide():
    assert smiles_to_iupac("O=S(=O)(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfonamide"


def test_n_methylcyclohexanesulfonamide():
    # PubChem structure match: "N-methylcyclohexanesulfonamide" (CID 23534457).
    assert smiles_to_iupac("O=S(=O)(NC)C1CCCCC1") == "N-methylcyclohexanesulfonamide"


def test_polycyclic_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(=O)(N)C1CCCC=C1")


def test_sulfonamide_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CC1CCCCC1")


def test_sulfonamide_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)(=O)CCO")


def test_n_methylmethanesulfonamide():
    # PubChem structure match: "N-methylmethanesulfonamide" (CID 97632).
    assert smiles_to_iupac("CS(=O)(=O)NC") == "N-methylmethanesulfonamide"


def test_n_n_dimethylmethanesulfonamide():
    # PubChem structure match: "N,N-dimethylmethanesulfonamide" (CID 70191).
    assert smiles_to_iupac("CS(=O)(=O)N(C)C") == "N,N-dimethylmethanesulfonamide"


def test_n_ethyl_n_methylethanesulfonamide():
    # PubChem structure match: "N-ethyl-N-methylethanesulfonamide" (CID 21102946).
    assert smiles_to_iupac("CCS(=O)(=O)N(C)CC") == "N-ethyl-N-methylethanesulfonamide"


def test_branched_n_substituted_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)(=O)NC(C)C")


def test_unsaturated_n_substituted_sulfonamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)(=O)NCC=C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_sulfonic_acid.py`. PubChem CID 93472539.
        ("CC[C@@H](C)S(=O)(=O)N", "(2R)-butane-2-sulfonamide"),
        ("CC[C@H](C)S(=O)(=O)N", "(2S)-butane-2-sulfonamide"),
    ],
)
def test_acyclic_sulfonamide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_acyclic_sulfonamide_stereocenter_with_coexisting_substituent():
    # A stereocenter that also bears a halogen substituent: the suffix's
    # own locant is still cited on the 2-carbon chain despite the
    # substituent sharing its position, same project-wide convention as
    # the identical pre-existing 'ClC(C)S(=O)(=O)N' case. PubChem's own
    # auto-generated name omits that locant ('(1R)-1-chloroethanesulfonamide',
    # CID 92264750) -- a non-PIN quirk already documented elsewhere in this
    # project -- so only the structure is cross-checked there.
    assert smiles_to_iupac("C[C@@H](Cl)S(=O)(=O)N") == "(1R)-1-chloroethane-1-sulfonamide"


def test_cyclic_sulfonamide_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_sulfonic_acid.py`'s `_name_cyclic_sulfonic_acid`. PubChem has no
    # registered CID for this exact stereoisomer (CID 0), so only the
    # achiral structure is cross-checked (CID 130649687,
    # '2-chlorocyclohexane-1-sulfonamide').
    assert (
        smiles_to_iupac("N[S](=O)(=O)[C@H]1CCCC[C@@H]1Cl")
        == "(1S,2S)-2-chlorocyclohexane-1-sulfonamide"
    )


def test_sulfonamide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(Cl)S(=O)(=O)N") == "1-chloropropane-1-sulfonamide"


def test_sulfonamide_partially_specified_stereocenters_raises():
    # Only one of the ring's two genuine stereocenters is marked -- must
    # raise rather than silently dropping the marker.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[S](=O)(=O)[C@H]1CCCCC1Cl")
