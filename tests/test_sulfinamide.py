import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # '-sulfinamide' mirrors '-sulfonamide' (see test_sulfonamide.py)
        # with one fewer oxygen; same locant rules as '-sulfinic acid'.
        ("CS(=O)N", "methanesulfinamide"),
        ("CCS(=O)N", "ethanesulfinamide"),
        ("CCCS(=O)N", "propane-1-sulfinamide"),
        ("CC(S(=O)N)C", "propane-2-sulfinamide"),
        ("CCCCS(=O)N", "butane-1-sulfinamide"),
    ],
)
def test_saturated_sulfinamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_sulfinamide():
    assert smiles_to_iupac("C=CCS(=O)N") == "prop-2-ene-1-sulfinamide"


def test_halogen_substituent():
    # PubChem's own generated name ("1-chloroethanesulfinamide") omits the
    # locant; this project cites it once a substituent is present on a
    # 2-carbon chain, the same accepted divergence as
    # '2-chloroethane-1-selenol' (see test_selenol.py).
    assert smiles_to_iupac("CC(Cl)S(=O)N") == "1-chloroethane-1-sulfinamide"


def test_ene_carbon_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=C(S(=O)N)C")


def test_two_sulfinamides_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CS(=O)N")


def test_sulfonamide_not_confused_with_sulfinamide():
    assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"


def test_sulfinic_acid_not_confused_with_sulfinamide():
    assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"


def test_cyclohexanesulfinamide():
    # PubChem structure match: "cyclohexanesulfinamide".
    assert smiles_to_iupac("O=S(N)C1CCCCC1") == "cyclohexanesulfinamide"


def test_2_methylcyclohexane_1_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1C") == "2-methylcyclohexane-1-sulfinamide"


def test_cyclopentanesulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCC1") == "cyclopentanesulfinamide"


def test_2_chlorocyclohexane_1_sulfinamide():
    assert smiles_to_iupac("O=S(N)C1CCCCC1Cl") == "2-chlorocyclohexane-1-sulfinamide"


def test_polycyclic_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CC2CCC1CC2")


def test_unsaturated_ring_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=S(N)C1CCCC=C1")


def test_sulfinamide_on_ring_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CC1CCCCC1")


def test_sulfinamide_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NS(=O)CCO")


def test_phenyl_chain_sulfinamide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#286's
    # carboxylic-acid/.../sulfonamide chains): the ring is cited as a
    # "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCS(=O)N") == "3-phenylpropane-1-sulfinamide"


def test_phenyl_directly_attached_sulfinamide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1S(=O)N")


def test_phenyl_chain_sulfinamide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CCCS(=O)NC")


def test_phenyl_substituted_benzene_ring_sulfinamide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CCS(=O)N")


def test_phenyl_chain_sulfinamide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS(=O)N")


def test_n_methylmethanesulfinamide():
    # PubChem structure match: "N-methylmethanesulfinamide" (CID 12734342).
    assert smiles_to_iupac("CS(=O)NC") == "N-methylmethanesulfinamide"


def test_n_n_dimethylmethanesulfinamide():
    assert smiles_to_iupac("CS(=O)N(C)C") == "N,N-dimethylmethanesulfinamide"


def test_n_methylcyclohexanesulfinamide():
    assert smiles_to_iupac("O=S(NC)C1CCCCC1") == "N-methylcyclohexanesulfinamide"


def test_branched_n_substituted_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)NC(C)C")


def test_unsaturated_n_substituted_sulfinamide_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CS(=O)NCC=C")


def test_sulfinamide_unspecified_stereocenter_unaffected():
    # A genuine chain stereocenter left unspecified (no @/@@) is named
    # exactly as before -- no error, matching this project's long-standing
    # convention for unspecified stereochemistry.
    assert smiles_to_iupac("CCC(C)S(=O)N") == "butane-2-sulfinamide"


def test_sulfinamide_specified_chain_stereocenter_raises():
    # The sulfinamide sulfur (-R, =O, -NH2) is itself a potential
    # stereocenter in this molecule too (unlike `_sulfonamide.py`'s sulfur,
    # whose two identical =O make it never stereogenic) -- so a specified
    # chain stereocenter here always coexists with an unspecified sulfur
    # one, and `specified_stereocenters` correctly rejects the combination
    # (P-92), same pattern as `_sulfinic_acid.py`.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC[C@@H](C)S(=O)N")


def test_sulfinamide_specified_sulfur_stereocenter_raises():
    # A specified sulfur configuration, with no chain stereocenter at all,
    # is explicitly out of scope too -- this project has no established
    # locant/prefix convention for a heteroatom-centered stereodescriptor.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCC[S@](=O)N")


def test_ring_sulfinamide_non_stereocenter_marker_unaffected():
    # The ring carbon bearing -S(=O)NH2 here isn't a genuine stereocenter
    # (both ring neighbors are identical -CH2- groups), so RDKit never
    # reports it as a stereo element at all -- this stray '@' marker is
    # silently ignored exactly as before, same policy as every other
    # module's "unspecified/non-stereogenic marker" handling.
    assert smiles_to_iupac("O=S(N)[C@H]1CCCCC1") == "cyclohexanesulfinamide"
