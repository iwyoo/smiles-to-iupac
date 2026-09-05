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


def test_phenyl_chain_aminide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_alkoxide.py`'s
    # phenyl-chain path): the ring is cited as a "phenyl" substituent
    # prefix. PubChem structure match only (this module already uses its
    # own '-aminide' convention rather than PubChem's 'azanide' naming).
    assert smiles_to_iupac("c1ccccc1CC[NH-]") == "2-phenylethan-1-aminide"
    assert smiles_to_iupac("c1ccccc1CCC[NH-]") == "3-phenylpropan-1-aminide"


def test_phenyl_directly_attached_aminide_raises():
    # Anilinide-type naming (-NH(-) directly on the ring) is a separate
    # construction, out of scope for this acyclic-chain-parent module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1[NH-]")


def test_phenyl_chain_aminide_ring_methyl():
    # A plain methyl ring substituent, like halogens, no longer forces
    # the benzene ring to win the ring-vs-chain parent competition (see
    # test_carboxylic_acid.py's identical methylated-ring cases; this
    # project's own 'aminide' suffix convention as in the
    # unsubstituted-ring test above).
    assert smiles_to_iupac("Cc1ccccc1CC[NH-]") == "2-(2-methylphenyl)ethan-1-aminide"


def test_phenyl_chain_aminide_ring_methyl_para():
    assert smiles_to_iupac("CC1=CC=C(CC[NH-])C=C1") == "2-(4-methylphenyl)ethan-1-aminide"


def test_phenyl_chain_aminide_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propylazanide" for this
    # structure -- a different naming convention (azanide parent +
    # substituent) than this project's own established 'aminide' suffix
    # (see the module's own unsubstituted-ring test), so this checks
    # self-consistency with that convention rather than a PubChem match.
    assert smiles_to_iupac("Clc1ccc(CCC[NH-])cc1") == "3-(4-chlorophenyl)propan-1-aminide"


def test_phenyl_chain_aminide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[NH-]")
