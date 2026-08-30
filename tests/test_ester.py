import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases (methyl formate, methyl/ethyl acetate, ethyl
        # propanoate), cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName).
        ("C(=O)OC", "methyl methanoate"),
        ("CC(=O)OC", "methyl ethanoate"),
        ("CC(=O)OCC", "ethyl ethanoate"),
        ("CCC(=O)OCC", "ethyl propanoate"),
        # A real positional choice for a substituent on the acyl chain: the
        # ester carbon still fixes C1 regardless.
        ("CC(C)C(=O)OC", "methyl 2-methylpropanoate"),
        # -oate combined with existing unsaturation support, cross-checked
        # against PubChem (methyl crotonate).
        ("CC=CC(=O)OC", "methyl but-2-enoate"),
        # -oate + halogen substituent prefix on the acyl chain, cross-checked
        # against PubChem.
        ("ClCC(=O)OC", "methyl 2-chloroethanoate"),
        # A longer, unbranched alcohol part.
        ("CC(=O)OCCC", "propyl ethanoate"),
    ],
)
def test_ester_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)CC(=O)OC")


def test_branched_alcohol_part_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OC(C)C")


def test_ring_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OC1CCCCC1")


def test_amine_coexisting_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(=O)OC")


def test_aryl_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)OC")


def test_acyl_carbon_off_longest_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)C(CCC)CCCC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter on the acyl chain
        # (P-92): the stereodescriptor is cited immediately ahead of the
        # acyl part specifically, not the whole two-word ester name (CIP
        # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID 7156991.
        ("CC[C@@H](C)C(=O)OCC", "ethyl (2R)-2-methylbutanoate"),
        ("CC[C@H](C)C(=O)OCC", "ethyl (2S)-2-methylbutanoate"),
    ],
)
def test_ester_acyl_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ester_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(=O)OCC") == "ethyl 2-methylbutanoate"
