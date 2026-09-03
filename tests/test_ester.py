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
        # A plain, unsubstituted benzene ring on the acyl chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/#271's
        # carboxylic-acid/ketone/alcohol chains): the ring is cited as a
        # "phenyl" substituent prefix. Cross-checked against PubChem CID
        # 7643 ("methyl 3-phenylpropanoate").
        ("c1ccccc1CCC(=O)OC", "methyl 3-phenylpropanoate"),
        # Two-carbon acyl chain: PubChem's own name for this SMILES
        # ('ethyl 2-phenylacetate', CID 7590) uses the retained 'acetate'
        # stem, but this module always uses the systematic 'ethanoate' stem
        # (see 'methyl ethanoate' above), so this is an accepted, reviewed
        # result rather than a PubChem-confirmed one -- same policy as the
        # carboxylic-acid module's '2-phenylethanoic acid' (PR #269).
        ("c1ccccc1CC(=O)OCC", "ethyl 2-phenylethanoate"),
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


def test_phenyl_ester_oxygen_side_raises():
    # The benzene ring attached via the ester oxygen (an aryl ester, e.g.
    # phenyl acetate) rather than the acyl chain is out of scope for this
    # first slice (module docstring: R' must be a plain unbranched alkyl).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)Oc1ccccc1")


def test_phenyl_substituted_benzene_ring_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)OC")


def test_phenyl_acyl_chain_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)OC")


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
