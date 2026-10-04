import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C1CCC(=[SnH2])CC1", "4-stannylidenecyclohexane-1-carboxylic acid"),
        ("OC(=O)C1CCC(=[AsH])CC1", "4-arsanylidenecyclohexane-1-carboxylic acid"),
        ("OC(=O)C#[SiH]", "2-silylidyneethanoic acid"),
        ("OC(=O)C#P", "2-phosphanylidyneethanoic acid"),
        ("NNNc1ccc(C(=O)O)cc1", "4-(triazan-1-yl)benzoic acid"),
        ("[SiH3]N[SiH2]c1ccc(C(=O)O)cc1", "4-[(silylamino)silyl]benzoic acid"),
        ("[SiH3]N[SiH](N[SiH3])c1ccc(C(=O)O)cc1", "4-[bis(silylamino)silyl]benzoic acid"),
        ("OC(=O)C1CCC(=NN)CC1", "4-hydrazinylidenecyclohexane-1-carboxylic acid"),
        ("[SiH3][SiH]([SiH3])[SiH2][SiH2]c1ccc(C(=O)O)cc1", "4-(3-silyltetrasilan-1-yl)benzoic acid"),
        ("OC(=O)Cc1cccc([SiH2]O[SiH3])n1", "2-(6-disiloxanylpyridin-2-yl)ethanoic acid"),
        ("C(OC)PC1CCCCC1(C#N)", "2-[(methoxymethyl)phosphanyl]cyclohexane-1-carbonitrile"),
    ],
)
def test_heteroatom_hydride_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles", ["C1CCCCC1=NN", "CO[Al](C)C", "NCB(C)O"])
def test_heteroatom_hydride_parent_without_principal_group_raises(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
