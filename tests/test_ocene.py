import pytest

from smiles_to_iupac import smiles_to_iupac

PRIME = "′"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Os+2].OCC[c-]1cccc1.[cH-]1cccc1", "2-(osmocen-1-yl)ethanol"),
        ("[Fe+2].CC(=O)[c-]1cccc1.CC(=O)[c-]1cccc1", f"1,1{PRIME}-(ferrocene-1,1{PRIME}-diyl)di(ethanone)"),
        ("[Fe+2].C[c-]1cccc1.[cH-]1cccc1", "1-methylferrocene"),
        ("[Fe+2].OC(=O)[c-]1cccc1.[cH-]1cccc1", "ferrocene-1-carboxylic acid"),
        ("[Fe+2].OC(=O)[c-]1cccc1.OC(=O)[c-]1cccc1", f"ferrocene-1,1{PRIME}-dicarboxylic acid"),
        ("[Fe+2].N[c-]1cccc1.[cH-]1cccc1", "ferrocen-1-amine"),
        ("[Fe+2].O[c-]1cccc1.[cH-]1cccc1", "ferrocen-1-ol"),
        ("[Fe+2].N#C[c-]1cccc1.[cH-]1cccc1", "ferrocene-1-carbonitrile"),
        ("[Fe+2].OC(=O)[c-]1cccc1.CC(=O)[c-]1cccc1", f"1{PRIME}-acetylferrocene-1-carboxylic acid"),
        ("[Fe+2].OC(=O)CC[c-]1cccc1.[cH-]1cccc1", "3-(ferrocen-1-yl)propanoic acid"),
        ("[V+2].CN(C)C(C)[c-]1cccc1.[cH-]1cccc1", "N,N-dimethyl-1-(vanadocen-1-yl)ethanamine"),
    ],
)
def test_substituted_metallocenes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ansa_metallocene_cited_as_divalent_group():
    name = smiles_to_iupac("OC(=O)CC([c-]1cccc1)CC[c-]1cccc1.[Fe+2]")
    assert name == f"3,5-(ferrocene-1,1{PRIME}-diyl)pentanoic acid"


def test_two_ruthenocenes_joined_by_a_chain():
    smiles = (
        "C[C]12->[Ru]3456789(<-[CH](=[CH]->3[CH2]->41)[CH]->5=2)<-[CH]1=[CH]->6[CH2]->7[C]->8(CC[C]23->"
        "[Ru]45678%10%11(<-[CH]%12=[CH]->4[CH2]->5[C]->6(C)=[CH]->7%12)<-[CH](=[CH]->8[CH2]->%102)[CH]->%11=3)=[CH]->91"
    )
    assert smiles_to_iupac(smiles) == f"1,1{PRIME}{PRIME}-(ethane-1,2-diyl)bis(1{PRIME}-methylruthenocene)"


def test_benzoferrocene():
    smiles = "c1cc[c]23->[Fe]456789%10(<-[CH]%11=[CH]->4[CH2]->5[CH]->6=[CH]->7%11)<-[CH](=[CH]->8[c]->92c1)[CH2]->%103"
    assert smiles_to_iupac(smiles) == "benzoferrocene"


def test_two_ferrocenes_in_a_ring_are_named_as_a_phane():
    smiles = (
        "C1[C]23->[Fe]456789%10(<-[CH](=[CH]->4[CH2]->52)[CH]->6=3)<-[CH]2=[CH]->7[CH2]->8[C]->9(C[C]34->"
        "[Fe]56789%11%12(<-[CH]%13=[CH]->5[CH2]->6[C]->71=[CH]->8%13)<-[CH](=[CH]->9[CH2]->%113)[CH]->%12=4)=[CH]->%102"
    )
    assert smiles_to_iupac(smiles) == f"1,3(1,1{PRIME})-diferrocenacyclotetraphane"


def test_three_ferrocenes_in_a_ring_are_named_as_a_phane():
    smiles = (
        "C1[C]23->[Fe]456789%10(<-[CH](=[CH]->4[CH2]->52)[CH]->6=3)<-[CH]2=[CH]->7[CH2]->8[C]->9(C[C]34->[Fe]"
        "56789%11%12(<-[CH](=[CH]->5[CH2]->63)[CH]->7=4)<-[CH]3=[CH]->8[CH2]->9[C]->%11(C[C]45->[Fe]6789%11%1"
        "3%14(<-[CH]%15=[CH]->6[CH2]->7[C]->81=[CH]->9%15)<-[CH](=[CH]->%11[CH2]->%134)[CH]->%14=5)=[CH]->%12"
        "3)=[CH]->%102"
    )
    assert smiles_to_iupac(smiles) == f"1,3,5(1,1{PRIME})-triferrocenacyclohexaphane"
