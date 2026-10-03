import pytest

from smiles_to_iupac import smiles_to_iupac

H = "\u03b7"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[cH]12->[Cr]3456<-[cH]1[cH]->3[cH]->4[cH]->5[cH]->62", f"({H}6-benzene)chromium"),
        (
            "[O+]#[C][Mn]1234([C]#[O+])([C]#[O+])[CH]5[CH]->1=[CH]->2[CH]->3=[CH]->45",
            f"tricarbonyl({H}5-cyclopenta-2,4-dien-1-yl)manganese",
        ),
        ("[CH2]1[CH]2=[CH2]->[Cr]<-21", f"({H}3-allyl)chromium"),
        (
            "[CH2]1[CH]2=[CH2]->[Cr]<-213456([CH2][CH]->3=[CH2]->4)[CH2][CH]->5=[CH2]->6",
            f"tris({H}3-allyl)chromium",
        ),
        ("C1C[CH]2->[Rh]34<-[CH]1=[CH]->3CC[CH]->4=2", f"[(1,2,5,6-{H})-cycloocta-1,5-diene]rhodium"),
        (
            "[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[CH]4=[CH]->1C1CC4[CH]->2=[CH]->31",
            f"[(2,3,5,6-{H})-bicyclo[2.2.1]hepta-2,5-diene]tricarbonyliron",
        ),
        (
            "[O+]#[C][Mo+]123456([C]#[O+])([C]#[O+])[CH]7[CH]->1=[CH]->2[CH]->3=[CH]->4[CH]->5=[CH]->67",
            f"tricarbonyl({H}7-cyclohepta-2,4,6-trien-1-yl)molybdenum(1+)",
        ),
        ("[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[CH2]=[CH]->1[CH]->2=[CH2]->3", f"({H}4-buta-1,3-diene)tricarbonyliron"),
    ],
)
def test_hapto_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "NCC[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->41)[cH]->5[cH]->62",
            f"tricarbonyl[2-({H}6-phenyl)ethanamine]chromium",
        ),
        (
            "C[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[c]->51CCN)[cH]->62",
            f"tricarbonyl[2-(2-methyl-{H}6-phenyl)ethanamine]chromium",
        ),
        (
            "CC(N(C)C)[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[c]->51P(c1ccccc1)c1ccccc1)[cH]->62",
            f"tricarbonyl{{1-[2-(diphenylphosphanyl)-{H}6-phenyl]-N,N-dimethylethanamine}}chromium",
        ),
        (
            "C[c]12->[Cr]3456([C]#[O+])([C]#[O+])([C]#[O+])<-[cH]([cH]->3[cH]->4[cH]->51)[cH]->62",
            f"tricarbonyl({H}6-methylbenzene)chromium",
        ),
        (
            "c1ccc([B-](c2ccccc2)(c2ccccc2)[c]23->[Rh+]456789%10(<-[CH]%11=[CH]->4CC[CH]->5=[CH]->6CC%11)"
            "<-[cH]([cH]->7[cH]->82)[cH]->9[cH]->%103)cc1",
            f"[(1,2,5,6-{H})-cycloocta-1,5-diene][triphenyl({H}6-phenyl)borato]rhodium",
        ),
    ],
)
def test_hapto_ring_inside_larger_ligand(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
