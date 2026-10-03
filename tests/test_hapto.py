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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('[O+]#[C][Mo]123456([C]#[O+])([CH]7C=CC=C[CH]->1=[CH]->27)[CH]1[CH]->3=[CH]->4[CH]->5=[CH]->61', 'dicarbonyl[(1–3-η)-cyclohepta-2,4,6-trien-1-yl](η5-cyclopenta-2,4-dien-1-yl)molybdenum'),
        ('[O+]#[C][Fe]1234([C]#[O+])[CH]5[CH]->1=[CH]->2[c]->31cccc[c]->415', 'dicarbonyl[(1–3,3a,7a-η)-1H-inden-1-yl]iron'),
        ('[O+]#[C][Fe]1234([C]#[O+])[CH]5[c]->16cccc[c]->26-[c]->31cccc[c]->415', 'dicarbonyl[(4a,4b,8a,9,9a-η)-9H-fluoren-9-yl]iron'),
        ('C[C]12->[Fe]345([C]#[O+])([C]#[O+])[CH]1[c]->31cccc[c]->41[CH]->5=2', 'dicarbonyl[2-methyl-(1–3,3a,7a-η)-1H-inden-1-yl]iron'),
        ('CO[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52', 'tricarbonyl(1-methoxy-η5-cyclopenta-2,4-dien-1-yl)manganese'),
        ('COC(=O)[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52', 'tricarbonyl[1-(methoxycarbonyl)-η5-cyclopenta-2,4-dien-1-yl]manganese'),
        ('[O+]#[C][Mn]1234([C]#[O+])([C]#[O+])<-[cH]5[cH]->1[cH]->2[cH-]->3[cH]->45', 'tricarbonyl(η5-cyclopenta-2,4-dien-1-ido)manganate(1-)'),
        ('C[Si](C)(C)[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52', 'tricarbonyl[1-(trimethylsilyl)-η5-cyclopenta-2,4-dien-1-yl]manganese'),
        ('CN(C)[C]12[Mn]345([C]#[O+])([C]#[O+])([C]#[O+])<-[CH](=[CH]->31)[CH]->4=[CH]->52', 'tricarbonyl[1-(dimethylamino)-η5-cyclopenta-2,4-dien-1-yl]manganese'),
        ('C[C]12->[Fe]345([C]#[O+])([C]#[O+])<-[CH](=[CH]->3[N]->41)[CH]->5=2', 'dicarbonyl(2-methyl-η5-1H-pyrrol-1-yl)iron'),
        ('[O+]#[C][Fe]1234([C]#[O+])<-[CH]5=[CH]->1[N]->2[CH]->3=[CH]->45', 'dicarbonyl(η5-1H-pyrrol-1-yl)iron'),
        ('[O+]#[C][Cr]1234([C]#[O+])([C]#[O+])<-[cH]5[cH]->1[cH]->2[s]->3[cH]->45', 'tricarbonyl(η5-thiophene)chromium'),
        ('CCCC[CH]1=[CH2]->[Pt]<-1([Cl])[Cl]', 'dichlorido(η2-hex-1-ene)platinum'),
        ('C=CCC[CH]1=[CH2]->[Pt]<-1([Cl])[Cl]', 'dichlorido[(1,2-η)-hexa-1,5-diene]platinum'),
        ('CC[CH]1=[CH]2[CH]3=[CH2]->[Fe]<-3<-1<-2([C]#[O+])([C]#[O+])[C]#[O+]', 'tricarbonyl(η4-hexa-1,3-diene)iron'),
        ('C[CH]1=[CH]2[CH2][Pd]<-1<-2[Cl]', '[(1–3-η)-but-2-en-1-yl]chloridopalladium'),
        ('[O+]#[C][Fe]1234([C]#[O+])([C]#[O+])[CH2][CH]->1=[CH]->2[CH]->3=[CH2]->4', 'tricarbonyl(η5-penta-2,4-dien-1-yl)iron'),
        ('[O+]#[C][Cr]12345([C]#[O+])([C]#[O+])<-[cH]6[cH]->1[cH]->2[c]->31cccc[c]->41[cH]->56', 'tricarbonyl[(1–4,4a,8a-η)-naphthalene]chromium'),
        ('[O+]#[C][Mn]1234([C]#[O+])([C]#[O+])[CH]5C[CH]->1=[CH]->2[CH]->3=[CH]->45', 'tricarbonyl(η5-cyclohexa-2,4-dien-1-yl)manganese'),
        ('[Cl][Pt]1([Cl])<-[CH]#[CH]->1', '(η2-acetylene)dichloridoplatinum'),
        ('[O+]#[C][Fe]123([C]#[O+])([C]#[O+])<-[cH]4cc[c]56->[Fe]1789([C]#[O+])([C]#[O+])<-[cH]([cH]->75)[cH]->8[c]->9-6[cH]->2[cH]->34', '{μ-[2(1–3,3a,8a-η):1(4–6-η)]azulene}-(pentacarbonyl-1κ3C,2κ2C)diiron(Fe—Fe)'),
        ('[O+]#[C][Fe]12([C]#[O+])([C]#[O+])<-[CH2]=[CH]->1[CH]1=[CH2]->[Fe]<-12([C]#[O+])([C]#[O+])[C]#[O+]', '{μ-[2(1,2-η):1(3,4-η)]buta-1,3-diene}-(hexacarbonyl-1κ3C,2κ3C)diiron(Fe—Fe)'),
    ],
)
def test_extended_hapto_ligands(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
