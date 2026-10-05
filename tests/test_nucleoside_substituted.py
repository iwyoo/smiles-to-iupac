import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-(propane-2,2-diyl)adenosine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(=O)nc(N)[nH]c21",
            "2′,3′-O-(propane-2,2-diyl)guanosine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1ccc(N)nc1=O",
            "2′,3′-O-(propane-2,2-diyl)cytidine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1ccc(=O)[nH]c1=O",
            "2′,3′-O-(propane-2,2-diyl)uridine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(=O)nc[nH]c21",
            "2′,3′-O-(propane-2,2-diyl)inosine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(=O)[nH]c(=O)[nH]c21",
            "2′,3′-O-(propane-2,2-diyl)xanthosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OCO[C@H]21",
            "2′,3′-O-methyleneadenosine",
        ),
        (
            "O=c1ccn([C@@H]2O[C@H](CO)[C@H]3OCO[C@H]32)c(=O)[nH]1",
            "2′,3′-O-methyleneuridine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OC(c3ccccc3)O[C@H]21",
            "2′,3′-O-(phenylmethylene)adenosine",
        ),
        (
            "CC1O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-(ethane-1,1-diyl)adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OC3(CCCCC3)O[C@H]21",
            "2′,3′-O-(cyclohexane-1,1-diyl)adenosine",
        ),
        (
            "CCC1(CC)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-(pentane-3,3-diyl)adenosine",
        ),
        (
            "Cc1cn([C@H]2C[C@@H]3OC(C)(C)OC[C@H]3O2)c(=O)[nH]c1=O",
            "3′,5′-O-(propane-2,2-diyl)thymidine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@@H]2COCO[C@H]2[C@H]1O",
            "3′,5′-O-methyleneadenosine",
        ),
        (
            "CC1(C)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O1)[C@@H]2O",
            "2′,5′-O-(propane-2,2-diyl)adenosine",
        ),
        (
            "C[C@@H]1O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-[(1S)-ethane-1,1-diyl]adenosine",
        ),
        (
            "C[C@H]1O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-[(1R)-ethane-1,1-diyl]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@H](c3ccccc3)O[C@H]21",
            "2′,3′-O-[(S)-phenylmethylene]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@@H](c3ccccc3)O[C@H]21",
            "2′,3′-O-[(R)-phenylmethylene]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@H](/C=C/c3ccccc3)O[C@H]21",
            "2′,3′-O-[(1S,2E)-3-phenylprop-2-ene-1,1-diyl]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@@H](/C=C/c3ccccc3)O[C@H]21",
            "2′,3′-O-[(1R,2E)-3-phenylprop-2-ene-1,1-diyl]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@H](/C=C\\c3ccccc3)O[C@H]21",
            "2′,3′-O-[(1S,2Z)-3-phenylprop-2-ene-1,1-diyl]adenosine",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2O[C@@H](/C=C\\c3ccccc3)O[C@H]21",
            "2′,3′-O-[(1R,2Z)-3-phenylprop-2-ene-1,1-diyl]adenosine",
        ),
        (
            "COC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@@H]2OC(C)(C)O[C@@H]21",
            "5′-O-methyl-2′,3′-O-(propane-2,2-diyl)adenosine",
        ),
        (
            "CC(=O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@@H]2OC(C)(C)O[C@@H]21",
            "2′,3′-O-(propane-2,2-diyl)adenosine 5′-acetate",
        ),
        (
            "CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OC(C)(C)O[C@H]21",
            "N6-methyl-2′,3′-O-(propane-2,2-diyl)adenosine",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cc(C(=O)O)c(=O)[nH]c1=O",
            "1-[2,3-O-(propane-2,2-diyl)-β-D-ribofuranosyl]-2,4-dioxo-1,2,3,4-tetrahydropyrimidine-5-carboxylic acid",
        ),
        (
            "COC(=O)c1nc2c(N)ncnc2n1[C@@H]1O[C@H](CO)[C@H]2OC(C)(C)O[C@H]21",
            "methyl 6-amino-9-[2,3-O-(propane-2,2-diyl)-β-D-ribofuranosyl]-9H-purine-8-carboxylate",
        ),
        (
            "CC1(C)O[C@@H]2[C@H](O1)[C@@H](COS(=O)(=O)O)O[C@H]2n1cnc2c(N)ncnc21",
            "2′,3′-O-(propane-2,2-diyl)adenosine 5′-(hydrogen sulfate)",
        ),
        (
            "O=C1O[C@@H]2[C@H](O1)[C@@H](CO)O[C@H]2n1cc(C(=O)O)c(=O)[nH]c1=O",
            "1-(2,3-O-carbonyl-β-D-ribofuranosyl)-2,4-dioxo-1,2,3,4-tetrahydropyrimidine-5-carboxylic acid",
        ),
    ],
)
def test_nucleoside_acetal_bridge(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OS(=O)(=O)O[C@H]21",
            "2′,3′-dideoxyadenosine-2′,3′-diyl sulfate",
        ),
        (
            "Nc1nc(=O)c2ncn([C@@H]3O[C@H](CO)[C@H]4OS(=O)(=O)O[C@H]43)c2[nH]1",
            "2′,3′-dideoxyguanosine-2′,3′-diyl sulfate",
        ),
        (
            "O=c1ccn([C@@H]2O[C@H](CO)[C@H]3OS(=O)(=O)O[C@H]32)c(=O)[nH]1",
            "2′,3′-dideoxyuridine-2′,3′-diyl sulfate",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OS(=O)O[C@H]21",
            "2′,3′-dideoxyadenosine-2′,3′-diyl sulfite",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OC(=O)O[C@H]21",
            "2′,3′-dideoxyadenosine-2′,3′-diyl carbonate",
        ),
    ],
)
def test_nucleoside_cyclic_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "Nc1ncnc2c1nc(S(=O)(=O)O)n2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "6-amino-9-β-D-ribofuranosyl-9H-purine-8-sulfonic acid",
        ),
        (
            "O=c1[nH]c(=O)n([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)cc1S(=O)(=O)O",
            "2,4-dioxo-1-β-D-ribofuranosyl-1,2,3,4-tetrahydropyrimidine-5-sulfonic acid",
        ),
        (
            "Nc1ncnc2c1nc(C(=O)Cl)n2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "6-amino-9-β-D-ribofuranosyl-9H-purine-8-carbonyl chloride",
        ),
        (
            "O=C(Br)c1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",
            "2,4-dioxo-1-β-D-ribofuranosyl-1,2,3,4-tetrahydropyrimidine-5-carbonyl bromide",
        ),
        (
            "CC(=O)OC(=O)c1nc2c(N)ncnc2n1[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "acetic 6-amino-9-β-D-ribofuranosyl-9H-purine-8-carboxylic anhydride",
        ),
        (
            "CC(=O)OC(=O)c1nc2c(=O)nc(N)[nH]c2n1[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "acetic 2-amino-6-oxo-9-β-D-ribofuranosyl-6,9-dihydro-3H-purine-8-carboxylic anhydride",
        ),
        (
            "CCC(=O)OC(=O)c1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)nc1N",
            "4-amino-2-oxo-1-β-D-ribofuranosyl-1,2-dihydropyrimidine-5-carboxylic propanoic anhydride",
        ),
        (
            "CNC(=O)c1nc2c(N)ncnc2n1[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "6-amino-N-methyl-9-β-D-ribofuranosyl-9H-purine-8-carboxamide",
        ),
        (
            "CN(C)C(=O)c1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",
            "N,N-dimethyl-2,4-dioxo-1-β-D-ribofuranosyl-1,2,3,4-tetrahydropyrimidine-5-carboxamide",
        ),
        (
            "CCNC(=O)c1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)nc1N",
            "4-amino-N-ethyl-2-oxo-1-β-D-ribofuranosyl-1,2-dihydropyrimidine-5-carboxamide",
        ),
        (
            "O=C(O)CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]ethanoic acid",
        ),
        (
            "COC(=O)CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "methyl 2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]ethanoate",
        ),
        (
            "NC(=O)CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]acetamide",
        ),
        (
            "N#CCNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]acetonitrile",
        ),
        (
            "O=CCNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]acetaldehyde",
        ),
        (
            "O=S(=O)(O)CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]methanesulfonic acid",
        ),
        (
            "O=C(Cl)CNc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]ethanoyl chloride",
        ),
        (
            "O=C(O)CNc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
            "2-[(2-oxo-1-β-D-ribofuranosyl-1,2-dihydropyrimidin-4-yl)amino]ethanoic acid",
        ),
        (
            "O=C(O)CCNc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
            "3-[(2-oxo-1-β-D-ribofuranosyl-1,2-dihydropyrimidin-4-yl)amino]propanoic acid",
        ),
        (
            "O=C(O)CNc1nc(=O)c2ncn([C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O)c2[nH]1",
            "2-[(6-oxo-9-β-D-ribofuranosyl-6,9-dihydro-3H-purin-2-yl)amino]ethanoic acid",
        ),
        (
            "Nc1nc(OCC(=O)O)c2ncn([C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O)c2n1",
            "2-[(2-amino-9-β-D-ribofuranosyl-9H-purin-6-yl)oxy]ethanoic acid",
        ),
        (
            "CN(CC(=O)O)c1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[methyl(9-β-D-ribofuranosyl-9H-purin-6-yl)amino]ethanoic acid",
        ),
        (
            "OC(=O)CNc1ncnc2c1nc(Br)n2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-[(8-bromo-9-β-D-ribofuranosyl-9H-purin-6-yl)amino]ethanoic acid",
        ),
    ],
)
def test_nucleoside_base_senior_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COS(=O)(=O)O)[C@@H](O)[C@H]1O",
            "adenosine 5′-(hydrogen sulfate)",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](OS(=O)(=O)O)[C@H]1OS(=O)(=O)O",
            "adenosine 2′,3′-bis(hydrogen sulfate)",
        ),
        (
            "Cc1cn([C@H]2C[C@H](O)[C@@H](COS(=O)(=O)O)O2)c(=O)[nH]c1=O",
            "thymidine 5′-(hydrogen sulfate)",
        ),
        (
            "Cc1ccc(S(=O)(=O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)cc1",
            "adenosine 5′-(4-methylbenzenesulfonate)",
        ),
        (
            "CS(=O)(=O)OC[C@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H]1O",
            "uridine 5′-methanesulfonate",
        ),
        (
            "CC(=O)OC[C@H]1O[C@@H](n2cnc3c(=O)nc(N)[nH]c32)[C@H](O)[C@@H]1OS(=O)(=O)O",
            "guanosine 5′-acetate 3′-(hydrogen sulfate)",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COC(=O)CC(=O)O)[C@@H](O)[C@H]1O",
            "adenosine 5′-(hydrogen propanedioate)",
        ),
        (
            "O=C(O)CCC(=O)OC[C@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H]1O",
            "uridine 5′-(hydrogen butanedioate)",
        ),
        (
            "O=C(O)c1cn([C@@H]2O[C@H](COS(=O)(=O)O)[C@@H](O)[C@H]2O)c(=O)[nH]c1=O",
            "1-(5-O-sulfo-β-D-ribofuranosyl)-2,4-dioxo-1,2,3,4-tetrahydropyrimidine-5-carboxylic acid",
        ),
    ],
)
def test_nucleoside_sugar_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C[n+]1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21",
            "7-methylguanosin-7-ium",
        ),
        ("Nc1[nH+]cnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "adenosin-1-ium"),
        (
            "C[n+]1cnc2c(ncn2[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c1N",
            "1-methyladenosin-1-ium",
        ),
        (
            "C[n+]1c(N)ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c1=O",
            "3-methylcytidin-3-ium",
        ),
        ("Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH+]1", "cytidin-3-ium"),
        (
            "C[n+]1cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21.[Cl-]",
            "7-methylguanosin-7-ium chloride",
        ),
        (
            "C[n+]1cnc2c(nc(Br)n2[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c1N",
            "8-bromo-1-methyladenosin-1-ium",
        ),
        (
            "C[n+]1c(C(=O)O)n([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21",
            "2-amino-7-methyl-6-oxo-9-β-D-ribofuranosyl-6,9-dihydro-1H-purin-7-ium-8-carboxylic acid",
        ),
        (
            "Nc1nc2c(c(=O)[nH]1)[n+](CC(=O)O)cn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "2-(2-amino-6-oxo-9-β-D-ribofuranosyl-6,9-dihydro-1H-purin-7-ium-7-yl)ethanoic acid",
        ),
        (
            "CC[n+]1cn([C@@H]2O[C@H](COC(C)=O)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21",
            "7-ethylguanosin-7-ium 5′-acetate",
        ),
        (
            "C[n+]1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c2nc(N)[nH]c(=O)c21",
            "2′-deoxy-7-methylguanosin-7-ium",
        ),
        ("Nc1[nH+]cnc2c1ncn2[C@H]1C[C@H](O)[C@@H](CO)O1", "2′-deoxyadenosin-1-ium"),
    ],
)
def test_nucleoside_cation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "C[n+]1c(N)c(C(=O)O)cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c1=O",
            "6-amino-1-methyl-2-oxo-3-β-D-ribofuranosyl-2,3-dihydropyrimidin-1-ium-5-carboxylic acid",
        ),
        (
            "Nc1c(C=O)cn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH+]1",
            "6-amino-2-oxo-3-β-D-ribofuranosyl-2,3-dihydropyrimidin-1-ium-5-carbaldehyde",
        ),
    ],
)
def test_nucleoside_cation_pyrimidine_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC1(C)O[C@H]2[C@@H](O1)[C@@H](CO)O[C@@H]2n1cnc2c(N)ncnc21",
    ],
)
def test_nucleoside_wrong_sugar_stereo_keeps_no_retained_name(smiles):
    assert "adenosine" not in smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COC=O)[C@@H](O)[C@H]1O",
            "adenosine 5′-formate",
        ),
        (
            "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COC(=O)C(=O)O)[C@@H](O)[C@H]1O",
            "adenosine 5′-(hydrogen oxalate)",
        ),
        (
            "Nc1ncnc2c1nc(C(=O)OC=O)n2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
            "6-amino-9-β-D-ribofuranosyl-9H-purine-8-carboxylic formic anhydride",
        ),
    ],
)
def test_nucleoside_retained_acid_words(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
