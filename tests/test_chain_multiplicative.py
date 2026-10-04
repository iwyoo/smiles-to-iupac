import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCOCCOCCO", "2,2'-[ethane-1,2-diylbis(oxy)]di(ethan-1-ol)"),
        ("OCCOCCOCCOCCO", "2,2'-[oxybis(ethane-2,1-diyloxy)]di(ethan-1-ol)"),
        ("OC(=O)COCOCC(=O)O", "2,2'-[methylenebis(oxy)]diethanoic acid"),
        ("OC(=O)CSCCSCC(O)=O", "2,2'-[ethane-1,2-diylbis(sulfanediyl)]diethanoic acid"),
        ("OC(=O)c1ccc(OCCOc2ccc(C(O)=O)cc2)cc1", "4,4'-[ethane-1,2-diylbis(oxy)]dibenzoic acid"),
        ("OC(=O)c1ccc(OCCOCCOc2ccc(C(O)=O)cc2)cc1", "4,4'-[oxybis(ethane-2,1-diyloxy)]dibenzoic acid"),
        ("OC(=O)CCOOCCC(O)=O", "3,3'-peroxydipropanoic acid"),
        ("OCCOP(C)OCCO", "2,2'-[(methylphosphanediyl)bis(oxy)]di(ethan-1-ol)"),
        ("OCCN(C)CCO", "2,2'-(methylazanediyl)di(ethan-1-ol)"),
        ("OCC[SiH2]CCO", "2,2'-silanediyldi(ethan-1-ol)"),
        ("OC[SiH2]CC[SiH2]CO", "[ethane-1,2-diylbis(silanediyl)]dimethanol"),
        ("C[Si](C)(CCC(O)=O)CCC(O)=O", "3,3'-(dimethylsilanediyl)dipropanoic acid"),
        ("OC(=O)CC[Si](C)(C)[Si](C)(C)CCC(O)=O", "3,3'-(1,1,2,2-tetramethyldisilane-1,2-diyl)dipropanoic acid"),
        ("OCCP(CCO)CCO", "2,2',2''-phosphanetriyltri(ethan-1-ol)"),
        ("OCCN(CCO)CCN(CCO)CCO", "2,2',2'',2'''-(ethane-1,2-diyldinitrilo)tetra(ethan-1-ol)"),
        ("OC(=O)CN(CC(O)=O)CCN(CC(O)=O)CC(O)=O", "2,2',2'',2'''-(ethane-1,2-diyldinitrilo)tetraethanoic acid"),
    ],
)
def test_heteroatom_and_concatenated_linkers_between_chain_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCOc1ccc(OCCO)cc1", "2,2'-[1,4-phenylenebis(oxy)]di(ethan-1-ol)"),
        ("OCCOc1ccc(Cl)c(OCCO)c1", "2,2'-[(4-chloro-1,3-phenylene)bis(oxy)]di(ethan-1-ol)"),
        ("OC(=O)COC1CCC(CC1)OCC(O)=O", "2,2'-[cyclohexane-1,4-diylbis(oxy)]diethanoic acid"),
        ("OC(=O)COc1ccc2ccccc2c1OCC(O)=O", "2,2'-[naphthalene-1,2-diylbis(oxy)]diethanoic acid"),
        ("OCCOc1ccc(cc1)Cc1ccc(OCCO)cc1", "2,2'-[methylenebis(4,1-phenyleneoxy)]di(ethan-1-ol)"),
        ("OCCc1ccc(cc1)c1ccc(cc1)CCO", "2,2'-([1,1'-biphenyl]-4,4'-diyl)di(ethan-1-ol)"),
        ("OC(=O)c1ccc(cc1)[SiH2]c1ccc(cc1)C(=O)O", "4,4'-silanediyldibenzoic acid"),
        ("Oc1ccc(cc1)P(C)c1ccc(cc1)O", "4,4'-(methylphosphanediyl)diphenol"),
        ("OCCCc1ccc(cc1)c1ccc(cc1)c1ccc(cc1)CCCO", "3,3'-([11,21:24,31-terphenyl]-14,34-diyl)di(propan-1-ol)"),
        ("OCCOCCOc1ccc(OCCOCCO)cc1", "2,2'-[1,4-phenylenebis(oxyethane-2,1-diyloxy)]di(ethan-1-ol)"),
        ("OCCOc1ccc(cc1)-c1ccc(OCCO)cc1", "2,2'-[[1,1'-biphenyl]-4,4'-diylbis(oxy)]di(ethan-1-ol)"),
        (
            "OC(=O)CCCCCCCCCOc1ccc(cc1)-c1ccc(OCCCCCCCCCC(O)=O)cc1",
            "10,10'-[[1,1'-biphenyl]-4,4'-diylbis(oxy)]di(decanoic acid)",
        ),
        (
            "OCCOc1ccc(cc1)-c1ccc(cc1)-c1ccc(OCCO)cc1",
            "2,2'-[[11,21:24,31-terphenyl]-14,34-diylbis(oxy)]di(ethan-1-ol)",
        ),
    ],
)
def test_ring_and_assembly_components_inside_a_linker(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H](O)COCCOC[C@@H](C)O", "(2R,2'S)-1,1'-[ethane-1,2-diylbis(oxy)]di(propan-2-ol)"),
        ("C[C@H](O)COCCOC[C@H](C)O", "(2S,2'S)-1,1'-[ethane-1,2-diylbis(oxy)]di(propan-2-ol)"),
        ("C[C@H](O)CN(C)C[C@H](C)O", "(2S,2'S)-1,1'-(methylazanediyl)di(propan-2-ol)"),
    ],
)
def test_stereodescriptors_of_multiplied_units(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NCCNCCNCCN", "N1-{2-[(2-aminoethyl)amino]ethyl}ethane-1,2-diamine"),
        ("OCCOc1cc(O)ccc1OCCO", None),
        ("OCCOCCOCCOCCOCCO", "3,6,9,12-tetraoxatetradecane-1,14-diol"),
    ],
)
def test_substitutive_or_skeletal_replacement_names_stay_when_multiplicative_does_not_apply(smiles, expected):
    if expected is None:
        with pytest.raises(Exception):
            smiles_to_iupac(smiles)
    else:
        assert smiles_to_iupac(smiles) == expected
