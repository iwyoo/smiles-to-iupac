import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCC(CC1)OC(C)=O", "cyclohexane-1,4-diyl diethanoate"),
        ("CC(=O)OC1CCCCC1OC(C)=O", "cyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)OC1CCCC(C1)OC(C)=O", "cyclohexane-1,3-diyl diethanoate"),
        ("CC(=O)OC1CCCC1OC(C)=O", "cyclopentane-1,2-diyl diethanoate"),
        ("CC(=O)OC1CC1OC(C)=O", "cyclopropane-1,2-diyl diethanoate"),
        ("CC(=O)OC1(OC(C)=O)CCCCC1", "cyclohexane-1,1-diyl diethanoate"),
        ("CCC(=O)OC1CCC(CC1)OC(=O)CC", "cyclohexane-1,4-diyl dipropanoate"),
        ("CC(=O)OC1CCC(C)CC1OC(C)=O", "4-methylcyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)OC1CCCC(OC(C)=O)C1Cl", "2-chlorocyclohexane-1,3-diyl diethanoate"),
        ("CC(=O)OC1CC(OC(=O)CC)CC(OC(C)=O)C1", "cyclohexane-1,3,5-triyl 1,3-diethanoate 5-propanoate"),
        ("CC(=O)OC1C=CC(OC(C)=O)CC1", "cyclohex-2-ene-1,4-diyl diethanoate"),
    ],
)
def test_carbocycle_diyl_diester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccccc1OC(C)=O", "1,2-phenylene diethanoate"),
        ("CC(=O)Oc1cccc(OC(C)=O)c1", "1,3-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(OC(C)=O)cc1", "1,4-phenylene diethanoate"),
        ("CC(=O)Oc1cc(OC(C)=O)cc(OC(C)=O)c1", "benzene-1,3,5-triyl triethanoate"),
        ("CC(=O)Oc1ccc(C)c(OC(C)=O)c1", "4-methyl-1,3-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(Br)c(Cl)c1OC(C)=O", "4-bromo-3-chloro-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(OC(C)=O)c(CC)c1", "2-ethyl-1,4-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(OC(C)=O)c(C(C)C)c1C", "2-methyl-3-(propan-2-yl)-1,4-phenylene diethanoate"),
    ],
)
def test_benzene_diyl_diester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1cccnc1OC(C)=O", "pyridine-2,3-diyl diethanoate"),
        ("CC(=O)Oc1ccc(OC(C)=O)cn1", "pyridine-2,5-diyl diethanoate"),
        ("CC(=O)Oc1ccoc1OC(C)=O", "furan-2,3-diyl diethanoate"),
        ("CC(=O)Oc1ccsc1OC(C)=O", "thiophene-2,3-diyl diethanoate"),
        ("CC(=O)Oc1cc[nH]c1OC(C)=O", "1H-pyrrole-2,3-diyl diethanoate"),
        ("CC(=O)Oc1ncnc(OC(C)=O)c1", "pyrimidine-4,6-diyl diethanoate"),
        ("CC(=O)Oc1nc(OC(C)=O)cs1", "1,3-thiazole-2,4-diyl diethanoate"),
        ("CC(=O)OC1CCOCC1OC(C)=O", "oxane-3,4-diyl diethanoate"),
        ("CC(=O)OC1CCCOC1OC(C)=O", "oxane-2,3-diyl diethanoate"),
        ("CC(=O)OC1CC(OC(C)=O)CO1", "oxolane-2,4-diyl diethanoate"),
        ("CC(=O)OC1CCNCC1OC(C)=O", "piperidine-3,4-diyl diethanoate"),
        ("CC(=O)OC1CCN(C)CC1OC(C)=O", "1-methylpiperidine-3,4-diyl diethanoate"),
        ("CC(=O)OC1OCC(OC(C)=O)O1", "1,3-dioxolane-2,4-diyl diethanoate"),
        ("CC(=O)OC1COC(CO1)OC(C)=O", "1,4-dioxane-2,5-diyl diethanoate"),
    ],
)
def test_heterocycle_diyl_diester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCC(CC1)OC(=O)CC", "cyclohexane-1,4-diyl ethanoate propanoate"),
        ("CC(=O)Oc1ccc(OC(=O)C(Cl)Cl)cc1", "1,4-phenylene 2,2-dichloroethanoate ethanoate"),
        ("CC(=O)Oc1ccc(C)c(OC(=O)CC)c1", "4-methyl-1,3-phenylene 1-ethanoate 3-propanoate"),
        ("CC(=O)Oc1ccc(Cl)cc1OC(=O)c1ccc(Cl)cc1", "4-chloro-1,2-phenylene 2-(4-chlorobenzoate) 1-ethanoate"),
        ("CC(=O)Oc1ccccc1OC(=O)c1ccccc1", "1,2-phenylene benzoate ethanoate"),
        ("CC(=O)OC1CCCCC1OC(=O)c1ccccc1", "cyclohexane-1,2-diyl benzoate ethanoate"),
    ],
)
def test_ring_diyl_differing_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("ClCC(=O)OC1CCC(CC1)OC(=O)CCl", "cyclohexane-1,4-diyl bis(2-chloroethanoate)"),
        ("O=C(OC1CCCCC1OC(=O)c1ccccc1)c1ccccc1", "cyclohexane-1,2-diyl dibenzoate"),
        ("O=C(Oc1ccc(OC(=O)c2ccc(Cl)cc2)cc1)c1ccc(Cl)cc1", "1,4-phenylene bis(4-chlorobenzoate)"),
        ("C=CC(=O)Oc1ccccc1OC(=O)C=C", "1,2-phenylene di(prop-2-enoate)"),
        ("O=C(OC1CCCCC1OC(=O)C1CCCCC1)C1CCCCC1", "cyclohexane-1,2-diyl di(cyclohexanecarboxylate)"),
    ],
)
def test_ring_diyl_substituted_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("O=C(OCCOC(=O)c1ccccc1)c1ccccc1", "ethane-1,2-diyl dibenzoate"),
        ("ClCC(=O)OCCCOC(=O)CCl", "propane-1,3-diyl bis(2-chloroethanoate)"),
        ("C=CC(=O)OCCOC(=O)C=C", "ethane-1,2-diyl di(prop-2-enoate)"),
        ("CC(=O)OCCOC(=O)c1ccc(Cl)cc1", "ethane-1,2-diyl 4-chlorobenzoate ethanoate"),
    ],
)
def test_acyclic_diyl_with_generalized_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)O[C@H]1CCCC[C@@H]1OC(C)=O", "(1S,2S)-cyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)O[C@@H]1CCCC[C@H]1OC(C)=O", "(1R,2R)-cyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)O[C@H]1CCCC[C@H]1OC(C)=O", "(1R,2S)-cyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)O[C@H]1CC[C@@H](CC1)OC(C)=O", "(1r,4r)-cyclohexane-1,4-diyl diethanoate"),
        ("CC(=O)O[C@H]1CC[C@H](CC1)OC(C)=O", "(1s,4s)-cyclohexane-1,4-diyl diethanoate"),
    ],
)
def test_ring_diyl_stereo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccc2ccccc2c1OC(C)=O", "naphthalene-1,2-diyl diethanoate"),
        ("CC(=O)Oc1cccc2c(OC(C)=O)cccc12", "naphthalene-1,5-diyl diethanoate"),
        ("CC(=O)Oc1ccc2cc(OC(C)=O)ccc2c1", "naphthalene-2,6-diyl diethanoate"),
        ("CC(=O)Oc1c2ccccc2c(OC(C)=O)c2ccccc12", "anthracene-9,10-diyl diethanoate"),
        ("CC(=O)Oc1ccc2[nH]ccc2c1OC(C)=O", "1H-indole-4,5-diyl diethanoate"),
        ("CC(=O)Oc1cc(OC(C)=O)c2ccccc2n1", "quinoline-2,4-diyl diethanoate"),
        ("CC(=O)OC1CCc2ccccc2C1OC(C)=O", "1,2,3,4-tetrahydronaphthalene-1,2-diyl diethanoate"),
        ("CC(=O)OC1CNc2ccccc2C1OC(C)=O", "1,2,3,4-tetrahydroquinoline-3,4-diyl diethanoate"),
        ("CC(=O)OC1C=CC2=CC=CC=C2C1OC(C)=O", "1,2-dihydronaphthalene-1,2-diyl diethanoate"),
        ("CC(=O)OC1CC2CCC1C2OC(C)=O", "bicyclo[2.2.1]heptane-2,7-diyl diethanoate"),
        ("CC(=O)OC1CCC2CCCCC2C1OC(C)=O", "decahydronaphthalene-1,2-diyl diethanoate"),
        ("CC(=O)OC1CCC2(CC1)CCC(OC(C)=O)CC2", "spiro[5.5]undecane-3,9-diyl diethanoate"),
        ("CC(=O)OC1CC2CC1C(OC(C)=O)O2", "2-oxabicyclo[2.2.1]heptane-3,5-diyl diethanoate"),
        ("CC(=O)OC1COC2(C1)CC(OC(C)=O)C2", "5-oxaspiro[3.4]octane-2,7-diyl diethanoate"),
    ],
)
def test_fused_bridged_spiro_diyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccc(O)c(OC(C)=O)c1", "4-hydroxy-1,3-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(C#N)cc1OC(C)=O", "4-cyano-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(N(=O)=O)cc1OC(C)=O", "4-nitro-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(OC)cc1OC(C)=O", "4-methoxy-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(OC(C)C)cc1OC(C)=O", "4-(propan-2-yloxy)-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(C=C)cc1OC(C)=O", "4-ethenyl-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(C(C)=O)cc1OC(C)=O", "4-ethanoyl-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(C=O)cc1OC(C)=O", "4-formyl-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(N(C)CC)cc1OC(C)=O", "4-[ethyl(methyl)amino]-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(SC)cc1OC(C)=O", "4-(methylsulfanyl)-1,2-phenylene diethanoate"),
        ("CC(=O)Oc1ccc(Oc2ccccc2)cc1OC(C)=O", "4-phenoxy-1,2-phenylene diethanoate"),
        ("CC(=O)OC1CCC(=O)CC1OC(C)=O", "4-oxocyclohexane-1,2-diyl diethanoate"),
        ("CC(=O)OC1CCNC(=O)C1OC(C)=O", "2-oxopiperidine-3,4-diyl diethanoate"),
        ("CC(=O)OCC(OC(C)=O)c1ccc(O)cc1", "1-(4-hydroxyphenyl)ethane-1,2-diyl diethanoate"),
        ("CC(=O)OCC(OC(C)=O)c1ccncc1", "1-(pyridin-4-yl)ethane-1,2-diyl diethanoate"),
    ],
)
def test_functional_substituents_on_diyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCOC=C1OC(C)=O", "3,4-dihydro-2H-pyran-4,5-diyl diethanoate"),
        ("CC(=O)OC1CN=CCC1OC(C)=O", "2,3,4,5-tetrahydropyridine-3,4-diyl diethanoate"),
        ("CC(=O)OC1COCOC1OC(C)=O", "1,3-dioxane-4,5-diyl diethanoate"),
        ("CC(=O)OC1CSCNC1OC(C)=O", "1,3-thiazinane-4,5-diyl diethanoate"),
        ("CC(=O)OC1CCCNCC1OC(C)=O", "azepane-3,4-diyl diethanoate"),
        ("CC(=O)Oc1ncnc(OC(C)=O)n1", "1,3,5-triazine-2,4-diyl diethanoate"),
        ("CC(=O)Oc1nnc(OC(C)=O)nn1", "1,2,4,5-tetraazine-3,6-diyl diethanoate"),
        ("CC(=O)Oc1nnc(OC(C)=O)o1", "1,3,4-oxadiazole-2,5-diyl diethanoate"),
        ("CC(=O)Oc1nc(OC(C)=O)n[nH]1", "1H-1,2,4-triazole-3,5-diyl diethanoate"),
        ("CC(=O)OC1OC1OC(C)=O", "oxirane-2,3-diyl diethanoate"),
    ],
)
def test_partly_saturated_and_polyheteroatom_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)Oc1ccccc1OC(=O)Cc1ccccc1", "1,2-phenylene ethanoate 2-phenylethanoate"),
        ("CC(=O)Oc1ccccc1OC(=O)c1ccncc1", "1,2-phenylene ethanoate pyridine-4-carboxylate"),
        ("CC(=O)Oc1ccccc1OC(=O)c1ccc(O)cc1", "1,2-phenylene ethanoate 4-hydroxybenzoate"),
        ("O=C(OCCOC(=O)c1ccncc1)c1ccncc1", "ethane-1,2-diyl di(pyridine-4-carboxylate)"),
        ("CC(=O)OCCOC(=O)CC#N", "ethane-1,2-diyl 2-cyanoethanoate ethanoate"),
        ("CC(=O)OCCOC(=O)CCC(=O)C", "ethane-1,2-diyl ethanoate 4-oxopentanoate"),
        ("CC(=O)OCCOC(=O)Cc1ccc(O)cc1", "ethane-1,2-diyl ethanoate 2-(4-hydroxyphenyl)ethanoate"),
        ("CC(=O)OC1CCCCC1OC(=O)C(N)C", "cyclohexane-1,2-diyl 2-aminopropanoate ethanoate"),
        ("CC(=O)OC1CCCCC1OC(=O)CC(C)=O", "cyclohexane-1,2-diyl ethanoate 3-oxobutanoate"),
        ("CC(=O)Oc1ccccc1OC(=O)/C=C/C", "1,2-phenylene (2E)-but-2-enoate ethanoate"),
        ("CC(=O)OC1CCCC1OC(=O)[C@H](C)Cl", "cyclopentane-1,2-diyl (2S)-2-chloropropanoate ethanoate"),
    ],
)
def test_generalized_acyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OC1CCC(COC(C)=O)CC1", "4-[(ethanoyloxy)methyl]cyclohexyl ethanoate"),
        ("CC(=O)Oc1ccccc1COC(C)=O", "2-[(ethanoyloxy)methyl]phenyl ethanoate"),
    ],
)
def test_ring_and_chain_polyol_uses_acyloxy_prefix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)OCCOCCOC(C)=O",
        "CC(=O)OCC(=C)COC(C)=O",
        "CC(=O)OC1CCCC[C@H]1OC(=O)[C@H](C)Cl",
        "CC(=O)OC1CCC(CC1)C1CCC(OC(C)=O)CC1",
    ],
)
def test_ring_diyl_out_of_scope(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
