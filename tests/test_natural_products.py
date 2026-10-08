"""Natural products named on Table 10.1 parents with the modifications of P-101.3 - P-101.8."""

import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCC[C@]4(C)[C@H]3CC[C@]12C', '4-nor-5β-pregnane'),
    ('C1=CC2=C(C1)CC1NCCC23CCCCC13', '1H-4-normorphinan'),
])
def test_nor(smiles, expected):
    # P-101.3.1 nor
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CCCC[C@]4(CC)[C@H]3CC[C@]12C', '19a-homo-5β-pregnane'),
    ('CC[C@H]1CCC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C', '16a-homo-5α-pregnane'),
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]2(C)C1', '13(17)a-homo-5α-pregnane'),
    pytest.param('CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]2(C)CC1', '13(17)a,13(17)b-dihomo-5α-pregnane', marks=pytest.mark.slow),
    pytest.param('CC[C@H]1CC[C@@H]2C[C@](C)(CC[C@H]3[C@H]2CC[C@H]2CCCC[C@@]23C)C1', '13(14)a,13(17)b-dihomo-5α-pregnane', marks=pytest.mark.slow),
    ('C1=CCC2=C(C=C1)C13CCCCC1C(C2)NCC3', '1H-4a-homomorphinan'),
    ('C1=CC2=NC1=CCc1ccc([nH]1)C=C1C=CC(=N1)C=c1ccc([nH]1)=C2', '20aH-20a-homoporphyrin'),
    ('C1=CC2=NC1=Cc1ccc([nH]1)C=C1C=CC(=N1)CC=c1ccc([nH]1)=C2', '20H-20a-homoporphyrin'),
])
def test_homo(smiles, expected):
    # P-101.3.2 homo
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@]45C[C@H]4CC[C@]5(C)[C@H]3CC[C@]12C', '3α,5-cyclo-5α-pregnane'),
    pytest.param('C[C@]12CCCC[C@@H]1CC[C@@H]1[C@@H]2CC[C@]2(C)[C@@H]3[C@@H]4C[C@]12C[C@@H]43', '(20S)-14,21:16β,20-dicyclo-5α,14β-pregnane', marks=pytest.mark.slow),
    ('CC[C@H]1CN2CCc3c4n(c5ccccc35)[C@@H](C)[C@H]1C[C@@H]42', '(16βH)-1,16-cyclocorynan'),
])
def test_cyclo(smiles, expected):
    # P-101.3.3 cyclo
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize("smiles, expected", [
    pytest.param('CC[C@]1(C)[C@H]2CC[C@@H]3[C@@]4(C)CC[C@H](C(C)C)[C@@H]4CC[C@@]3(C)[C@]2(C)CC[C@H]1C(C)(C)C', '2,3-secohopane', marks=pytest.mark.slow),
    ('CC[C@@H]1CNCC[C@@]23CC[C@@H]1[C@@H](C)[C@@H]2Nc1ccccc13', '3,4-secocuran'),
])
def test_seco(smiles, expected):
    # P-101.3.4.1 seco
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1CC[C@H]2[C@H](CC[C@H]3CCCCC[C@H]3C)CCC[C@]12C', '9,10-seco-4a-homo-5α-pregnane'),
    ('CCCC[C@@]1(C)CC[C@@H]2[C@@H]1CC[C@]1(C)[C@@H](CC)CC[C@@H]21', '4,5-seco-7-norpregnane'),
    ('CC[C@]1(C)CC[C@H]([C@@]2(C)CCCCC(C)(C)[C@@H]2C)[C@@H](C)C1', '6,7-seco-3a-homopimarane'),
    ('c1cc2c3c(c[nH]c3c1)C[C@H]1NCC[C@@H]1C2', '10(11)a-homo-9-norergoline'),
    ('C[C@@H]1CC[C@@H]2C[C@@]21CC[C@@H]1CCC[C@]2(C)CCC[C@@H]12', '3α,5-cyclo-9,10-seco-5α-androstane'),
])
def test_combined(smiles, expected):
    # P-101.3.7 combinations of modifying prefixes
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('C[C@@H]1CCC[C@@H]2[C@@H]1CC[C@]1(C)CCC[C@@H]21', 'des-A-androstane'),
    ('C[C@H]1CCC[C@@H]2[C@@H]1CC[C@]1(C)CCC[C@@H]21', 'des-A-10α-androstane'),
])
def test_des(smiles, expected):
    # P-101.3.6 des
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C', '6′-apo-β-carotene'),
])
def test_apo(smiles, expected):
    # P-101.3.4.2 apo
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC(C)=CCC=C(C)C=CC=C(C)C=CC=C(C)C=CC=CC(C)=CC/C=C(C)/C=C/C1=C(C)CCCC1(C)C', '4′,11-retro-β,ψ-carotene'),
])
def test_retro(smiles, expected):
    # P-101.3.5.2 retro
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC(C)[C@@H]1CC[C@H](C)[C@@H]2CNC[C@@]2(C)C1', '3-azaambrosane'),
    pytest.param('C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CC[Te]CC[C@]3(C)[C@H]1CC2', '3-tellura-4a-homo-5α-androstane', marks=pytest.mark.slow),
    ('c1cc2c3c(csc3c1)C[C@H]1NCCC[C@H]21', '1-thiaergoline'),
    ('C1=CC2=C(CC3NCCC24CCCCC34)OC1', '2H-1-oxamorphinan'),
    ('c1ccc2c3c([nH]c2c1)[C@@H]1C[C@@H]2CCCC[C@H]2C[C@H]1CC3', '(4βH)-4-carbayohimban'),
    pytest.param('C[C@H]1[C@H]2[C@@H](C[C@H]3[C@@H]4CC[C@H]5CCCC[C@]5(C)[C@H]4CC[C@]23C)C[C@]12CC[C@@H](C)CC2', '(22r)-16a,22a-dicarba-5α,16β,25α-spirostan', marks=pytest.mark.slow),
])
def test_replacement(smiles, expected):
    # P-101.4 skeletal replacement
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3Cc4ccccc4C[C@]3(C)[C@H]1CC2', 'benzo[2,3]-5α-androstane'),
    ('C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CC4=C(CC=C4)C[C@]3(C)[C@H]1CC2', '5′H-cyclopenta[2,3]-5α-androstane'),
    ('C1=N[C@@]23[C@H](CCC[C@@]24CCN[C@@H]3Cc2ccccc24)O1', '(8αH)-[1,3]oxazolo[5′,4′:8,14]morphinan'),
    pytest.param('C[C@]12CCCC[C@@H]1c1conc1[C@@H]1[C@@H]2CC[C@]2(C)c3cnoc3C[C@@H]12', 'bis[1,2]oxazolo[4′,3′:6,7;5″,4″:16,17]-5α-androstane', marks=pytest.mark.slow),
    pytest.param('CC(C)CCCC(C)C1CCC2C3Cc4nsc5c4C(C)(CCC5)C3CCC12C', '[1,2]thiazolo[5′,4′,3′:4,5,6]cholestane', marks=pytest.mark.slow),
    ('C[C@]12C3=COCN=C1C[C@H]1[C@@H](CC[C@H]4CCCC[C@@]41C)[C@@H]2CC3', '2′H-[1,3]oxazepino[4′,5′,6′:12,13,17]-5α-androstane'),
    ('c1ccc2c(c1)C[C@H]1NCC[C@@]23CCC2=C(CCCC2)[C@@H]13', '3′,4′,5′,6′-tetrahydrobenzo[7,8]morphinan'),
    ('C[C@@]12CCC[C@H]1[C@@H]1C[C@H]3N[C@]34CCCC[C@]4(C)[C@H]1CC2', '(6αH)-1′,6-dihydroazirino[2′,3′:5,6]-5β-androstane'),
])
def test_fusion(smiles, expected):
    # P-101.5.1 fused rings and P-101.6.5 hydro prefixes on them
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('C[C@@]12CCC[C@H]1[C@@]13CC[C@H]4C[C@@H](CC[C@]4(C)[C@H]1CC2)OO3', '3α,8-epidioxy-5α,8α-androstane'),
    ('CC[C@@]12S[C@@H]1C[C@H]1[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@@]12C', '16α,17-epithio-5α-pregnane'),
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3[C@@H]3CCC[C@]12C3', '11α,18-ethano-5α,13α-pregnane'),
    pytest.param('CC(/C=C/C=C(\\C)C1C=C2C(C)(C)CCCC2(C)O1)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)C1C=C2C(C)(C)CCCC2(C)O1', '5,8:5′,8′-diepoxy-5,8,5′,8′-tetrahydro-β,β-carotene', marks=pytest.mark.slow),
    pytest.param('C1(Oc2ccccc2)OC[C@@H]3[C@H]1CO[C@H]3c4ccccc4', '(7R,8S,9′S)-7,9a′:8′,9-diepoxy-7′-oxa-9a′-homo-8,9′-neolignane', marks=pytest.mark.slow),
])
def test_bridge(smiles, expected):
    # P-101.5.2 bridges
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize("smiles, expected", [
    ('CN1C[C@]2(C[C@@H]3c4cccc5[nH]cc(c45)C[C@H]31)OCC(C)(C)O2', '(2R)-4,4,6′-trimethylspiro[1,3-dioxolane-2,8′-ergoline]'),
])
def test_spiro(smiles, expected):
    # P-101.5.3 spiro rings
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('C[C@@]12CCC[C@H]1C1=CC=C3CCCC[C@]3(C)[C@H]1CC2', 'androsta-5,7-diene'),
    ('CC(C)/C=C/C[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C', '(23E)-5α-cholest-23-ene'),
    ('C=C1CCCC/C1=C/C=C1\\CCC[C@]2(C)[C@@H]([C@H](C)CCCC(C)C)CC[C@@H]12', '(5Z,7E)-9,10-secocholesta-5,7,10(19)-triene'),
    ('c1ccc2c(c1)NC[C@]21CCN2C[C@@H]3COCC[C@H]3C[C@H]21', '16,17-dihydroformosanan'),
    ('O=C1C[C@H]2SC=CN12', '2,3-didehydropenam'),
    ('CC1=CCCC(C)(C)C1C#C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1C(C)=CCCC1(C)C', '7,8-didehydro-ε,ε-carotene'),
    pytest.param('C=C1CCO[C@H]2CCN3c4ccccc4[C@]45CCN[C@H]4C[C@@H]1[C@@H]2[C@H]35', '20,21-didehydro-21,22-dihydro-19,20-secostrychnidine', marks=pytest.mark.slow),
    ('C1=CC2OCc3ccccc3C2C2NCC=C12', '3,5-didehydrolycorenan'),
    ('C[C@H]1[C@H]2CC[C@H]3[C@@H]4CC=C5CCCC[C@]5(C)[C@H]4CC[C@]23CN1C', 'con-5-enine'),
    ('CC(/C=C/C=C(\\C)CCC1CCC(C)C(C)C1C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CCC1C(C)CCCC1(C)C', '5,6,7,8,1′,2′,3′,4′,5′,6′,7′,8′-dodecahydro-β,χ-carotene'),
    ('CC1=C(CCC(C)C)O[C@H]2C[C@H]3[C@@H]4CC[C@@H]5CCCC[C@]5(C)[C@H]4CC[C@]3(C)[C@@H]12', '5β-furost-20(22)-ene'),
    ('C/C(=C\\c1ccccc1)[C@@H](C)Cc1ccccc1', '(7E,8′S)-lign-7-ene'),
])
def test_unsaturation(smiles, expected):
    # P-101.6 degree of hydrogenation
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1C[C@@]2(C)CCC[C@H](C)[C@@H]2C[C@@H]1C(C)C', '8α-ethyleudesmane'),
    ('C[C@]12CCCC[C@@H]1CC[C@@H]1[C@@H]2CC[C@@]2(C)[C@H]1CC[C@@]2(C)CO', '(17β-methyl-5α-androstan-17α-yl)methanol'),
    ('CCCCCCCC[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(=O)O', '11α-hydroxy-9-oxoprostan-1-oic acid'),
    ('CCCC(=O)O[C@@H]1CCC[C@]23c4ccccc4CCN2CC[C@H]13', 'erythrinan-1β-yl butanoate'),
    ('CC(C)[C@H](C)CC[C@@H](C)[C@@]1(C)CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@@]21C', '17-methyl-5α-campestane'),
    pytest.param('CC1(C)CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@@]21C', '17,17-dimethyl-5α-androstane', marks=pytest.mark.slow),
    ('C[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C', '17β-methyl-5α-androstane'),
])
def test_substituents(smiles, expected):
    # P-101.7.1 substituents
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    pytest.param('O=C1C[C@@]23CCCN4CC[C@]5(c6ccccc6N[C@@H]5CC2)[C@]43O1', '19,21-epoxyaspidospermidin-21-one', marks=pytest.mark.slow),
    ('CC[C@@]12CCCN3CC[C@]4(c5ccccc5N[C@@H]4[C@H]4OC(=O)O[C@H]41)[C@H]32', 'aspidospermidine-3α,4α-diyl carbonate'),
    ('CC1(C)O[C@H]2[C@H]3CN4CCCC[C@@H]4[C@H]4CCCN(C[C@H]2O1)[C@@H]34', 'propan-2-one matridine-3β,4β-diyl ketal'),
    ('Oc1ccc2cc1-c1cccc(c1)CCC(O)CC1CCCC(CC(O)CC2)N1', '2″-demethoxylythranidine'),
    pytest.param('COc1ccc2cc1-c1cccc(c1)CCC(O)CC1CCCC(CC(O)CC2)N1', '6′-deoxylythranidine', marks=pytest.mark.slow),
])
def test_groups(smiles, expected):
    # P-101.7.4 and P-101.7.5 ring groups and de
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CCCC[C@@]4(C)[C@@H]3CC[C@]12C', '5β,9β,10α-pregnane'),
    ('CC(C)[C@@H]1CC[C@H]2[C@@H](CC[C@H]3C(C)(C)CCC[C@]23C)C1', '13β-abietane'),
    ('CC[C@]12C=Cn3c4c(c5ccccc53)CCN(CCC1)[C@H]42', '3α-eburnamenine'),
    ('C=C1CC[C@H](O)C/C1=C/C=C1\\CCC[C@]2(C)[C@@H]([C@H](C)CCCC(C)C)CC[C@@H]12', '(3S,5Z,7E)-9,10-secocholesta-5,7,10(19)-trien-3-ol'),
])
def test_configuration(smiles, expected):
    # P-101.2.6 and P-101.8 configuration
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles, expected", [
    ('CC(CCCc1ccccc1)Cc1ccccc1', '8,9′-neolignane'),
])
def test_parents(smiles, expected):
    # P-101.2.7 parent structures
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC12CC[C@@H]3[C@H](CCC4CC(=O)CC[C@]43C)[C@H]1CCC2O |&1:4,5,14,16|", "rac-17ξ-hydroxy-5ξ,8α,9β,10α,13ξ,14β-androstan-3-one"),
        ("CC12CC[C@@H]3[C@H](CCC4CC(=O)CC[C@]43C)[C@H]1CCC2O |o1:4,5,14,16|", "rel-17ξ-hydroxy-5ξ,8α,9β,10α,13ξ,14β-androstan-3-one"),
    ],
)
def test_racemates_and_relative_configuration(smiles, expected):
    # P-101.8.2 and P-101.8.3
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@@H]1C[C@]23CC[C@H]4C(C)(C)CCC[C@]4(C)[C@H]2CC[C@H]1C3", "kaurane"),
        ("C[C@H]1C[C@@]23CC[C@@H]4C(C)(C)CCC[C@@]4(C)[C@@H]2CC[C@@H]1C3", "ent-kaurane"),
    ],
)
def test_inversion_of_all_chirality_centres(smiles, expected):
    # P-101.8.1
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, operations",
    [
        ("CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]2(C)CC1", 2),
        ("C1(Oc2ccccc2)OC[C@@H]3[C@H]1CO[C@H]3c4ccccc4", 2),
        ("C[C@@]12CCC[C@H]1[C@@H]1CC[C@H]3CC[Te]CC[C@]3(C)[C@H]1CC2", 2),
    ],
)
def test_operations_counted_for_preferred_semisystematic_names(smiles, operations):
    # P-101.3.7.2: modification prefixes and replacements are counted; two is the most for a preferred name
    from rdkit import Chem

    from smiles_to_iupac._np import name_natural_product_ranked

    assert name_natural_product_ranked(Chem.MolFromSmiles(smiles))[1] == operations


@pytest.mark.slow
def test_spiro_on_a_skeleton_that_is_not_a_ring_system_is_named_substitutively():
    smiles = "COc1ccc(-c2ccc(C[C@@]3(C(N)=O)CCCN(C(=O)c4cccnc4)C3)cc2)cc1"
    assert smiles_to_iupac(smiles) == "(3S)-3-[(4'-methoxy[1,1'-biphenyl]-4-yl)methyl]-1-(pyridine-3-carbonyl)piperidine-3-carboxamide"
