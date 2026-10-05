import pytest
from fractions import Fraction
from rdkit import Chem
from rdkit.Chem import Atom, BondType, RWMol
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure, adjacency, ring_cycle
from smiles_to_iupac._pyrimidinedione import has_pyrimidinedione_shape


def test_benzo_d_aceanthrylene():
    assert smiles_to_iupac("C1=Cc2c3ccccc3cc3cc4ccccc4c1c23") == "cyclopenta[fg]tetracene"


def test_benzo_a_acephenanthrylene():
    assert smiles_to_iupac("C1=Cc2cc3ccccc3c3c2c1cc1ccccc13") == "benzo[a]acephenanthrylene"


def test_anthanthrene():
    assert (
        smiles_to_iupac("C1=CC2=C3C(=C1)C=C4C=CC5=C6C4=C3C(=CC6=CC=C5)C=C2")
        == "naphtho[7,8,1,2,3-nopqr]tetraphene"
    )


def test_1h_cyclopenta_a_anthracene():
    assert smiles_to_iupac("C1C=CC2=C1C3=CC4=CC=CC=C4C=C3C=C2") == "1H-cyclopenta[a]anthracene"


def test_benzo_a_anthracene():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32") == "tetraphene"


def test_benzo_a_azulene():
    assert smiles_to_iupac("C1=CC=C2C=C3C=CC=CC3=C2C=C1") == "benzo[a]azulene"


def test_azulene():
    assert smiles_to_iupac("C1=CC2=CC=CC=CC2=C1") == "azulene"


def test_benzo_cd_indole():
    assert smiles_to_iupac("C1=CC2=C3C(=C1)C=NC3=CC=C2") == "benzo[cd]indole"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2[nH]cnc2c1", "1H-1,3-benzimidazole"),
        ("c1ccc2c(c1)OCCO2", "2,3-dihydro-1,4-benzodioxine"),
        ("c1ccc2[nH]ncc2c1", "1H-indazole"),
        ("O1C=CC=Cc2ccccc12", "1-benzoxepine"),
    ],
)
def test_benzo_heterocycle_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_triphenylene_is_named():
    assert smiles_to_iupac("Cc1ccc2c(c1)c1ccccc1c1ccccc21") == "2-methyltriphenylene"


def test_dihydronaphthalene_itself_is_unaffected():
    assert smiles_to_iupac("C1CC=Cc2ccccc12") == "1,2-dihydronaphthalene"


def test_fully_saturated_bridge_is_von_baeyer_not_fusion_name():
    assert smiles_to_iupac("C1CC2CC1c1ccccc12") == "1,2,3,4-tetrahydro-1,4-methanonaphthalene"


def test_substituted_epoxy_aromatic_ring_is_named():
    assert smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2O1") == "6-methyl-1,4-dihydro-1,4-epoxynaphthalene"


@pytest.mark.parametrize(
    "smiles",
    [
        "c12ccccc1C1c3ccccc3C2C1",
    ],
)
def test_bridged_anthracene(smiles):
    assert smiles_to_iupac(smiles) == "9,10-dihydro-9,10-methanoanthracene"


def test_anthracene_itself_is_unaffected():
    assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1") == "anthracene"


def test_halogen_on_bridged_anthracene_terminal_ring():
    assert (
        smiles_to_iupac("Clc1ccc2cc3c(cc2c1)C1C=CC3C1")
        == "6-chloro-1,4-dihydro-1,4-methanoanthracene"
    )


def test_bridged_phenanthrene():
    assert smiles_to_iupac("C1=CC2CC1c1ccc3ccccc3c12") == "1,4-dihydro-1,4-methanophenanthrene"


def test_bridged_tetracene():
    assert smiles_to_iupac("C1=CC2CC1c1cc3cc4ccccc4cc3cc12") == "1,4-dihydro-1,4-methanotetracene"


def test_substituted_ethano_aromatic_ring_is_named():
    assert smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2CC1") == "6-methyl-1,4-dihydro-1,4-ethanonaphthalene"


@pytest.mark.parametrize(
    "smiles",
    [
        "c1ccc2c(c1)C1c3ccccc3C2c2ccccc21",
    ],
)
def test_bridged_anthracene_benzo(smiles):
    assert smiles_to_iupac(smiles) == "9,10-dihydro-9,10-[1,2]benzenoanthracene"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)CN1C=CN2C=CC=CC12", "2-[imidazo[1,2-a]pyridin-1(8aH)-yl]ethanoic acid"),
    ],
)
def test_bridgehead_heteroatom_fused_yl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_chrysene_fusion_is_named():
    assert smiles_to_iupac("Cc1ccc2cc3c(ccc4c5ccccc5ccc34)cc2c1") == "9-methylbenzo[c]tetraphene"


def test_cyclopenta_a_naphthalene_1h():
    assert smiles_to_iupac("C1C=CC2=C1C3=CC=CC=C3C=C2") == "1H-cyclopenta[a]naphthalene"


def test_cyclopenta_b_naphthalene_1h():
    assert smiles_to_iupac("C1C=CC2=CC3=CC=CC=C3C=C21") == "1H-cyclopenta[b]naphthalene"


def test_didehydropiperidine_nitrogen_ring():
    assert smiles_to_iupac("C1CC=CCN1") == "3,4-didehydropiperidine"


def test_benzo_a_fluoranthene():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC=C4C3=C2C5=CC=CC=C54") == "benzo[a]aceanthrylene"


def test_benzo_a_fluorene():
    assert smiles_to_iupac("C1C2=CC=CC=C2C3=C1C4=CC=CC=C4C=C3") == "11H-benzo[a]fluorene"


def test_fluorene():
    assert smiles_to_iupac("C1C2=CC=CC=C2C3=CC=CC=C31") == "9H-fluorene"


def test_furo_3_2_b_pyran():
    assert smiles_to_iupac("C1=COC2=CCOC2=C1") == "2H-furo[3,2-b]pyran"


def _furanobenzoquinoline_smiles():
    base = Chem.MolFromSmiles("C1=CC=C2C=C3C(=CC2=C1)C=CC=N3")
    furan = Chem.MolFromSmiles("c1ccoc1")
    combo = Chem.CombineMols(base, furan)
    rw = RWMol(combo)
    n_base = base.GetNumAtoms()
    rw.AddBond(4, n_base + 4, BondType.SINGLE)
    rw.AddBond(7, n_base + 0, BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def test_furano_bridge_benzo_g_quinoline():
    assert smiles_to_iupac(_furanobenzoquinoline_smiles()) == "10,5-[2,3]furanobenzo[g]quinoline"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2cs(=O)cc2c1", "2-benzothiophene 2-oxide"),
    ],
)
def test_fused_hetero_ring_oxide_resolves(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Oc1ccc2ccc3cccc4ccc1c2c34", "pyren-1-ol"),
        ("Oc1cc2c3c(N)cccc3cc3ccc4cccc1c4c32", "10-aminobenzo[pqr]tetraphen-12-ol"),
    ],
)
def test_substituted_larger_fused_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CN(C)CCC=C1c2ccccc2CCc2ccccc12",
            "3-(10,11-dihydro-5H-dibenzo[a,d][7]annulen-5-ylidene)-N,N-dimethylpropan-1-amine",
        ),
    ],
)
def test_seven_membered_fused_parents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("NC(=O)N1CCCC1", "pyrrolidine-1-carboxamide"),
        ("OC(=O)N1CCOCC1", "morpholine-4-carboxylic acid"),
        ("NC(=O)N1c2ccccc2Sc2ccccc21", "10H-phenothiazine-10-carboxamide"),
    ],
)
def test_acyl_groups_on_a_ring_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


_PHENANTHRENE_REF = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
_PHEN_1_METHYL = Chem.MolFromSmiles("CC1=C2C=CC3=CC=CC=C3C2=CC=C1")
_PHEN_9_METHYL = Chem.MolFromSmiles("CC1=CC2=CC=CC=C2C3=CC=CC=C13")
_PHEN_4_METHYL = Chem.MolFromSmiles("CC1=C2C(=CC=C1)C=CC3=CC=CC=C32")

_TET_1_METHYL = Chem.MolFromSmiles("CC1=CC=CC2=CC3=CC4=CC=CC=C4C=C3C=C12")
_TET_2_METHYL = Chem.MolFromSmiles("CC1=CC2=CC3=CC4=CC=CC=C4C=C3C=C2C=C1")
_TET_5_METHYL = Chem.MolFromSmiles("CC1=C2C=CC=CC2=CC3=CC4=CC=CC=C4C=C13")








_CHRYSENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")
_BENZO_C_PHENANTHRENE = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C=C3")






def _sp3_locant(mol, locants):
    (sp3,) = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 2
    ]
    return locants[sp3]












_PUBCHEM_SAMPLE = [
    (7080, "naphthacene", "C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1", 4),
    (8671, "pentacene", "C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1", 5),
    (9162, "picene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC=CC=C54", 4),
    (519935, "pentaphene", "C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC=CC=C5C=C43", 3),
    (9135, "benzo[c]chrysene", "C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C5=CC=CC=C5C=C4", 4),
    (9140, "benzo[g]chrysene", "C1=CC=C2C(=C1)C=CC3=C2C4=CC=CC=C4C5=CC=CC=C35", 4),
    (9164, "dibenzo[a,c]anthracene", "C1=CC=C2C=C3C4=CC=CC=C4C5=CC=CC=C5C3=CC2=C1", 3),
    (5889, "dibenz[a,h]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=C(C=CC5=CC=CC=C54)C=C32", 3),
    (9176, "dibenz[a,j]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=C(C=C32)C5=CC=CC=C5C=C4", 4),
    (5954, "benzo[a]anthracene", "C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32", 3),
    (9167, "dibenzo[a,c]naphthacene", "C1=CC=C2C=C3C=C4C5=CC=CC=C5C6=CC=CC=C6C4=CC3=CC2=C1", 4),
    (98863, "hexahelicene", "C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=C4C6=CC=CC=C6C=C5", 6),
]








def _synthetic_chain(sizes, edge_pair_at_middle):
    """A 3-ring linear chain ring0-ring1-ring2 where ring1 (the middle
    ring, of size `sizes[1]`) connects to ring0 via edge index
    `edge_pair_at_middle[0]` and to ring2 via `edge_pair_at_middle[1]` (in
    ring1's own cyclic bond order) -- ring0/ring2's own edge indices don't
    affect any result here (a leaf ring's `ref_edge` is never used to
    compute anything further)."""
    adj = {0: {1}, 1: {0, 2}, 2: {1}}
    ring_sizes = {0: sizes[0], 1: sizes[1], 2: sizes[2]}
    edge_index_of = {
        (0, 1): 0,
        (1, 0): edge_pair_at_middle[0],
        (1, 2): edge_pair_at_middle[1],
        (2, 1): 0,
    }
    return adj, edge_index_of, ring_sizes, 3










def test_benzo_a_heptacene():
    assert (
        smiles_to_iupac("c1ccc2cc3cc4cc5cc6cc7c(ccc8ccccc87)cc6cc5cc4cc3cc2c1")
        == "benzo[a]heptacene"
    )


def test_3h_indole_tautomer_is_distinct_from_1h_indole():
    assert smiles_to_iupac("C1(C=Nc2ccccc12)") == "3H-indole"




def test_benzo_a_hexacene():
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6c(ccc7ccccc76)cc5cc4cc3cc2c1") == "benzo[a]hexacene"


def test_plain_hexacene_still_resolves():
    assert smiles_to_iupac("c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1") == "hexacene"


def test_substituent_on_a_naphthalene_amplificant_gets_a_composite_locant():
    assert (
        smiles_to_iupac("Cc1cc2ccc3cc2cc1CCc1ccc(cc1)CCC3")
        == "13-methyl-1(2,7)-naphthalena-4(1,4)-benzenacycloheptaphane"
    )


def test_pyridine_benzene_equal_bridge_phane():
    assert (
        smiles_to_iupac("c1cc2ccc1CCc1ccc(nc1)CC2")
        == "1(2,5)-pyridina-4(1,4)-benzenacyclohexaphane"
    )


def test_pyridine_benzene_unequal_bridge_phane():
    assert (
        smiles_to_iupac("c1cc2ccc1CCCc1ccc(nc1)CC2")
        == "1(2,5)-pyridina-4(1,4)-benzenacycloheptaphane"
    )


def test_pyridinium_nitrogen_attachment_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1Cc2ccc(cc2)CC[n+]2ccccc21")


def test_benzo_a_pentacene():
    assert smiles_to_iupac("c1ccc2cc3cc4cc5c(ccc6ccccc65)cc4cc3cc2c1") == "benzo[a]pentacene"


def test_hexaphene():
    assert smiles_to_iupac("C1=CC=C2C=C3C(=CC2=C1)C=CC4=CC5=CC6=CC=CC=C6C=C5C=C43") == "hexaphene"


def test_substituted_acenaphthylene_is_named():
    assert smiles_to_iupac("Cc1ccc2cccc3C=Cc1c23") == "3-methylacenaphthylene"


def test_two_peri_fused_cyclopenta_rings_on_perylene():
    assert smiles_to_iupac("c1cc2ccc3cc4ccc5ccc6cc1c1c2c3c4c5c61") == "dicyclopenta[ghi,pqr]perylene"


def test_benzo_a_perylene():
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=C3C=CC=C4C3=C2C5=CC=CC6=C5C4=CC=C6") == "benzo[a]perylene"
    )


def test_benzo_ghi_perylene():
    assert (
        smiles_to_iupac("C1=CC2=C3C(=C1)C4=CC=CC5=C4C6=C(C=C5)C=CC(=C36)C=C2")
        == "benzo[ghi]perylene"
    )


def test_chrysene():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43") == "chrysene"


def test_4_7_phenanthroline():
    assert smiles_to_iupac("C1=CC2=C(C=CC3=C2C=CC=N3)N=C1") == "4,7-phenanthroline"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccc2nnccc2c1", "cinnoline"),
    ],
)
def test_diazanaphthalenes_with_both_nitrogens_in_one_ring_keep_their_retained_names(
    smiles, expected
):
    assert smiles_to_iupac(smiles) == expected


def test_benzo_b_picene():
    assert (
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=C3C=CC5=CC6=CC=CC=C6C=C54") == "benzo[b]picene"
    )


def test_benzo_e_1_benzofuran():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CO3") == "naphtho[2,1-b]furan"


def _fuse_new_ring(base_smiles, atom_a, atom_b):
    mol = Chem.RWMol(Chem.MolFromSmiles(base_smiles))
    new_idxs = []
    for _ in range(4):
        atom = Chem.Atom(6)
        atom.SetIsAromatic(True)
        new_idxs.append(mol.AddAtom(atom))
    mol.AddBond(atom_a, new_idxs[0], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[0], new_idxs[1], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[1], new_idxs[2], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[2], new_idxs[3], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[3], atom_b, Chem.BondType.AROMATIC)
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return Chem.MolToSmiles(result)


def test_substituted_variant_names_with_locant():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2NC(C)=C3") == "2-methyl-1H-benzo[g]indole"


def _epipyranobenzoquinoline_smiles():
    base = Chem.MolFromSmiles("c1ccc2cc3ncccc3cc2c1")
    rw = RWMol(base)
    o1 = rw.AddAtom(Atom(8))
    c2 = rw.AddAtom(Atom(6))
    c3 = rw.AddAtom(Atom(6))
    c4 = rw.AddAtom(Atom(6))
    c5 = rw.AddAtom(Atom(6))
    c6 = rw.AddAtom(Atom(6))
    rw.AddBond(o1, c2, BondType.SINGLE)
    rw.AddBond(c2, c3, BondType.DOUBLE)
    rw.AddBond(c3, c4, BondType.SINGLE)
    rw.AddBond(c4, c5, BondType.DOUBLE)
    rw.AddBond(c5, c6, BondType.SINGLE)
    rw.AddBond(c6, o1, BondType.SINGLE)
    rw.AddBond(c2, 11, BondType.SINGLE)
    rw.AddBond(c5, 4, BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def test_epipyrano_bridge_benzo_g_quinoline():
    assert (
        smiles_to_iupac(_epipyranobenzoquinoline_smiles())
        == "12H-5,10-[2,5]epipyranobenzo[g]quinoline"
    )


def test_benzo_a_pyrene():
    assert smiles_to_iupac("C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC=C5)C=C3") == "benzo[pqr]tetraphene"


def test_benzo_f_isoquinoline():
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CN=C3") == "benzo[f]isoquinoline"








def test_non_aromatic_ketone_unaffected():
    assert smiles_to_iupac("O=C1C=CC=CC1") == "cyclohexa-2,4-dien-1-one"


def test_barbituric_acid_names_as_saturated_trione():
    assert smiles_to_iupac("O=C1CC(=O)NC(=O)N1") == "1,3-diazinane-2,4,6-trione"


def test_non_1_3_diazine_dione_not_matched():
    mol = Chem.MolFromSmiles("O=C1C=CC(=O)NN1")
    assert not has_pyrimidinedione_shape(mol)


def test_plain_pyrimidine_unaffected():
    assert smiles_to_iupac("c1ccncn1") == "pyrimidine"


def test_pyrrolo_ij_quinoline():
    assert smiles_to_iupac("C1C=CC2=CC=CC3=C2N1C=C3") == "4H-pyrrolo[3,2,1-ij]quinoline"








def test_benzo_a_tetracene():
    assert smiles_to_iupac("c1ccc2cc3cc4c(ccc5ccccc54)cc3cc2c1") == "benzo[a]tetracene"


def test_benzo_b_triphenylene():
    assert smiles_to_iupac("c1ccc2cc3c4ccccc4c4ccccc4c3cc2c1") == "benzo[f]tetraphene"


def test_thieno_3_4_b_thiophene():
    assert smiles_to_iupac("C1=CSC2=CSC=C21") == "thieno[3,4-b]thiophene"


def test_thieno_3_2_b_furan():
    assert smiles_to_iupac("C1=COC2=C1SC=C2") == "thieno[3,2-b]furan"


def test_furo_2_3_b_pyrrole():
    assert smiles_to_iupac("C1=CNC2=C1C=CO2") == "6H-furo[2,3-b]pyrrole"


def test_selenopheno_furan():
    assert smiles_to_iupac("C1=C[Se]C2=C1C=CO2") == "selenopheno[2,3-b]furan"


@pytest.mark.parametrize(
    "smiles, indicated_hydrogen_locant",
    [
        ("c1ccc2c(c1)Cc1ccccc12", "9"),
        ("c1ccc2c(c1)Cc1cc3ccccc3cc12", "11"),
        ("c1ccc2c(c1)Cc1ccc3ccccc3c12", "7"),
        ("c1ccc2c(c1)Cc1c2ccc2ccccc12", "11"),
        ("C1C=Cc2ccccc12", "1"),
        ("C1C=CC=Cc2ccccc12", "5"),
        ("C1C=Cc2ccccc2-c2ccccc12", "5"),
    ],
)
def test_numbering_of_five_six_and_seven_membered_fused_systems(smiles, indicated_hydrogen_locant):
    from smiles_to_iupac._fused_numbering import fused_numberings

    mol = Chem.MolFromSmiles(smiles)
    locants = {
        numbering[atom.GetIdx()]
        for numbering in fused_numberings(mol)
        for atom in mol.GetAtoms()
        if atom.GetTotalNumHs() == 2
    }
    assert min(locants, key=int) == indicated_hydrogen_locant


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)Cc1ccc2c(c1)Cc1cc3ccccc3cc12", "2-(11H-benzo[b]fluoren-2-yl)ethanoic acid"),
        ("OC(=O)Cc1ccc2c(c1)Cc1ccc3ccccc3c12", "2-(7H-benzo[c]fluoren-9-yl)ethanoic acid"),
        ("OC(=O)Cc1ccc2c(c1)Cc1c2ccc2ccccc12", "2-(11H-benzo[a]fluoren-9-yl)ethanoic acid"),
    ],
)
def test_yl_groups_of_fused_systems_with_a_five_membered_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)Cc1cn2ccccc2n1", "2-(imidazo[1,2-a]pyridin-2-yl)ethanoic acid"),
        ("OC(=O)CN1CC=NC2=NC=CN12", "2-[imidazo[1,2-b][1,2,4]triazin-1(2H)-yl]ethanoic acid"),
    ],
)
def test_bridgehead_heteroatom_fused_yl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1C=CC=Cc2ccccc12", "5H-benzo[7]annulene"),
        ("C1C=Cc2ccccc2-c2ccccc12", "5H-dibenzo[a,c][7]annulene"),
        ("C1C=Cc2c1ccc1c2ccc2ccccc21", "17H-cyclopenta[a]phenanthrene"),
    ],
)
def test_fusion_names_of_a_parent_component_with_attached_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)CC1C=CC=Cc2ccccc12", "2-(5H-benzo[7]annulen-5-yl)ethanoic acid"),
        ("OC(=O)CC1C=Cc2ccccc2-c2ccccc12", "2-(5H-dibenzo[a,c][7]annulen-5-yl)ethanoic acid"),
        ("OC(=O)CC1C=CC2=C1C=CC1=C2C=Cc2ccccc21", "2-(17H-cyclopenta[a]phenanthren-17-yl)ethanoic acid"),
    ],
)
def test_yl_groups_of_parent_component_fusion_systems(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("c1cc2ccc1CCc1ccc3ccccc3c1CC2", "1(1,2)-naphthalena-4(1,4)-benzenacyclohexaphane"),
        ("C1CCc2ccc3cc(ccc3c2)CCc2ccc(cc2)C1", "1(2,6)-naphthalena-4(1,4)-benzenacyclooctaphane"),
        (
            "OC(=O)Cc1cc2ccc1CCc1ccc3ccccc3c1CC2",
            "2-[1(1,2)-naphthalena-4(1,4)-benzenacyclohexaphan-42-yl]ethanoic acid",
        ),
        (
            "OC(=O)CC1Cc2ccc(cc2)CCc2ccc3ccccc3c2C1",
            "2-[1(1,2)-naphthalena-4(1,4)-benzenacycloheptaphan-6-yl]ethanoic acid",
        ),
    ],
)
def test_phane_with_a_naphthalene_amplificant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1cc2c3cc4c5cocc5c4cc3c2o1", "benzo[1'',2'':3,4;4'',5'':3',4']dicyclobuta[1,2-b:1',2'-c']difuran"),
        ("C1=CC2=CC=NC3=NC=CC(=N1)N23", "1,3a1,4,9-tetraazaphenalene"),
        ("C1=CC2=C3C=CC=CN3C=CN2C=C1", "dipyrido[1,2-a:2',1'-c]pyrazine"),
        ("c1cnc2cc3cc4cnoc4cc3cc2c1", "[1,2]benzoxazolo[6,5-g]quinoline"),
        ("C1=Cc2cc3c(cc2=Cc2ccccc21)C=c1ccccc1=c1ccccc1=3", "tribenzo[c,d',e]benzo[1,2-a:4,5-a']di[7]annulene"),
        ("C1=Cc2cc3c(cc2=C1)-c1cc2c(nc1C=3)C=c1c-2ccc2c1=CC=C2", "as-indaceno[2,3-b]-s-indaceno[1,2-e]pyridine"),
        ("C1=S=CC2=C1C=S=C2", "2λ4δ2,5λ4δ2-thieno[3,4-c]thiophene"),
        ("C1=CSC23OC=CSC2=CC=C3O1", "cyclopenta[1,2-b:5,1-b']bis([1,4]oxathiine)"),
        (
            "C1=CC2OC1C1=C2C2C3=C(C4C=CC3O4)C1C1=C2C2C=CC1O2",
            "1,4,5,8,9,10,13,16-octahydro-13,16-epoxy-9,10-[1,2]benzeno-1,4:5,8-diepoxyanthracene",
        ),
        ("c1cc2c3c(cccc3c1)C13c4cccc5cccc(c45)C21c1cccc2cccc3c12", "6b,12b-[1,8]naphthalenoacenaphthyleno[1,2-a]acenaphthylene"),
        ("C1C=C2c3ccccc3C1c1ccccc12", "9H-9,10-(epiethanylylidene)anthracene"),
        ("C1%10c2cc3ccccc3cc2C(c2cc3ccccc3cc12)Cc1ccccc1C%10", "6,13-dihydro-6,13-(methano[1,2]benzenomethano)pentacene"),
        ("C1=CC=C2C(=C1)C3C4=CC=CC=C4C2[SH2]3", "9,10-dihydro-9,10-λ4-sulfanoanthracene"),
        ("CC1=CC=CC=C2C(C)=CC=CC=C12", "1,6-dimethyl-Δ1-heptalene"),
    ],
)
def test_fusion_engine_decompositions_and_bridges(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
