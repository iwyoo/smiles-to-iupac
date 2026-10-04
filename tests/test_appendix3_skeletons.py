"""Appendix 3 parent structures (P-101.2.7) and the compounds named on them (P-101.2.6, P-101.6, P-101.7)."""

import pytest
from rdkit import Chem

from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac
from smiles_to_iupac._appendix3_stereo_data import STEREO
from smiles_to_iupac._appendix3_table import SKELETONS
from smiles_to_iupac._common import UnsupportedStructure

YOHIMBAN = "C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2"


def _bare(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def _graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return sorted(
        tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()))) for b in mol.GetBonds()
    )


def test_every_atom_of_a_skeleton_has_its_own_locant():
    for name, smiles in SKELETONS.items():
        maps = [a.GetAtomMapNum() for a in Chem.MolFromSmiles(smiles).GetAtoms()]
        assert all(maps) and len(set(maps)) == len(maps), name


def test_stereo_parents_are_drawn_on_the_same_graph_as_their_skeleton():
    for name, (smiles, ref, _) in STEREO.items():
        assert _graph(smiles) == _graph(SKELETONS[name]), name
        assert ref, name


def test_every_bare_skeleton_is_named_by_its_retained_name():
    # a parent that has an implied configuration is given in it, which also tells dammarane from protostane
    wrong = {}
    for name, smiles in SKELETONS.items():
        got = smiles_to_iupac(_bare(STEREO[name][0] if name in STEREO else smiles))
        if got != name:
            wrong[name] = got
    assert wrong == {}


def test_retained_names_are_valid_but_not_preferred():
    with pytest.warns(NonPreferredNameWarning, match="P-101"):
        smiles_to_iupac(YOHIMBAN)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.6.1: ene endings in a saturated portion, 'enine' for -anine, aromatic ring with lowest locants
        ("CC(C)C1CC2=CCC3C(C)(C)CCCC3(C)C2CC1", "abiet-7-ene"),
        ("CC(C)C1=CC2=CCC3C(C)(C)CCCC3(C)C2CC1", "abieta-7,13-diene"),
        ("CC(C)c1ccc2c(c1)CCC1C(C)(C)CCCC21C", "abieta-8,11,13-triene"),
        ("C=C(C)C1CCC2(C)CCC3(C)C(CCC4C5(C)CCCC(C)(C)C5CCC43C)C12", "lup-20(29)-ene"),
        ("CC1C2CCC3C4CC=C5CCCCC5(C)C4CCC32CN1C", "con-5-enine"),
        ("C1=C2CCCCC2CN2CCc3c([nH]c4ccccc34)C12", "yohimb-14-ene"),
        # P-101.6.4 and P-101.6.6: hydro and dehydro prefixes on a parent that is not saturated
        ("CCC1CN2CCC1CC2Cc1ccnc2ccccc12", "10,11-dihydrocinchonan"),
        ("C1=C2c3cccc4[nH]cc(c34)CC2NCC1", "9,10-didehydroergoline"),
        ("c1ccc2c(c1)NC1C2CCN2CC3CCCCC3CC12", "2,7-dihydroyohimban"),
        ("CC(C=CC=C(C)C=CC1=C(C)CCCC1(C)C)=CC=CC=C(C)C=CC=C(C)CCC1=C(C)CCCC1(C)C", "7,8-dihydro-β,β-carotene"),
        # the skeleton that needs the fewest changes is the parent
        ("CCC12CCCN3CCc4c(n(c5ccccc45)CC1)C32", "vincane"),
        # P-31.1.4.3.4: a ketone on a doubly bonded atom takes added hydrogen
        ("C=CC1CN2CCC1CC2Cc1cc(=O)[nH]c2ccccc12", "cinchonan-2′(1′H)-one"),
        ("CC(C)c1ccc2c(c1)C(=O)CC1C(C)(C)CCCC21C", "abieta-8,11,13-trien-7-one"),
        # P-31.1.4.2.4: the indicated hydrogen of the parent moves to a ketone atom that has no partner
        (
            "CC1=C(C=C)C(=O)NC1=Cc1[nH]c(Cc2[nH]c(C=C3NC(=O)C(C=C)=C3C)c(C)c2C)c(C)c1C",
            "2,18-diethenyl-3,7,8,12,13,17-hexamethyl-10,21,22,23-tetrahydro-1H-biline-1,19(24H)-dione",
        ),
    ],
)
def test_degree_of_hydrogenation(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-ol"),
        ("O=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-one"),
        ("OC1CC(O)C2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-16,18-diol"),
        ("SC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-thiol"),
        ("CN(C)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N,N-dimethylyohimban-18-amine"),
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carboxylic acid"),
        ("COC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "methyl yohimban-18-carboxylate"),
        ("O=C(N(C)C)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N,N-dimethylyohimban-18-carboxamide"),
        ("N#CC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carbonitrile"),
        ("O=CC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carbaldehyde"),
        # P-41: the senior class is the suffix, the others are prefixes
        ("OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(=O)O", "18-hydroxyyohimban-14-carboxylic acid"),
        ("OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(=O)OC", "methyl 18-hydroxyyohimban-14-carboxylate"),
        ("O=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2N", "14-aminoyohimban-18-one"),
        ("COC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2Cl", "14-chloro-18-methoxyyohimban"),
        ("OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2c1ccccc1", "14-phenylyohimban-18-ol"),
        ("C=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "18-methylideneyohimban"),
        ("Cn1c2ccccc2c2CCN3CC4CCCCC4CC3c12", "1-methylyohimban"),
        ("Oc1ccc2c(c1)[nH]c1c2CCN2CC3CCCCC3CC12", "yohimban-11-ol"),
        # P-65.5.1, P-66.1, P-66.3, P-66.4, P-62.3.1, P-68.3.1.1.2: the other suffix classes of P-41 and their N locants
        ("CC(=O)OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "ethanoic yohimban-18-carboxylic anhydride"),
        ("ClC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carbonyl chloride"),
        ("NNC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carbohydrazide"),
        ("CCNNC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N′-ethylyohimban-18-carbohydrazide"),
        ("NC(=N)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-carboximidamide"),
        ("CNC(=NC)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N,N′-dimethylyohimban-18-carboximidamide"),
        ("CNS(=O)(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N-methylyohimban-18-sulfonamide"),
        ("N=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-imine"),
        ("ON=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N-hydroxyyohimban-18-imine"),
        ("CNC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2NC", "N14,N18-dimethylyohimban-14,18-diamine"),
        # P-65.5.4, P-66.3, P-66.4.1.3, P-65.3.2.3: prefixes of the same groups below a carboxylic acid
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(=O)Cl", "14-carbonochloridoylyohimban-18-carboxylic acid"),
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(=O)NN", "14-(hydrazinecarbonyl)yohimban-18-carboxylic acid"),
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(N)=N", "14-carbamimidoylyohimban-18-carboxylic acid"),
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2S(N)(=O)=O", "14-sulfamoylyohimban-18-carboxylic acid"),
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C=NO", "14-[(hydroxyimino)methyl]yohimban-18-carboxylic acid"),
        # P-65.6.3.3.2.1: unlike alkyl groups in alphanumerical order, identical ones grouped under their locants
        ("COC(=O)C1CC(C(=O)OCC)C2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "16-ethyl 18-methyl yohimban-16,18-dicarboxylate"),
        (
            "CCOC(=O)C1CC(C(=O)OCC)C2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C(=O)OC",
            "16,18-diethyl 14-methyl yohimban-14,16,18-tricarboxylate",
        ),
        # P-65.6.3: an O-acyl group makes the parent the alcohol part of an ester
        ("CC(=O)OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "yohimban-18-yl ethanoate"),
        ("CC(=O)OC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2OC(=O)CC", "yohimban-14,18-diyl 18-ethanoate 14-propanoate"),
        ("COC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2OC(C)=O", "methyl 14-(ethanoyloxy)yohimban-18-carboxylate"),
    ],
)
def test_characteristic_groups_on_the_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-44.1.1: a chain that carries more of the senior class than the parent becomes the parent, and the
        # Appendix 3 structure a substituent group (P-101.7.3)
        ("OCC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "(yohimban-18-yl)methanol"),
        ("OC(=O)CC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "2-(yohimban-18-yl)ethanoic acid"),
        ("CC(O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "1-(yohimban-18-yl)ethanol"),
        ("CC(C)Oc1ccc2c(c1)[nH]c1c2CCN2CC3CCCCC3CC12", "11-(propan-2-yloxy)yohimban"),
        # the parent carries as many acids as the chain does, so it stays the parent
        ("OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2CCC(=O)O", "14-(2-carboxyethyl)yohimban-18-carboxylic acid"),
        # an acid in the substituent of an amine group of the parent makes the chain the parent
        ("OC(=O)CCNC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "3-[(yohimban-18-yl)amino]propanoic acid"),
        # an oxime hydroxy group is no alcohol group, so the amine and oxime classes decide as before
        ("ON=C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2", "N-hydroxyyohimban-18-imine"),
    ],
)
def test_parent_selection_and_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


BILIRUBIN = "CC1=C(C(=O)N/C1=C\\C2=C(C(=C(N2)CC3=C(C(=C(N3)/C=C\\4/C(=C(C(=O)N4)C=C)C)C)CCC(=O)O)CCC(=O)O)C)C=C"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-15.3.2: identical chain parents on the multivalent group of the parent (P-101.7.3), numbered as the parent
        ("OC(=O)CCC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2CCC(=O)O", "3,3'-(yohimban-14,18-diyl)dipropanoic acid"),
        ("OCCC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2CCO", "2,2'-(yohimban-14,18-diyl)di(ethan-1-ol)"),
        ("NC(=O)CCC1CC2=CCC3C(C)(C)C(CCC(N)=O)CCC3(C)C2CC1", "3,3'-(podocarp-7-ene-3,13-diyl)dipropanamide"),
        # the other substituents of the parent stay in the group
        ("OC(=O)CCC1CCC2C(C1)CN1CCc3c([nH]c4ccc(C)cc34)C1C2CCC(=O)O", "3,3'-(10-methylyohimban-14,18-diyl)dipropanoic acid"),
        # P-101.6.1: ene ending of the group; P-101.6.4: hydro prefixes of the group (bilirubin); P-16.5.1: nesting
        ("OC(=O)CCC1CC2=CCC3C(C)(C)CCCC3(C)C2CC1CCC(=O)O", "3,3'-(podocarp-7-ene-12,13-diyl)dipropanoic acid"),
        (
            "C=CC1=C(C)C(=O)NC1=Cc1[nH]c(Cc2[nH]c(C=C3NC(=O)C(C=C)=C3C)c(C)c2CCC(=O)O)c(CCC(=O)O)c1C",
            "3,3'-(3,18-diethenyl-2,7,13,17-tetramethyl-1,19-dioxo-1,10,19,22,23,24-hexahydro-21H-biline-8,12-diyl)dipropanoic acid",
        ),
        # P-101.8.2: double bond configuration of the group; the descriptors make the enclosing marks square
        (
            BILIRUBIN,
            "3,3'-[(4Z,15Z)-2,18-diethenyl-3,7,13,17-tetramethyl-1,19-dioxo-1,10,19,22,23,24-hexahydro-21H-biline-8,12-diyl]"
            "dipropanoic acid",
        ),
        # P-101.2.6, P-101.7.1: the face of a free valence on a ring atom of the parent precedes its locant
        (
            "O=C(O)CCC1C[C@H](CCC(=O)O)[C@H]2C[C@H]3c4[nH]c5ccccc5c4CCN3C[C@@H]2C1",
            "3,3'-(yohimban-16α,18-diyl)dipropanoic acid",
        ),
        (
            "O=C(O)CCC1C[C@@H](CCC(=O)O)[C@H]2C[C@H]3c4[nH]c5ccccc5c4CCN3C[C@@H]2C1",
            "3,3'-(yohimban-16β,18-diyl)dipropanoic acid",
        ),
        ("O=C(O)CCC1CC(CCC(=O)O)[C@H]2C[C@H]3c4[nH]c5ccccc5c4CCN3C[C@@H]2C1", "3,3'-(yohimban-16,18-diyl)dipropanoic acid"),
        # P-44.1.1: the acids of the chains outrank the amine, ether and ester groups that join them to the parent
        ("OC(=O)CCNC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2NCCC(=O)O", "3,3'-[yohimban-14,18-diylbis(azanediyl)]dipropanoic acid"),
        ("OC(=O)CCOC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2OCCC(=O)O", "3,3'-[yohimban-14,18-diylbis(oxy)]dipropanoic acid"),
        # three identical chains
        (
            "OC(=O)CCC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2(CCC(O)=O)CCC(O)=O",
            "3,3',3''-(yohimban-14,14,18-triyl)tripropanoic acid",
        ),
    ],
)
def test_multiplicative_names_with_a_parent_as_the_central_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1C[C@H](CCC(=O)O)[C@H]2C[C@H]3c4[nH]c5ccccc5c4CCN3C[C@@H]2C1", "3-(yohimban-16α-yl)propanoic acid"),
        ("C1C[C@@H](CCC(=O)O)[C@H]2C[C@H]3c4[nH]c5ccccc5c4CCN3C[C@@H]2C1", "3-(yohimban-16β-yl)propanoic acid"),
    ],
)
def test_face_of_the_free_valence_of_a_substituent_group(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_arms_with_stereo_outside_the_group_are_not_multiplied():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)[C@@H](C)CC1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2C[C@H](C)C(O)=O")


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # P-101.2.6: the name implies the configuration of the drawing
        ("CC(C)[C@H]1CC[C@H]2[C@@H](CC[C@H]3C(C)(C)CCC[C@]23C)C1", "abietane"),
        ("CC[C@]12C=Cn3c4c(c5ccccc53)CCN(CCC1)[C@@H]42", "eburnamenine"),
        # P-101.2.6.1.1: the two worked examples of a changed centre
        ("CC(C)[C@@H]1CC[C@H]2[C@@H](CC[C@H]3C(C)(C)CCC[C@]23C)C1", "13β-abietane"),
        ("CC[C@]12C=Cn3c4c(c5ccccc53)CCN(CCC1)[C@H]42", "3α-eburnamenine"),
        # P-101.7.1: α/β before the locant of a substituent
        (
            "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C",
            "lup-20(29)-en-3β-ol",
        ),
        (
            "COC(=O)[C@H]1[C@@H](O)CC[C@H]2CN3CCc4c([nH]c5ccccc45)[C@@H]3C[C@H]12",
            "methyl 17α-hydroxyyohimban-16α-carboxylate",
        ),
        # P-101.2.6.1.3: R and S for a parent drawn in perspective
        (
            "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O",
            "(3R,4S,8S,9R)-6′-methoxycinchonan-9-ol",
        ),
        (
            "C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75",
            "(7R,8S,12S,13R,14R,16S)-strychnidin-10-one",
        ),
        ("CN1CC[C@]23CCCC[C@H]2[C@H]1Cc1ccc(OC)cc13", "(9R,13R,14R)-3-methoxy-17-methylmorphinan"),
        # P-91.2: centres of a side branch are cited inside the substituent prefix, cepham numbering as in P-101
        (
            "CC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)[C@@H](C3=CC=CC=C3)N)SC1)C(=O)O",
            "7β-{[(2R)-2-amino-2-phenylethanoyl]amino}-3-methyl-3,4-didehydrocepham-4-carboxylic acid",
        ),
        # P-101.2.6: α and β only for rings fused to the drawn plane; a centre the drawing leaves open (C-22 of
        # spirostan) or one beyond a spiro atom is cited with R and S
        ("CC[C@@H]1CN2CC[C@@]3(CNc4ccccc43)[C@H]2C[C@@H]1CC", "3β-corynoxan"),
        (
            "C[C@@H]1CC[C@@]2([C@H]([C@H]3[C@@H](O2)C[C@@H]4[C@@]3(CC[C@H]5[C@H]4CC=C6[C@@]5(CC[C@@H](C6)O)C)C)C)OC1",
            "(22R,25R)-spirost-5-en-3β-ol",
        ),
        # P-101.6.2: the geometry of a skeleton double bond is implied by the parent
        ("C/C=C1\\CO[C@@H]2CCN3CC=C(COC[C@@H](C)[C@H](C)C1)[C@H]23", "senecionan"),
        ("C/C=C1/CO[C@@H]2CCN3CC=C(COC[C@@H](C)[C@H](C)C1)[C@H]23", "(15E)-senecionan"),
        # P-101.6.2, P-101.6.3: cis and trans in retinoids and carotenoids
        ("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C", "all-trans-retinal"),
        ("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C\\C(=C\\C=O)\\C)/C", "11-cis-retinal"),
        (
            "CC1=C(C(C[C@@H](C1)O)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C2=C(C[C@H](CC2(C)C)O)C)\\C)\\C)/C)/C",
            "all-trans-(3R,3′R)-β,β-carotene-3,3′-diol",
        ),
    ],
)
def test_configuration(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        # a further ring fused to the parent is a different ring system
        "c1ccc2c(c1)ccc1c2CCN2CC3CCCCC3CC12",
        # a stereocentre in a side branch that no substituent name cites (an acyl group of a retained amino acid)
        "OC(=O)C1CCC2C(C1)CN1CCc3c([nH]c4ccccc34)C1C2NC(=O)[C@H]1CCCN1",
    ],
)
def test_structures_the_parent_cannot_carry_are_left_to_the_other_engines(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC1CCC(C(C)C)CC1", "1-methyl-4-(propan-2-yl)cyclohexane"),
        ("CC1(C)C2CCC1(C)C(=O)C2", "1,7,7-trimethylbicyclo[2.2.1]heptan-2-one"),
    ],
)
def test_monoterpenes_keep_their_systematic_names(smiles, expected):
    # P-13.4.3.2 gives the bicyclo[2.2.1]heptane name as the PIN for bornane-based structures
    assert smiles_to_iupac(smiles) == expected
