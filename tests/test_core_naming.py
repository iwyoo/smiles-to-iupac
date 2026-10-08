import pathlib
import pytest
import re
import tokenize
from rdkit import Chem
from rdkit.Chem import RWMol
from smiles_to_iupac import NonPreferredNameWarning, __version__, smiles_to_iupac
from smiles_to_iupac.__main__ import main
from smiles_to_iupac._common import UnsupportedStructure, adjacency
from smiles_to_iupac._numerals import numerical_term
from smiles_to_iupac._parent_hydride_stripping import strip_substituents
from smiles_to_iupac._substituents import name_branch


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C[C@H](Cl)CC", "(2S)-2-chlorobutane"),
    ],
)
def test_acyclic_alkane_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("C1(=CC=CC=C1)/C=C/C(=O)C1=CC=CC=C1", "chalcone", id="retained_name"),
        pytest.param("OC1=C(C(/C=C/C2=CC(=CC=C2)OC)=O)C=CC(=C1OC)O", "2′,4′-dihydroxy-3,3′-dimethoxychalcone", id="primed_locants_on_the_benzoyl_ring"),
        pytest.param("O=C(/C=C/c1ccc(Cl)cc1)c1ccccc1", "4-chlorochalcone", id="unprimed_locants_on_the_styryl_ring"),
        pytest.param("O=C(/C=C/c1ccccc1Br)c1cccc([N+](=O)[O-])c1", "2-bromo-3′-nitrochalcone", id="nitro_group"),
        pytest.param("OC1=C(C=CC(=C1)O)C(/C=C/C1=CC=C(C(=O)N)C=C1)=O", "4-[(1E)-3-(2,4-dihydroxyphenyl)-3-oxoprop-1-en-1-yl]benzamide", id="senior_group_keeps_the_substitutive_name"),
        pytest.param("C1(=CC=CC=C1)/C=C\\C(=O)C1=CC=CC=C1", "(2Z)-1,3-diphenylprop-2-en-1-one", id="only_the_e_isomer_is_chalcone"),
    ],
)
def test_chalcone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("ClC=[C@AL1]=CCl", "(1M)-1,3-dichloropropa-1,2-diene", id="allene_anticlockwise"),
        pytest.param("ClC=[C@AL2]=CCl", "(1P)-1,3-dichloropropa-1,2-diene", id="allene_clockwise_is_the_enantiomer"),
        pytest.param("NC(Br)=[C@AL1]=C(F)Cl", "(1P)-1-bromo-3-chloro-3-fluoropropa-1,2-dien-1-amine", id="priority_order_swaps_the_sense"),
        pytest.param("OC(=O)C=[C@AL2]=C=C=CCl", "(2P)-6-chlorohexa-2,3,4,5-tetraenoic acid", id="even_cumulene_with_four_double_bonds"),
        pytest.param("O[C@H](C)C=[C@AL1]=CCl", "(2R,3M)-5-chloropenta-3,4-dien-2-ol", id="axis_merged_into_the_stereo_set"),
        pytest.param("C/C=C=C=C/C", "(2E)-hexa-2,3,4-triene", id="odd_cumulene_trans"),
        pytest.param("C/C=C=C=C\\C", "(2Z)-hexa-2,3,4-triene", id="odd_cumulene_cis"),
        pytest.param("CC=C=CC", "penta-2,3-diene", id="unspecified_allene_stays_bare"),
        pytest.param("Cl[C@H]1C[C@]2(C1)C[C@@H](Cl)C2", "(2R,4S,6R)-2,6-dichlorospiro[3.3]heptane", id="spiro_atom_descriptor"),
    ],
)
def test_axial_and_spiro_stereodescriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("CC(C)(C)C1(C)CCCC(=O)C1", "3-tert-butyl-3-methylcyclohexan-1-one", id="tert_butyl_cited_before_methyl_on_ring"),
        pytest.param("Cc1ccc2c(c1)c1ccccc1c1ccccc21", "2-methyltriphenylene", id="substituted_branched_fusion_is_named"),
        pytest.param("c1nc[nH]n1", "1H-1,2,4-triazole", id="heteroaromatic"),
    ],
)
def test_tert_butyl_cited_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_saturated_rings_still_resolve_unaffected():
    assert smiles_to_iupac("C1CCCCC1") == "cyclohexane"
    assert smiles_to_iupac("C1CC2CCC1C2") == "bicyclo[2.2.1]heptane"
    assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "adamantane"
    assert smiles_to_iupac("C1C2C3C2C4C1C34") == "tetracyclo[3.2.0.0^2,7.0^4,6]heptane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C#Cc1ccccc1", "ethynylbenzene"),
    ],
)
def test_exocyclic_unsaturated_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1ccccc1C[C@@H](Cl)C", "[(2S)-2-chloropropyl]benzene"),
        ("Cl[C@H](CC)c1ccc2ccccc2c1", "2-[(1R)-1-chloropropyl]naphthalene"),
    ],
)
def test_substituent_branch_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("Clc1ccccc1[C@@H](Cl)CC", "1-chloro-2-[(1S)-1-chloropropyl]benzene", id="stereo_marker_no_longer_silently_dropped_without_scope_is_named"),
        pytest.param("c1ccccc1[C@@H](Cl)[C@@H](Cl)C", "[(1R,2S)-1,2-dichloropropyl]benzene", id="two_stereocenters_on_one_substituent"),
    ],
)
def test_stereo_marker_no_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
def test_stereocenter_on_three_ring_fused_substituent():
    assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1[C@@H](Cl)CC") == "1-[(1S)-1-chloropropyl]anthracene"


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    captured = capsys.readouterr()
    assert exc_info.value.code == 0
    assert captured.out.strip() == __version__


def test_continues_after_error(capsys):
    exit_code = main(["not_a_smiles", "CC"])
    captured = capsys.readouterr()
    assert exit_code != 0
    assert "CC\tethane" in captured.out
    assert "not_a_smiles" in captured.err


P = "P(c1ccccc1)(c1ccccc1)c1ccccc1"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[NH3][Pt@SP1](Cl)(Cl)[NH3]", "(SP-4-2)-diamminedichloridoplatinum"),
        ("[NH3][Pt@SP1](Cl)([NH3])Cl", "(SP-4-1)-diamminedichloridoplatinum"),
        ("[NH3][Co@OH1](Cl)([NH3])(Cl)([NH3])Cl", "(OC-6-21)-triamminetrichloridocobalt"),
        ("[NH3][Co@OH4](Cl)([NH3])(Cl)([NH3])Cl", "(OC-6-22)-triamminetrichloridocobalt"),
        (
            f"[O+]#[C][Fe@TB17]([C]#[O+])([C]#[O+])({P}){P}",
            "(TBPY-5-11)-tricarbonylbis(triphenylphosphane)iron",
        ),
        (
            "[O+]#[C][Fe@TB1]([C]#[O+])(I)(" + P + ")Cl",
            "(TBPY-5-24-C)-dicarbonylchloridoiodido(triphenylphosphane)iron",
        ),
        ("I[Fe@](Cl)(Br)F", "(T-4-R)-bromidochloridofluoridoiodidoiron"),
        ("I[Fe@@](Cl)(Br)F", "(T-4-S)-bromidochloridofluoridoiodidoiron"),
    ],
)
def test_configuration_index_and_chirality(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def _co_smiles(cls):
    return f"[Co@OH{cls}]%10%11%12%13%14%15.Br%10.Br%11.[NH2]%12CC[NH2]%13.[NH3]%14.[NH3]%15"


@pytest.mark.parametrize(
    "cls,expected",
    [
        (12, "(OC-6-22)-diamminedibromido(ethane-1,2-diamine-κ2N,N')cobalt"),
    ],
)
def test_octahedral_chelate_configuration_index(cls, expected):
    assert smiles_to_iupac(_co_smiles(cls)) == expected


def _bis_tridentate(cls):
    a = "[NH2]%10CC[NH]%11CC[NH2]%12"
    b = "[NH2]%13CC[NH]%14CC[NH2]%15"
    return f"[Co@OH{cls}]%10%11%12%13%14%15.{a}.{b}"


@pytest.mark.parametrize(
    "cls,prefix",
    [(4, "(OC-6-1′2′)-"), (19, "(OC-6-1′2)-"), (1, "(OC-6-2′2-C)-"), (2, "(OC-6-2′2-A)-")],
)
def test_bis_tridentate_uses_the_priming_convention(cls, prefix):
    assert smiles_to_iupac(_bis_tridentate(cls)).startswith(prefix)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OCCc4ccc(c7ccc8ccccc8c7)cc4", "2-[4-(naphthalen-2-yl)phenyl]ethan-1-ol"),
        ("OCCc4ccc(c7ccc8ccccc8c7)nc4", "2-[6-(naphthalen-2-yl)pyridin-3-yl]ethan-1-ol"),
        ("OCCc4ccc(C7CCCCO7)nc4", "2-[6-(oxan-2-yl)pyridin-3-yl]ethan-1-ol"),
    ],
)
def test_complex_substituted_substituent_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_smiles_to_iupac_not_implemented():
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("[Li]Cl")


def test_tetrahydronaphthalene():
    assert smiles_to_iupac("C1CCCc2ccccc12") == "1,2,3,4-tetrahydronaphthalene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=CC2CCCCC2C=C1", "1,2,3,4,4a,8a-hexahydronaphthalene"),
    ],
)
def test_partially_unsaturated_naphthalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_arbitrary_double_bond_placement_not_matched():
    assert smiles_to_iupac("C1CC2CC=CCC2C=C1") == "1,2,4a,5,8,8a-hexahydronaphthalene"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C#CF", "fluoroethyne"),
    ],
)
def test_halogen_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_lone_halogen_atom_raises():
    # No carbon atom at all: no hydrocarbon parent hydride to substitute.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("FCl")


def test_non_halogen_heteroatom():
    assert smiles_to_iupac("c1ccccc1[SiH3]") == "phenylsilane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C1CCC(=[AsH])CC1", "4-arsanylidenecyclohexane-1-carboxylic acid"),
        ("NNNc1ccc(C(=O)O)cc1", "4-(triazan-1-yl)benzoic acid"),
        ("N(N=N)C1=CC=C(C(=O)O)C=C1", "4-(triaz-2-en-1-yl)benzoic acid"),
        ("[SiH3]N[SiH2]c1ccc(C(=O)O)cc1", "4-[(silylamino)silyl]benzoic acid"),
        ("OC(=O)C1CCC(=NN)CC1", "4-hydrazinylidenecyclohexane-1-carboxylic acid"),
        ("CN(C)C(=O)NN=CCC(=O)O", "3-[(dimethylcarbamoyl)hydrazinylidene]propanoic acid"),
        (
            "[SiH3][SiH]([SiH3])[SiH2][SiH2]c1ccc(C(=O)O)cc1",
            "4-(3-silyltetrasilan-1-yl)benzoic acid",
        ),
        ("OC(=O)Cc1cccc([SiH2]O[SiH3])n1", "(6-disiloxanylpyridin-2-yl)acetic acid"),
    ],
)
def test_heteroatom_hydride_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("c1cc[se]c1", "selenophene"),
    ],
)
def test_smiles_to_iupac_hetero_monocyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Cc1cnccn1", "2-methylpyrazine"),  # CID 7976
        ("Cc1c(C)coc1", "3,4-dimethylfuran"),  # CID 34338
    ],
)
def test_smiles_to_iupac_hetero_monocyclic_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        pytest.param("[Te]1CC[Se]CC1", "1,4-selenatellurane", id="two_heteroatoms"),
        pytest.param("Cc1cc[nH]n1", "3-methyl-1H-pyrazole", id="pyrazole_ring_carbon_substituent"),
    ],
)
def test_two_heteroatoms_and_related(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "n,expected",
    [
        (100, "hecta"),
        (111, "undecahecta"),
    ],
)
def test_numerical_term(n, expected):
    assert numerical_term(n) == expected


def _androstane_with_two_hydroxyls():
    """A synthetic androstane-3,17-diol-shaped molecule: two -OH oxygens
    added at the C3 and C17 ring carbons (one -CH2- position each, valence
    reduced by one to make room), built programmatically off the already-
    verified plain androstane parent so its own ring-atom indices are known
    -- proves `strip_substituents` handles N=2 simultaneous removals, not
    just #1144's original single-substituent case."""
    mol = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
    rw = RWMol(mol)
    c1, c2 = 4, 15
    for c in (c1, c2):
        atom = rw.GetAtomWithIdx(c)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
    o1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c1, o1, Chem.BondType.SINGLE)
    o2 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c2, o2, Chem.BondType.SINGLE)
    diol = rw.GetMol()
    Chem.SanitizeMol(diol)
    return diol, c1, c2, o1, o2


def test_strip_substituents_removes_two_atoms_and_reindexes():
    diol, c1, c2, o1, o2 = _androstane_with_two_hydroxyls()
    stripped, old_to_new = strip_substituents(diol, {o1, o2})
    assert stripped is not None
    assert stripped.GetNumAtoms() == diol.GetNumAtoms() - 2
    # Both original ring carbons survive, correctly reindexed.
    assert c1 in old_to_new and c2 in old_to_new
    assert old_to_new[c1] != old_to_new[c2]


_STEM = r"(?:meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
_BAD_COMPOUND_SUBSTITUENT_RE = re.compile(rf"\b1-{_STEM}yl{_STEM}yl\b|\b1,1-di{_STEM}yl{_STEM}yl\b")

_TESTS_DIR = pathlib.Path(__file__).parent


def _string_literals(path):
    with open(path, "rb") as f:
        tokens = tokenize.tokenize(f.readline)
        for tok in tokens:
            if tok.type == tokenize.STRING:
                yield tok.start[0], tok.string


def _violations():
    for path in sorted(_TESTS_DIR.glob("*.py")):
        if path == pathlib.Path(__file__):
            continue
        for lineno, literal in _string_literals(path):
            if _BAD_COMPOUND_SUBSTITUENT_RE.search(literal):
                yield f"{path.name}:{lineno}: {literal}"


def test_no_pre_pin_branch_point_substituent_names():
    violations = list(_violations())
    assert not violations, (
        "found pre-PIN CAS-style compound substituent name(s) in test "
        "expected values -- renumber through the branch point per "
        "P-29.3.2.2 instead:\n" + "\n".join(violations)
    )


@pytest.mark.parametrize(
    "smiles",
    ["C[Ti](Cl)(Cl)Cl", "C[Li]", "C[Zn]C", "CC1C[Pt](Cl)(Cl)C1"],
)
def test_organometallic_without_a_pin_warns(smiles):
    with pytest.warns(NonPreferredNameWarning, match="no PIN"):
        smiles_to_iupac(smiles)


def test_warning_is_emitted_once_per_call():
    with pytest.warns(NonPreferredNameWarning) as record:
        smiles_to_iupac("COC(=O)C(C)CC1C[Pt](C)(I)(P(CC)(CC)CC)(P(CC)(CC)CC)C1")
    assert len(record) == 1


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("OC(=O)C1CCC(=[Se])CC1", "4-selanylidenecyclohexane-1-carboxylic acid"),
    ],
)
def test_selanyl_and_tellanyl_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("COCCOCCOCCOCCCc1ccc(C(=O)O)cc1", "4-(2,5,8,11-tetraoxatetradecan-14-yl)benzoic acid"),
        ("OC(=O)CC1CCCCCOCCCCC1C", "(6-methyl-1-oxacyclododecan-7-yl)acetic acid"),
    ],
)
def test_skeletal_replacement_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("N[C@H]1CCC[C@@H](N)C1", "(1R,3S)-cyclohexane-1,3-diamine"),
        ("OC(=O)[C@H]1CCC[C@@H](C(O)=O)C1", "(1R,3S)-cyclohexane-1,3-dicarboxylic acid"),
    ],
)
def test_r_is_cited_at_the_lower_locant_of_a_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.slow
@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H]1CCCC[C@@H]1C1CCCCC1C1CCCCC1", "(1¹S,1²S)-1²-methyl-1¹,2¹:2²,3¹-tercyclohexane"),
    ],
)
def test_stereodescriptors_in_an_assembly_of_three_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_substituent_is_named_by_the_general_ring_namer():
    plain = Chem.MolFromSmiles("CC1CC1")
    assert name_branch(adjacency(plain), 1, 0, {}, mol=plain) == ("cyclopropyl", False)
    substituted = Chem.MolFromSmiles("CC1CC1C")
    assert name_branch(adjacency(substituted), 1, 0, {}, mol=substituted) == ("2-methylcyclopropyl", True)


def test_substituted_ring_substituent_raises():
    graph = {0: [1], 1: [0, 2, 3], 2: [1, 3, 4], 3: [1, 2], 4: [2]}
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0)


def test_mol_names_unsaturated_branch():
    mol = Chem.MolFromSmiles("CCC=C")
    graph = adjacency(mol)
    assert name_branch(graph, 1, 0, {}) == ("propyl", False)
    assert name_branch(graph, 1, 0, {}, mol=mol) == ("prop-2-en-1-yl", True)
    with pytest.raises(UnsupportedStructure):
        name_branch(graph, 1, 0, {}, mol=mol, unsaturated=False)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C=C", "ethene"),
        ("FC=CF", "1,2-difluoroethene"),
        ("FC=CCl", "1-chloro-2-fluoroethene"),
    ],
)
def test_smiles_to_iupac_unsaturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C/C=C/C", "(2E)-but-2-ene"),
    ],
)
def test_smiles_to_iupac_ez_double_bond(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("C/C=C/CC#C", id="specified_double_bond_with_triple_bond_raises"),
        pytest.param("C/C=C/C=CC", id="partially_specified_diene_raises"),
    ],
)
def test_specified_double_bond_and_related_raise(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize("smiles", ["", "   ", "xyz", "C("])
def test_invalid_smiles_raises_value_error(smiles):
    with pytest.raises(ValueError):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize("value", [None, 5, ["CCO"]])
def test_non_string_input_raises_type_error(value):
    with pytest.raises(TypeError):
        smiles_to_iupac(value)


@pytest.mark.parametrize(
    "smiles",
    [
        "CN1CCC[C@H]1c1cccnc1",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "c1ccc2[nH]ccc2c1",
        "OC(=O)Cn1c2ccccc2c2ccccc21",
        pytest.param("C1=C/c2cccc(c2)CCCCCCCc2cccc(c2)CCCCC/1", marks=pytest.mark.slow),
    ],
)
def test_valid_smiles_do_not_write_rdkit_logs_to_stderr(smiles, capfd):
    smiles_to_iupac(smiles)
    assert capfd.readouterr().err == ""


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-93.5.1: a centre bearing one R and one S ligand is pseudoasymmetric (lowercase s) in the group alone
        ("c1ccccc1[C@@H]([C@H](C)O)[C@@H](C)O", {6: "s"}),
        ("c1ccncc1[C@@H]([C@H](C)O)[C@@H](C)O", {6: "s"}),
        # two like ligands leave the centre achiral, so the group has no pseudoasymmetric centre
        ("c1ccccc1[C@@H]([C@H](C)O)[C@H](C)O", {}),
    ],
)
def test_pseudoasymmetric_centre_in_a_group_on_an_aromatic_atom(smiles, expected):
    from rdkit.Chem import rdCIPLabeler

    from smiles_to_iupac._diester_anions import _pseudoasymmetric_in_group

    mol = Chem.MolFromSmiles(smiles)
    side = {a.GetIdx() for a in mol.GetAtoms() if not a.GetIsAromatic()}
    rdCIPLabeler.AssignCIPLabels(mol)
    assert _pseudoasymmetric_in_group(mol, side) == expected


@pytest.mark.slow
def test_nested_fragment_is_named_once_per_top_level_call(monkeypatch):
    from smiles_to_iupac import core

    seen = []
    original = core._name_unabridged
    monkeypatch.setattr(core, "_name_unabridged", lambda smiles: seen.append(smiles) or original(smiles))
    core._NESTED_NAMES.clear()
    smiles = "CS[C@H](C(=O)N[C@H](C(=O)O)C(C)(C)C)C(C)(C)c1ccccc1"
    name = smiles_to_iupac(smiles)
    assert len(seen) == len(set(seen))
    assert smiles_to_iupac(smiles) == name


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Nc1ccccc1S(=O)(=O)N1CCC1", "2-[(azetidin-1-yl)sulfonyl]aniline"),
        ("OC(=O)CS(=O)(=O)N1CCCC1", "[(pyrrolidin-1-yl)sulfonyl]acetic acid"),
        ("OC(=O)c1ccccc1S(=O)(=O)N1CCOCC1", "2-[(morpholin-4-yl)sulfonyl]benzoic acid"),
        ("OC(=O)CNC(=O)N1CCCC1", "(pyrrolidine-1-carboxamido)acetic acid"),
    ],
)
def test_ring_nitrogen_acyl_prefix_keeps_the_ring_intact(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CN=C(NCC)NCCC(=O)O", "3-[(N-ethyl-N'-methylcarbamimidoyl)amino]propanoic acid"),
        ("CN=C(NC)N(C)CCC(=O)O", "3-[methyl(N,N'-dimethylcarbamimidoyl)amino]propanoic acid"),
    ],
)
def test_guanidine_prefix_cites_the_substituents_of_every_nitrogen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CCOC(=O)N[C@@H](C)CC",
        "C/C=C/c1ccc(C(=O)OC)cc1",
    ],
)
def test_stereo_is_never_silently_dropped_from_a_name(smiles):
    try:
        name = smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return
    assert re.search(r"\((?:\d+[A-Za-z]?,?)*[RSEZ]\)|[RSEZ]\)", name)


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=[N+]([O-])c1ccc(NC)cc1", "N-methyl-4-nitroaniline"),
        ("N#Cc1ccc(NC)cc1", "4-(methylamino)benzonitrile"),
        ("CCN1CCC(CNc2ccc([N+](=O)[O-])cc2I)CC1", "N-[(1-ethylpiperidin-4-yl)methyl]-2-iodo-4-nitroaniline"),
    ],
)
def test_amine_parent_leaves_other_nitrogens_of_the_arm_untouched(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_interior_fusion_locants_order_between_their_peripheral_neighbours():
    from smiles_to_iupac._phane_amplificant import AmpLoc

    assert [str(x) for x in sorted(map(AmpLoc, ["4b", "5", "4a1", "4a", "4"]))] == ["4", "4a", "4a1", "4b", "5"]


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O=C(Cl)CCC(=O)c1ccc(S(=O)(=O)Cl)cc1", "4-[4-(chlorosulfonyl)phenyl]-4-oxobutanoyl chloride"),
        (
            "CCCCOC(=O)N1CCC(NC(=O)c2ccc(C(=O)OC)cc2)CC1",
            "butyl 4-{[4-(methoxycarbonyl)benzoyl]amino}piperidine-1-carboxylate",
        ),
        (
            "CC(C)(C)OC(=O)CC(N)C1COCC(c2ccc(Br)cc2)N1C(=O)OC(C)(C)C",
            "tert-butyl 3-[1-amino-2-(tert-butoxycarbonyl)ethyl]-5-(4-bromophenyl)morpholine-4-carboxylate",
        ),
        ("CCOC(=O)c1ccc(C(=O)OC)cc1", "ethyl methyl benzene-1,4-dicarboxylate"),
        ("ClC(=O)CCCC(Cl)=O", "pentanedioyl dichloride"),
    ],
)
def test_ester_and_halide_groups_the_parent_cannot_carry_are_cited_as_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("N[C@@H](C)C(=O)N[C@@H](C)C(=O)O", "N-[(2S)-2-aminopropanoyl]-L-alanine"),
        ("N[C@@H](CS)C(=O)N[C@@H](C)C(=O)O", "N-[(2R)-2-amino-3-sulfanylpropanoyl]-L-alanine"),
        (
            "CC[C@H](C)[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O",
            "N-((2S)-2-{[(2S,3S)-2-amino-3-methylpentanoyl]amino}propanoyl)-L-alanine",
        ),
    ],
)
def test_peptide_acyl_residues_cite_their_own_stereodescriptors(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCOP(=O)(O)OP(=O)(O)OC", "1-ethyl 3-methyl dihydrogen diphosphate"),
        ("CCOP(=O)(OC)OP(=O)(O)O", "1-ethyl 1-methyl dihydrogen diphosphate"),
        ("CCOP(=O)(O)OP(=O)(O)OCC", "1,3-diethyl dihydrogen diphosphate"),
        ("CCOP(=O)(OCC)OP(=O)(O)O", "1,1-diethyl dihydrogen diphosphate"),
        ("CCOP(=O)(OCC)OP(=O)(O)OC", "1,1-diethyl 3-methyl hydrogen diphosphate"),
        ("COP(=O)(OC)OP(=O)(O)OCC", "3-ethyl 1,1-dimethyl hydrogen diphosphate"),
        ("ClCCOP(=O)(O)OP(=O)(O)OC", "1-(2-chloroethyl) 3-methyl dihydrogen diphosphate"),
        ("CCCCOP(=O)(O)OP(=O)(O)O", "butyl trihydrogen diphosphate"),
        ("CCOP(=O)(OC)OP(=O)(OC)OC", "ethyl trimethyl diphosphate"),
        ("CCOP(=O)(OCC)OP(=O)(OCC)OCC", "tetraethyl diphosphate"),
    ],
)
def test_diphosphate_esters_cite_locants_when_the_arrangement_is_otherwise_ambiguous(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("[O-][Br+]c1ccccc1", "bromosylbenzene", id="bromosyl_charge_separated"),
        pytest.param("OCC[Cl+3]([O-])([O-])[O-]", "2-perchlorylethan-1-ol", id="perhalyl_on_chain_with_alcohol"),
    ],
)
def test_halogen_oxo_prefixes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        pytest.param("ClC(C(=O)O)CC(=O)O", "chlorobutanedioic acid", id="sole_substituent_on_symmetric_chain"),
        pytest.param("ClC(C(=O)O)C(=O)O", "chloropropanedioic acid", id="sole_substituent_between_two_acid_groups"),
        pytest.param("CC(Cl)C(=O)O", "2-chloropropanoic acid", id="asymmetric_chain_keeps_locant"),
        pytest.param("ClC(C(=O)O)=CC(=O)O", "2-chlorobut-2-enedioic acid", id="cited_ene_locant_keeps_locants"),
        pytest.param("FC(C(C(C(=O)O)(F)F)(F)F)(F)F", "heptafluorobutanoic acid", id="fully_substituted_chain_acid"),
        pytest.param("FC(F)(F)C(F)(F)F", "hexafluoroethane", id="fully_substituted_chain"),
        pytest.param("FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F", "dodecafluorocyclohexane", id="fully_substituted_ring"),
        pytest.param("FC(F)C(F)(F)C(F)(F)F", "1,1,1,2,2,3,3-heptafluoropropane", id="partial_substitution_keeps_locants"),
        pytest.param("ClC1=CC2=CC=C3C=CC4=CC=C5C=CC6=CC=C1C1=C6C5=C4C3=C21", "chlorocoronene", id="symmetric_fused_system"),
        pytest.param("Cc1cccc2ccccc12", "1-methylnaphthalene", id="asymmetric_fused_system_keeps_locant"),
        pytest.param("ClC1(OOO1)Cl", "dichlorotrioxetane", id="fully_substituted_hetero_ring"),
        pytest.param("CNC(=O)N", "methylurea", id="sole_substituent_on_urea"),
        pytest.param("FN=C(C(C(F)(F)F)(F)F)N(F)F", "octafluoropropanimidamide", id="fully_substituted_imidamide"),
        pytest.param("FN(F)C(=O)C(F)(F)F", "pentafluoroacetamide", id="fully_substituted_amide"),
        pytest.param("CNC=O", "N-methylformamide", id="single_nitrogen_position_keeps_locant"),
        pytest.param("ClC1=C(C=CC=C1)C(C(F)(F)F)(F)F", "1-chloro-2-(pentafluoroethyl)benzene", id="fully_substituted_substituent_group"),
        pytest.param("C1(=C(C(=C(C(=C1O)O)O)O)O)O", "benzenehexol", id="every_ring_position_carries_the_suffix"),
        pytest.param("OC1C(O)C(O)C(O)C(O)C1O", "cyclohexane-1,2,3,4,5,6-hexol", id="partly_modified_ring_keeps_suffix_locants"),
        pytest.param("CC(=O)CCl", "1-chloropropan-2-one", id="cited_suffix_locant_keeps_prefix_locant"),
    ],
)
def test_locants_without_information_are_omitted(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C(CCCCCCCCC)C1CC(CC(C1)CCCCCCCCCC)CCCCCCCCCC", "1,3,5-tri(decyl)cyclohexane"),
        ("C(CCCCCCCCCCCC)C1=CC=C(C=C1)CCCCCCCCCCCCC", "1,4-di(tridecyl)benzene"),
        ("C(C)(C)(C)C1=C(C=CC=C1)C(C)(C)C", "1,2-di-tert-butylbenzene"),
        ("SC(CC(=O)O)CS", "3,4-bis(sulfanyl)butanoic acid"),
        ("S(S)C=1C=C(C(=O)N)C=CC1SS", "3,4-bis(disulfanyl)benzamide"),
        ("CC(CC)N(C(C)(CC)O)C(C)CC", "2-[di(butan-2-yl)amino]butan-2-ol"),
        (
            "C12CC(CC(CC1)C2)C2=CC=CC1=CC3=CC=CC(=C3C=C21)C2CC1CCC(C2)C1",
            "1,8-di(bicyclo[3.2.1]octan-3-yl)anthracene",
        ),
        ("C1(=CC=CC=C1)S(=O)C(C(=O)O)S(=O)C1=CC=CC=C1", "di(benzenesulfinyl)acetic acid"),
        ("N#CC(=S)c1ccccc1C(=O)Cl", "2-(carbonocyanidothioyl)benzoyl chloride"),
        ("C1=C(C=CC2=CC=CC=C12)[Se](=NN)(=N)S", "naphthalene-2-selenonohydrazonimidothioic acid"),
        ("C=1SC=CN=CC=COC=C2C1C=CC=C2", "9,2,5-benzoxathiaazacyclododecine"),
    ],
)
def test_multiplying_prefixes_parentheses_and_elision(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("Br[C@H](F)Cl", "(R)-bromo(chloro)(fluoro)methane"),
        ("C1(CCC1)[C@@H](O)C1CC1", "(S)-cyclobutyl(cyclopropyl)methanol"),
        ("CCC(C)(C(=O)OCC)C(=O)OCC", "diethyl ethyl(methyl)propanedioate"),
        ("Br[C@H](F)C", "(1S)-1-bromo-1-fluoroethane"),
        ("ClC(Cl)Cl", "trichloromethane"),
    ],
)
def test_mononuclear_and_single_site_parents_enclose_later_prefixes_and_drop_stereo_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C12C3C4C5C3C1C5C24", "cubane"),
        ("OC(=O)C12C3C4C1C5C2C3C45", "cubane-1-carboxylic acid"),
        ("CC12C3C4C1C5C2C3C45", "1-methylcubane"),
        ("OC12C3C4C1C5C2C3C45", "cuban-1-ol"),
        ("Cl[C]12C3C4C1C5C2C3C45", "1-chlorocubane"),
        ("C12CC3CC(CC(C1)C3)C2", "adamantane"),
    ],
)
def test_retained_cage_names_cubane_and_adamantane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("Cc1ccccc1", "toluene"),
        ("Cc1ccccc1C", "1,2-xylene"),
        ("Cc1cccc(C)c1", "1,3-xylene"),
        ("Cc1ccc(C)cc1", "1,4-xylene"),
        ("Cc1cc(C)cc(C)c1", "1,3,5-trimethylbenzene"),
        ("Cc1ccccc1Cl", "1-chloro-2-methylbenzene"),
        ("Cc1ccc(C(=O)O)cc1", "4-methylbenzoic acid"),
    ],
)
def test_toluene_and_xylenes_only_when_unsubstituted(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("O1C=CC=CC=CC=CC=C1", "1-oxacycloundeca-2,4,6,8,10-pentaene"),
        ("N1C=CC=CC=CC=CC=CC=C1", "1-azacyclotrideca-2,4,6,8,10,12-hexaene"),
        ("O1NC=CC=CC=CC=CC=C1", "1-oxa-2-azacyclododeca-3,5,7,9,11-pentaene"),
        ("O1C=C[Se]C=CC=CC=CNC=C1", "1-oxa-4-selena-11-azacyclotrideca-2,5,7,9,12-pentaene"),
        ("O1C=CC=CC=COC=CC=CC=CC=CC=C1", "1,8-dioxacyclooctadeca-2,4,6,9,11,13,15,17-octaene"),
        ("O1CC=NC=CC=NC=CN=CC=C1", "1-oxa-4,8,11-triazacyclotetradeca-3,5,7,9,11,13-hexaene"),
        ("O1\\C=C/OCCOCCOCC1", "(2Z)-1,4,7,10-tetraoxacyclododec-2-ene"),
        ("O1C\\C=N\\CCCCCCCC1", "(3E)-1-oxa-4-azacyclododec-3-ene"),
    ],
)
def test_unsaturated_heteromacrocycles_by_skeletal_replacement(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_double_bond_to_an_ylidene_group_takes_the_parent_locant():
    assert (
        smiles_to_iupac("Cl/C(/C(/C=C/C(=O)O)=C/S(=O)(=O)O)=C\\C")
        == "(2E,4E,5Z)-5-chloro-4-(sulfomethylidene)hepta-2,5-dienoic acid"
    )
