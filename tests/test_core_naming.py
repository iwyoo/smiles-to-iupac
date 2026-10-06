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


def test_tert_butyl_cited_before_methyl_on_ring():
    assert smiles_to_iupac("CC(C)(C)C1(C)CCCC(=O)C1") == "3-tert-butyl-3-methylcyclohexan-1-one"


def test_substituted_branched_fusion_is_named():
    assert smiles_to_iupac("Cc1ccc2c(c1)c1ccccc1c1ccccc21") == "2-methyltriphenylene"


def test_heteroaromatic():
    assert smiles_to_iupac("c1nc[nH]n1") == "1H-1,2,4-triazole"


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


def test_stereo_marker_no_longer_silently_dropped_without_scope_is_named():
    assert smiles_to_iupac("Clc1ccccc1[C@@H](Cl)CC") == "1-chloro-2-[(1S)-1-chloropropyl]benzene"


def test_two_stereocenters_on_one_substituent():
    assert smiles_to_iupac("c1ccccc1[C@@H](Cl)[C@@H](Cl)C") == "[(1R,2S)-1,2-dichloropropyl]benzene"


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
        ("[SiH3]N[SiH2]c1ccc(C(=O)O)cc1", "4-[(silylamino)silyl]benzoic acid"),
        ("OC(=O)C1CCC(=NN)CC1", "4-hydrazinylidenecyclohexane-1-carboxylic acid"),
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


def test_two_heteroatoms():
    assert smiles_to_iupac("[Te]1CC[Se]CC1") == "1,4-selenatellurane"


def test_pyrazole_ring_carbon_substituent():
    assert smiles_to_iupac("Cc1cc[nH]n1") == "3-methyl-1H-pyrazole"


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


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C[C@H]1CCCC[C@@H]1C1CCCCC1C1CCCCC1", "(11S,12S)-12-methyl-11,21:22,31-tercyclohexane"),
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


def test_specified_double_bond_with_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=C/CC#C")


def test_partially_specified_diene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C/C=C/C=CC")


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
        "C1=C/c2cccc(c2)CCCCCCCc2cccc(c2)CCCCC/1",
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
        ("OC(=O)CNC(=O)N1CCCC1", "[(pyrrolidine-1-carbonyl)amino]acetic acid"),
    ],
)
def test_ring_nitrogen_acyl_prefix_keeps_the_ring_intact(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CN=C(NCC)NCCC(=O)O", "3-[(N'-ethyl-N-methylcarbamimidoyl)amino]propanoic acid"),
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
            "methyl 4-{[1-(butoxycarbonyl)piperidin-4-yl]carbamoyl}benzoate",
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
