import pytest
from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._cyclo_steroid import _CYCLO_LOOKUP
from smiles_to_iupac._dinor_steroid import _DINOR_LOOKUP
from smiles_to_iupac._homo_steroid import _HOMO_LOOKUP
from smiles_to_iupac._nor_steroid import _NOR_LOOKUP
from smiles_to_iupac._seco_steroid import _SECO_LOOKUP


def _smiles_for_1(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _CYCLO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no cyclo-lookup entry for {label}-cyclo-{parent}")


def test_lookup_spans_multiple_parent_families():
    names = {name for _, name in _CYCLO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_gonane_pair():
    assert smiles_to_iupac(_smiles_for_1(2, 7, "gonane")) == "2,7-cyclo-gonane"


def test_no_bond_pair_excluded():
    labels_names = set(_CYCLO_LOOKUP.values())
    assert ("1,2", "androstane") not in labels_names


def test_methyl_locants_excluded():
    for label, _ in _CYCLO_LOOKUP.values():
        a, b = (int(x) for x in label.split(","))
        assert a <= 17 and b <= 17


def test_plain_cholane_still_resolves():
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"


def _smiles_for_2(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _DINOR_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no dinor-lookup entry for {label}-dinor-{parent}")


def test_lookup_spans_multiple_parent_families__dinor_steroid():
    names = {name for _, name in _DINOR_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_both_methyls_excludes_as_gonane_collision():
    labels_names = {(lbl, name) for lbl, name in _DINOR_LOOKUP.values()}
    assert ("18,19", "androstane") not in labels_names


def test_gonane_pair__dinor_steroid():
    assert smiles_to_iupac(_smiles_for_2(1, 15, "gonane")) == "1,15-dinor-gonane"


def test_dinor_lookup_uses_lowest_locant_pair():
    for lbl, name in _DINOR_LOOKUP.values():
        a, b = (int(x) for x in lbl.split(","))
        assert a < b


def _smiles_for_3(label, parent):
    for smi, (lbl, name) in _HOMO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no homo-lookup entry for {label}-homo-{parent}")


def test_lookup_spans_multiple_parent_families__homo_steroid():
    names = {name for _, name in _HOMO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_ring_b_bond_cholane():
    assert smiles_to_iupac(_smiles_for_3("5(6)a", "cholane")) == "5(6)a-homo-cholane"


def _smiles_for_4(locant, parent):
    for smi, (loc, name) in _NOR_LOOKUP.items():
        if loc == locant and name == parent:
            return smi
    raise AssertionError(f"no nor-lookup entry for {locant}-nor-{parent}")


def test_lookup_spans_multiple_parent_families__nor_steroid():
    names = {name for _, name in _NOR_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_17_nor_gonane_ring_position():
    assert smiles_to_iupac(_smiles_for_4(17, "gonane")) == "17-nor-gonane"


def test_ring_fusion_atom_not_eligible_for_nor():
    names = {(loc, name) for loc, name in _NOR_LOOKUP.values()}
    for parent in (
        "gonane",
        "androstane",
        "estrane",
        "pregnane",
        "cholane",
        "cholestane",
        "ergostane",
    ):
        for fusion_locant in (5, 8, 9, 10, 13, 14):
            assert (fusion_locant, parent) not in names


def _smiles_for_5(lo, hi, parent):
    label = f"{lo},{hi}"
    for smi, (lbl, name) in _SECO_LOOKUP.items():
        if lbl == label and name == parent:
            return smi
    raise AssertionError(f"no seco-lookup entry for {label}-seco-{parent}")


def test_lookup_spans_multiple_parent_families__seco_steroid():
    names = {name for _, name in _SECO_LOOKUP.values()}
    assert {"androstane", "pregnane", "cholestane", "estrane", "gonane"} <= names


def test_gonane_ring_a_bond():
    assert smiles_to_iupac(_smiles_for_5(3, 4, "gonane")) == "3,4-seco-gonane"


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C",
            "3α,7α,12α-trihydroxy-5β-cholan-24-oic acid",
        ),
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)O",
            "17α-ethynylestra-1,3,5(10)-triene-3,17β-diol",
        ),
        (
            "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2N)CCC4=CC(=O)CC[C@]34C",
            "17β-aminoandrost-4-en-3-one",
        ),
    ],
)
def test_natural_steroids_with_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        (
            "CC(=O)[C@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@@]21C",
            "17-hydroxy-17α-pregn-4-ene-3,20-dione",
        ),
        ("C[C@]12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", "17ξ-hydroxy-8ξ,9ξ,10ξ,14ξ-estr-4-en-3-one"),
    ],
)
def test_configuration_that_differs_from_the_parent_or_is_unspecified(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_non_natural_stereo_specified_parent_hydrides_are_not_named_without_their_stereo():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC2[C@H](C1)CCC3C2CCC4C3CCC4")


def test_extra_methyl_beyond_each_parent_hydride_falls_through_to_von_baeyer():
    assert (
        smiles_to_iupac("CC12CCCC1(C)CCC1C2CCC2C1CCCC2")
        == "11,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C3CC(C)C4CCCCC4(C3CC2)C")
        == "2,8,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C1CC(C)C3CCCCC3C1CC2")
        == "8,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCC1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-ethyl-2,12,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(pentan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)C(C)CC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-(4,5-dimethylhexan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_androstane_other_diastereomer_gets_its_own_alpha_beta_citation():
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "5alpha,8alpha-androstane"
    )


def test_androst_4_ene():
    assert smiles_to_iupac("CC12CCCC1C1CCC3=CCCCC3(C)C1CC2") == "androst-4-ene"
