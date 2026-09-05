import pytest

from smiles_to_iupac._coexisting_groups import name_via_senior_acyclic


def test_name_via_senior_acyclic_rejects_reversed_seniority():
    with pytest.raises(AssertionError):
        name_via_senior_acyclic(lambda *a, **k: "unused", "alcohol", "sulfonic_acid", (), {})


def test_name_via_senior_acyclic_calls_through_with_extra_names():
    def fake_namer(x, extra_names=None, required_atoms=frozenset()):
        return (x, extra_names, required_atoms)

    result = name_via_senior_acyclic(
        fake_namer, "sulfonic_acid", "alcohol", (42,), {1: "sulfanyl"}, required_atoms={1}
    )
    assert result == (42, {1: "sulfanyl"}, {1})
