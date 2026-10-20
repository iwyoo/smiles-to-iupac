import pytest
from rdkit import Chem

import smiles_to_iupac
from smiles_to_iupac import InternalInvariantError, InvalidSmilesError, RuleUndefinedStructure, SmilesToIupacError, UnsupportedStructure, smiles_to_alternative_names, smiles_to_iupac as name
from smiles_to_iupac._conjunctive import conjunctive_names


def _arylalkanoic_acid(chain_length, branch_at):
    chain = ["C(C)" if i == branch_at else "C" for i in range(chain_length)]
    return "c1ccccc1" + "".join(chain) + "C(=O)O"


def test_error_hierarchy():
    assert issubclass(UnsupportedStructure, NotImplementedError)
    assert issubclass(InternalInvariantError, RuntimeError)
    assert issubclass(InvalidSmilesError, ValueError)
    for kind in (UnsupportedStructure, RuleUndefinedStructure, InternalInvariantError, InvalidSmilesError):
        assert issubclass(kind, SmilesToIupacError)
    assert not issubclass(RuleUndefinedStructure, UnsupportedStructure)
    assert not issubclass(RuleUndefinedStructure, NotImplementedError)
    assert {"SmilesToIupacError", "UnsupportedStructure", "RuleUndefinedStructure", "InternalInvariantError", "InvalidSmilesError"} <= set(smiles_to_iupac.__all__)


@pytest.mark.parametrize("chain_length", [25, 26])
def test_side_chain_beyond_the_greek_alphabet_is_rule_undefined(chain_length):
    mol = Chem.MolFromSmiles(_arylalkanoic_acid(chain_length, 0))
    with pytest.raises(RuleUndefinedStructure, match="Greek alphabet"):
        conjunctive_names(mol, name)


def test_side_chain_within_the_greek_alphabet_is_named():
    mol = Chem.MolFromSmiles(_arylalkanoic_acid(23, 22))
    assert list(conjunctive_names(mol, name))


def test_rule_undefined_conjunctive_name_does_not_stop_the_other_alternatives():
    alternatives = smiles_to_alternative_names(_arylalkanoic_acid(25, 0))
    assert all(alternative.kind != "conjunctive" for alternative in alternatives)


@pytest.mark.parametrize("smiles", ["FCl", "CC.CC"])
def test_implementation_gap_is_not_rule_undefined(smiles):
    with pytest.raises(UnsupportedStructure) as raised:
        name(smiles)
    assert not isinstance(raised.value, RuleUndefinedStructure)


def test_substituted_fullerene_without_defined_numbering_is_rule_undefined():
    with pytest.raises(RuleUndefinedStructure, match="P-27.3"):
        name("CC12C3C4C5C1C1C6C2C2C3C3C4C4C5C1C1C6C2C3C41")


@pytest.mark.parametrize("function", [name, smiles_to_alternative_names])
@pytest.mark.parametrize("smiles", ["not_a_smiles", "", "C(("])
def test_invalid_smiles_is_an_invalid_smiles_error(function, smiles):
    with pytest.raises(InvalidSmilesError):
        function(smiles)


def test_runtime_error_handlers_let_internal_invariant_errors_through():
    import pathlib
    import re

    source_dir = pathlib.Path(smiles_to_iupac.__file__).parent
    for path in source_dir.glob("*.py"):
        lines = path.read_text().split("\n")
        for i, line in enumerate(lines):
            if re.match(r"\s*except .*\bRuntimeError\b", line):
                assert lines[i - 2].strip() == "except InternalInvariantError:", f"{path.name}:{i + 1}"


def test_carbenium_beside_a_carbanion_and_a_cationic_heteroatom_is_refused_not_unpacked():
    with pytest.raises(UnsupportedStructure):
        name("C=[O+]C[CH+]C[CH2-]")


def test_phosphoryl_cation_in_a_branch_is_refused_not_an_internal_error():
    with pytest.raises(UnsupportedStructure):
        name("C[P+](=O)OP(=O)([O-])[O-]")


@pytest.mark.parametrize("smiles", ["C1=CC=C2C=C[Po]C2=C1", "C1=CC2=CC=C[At]2C=C1"])
def test_fused_ring_with_an_element_that_has_no_replacement_prefix_is_refused(smiles):
    with pytest.raises(UnsupportedStructure):
        name(smiles)
