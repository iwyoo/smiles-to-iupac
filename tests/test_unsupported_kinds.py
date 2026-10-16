import pytest
from rdkit import Chem

import smiles_to_iupac
from smiles_to_iupac import RuleUndefinedStructure, UnsupportedStructure, smiles_to_alternative_names, smiles_to_iupac as name
from smiles_to_iupac._conjunctive import conjunctive_names


def _arylalkanoic_acid(chain_length, branch_at):
    chain = ["C(C)" if i == branch_at else "C" for i in range(chain_length)]
    return "c1ccccc1" + "".join(chain) + "C(=O)O"


def test_exports_are_subclasses_of_not_implemented_error():
    assert issubclass(UnsupportedStructure, NotImplementedError)
    assert issubclass(RuleUndefinedStructure, UnsupportedStructure)
    assert {"UnsupportedStructure", "RuleUndefinedStructure"} <= set(smiles_to_iupac.__all__)


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
