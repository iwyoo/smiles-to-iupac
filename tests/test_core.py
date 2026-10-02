import pytest

from smiles_to_iupac import smiles_to_iupac


def test_smiles_to_iupac_not_implemented():
    # `UnsupportedStructure` is a `NotImplementedError` subclass
    # (`_common.py`); a metal carbonyl complex (P-69.2 ligand naming beyond
    # halido/hydrido/alkyl) is a genuinely unimplemented smoke-test case.
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("[Fe](C#[O+])(Cl)Cl")
