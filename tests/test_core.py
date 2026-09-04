import pytest

from smiles_to_iupac import smiles_to_iupac


def test_smiles_to_iupac_not_implemented():
    # `UnsupportedStructure` is a `NotImplementedError` subclass
    # (`_common.py`); an organometallic compound (P-69, no metal-atom
    # handling anywhere in this project) is a genuinely unimplemented
    # smoke-test case -- unlike the aromatic ammonium example this test
    # used to use, which `_amine.py`'s aniline support (and `_ammonium.py`'s
    # existing neutralize-and-derive path) now handles.
    with pytest.raises(NotImplementedError):
        smiles_to_iupac("[Fe](Cl)(Cl)(Cl)")
