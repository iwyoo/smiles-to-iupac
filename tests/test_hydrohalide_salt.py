import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_aniline_hydrochloride():
    # PubChem PUG REST IUPACName for this exact SMILES.
    assert smiles_to_iupac("Nc1ccccc1.Cl") == "aniline;hydrochloride"


def test_methanamine_hydrobromide():
    assert smiles_to_iupac("CN.Br") == "methanamine;hydrobromide"


def test_methanamine_hydroiodide():
    assert smiles_to_iupac("CN.I") == "methanamine;hydroiodide"


def test_methanamine_hydrofluoride():
    assert smiles_to_iupac("CN.F") == "methanamine;hydrofluoride"


def test_secondary_amine_hydrochloride():
    assert smiles_to_iupac("CNC.Cl") == "N-methylmethanamine;hydrochloride"


def test_non_base_fragment_still_names():
    # PubChem applies the same mechanical ';hydrochloride' pattern even to
    # a non-basic fragment -- not restricted to organic bases.
    assert smiles_to_iupac("CCO.Cl") == "ethanol;hydrochloride"


def test_halide_fragment_order_independent():
    assert smiles_to_iupac("Cl.CN") == "methanamine;hydrochloride"


def test_ionized_form_routes_to_salt_module():
    # [NH3+]/[Cl-] (rather than the neutral N/Cl atoms this module itself
    # handles) is a different shape entirely -- an ammonium cation plus a
    # halide anion, now handled by `_salt.py` (see test_salt.py's own
    # anilinium chloride case), not this module.
    assert smiles_to_iupac("[NH3+]c1ccccc1.[Cl-]") == "anilinium chloride"


def test_two_halide_fragments_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCCN.Cl.Cl")


def test_plain_mixture_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCO.CCO")


def test_unsupported_base_fragment_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)Nc1ccccc1C(=O)OCCCC.Cl")
