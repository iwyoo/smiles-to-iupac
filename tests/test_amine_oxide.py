import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_tertiary_amine_oxide():
    # Blue Book P-62.5 worked example: (CH3)3N+-O- ->
    # "N,N-dimethylmethanamine N-oxide (PIN)" (trimethylamine N-oxide).
    assert smiles_to_iupac("C[N+](C)(C)[O-]") == "N,N-dimethylmethanamine N-oxide"


def test_secondary_amine_oxide():
    # Structure/base-name cross-check against _amine.py's own
    # "N-ethylethanamine" for the un-oxidized diethylamine.
    assert smiles_to_iupac("CC[NH+](CC)[O-]") == "N-ethylethanamine N-oxide"


def test_tertiary_amine_oxide_asymmetric():
    assert smiles_to_iupac("CC[N+](C)(CCC)[O-]") == "N-ethyl-N-methylpropan-1-amine N-oxide"


def test_primary_amine_oxide_not_matched():
    # PubChem itself doesn't name a primary amine N-oxide with the same
    # "<amine> N-oxide" pattern (it computes an ammonium-salt-style name
    # instead), so this module deliberately doesn't claim this shape --
    # it falls through to _amine.py's own charged-atom rejection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+][O-]")


def test_amine_oxide_ring_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[N+]1([O-])CCCCC1")


def test_amine_oxide_halogen_on_parent_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("ClCC[N+](CC)(CC)[O-]")


def test_plain_tertiary_amine_still_works():
    assert smiles_to_iupac("CN(C)C") == "N,N-dimethylmethanamine"


def test_amine_oxide_stereocenter():
    # This module always delegates to `name_amine` on the reduced (neutral)
    # molecule, so it inherited full R/S support automatically once
    # `_amine.py` gained it (P-91.3/P-92) -- no separate wiring needed here,
    # this is a regression test locking that in.
    assert smiles_to_iupac("CC[C@@H](C)[N+](C)(C)[O-]") == "(2R)-N,N-dimethylbutan-2-amine N-oxide"
    assert smiles_to_iupac("CC[C@H](C)[N+](C)(C)[O-]") == "(2S)-N,N-dimethylbutan-2-amine N-oxide"
