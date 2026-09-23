import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # No real PubChem-registered polycyclic/spiro tellurol structure
        # was found (sparse tellurol coverage in general, same caveat
        # `_tellurol.py`'s own module docstring and existing monocyclic
        # tests already document) -- verified instead by structural
        # analogy against the identical, already-independently-verified
        # mechanism on the same ring shapes for the sulfur/selenium
        # analogues (`test_von_baeyer_thiol.py`'s `SC1CC2CCC1C2`/CID
        # 164174492, `test_von_baeyer_spiro_selenol.py`'s matching
        # bicyclic case), reviewed rather than independently
        # PubChem-confirmed (see `implementation-notes.md`'s test-writing
        # policy) -- the same method `_tellurol.py`'s sibling
        # `_tellone.py`'s own polycyclic step (#887/#892) already used.
        ("[TeH]C1CC2CCC1C2", "bicyclo[2.2.1]heptane-2-tellurol"),
        ("[TeH]C1CCCC2(C1)CCCCC2", "spiro[5.5]undecane-2-tellurol"),
        ("[TeH]C1CCC2CCC1C2", "bicyclo[3.2.1]octane-2-tellurol"),
    ],
)
def test_von_baeyer_spiro_tellurol_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_multiple_ring_tellurols_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CC2CCC1C2[TeH]")


def test_unsaturated_von_baeyer_tellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CC2C=CC1C2")
