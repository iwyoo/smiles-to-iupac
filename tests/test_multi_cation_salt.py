import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-65.6.2's own worked example shape: two different +1 cations
        # balancing one -2 carbonate, cited alphabetically, no
        # multiplying prefix on either cation.
        ("[K+].[Na+].[O-]C(=O)[O-]", "potassium sodium carbonate"),
        # Ammonium (this project's own 'azanium' PIN convention) plus a
        # metal cation balancing sulfate.
        ("[NH4+].[K+].[O-]S(=O)(=O)[O-]", "azanium potassium sulfate"),
        # A +1 and a +2 cation together balancing phosphate's -3 charge.
        ("[NH4+].[Ca+2].[O-]P(=O)([O-])[O-]", "azanium calcium phosphate"),
        # A different +2/+1 combination, same mechanism.
        ("[Ca+2].[Na+].[O-]P(=O)([O-])[O-]", "calcium sodium phosphate"),
    ],
)
def test_multi_cation_salt_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_same_cation_type_prefers_multiplying_prefix():
    # Two copies of the *same* cation still use the multiplying-prefix
    # shape (#788), not this module's different-cation-types shape.
    assert smiles_to_iupac("[Na+].[Na+].[O-]C(=O)[O-]") == "disodium carbonate"


def test_multi_cation_mismatched_charge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Al+3].[K+].[O-]S(=O)(=O)[O-]")
