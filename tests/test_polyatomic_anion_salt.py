import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Potassium sulfate (PubChem CID 24507's own structure): 2 K+ for
        # 1 SO4^2- puts the multiplying prefix on the cation side.
        ("[K+].[K+].[O-]S(=O)(=O)[O-]", "dipotassium sulfate"),
        # Sodium carbonate (PubChem CID 10340), matching the Blue Book's
        # own P-65.6.2 worked example 'disodium carbonate (PIN)'.
        ("[Na+].[Na+].[O-]C(=O)[O-]", "disodium carbonate"),
        # Potassium nitrate (PubChem CID 24434): nitrate's charge
        # magnitude 1 needs no cation multiplying prefix.
        ("[K+].[O-][N+](=O)[O-]", "potassium nitrate"),
        # Trisodium phosphate (PubChem CID 24243): magnitude 3, three Na+.
        ("[Na+].[Na+].[Na+].[O-]P(=O)([O-])[O-]", "trisodium phosphate"),
        # Ammonium sulfate (PubChem CID 62648): the ammonium-shaped cation
        # path (this project's own 'azanium' PIN, not PubChem's own
        # 'diazanium sulfate' divergence-free match here).
        ("[NH4+].[NH4+].[O-]S(=O)(=O)[O-]", "diazanium sulfate"),
        # Calcium carbonate (PubChem CID 10112): a divalent cation and a
        # -2 anion in a 1:1 ratio needs no multiplying prefix on either
        # side (magnitude 2 / charge 2 = 1 copy).
        ("[Ca+2].[O-]C(=O)[O-]", "calcium carbonate"),
        # Aluminium phosphate (PubChem CID 24418): a trivalent cation and
        # a -3 anion, same 1:1 ratio reasoning.
        ("[Al+3].[O-]P(=O)([O-])[O-]", "aluminium phosphate"),
        # Sodium nitrate (PubChem CID 24268).
        ("[Na+].[O-][N+](=O)[O-]", "sodium nitrate"),
    ],
)
def test_polyatomic_anion_salt_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_polyatomic_anion_mismatched_charge_ratio_raises():
    # Al3+ paired with SO4^2- needs a mixed 2:3 cation:anion count, out of
    # scope for this single-anion-fragment mechanism -- falls through to
    # an unrelated later module's own (unhelpful) rejection, which is
    # expected and accepted per this feature's own scope note.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Al+3].[O-]S(=O)(=O)[O-]")


def test_plain_carboxylate_salt_still_works():
    # Regression check: the new polyatomic-anion path sits alongside, not
    # in place of, the existing organic-anion path.
    assert smiles_to_iupac("[Na+].CC(=O)[O-]") == "sodium ethanoate"
