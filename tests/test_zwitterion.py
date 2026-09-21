import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Glycine zwitterion (PubChem CID 750's own structure, InChI
        # NH2CH2COOH as the H3N+-CH2-COO- zwitterion). PubChem's own
        # IUPACName is "2-amino acetic acid" (drawn/named as the neutral
        # tautomer, not the zwitterion) -- this project instead follows
        # P-74.1.3's own zwitterion citation order directly (ammonium
        # nitrogen prefix on the carboxylate parent), and, like
        # `_carboxylic_acid_amine.py`'s own glycine divergence, keeps the
        # systematic "ethanoate" stem rather than the retained "acetate"
        # one for consistency with `_carboxylate.py`'s established
        # convention.
        ("C(C(=O)[O-])[NH3+]", "2-azaniumylethanoate"),
        # Alanine zwitterion.
        ("CC(C(=O)[O-])[NH3+]", "2-azaniumylpropanoate"),
        # Valine zwitterion (branched chain).
        ("CC(C)C(C(=O)[O-])[NH3+]", "2-azaniumyl-3-methylbutanoate"),
        # Leucine zwitterion.
        ("CC(C)CC(C(=O)[O-])[NH3+]", "2-azaniumyl-4-methylpentanoate"),
        # Beta-alanine zwitterion (ammonium not on the alpha carbon).
        ("C(CC(=O)[O-])[NH3+]", "3-azaniumylpropanoate"),
        # Betaine (trimethylglycine) zwitterion, PubChem CID 247's own
        # structure -- the quaternary nitrogen's own name
        # ("N,N-dimethylmethanaminium", `_ammonium.py`'s established PIN
        # derivation, not the alternative direct 'azanium'-multiplicative
        # style) gets 'yl' appended and is parenthesized since it carries
        # its own internal locant.
        ("C[N+](C)(C)CC(=O)[O-]", "2-(N,N-dimethylmethanaminiumyl)ethanoate"),
        # Taurine zwitterion (PubChem CID 1123's own structure, neutral
        # tautomer NH2CH2CH2SO3H) -- the ammonium+sulfonate family (WS2).
        ("[NH3+]CCS(=O)(=O)[O-]", "2-azaniumylethanesulfonate"),
        # Homotaurine zwitterion (PubChem CID 1646's own structure).
        ("[NH3+]CCCS(=O)(=O)[O-]", "3-azaniumylpropane-1-sulfonate"),
        # A quaternary-nitrogen sulfobetaine, same shape as the
        # trimethylglycine case above but on a sulfonate parent.
        ("C[N+](C)(C)CCS(=O)(=O)[O-]", "2-(N,N-dimethylmethanaminiumyl)ethanesulfonate"),
    ],
)
def test_zwitterion_name(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_zwitterion_ammonium_bonded_directly_to_carboxylate_carbon_raises():
    # P-74.1.1/1.2's same-parent case -- a different mechanism, not
    # handled by this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3+]C(=O)[O-]")


def test_zwitterion_ionic_center_in_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC([NH3+])C1C(=O)[O-]")


def test_zwitterion_multiple_carboxylates_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C(C(=O)[O-])(C(=O)[O-])[NH3+]")


def test_zwitterion_ammonium_bonded_directly_to_sulfonate_carbon_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3+]C(S(=O)(=O)[O-])")


def test_zwitterion_sulfonate_ionic_center_in_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CC([NH3+])C1S(=O)(=O)[O-]")


def test_zwitterion_multiple_sulfonates_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C(S(=O)(=O)[O-])(S(=O)(=O)[O-])[NH3+]")


def test_plain_metal_carboxylate_salt_still_works():
    # Regression check: the new zwitterion routing sits ahead of
    # `has_salt_shape` in `core.py`, but a genuine two-fragment salt (no
    # single-fragment ammonium+carboxylate shape) must still fall through
    # to `_salt.py` unaffected.
    assert smiles_to_iupac("[Na+].CC(=O)[O-]") == "sodium ethanoate"


def test_plain_ammonium_carboxylate_salt_still_works():
    assert smiles_to_iupac("C[NH3+].CC(=O)[O-]") == "methanaminium ethanoate"
