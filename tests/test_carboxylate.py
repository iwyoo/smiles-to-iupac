import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Acetate (PubChem CID 175, structure-verified via InChIKey
        # QTBSBXVTEAMEQO-UHFFFAOYSA-M; PubChem's own PIN is "acetate", but
        # this project's `_carboxylic_acid.py` always uses the systematic
        # "ethanoic acid" stem rather than the retained "acetic acid" one,
        # so this module follows the same "ethanoate" convention).
        ("CC(=O)[O-]", "ethanoate"),
        # Formate: same terminal-carbon pattern with no chain carbon at all.
        ("C(=O)[O-]", "methanoate"),
        # A real positional choice for a substituent: the carboxylate carbon
        # still fixes C1 regardless (same mechanic as `_carboxylic_acid.py`).
        ("CC(C)C(=O)[O-]", "2-methylpropanoate"),
        # -oate combined with existing unsaturation support.
        ("C=CC(=O)[O-]", "prop-2-enoate"),
        # -oate + halogen substituent prefix, same mechanic as
        # `_carboxylic_acid.py`'s own halogen coexistence.
        ("ClCC(=O)[O-]", "2-chloroethanoate"),
    ],
)
def test_carboxylate_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_neutral_carboxylic_acid_still_routes_to_acid():
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_dicarboxylate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]C(=O)CC(=O)[O-]")


def test_ring_carboxylate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]C(=O)C1CCCCC1")


def test_extra_charged_atom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH4+].CC(=O)[O-]")


def test_amine_coexisting_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NCC(=O)[O-]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): the
        # carboxylate carbon is always chain-terminal (fixed C1), same
        # pattern as `_carboxylic_acid.py`. PubChem CID 6950478.
        ("CC[C@@H](C)C(=O)[O-]", "(2R)-2-methylbutanoate"),
        ("CC[C@H](C)C(=O)[O-]", "(2S)-2-methylbutanoate"),
    ],
)
def test_carboxylate_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_carboxylate_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)C(=O)[O-]") == "2-methylbutanoate"


def test_phenyl_chain_carboxylate():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring
    # `_carboxylic_acid.py`'s `test_phenyl_chain_carboxylic_acid`): the
    # ring is cited as a "phenyl" substituent prefix. PubChem PUG REST
    # (structure match only; this project keeps its own "ethanoate"
    # naming convention): "2-phenylacetate"/"3-phenylpropanoate".
    assert smiles_to_iupac("c1ccccc1CC(=O)[O-]") == "2-phenylethanoate"
    assert smiles_to_iupac("c1ccccc1CCC(=O)[O-]") == "3-phenylpropanoate"


def test_phenyl_directly_attached_carboxylate_raises():
    # Benzoate-type naming (the carboxylate carbon directly on the ring)
    # is a separate construction, out of scope for this acyclic-
    # chain-parent module, mirroring `_carboxylic_acid.py`'s benzoic-acid
    # rejection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)[O-]")


def test_phenyl_substituted_benzene_ring_carboxylate_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)[O-]")


def test_phenyl_chain_carboxylate_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)[O-]")
