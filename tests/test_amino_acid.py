import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._amino_acid import has_amino_acid_shape
from smiles_to_iupac._common import UnsupportedStructure


def test_glycine():
    # No stereocenter (C-2 bears two hydrogens) -- PubChem CID 750.
    assert smiles_to_iupac("NCC(=O)O") == "glycine"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # L-/D-alanine: PubChem CID 5950 / CID 71080.
        ("C[C@@H](C(=O)O)N", "L-alanine"),
        ("C[C@H](C(=O)O)N", "D-alanine"),
        # L-/D-valine: PubChem CID 6287 / CID 439610.
        ("CC(C)[C@@H](C(=O)O)N", "L-valine"),
        ("CC(C)[C@H](C(=O)O)N", "D-valine"),
        # L-/D-leucine: PubChem CID 6106 / CID 439524.
        ("CC(C)C[C@@H](C(=O)O)N", "L-leucine"),
        ("CC(C)C[C@H](C(=O)O)N", "D-leucine"),
    ],
)
def test_common_amino_acid_ld(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_no_ld_prefix():
    # No stereo marker at all -- no L-/D- prefix, matching this project's
    # long-standing convention for an unspecified stereocenter.
    assert smiles_to_iupac("CC(N)C(=O)O") == "alanine"


def test_isoleucine_still_generic():
    # A second side-chain stereocenter (and the 'allo' prefix it implies)
    # is out of scope for this step -- still falls through to the generic
    # CIP-labeled path unchanged.
    assert smiles_to_iupac("CC[C@H](C)[C@@H](N)C(=O)O") == "(2R,3S)-2-amino-3-methylpentanoic acid"


def test_serine_still_unsupported():
    # A heteroatom-bearing side chain (hydroxyl) trips an unrelated
    # coexisting-oxygen guard, unaffected by this module -- deferred to a
    # batch-rollout step.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC[C@@H](N)C(=O)O")


def test_phenylalanine_still_unsupported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("N[C@@H](Cc1ccccc1)C(=O)O")


def test_two_amino_acid_shaped_stereocenters_not_matched():
    # A molecule with the alpha-carbon's exact H-count/degree shape but a
    # second, deeper stereocenter beyond what a plain alanine/valine/
    # leucine side chain has must not be misrecognized -- covered by
    # `test_isoleucine_still_generic` above; this asserts the negative
    # detector itself doesn't fire for it.
    mol = Chem.MolFromSmiles("CC[C@H](C)[C@@H](N)C(=O)O")
    assert has_amino_acid_shape(mol) is False
