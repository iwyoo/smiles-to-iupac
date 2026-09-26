import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._amino_acid import _SIDE_CHAIN_TABLE, has_amino_acid_shape
from smiles_to_iupac._common import UnsupportedStructure


def test_side_chain_table_has_no_collisions():
    # 7 non-glycine side chains (glycine has no fragment, handled as a
    # special zero-neighbor case) -- each canonical fragment SMILES must
    # map to exactly one retained name, or two different amino acids
    # would silently get the same name.
    assert len(_SIDE_CHAIN_TABLE) == 7
    assert len(set(_SIDE_CHAIN_TABLE.values())) == 7


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
        # L-/D-serine: PubChem CID 71077 / CID 5951.
        ("C([C@@H](C(=O)O)N)O", "L-serine"),
        ("C([C@H](C(=O)O)N)O", "D-serine"),
        # L-/D-cysteine: PubChem CID 92851 / CID 5862. P-103.1.3.1's
        # exception -- L is CIP 'R' here, not 'S' like every other case
        # in this table.
        ("C([C@@H](C(=O)O)N)S", "L-cysteine"),
        ("C([C@H](C(=O)O)N)S", "D-cysteine"),
        # L-/D-aspartic acid: PubChem CID 5960 / CID 83887. No exception --
        # ordinary S=L/R=D rule, confirmed against PubChem's own IUPACName
        # (an earlier guess at which SMILES tag was L got this backwards).
        ("OC(=O)C[C@H](N)C(=O)O", "L-aspartic acid"),
        ("OC(=O)C[C@@H](N)C(=O)O", "D-aspartic acid"),
        # L-/D-glutamic acid: PubChem CID 33032 / CID 23327.
        ("OC(=O)CC[C@H](N)C(=O)O", "L-glutamic acid"),
        ("OC(=O)CC[C@@H](N)C(=O)O", "D-glutamic acid"),
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


def test_asparagine_still_unsupported():
    # An amide terminus (-CO-NH2), not a second carboxylic acid -- a
    # different guard/module entirely from aspartic acid's shape.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)C[C@@H](N)C(=O)O")


def test_methionine_still_unsupported():
    # A two-carbon side chain to a thioether (not a terminal thiol) is a
    # different shape from cysteine's -CH2-SH -- still falls through to
    # the coexisting-heteroatom guard, deferred to a batch-rollout step.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CSCC[C@@H](N)C(=O)O")


def test_threonine_still_unsupported():
    # A second side-chain stereocenter (and 'allo' complexity), same as
    # isoleucine -- out of scope for this step even though its side chain
    # also carries a hydroxyl.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@@H](O)[C@H](N)C(=O)O")


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
