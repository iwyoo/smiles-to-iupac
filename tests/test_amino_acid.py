import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._amino_acid import _SIDE_CHAIN_TABLE, has_amino_acid_shape
from smiles_to_iupac._common import UnsupportedStructure


def test_side_chain_table_has_no_collisions():
    # 15 non-glycine side chains (glycine has no fragment, handled as a
    # special zero-neighbor case) -- each canonical fragment SMILES must
    # map to exactly one retained name, or two different amino acids
    # would silently get the same name.
    assert len(_SIDE_CHAIN_TABLE) == 15
    assert len(set(_SIDE_CHAIN_TABLE.values())) == 15


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



@pytest.mark.parametrize(
    "smiles,expected",
    [
        # L-/D-asparagine: PubChem CID 6267 / CID 439600.
        ("NC(=O)C[C@H](N)C(=O)O", "L-asparagine"),
        ("NC(=O)C[C@@H](N)C(=O)O", "D-asparagine"),
        # L-/D-glutamine: PubChem CID 5961 / CID 145815.
        ("NC(=O)CC[C@H](N)C(=O)O", "L-glutamine"),
        ("NC(=O)CC[C@@H](N)C(=O)O", "D-glutamine"),
        # L-/D-methionine: PubChem CID 6137 / CID 84815.
        ("CSCC[C@H](N)C(=O)O", "L-methionine"),
        ("CSCC[C@@H](N)C(=O)O", "D-methionine"),
        # L-/D-lysine: PubChem CID 5962 / CID 57449. Its own extra side-chain
        # amine no longer trips the whole-molecule "exactly one amine" guard
        # -- `_match` anchors the backbone amine directly instead.
        ("NCCCC[C@H](N)C(=O)O", "L-lysine"),
        ("NCCCC[C@@H](N)C(=O)O", "D-lysine"),
        # L-/D-phenylalanine: PubChem CID 6140 / CID 71567. Its ring-bearing
        # side chain no longer trips the whole-molecule "no rings" guard.
        ("N[C@@H](Cc1ccccc1)C(=O)O", "L-phenylalanine"),
        ("N[C@H](Cc1ccccc1)C(=O)O", "D-phenylalanine"),
        # L-/D-tyrosine: PubChem CID 6057 / CID 71098.
        ("N[C@@H](Cc1ccc(O)cc1)C(=O)O", "L-tyrosine"),
        ("N[C@H](Cc1ccc(O)cc1)C(=O)O", "D-tyrosine"),
        # L-/D-tryptophan: PubChem CID 6305 / CID 9060.
        ("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "L-tryptophan"),
        ("N[C@H](Cc1c[nH]c2ccccc12)C(=O)O", "D-tryptophan"),
    ],
)
def test_remaining_single_stereocenter_amino_acid_ld(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_arginine_unspecified_stereo_resolves():
    # The retained name resolves for the unspecified-stereo case (correct,
    # no L-/D- prefix per this project's convention), but a *specified*
    # stereocenter currently still fails -- see
    # test_arginine_specified_stereo_still_unsupported below, a separate,
    # shared limitation unrelated to this module's own side-chain table.
    assert smiles_to_iupac("NC(=N)NCCCC(N)C(=O)O") == "arginine"


def test_arginine_specified_stereo_still_unsupported():
    # `_common.specified_stereocenters` rejects any C=N double bond in the
    # whole molecule, including a non-stereogenic symmetric one like the
    # guanidino group's own C(=N)(N)N -- a separate, shared limitation
    # (affects every module that calls this helper, not just amino acids),
    # out of scope for this step. Filed as a follow-up rather than fixed
    # here to avoid widening this PR into `_common.py`'s general
    # double-bond-stereo detection.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=N)NCCC[C@H](N)C(=O)O")


def test_threonine_still_unsupported():
    # A second side-chain stereocenter (and 'allo' complexity), same as
    # isoleucine -- out of scope for this step even though its side chain
    # also carries a hydroxyl.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@@H](O)[C@H](N)C(=O)O")


def test_two_amino_acid_shaped_stereocenters_not_matched():
    # A molecule with the alpha-carbon's exact H-count/degree shape but a
    # second, deeper stereocenter beyond what a plain alanine/valine/
    # leucine side chain has must not be misrecognized -- covered by
    # `test_isoleucine_still_generic` above; this asserts the negative
    # detector itself doesn't fire for it.
    mol = Chem.MolFromSmiles("CC[C@H](C)[C@@H](N)C(=O)O")
    assert has_amino_acid_shape(mol) is False
