from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._bridged_alicyclic_parent import has_bridged_steroid_name


def test_epoxycholestane():
    # PubChem CID 281912; CAS-style synonym "5,6-epoxycholestane" (also
    # listed as the fully stereo-specified "5.alpha.,6.alpha.-Epoxycholestane").
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC4C5(C3(CCCC5)C)O4)C")
        == "5,6-epoxycholestane"
    )


def test_plain_cholestane_unaffected():
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"


def test_transannular_bridge():
    # Synthetic: androstane with an -O- bridge whose two neighbors are NOT
    # already directly bonded in the bare skeleton (a genuine transannular
    # span, not an ortho-fused three-membered ring). No real registered
    # structure with this exact shape was found on a steroid skeleton, so
    # this proves the generalized mechanism itself with a constructed case
    # rather than a named real compound.
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCC4OC(C2)C1C34C") == "1,11-epoxyandrostane"


def test_more_than_one_bridge_not_claimed():
    # Same skeleton as above plus a second synthetic -O- bridge elsewhere.
    mol = Chem.MolFromSmiles("CC12CC3OC4CCCC5CCC(C1C1CC2O1)C3C54C")
    assert not has_bridged_steroid_name(mol)
