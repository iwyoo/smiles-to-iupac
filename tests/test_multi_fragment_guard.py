import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles",
    [
        # Modules that allow a ring/chain heteroatom (O/N/S/...) can't reuse
        # `_common.validate_atoms_and_bonds`'s carbon-plus-halogen-only
        # multi-fragment check, so each wrote its own atom-validation loop --
        # every one of those loops originally forgot the multi-fragment
        # check itself, silently dropping the extra fragment instead of
        # rejecting (found via `smiles-to-iupac-realdata-test`'s pubchem
        # diff: "CC.CC1CCC12CNC2" used to come out as
        # "5-methyl-2-azaspiro[3.3]heptane", the ethane fragment vanishing).
        "CC.CC1CCC12CNC2",  # _spiro_heteroatom.py
        "CC.CCOCC",  # _ether.py
        "CC.CCC[N+](=O)[O-]",  # _nitro.py
        "CC.CCS(C)(=O)=O",  # _sulfone.py
        "CC.O=[N+]([O-])c1ccccc1",  # _nitro.py, aromatic-ring path
        "CC.Cc1ccncc1",  # _hetero_monocyclic.py substituent path
    ],
)
def test_multi_fragment_rejected_instead_of_silently_dropped(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
