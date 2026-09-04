import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A diester on a shared diol chain is not a multiplied 'oate'
        # suffix -- one ester stays the suffix parent (the longer acyl
        # chain) and the other is demoted to an 'acyloxy' prefix.
        # PubChem PUG REST confirms the structure/locants for all of these
        # (in its own retained-name style, e.g. 'acetyloxy'/'acetate');
        # this project's systematic-name convention (methanoate/ethanoate,
        # not formate/acetate -- see test_ester.py) carries over here too.
        ("CC(=O)OCCOC(=O)C", "2-ethanoyloxyethyl ethanoate"),
        ("CC(=O)OCCOC(=O)CC", "2-ethanoyloxyethyl propanoate"),
        ("CCC(=O)OCCOC(=O)C", "2-ethanoyloxyethyl propanoate"),
        ("CCCC(=O)OCCOC(=O)C", "2-ethanoyloxyethyl butanoate"),
        ("CC(=O)OCCOC(=O)CCC", "2-ethanoyloxyethyl butanoate"),
        ("O=COCCOC(=O)C", "2-methanoyloxyethyl ethanoate"),
        ("CC(=O)OCCCOC(=O)C", "3-ethanoyloxypropyl ethanoate"),
    ],
)
def test_diester_acyloxy_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        # A branched backbone (a real alkyl branch alongside the second
        # ester attachment) needs the full name_branch-style longest-chain
        # machinery, deferred -- see the module docstring.
        "CC(=O)OC(C)COC(=O)C",
        # A halogen on the backbone is out of scope for this narrow first
        # pass.
        "CC(=O)OCC(Cl)OC(=O)C",
    ],
)
def test_diester_acyloxy_out_of_scope(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)
