import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A *differing*-acid diester on a shared diol chain is not a
        # multiplied 'oate' suffix -- one ester stays the suffix parent
        # (the longer acyl chain) and the other is demoted to an
        # 'acyloxy' prefix (P-65.6.3.3.2 method 2, general nomenclature
        # only). PubChem PUG REST confirms the structure/locants for all
        # of these (in its own retained-name style, e.g. 'acetyloxy'/
        # 'acetate'); this project's systematic-name convention
        # (methanoate/ethanoate, not formate/acetate -- see
        # test_ester.py) carries over here too.
        ("CC(=O)OCCOC(=O)CC", "2-ethanoyloxyethyl propanoate"),
        ("CCC(=O)OCCOC(=O)C", "2-ethanoyloxyethyl propanoate"),
        ("CCCC(=O)OCCOC(=O)C", "2-ethanoyloxyethyl butanoate"),
        ("CC(=O)OCCOC(=O)CCC", "2-ethanoyloxyethyl butanoate"),
        ("O=COCCOC(=O)C", "2-methanoyloxyethyl ethanoate"),
        # An *identical*-acid diester instead uses the PIN's diyl+
        # multiplicative form (P-65.6.3.3.3.1, see module docstring for
        # the Blue Book's own 'ethane-1,2-diyl diacetate (PIN)' example).
        ("CC(=O)OCCOC(=O)C", "ethane-1,2-diyl diethanoate"),
        ("CC(=O)OCCCOC(=O)C", "propane-1,3-diyl diethanoate"),
        ("CCC(=O)OCCCCOC(=O)CC", "butane-1,4-diyl dipropanoate"),
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
