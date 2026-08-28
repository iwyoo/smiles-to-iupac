import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName) and matching the
        # Blue Book's own worked examples (P6a.pdf, P-66.3.1.1/.1.2.3):
        # 'butanehydrazide (PIN)', 'pentanehydrazide (PIN)'. The
        # hydrazide carbon is always C1, and its own locant is never
        # cited (P-14.3.3), same as `_amide.py`.
        ("CCCC(=O)NN", "butanehydrazide"),
        ("CCCCC(=O)NN", "pentanehydrazide"),
        ("CCC(=O)NN", "propanehydrazide"),
        # -hydrazide + halogen substituent prefix.
        ("ClCCC(=O)NN", "3-chloropropanehydrazide"),
        # -hydrazide combined with existing unsaturation support.
        ("CC=CC(=O)NN", "but-2-enehydrazide"),
        # 'hydrazide' outranks 'ol' in Table 3.3, so a coexisting
        # standalone -OH is cited as the 'hydroxy' prefix, same as
        # `_amide.py`.
        ("OCCC(=O)NN", "3-hydroxypropanehydrazide"),
    ],
)
def test_hydrazide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_methane_hydrazide_raises():
    # P-66.3.1.2.1: 'formohydrazide' (a retained name), not the
    # systematic 'methanehydrazide', is the actual PIN here -- out of
    # scope for this module's systematic-only first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C(=O)NN")


def test_ethane_hydrazide_raises():
    # Same reasoning: 'acetohydrazide' is the PIN, not 'ethanehydrazide'.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NN")


def test_n_substituted_hydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCC(=O)NNC")


def test_diacylhydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCC(=O)NNC(=O)C")


def test_amide_not_misnamed_as_hydrazide():
    # A plain primary amide (`_amide.py`) has only one nitrogen and must
    # not be routed here.
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"
