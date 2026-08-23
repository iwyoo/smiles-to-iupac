import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # The Blue Book's own three P-24.2.2 worked examples (built from
        # scratch to match the ring shapes shown in the PDF, not copied from
        # a database). Independently cross-checked against PubChem: each
        # SMILES below has a CID whose recorded IUPACName property matches
        # (modulo PubChem's plain-text rendering dropping the superscript
        # formatting, e.g. "dispiro[3.2.37.24]dodecane" for
        # "dispiro[3.2.3^7.2^4]dodecane") -- CID 147569901, 123755929, and
        # 58501616 respectively.
        ("C1CCC12CCC3(CC2)CCC3", "dispiro[3.2.3^7.2^4]dodecane"),
        ("C1CC12CCC3(CC2)CCC3", "dispiro[2.2.3^6.2^3]undecane"),
        ("C1CC12CCC4(CC2)CCC5(CC4)CC5", "trispiro[2.2.2.2^9.2^6.2^3]pentadecane"),
        # P-24.2.2.1's first worked example: a genuine (non-tied) case where
        # the internal ring's two arcs differ in length (1 vs 2) -- proves
        # 'always by the shorter path' is enforced as a hard structural rule,
        # not merely preferred by the substituent-locant tie-break. Also
        # cross-checked against PubChem (CID 147810604).
        ("C1CCCC12CC3(CC2)CCCC3", "dispiro[4.1.4^7.2^5]tridecane"),
        # P-24.2.2.1's second worked example: same rule, larger rings, only
        # verified by hand-derivation from the Blue Book text (no PubChem CID
        # found for this exact SMILES).
        ("C1CCCCC12CC3(CC2)CCCCCCCC3", "dispiro[5.1.8^8.2^6]octadecane"),
        # P-24.2.2.2's worked example: both terminal rings tie in size, so
        # both overall chain directions are structurally valid; the
        # descriptor-number-sequence tie-break (lowest at first point of
        # difference) must reject the reversed-direction reading. Only
        # verified by hand-derivation (no PubChem CID found).
        ("C1CC12CCC4(CCC2)CCC5(CC4)CC5", "trispiro[2.2.2.2^9.2^6.3^3]hexadecane"),
        # Five rings (tetraspiro, four spiro atoms): proves the walk/
        # numbering algorithm is genuinely parameterized by ring count, not
        # copy-pasted per dispiro/trispiro special case. Only verified by
        # hand-derivation (no PubChem CID found for this constructed SMILES).
        (
            "C1CC12CC3(C2)CC4(C3)CC5(C4)CC5",
            "tetraspiro[2.1.1.1.2^9.1^7.1^5.1^3]tetradecane",
        ),
        # halogen substituent on the starting terminal ring, threaded through
        # the same way _spiro.py/_cyclic.py do it.
        ("C1(Cl)CCC12CCC3(CC2)CCC3", "1-chlorodispiro[3.2.3^7.2^4]dodecane"),
    ],
)
def test_smiles_to_iupac_linear_polyspiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsaturated_linear_polyspiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC12CCC3(CC2)CCC3")


def test_heteroatom_linear_polyspiro_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC12CCC3(CC2)CCO3")


def test_branched_polyspiro_raises():
    # a single carbon shared by three rings at once (P-24.2.3, branched
    # polyspiro): out of scope, so this must not be mistaken for a linear
    # chain and must still raise.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C123CCC1CCC2CCC3")


def test_monospiro_two_rings_is_unaffected():
    # exactly two rings still goes through _spiro.py's monospiro path, not
    # this module (find_linear_polyspiro_chain requires >= 3 rings).
    assert smiles_to_iupac("C1CCCC12CCCCC2") == "spiro[4.5]decane"
