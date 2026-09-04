import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Symmetric ethers: P-14.3.4.2(b)'s omitted-locant convention for a
        # homogeneous two-carbon chain with exactly one substituent applies
        # here too, cross-checked against PubChem.
        ("COC", "methoxymethane"),
        ("CCOCC", "ethoxyethane"),
        # Asymmetric ethers: the longer chain is the parent (P-44.3), the
        # shorter side becomes the 'oxy' prefix (P-63.2.1). Cross-checked
        # against PubChem.
        ("COCC", "methoxyethane"),
        ("COCCC", "1-methoxypropane"),
        ("CCOCCC", "1-ethoxypropane"),
        # Longer unbranched alkoxy: no contracted retained name past
        # 'butoxy' (P-63.2.2.1.1), so the plain alkyl name plus 'oxy' is
        # used unchanged.
        ("COCCCCC", "1-methoxypentane"),
        # A branched parent chain is unaffected by this module's R'-only
        # restriction; only the (unbranched) shorter side is the prefix.
        ("CC(C)OCC", "2-ethoxypropane"),
        # A branched R' (the shorter, prefix side): P-63.2.2.1.1 encloses
        # only R' in parentheses, with 'oxy' outside them. Structure
        # verified against PubChem: CCCCOC(C)C -> CID 137240
        # ("1-propan-2-yloxybutane" -- PubChem elides the parentheses this
        # project's usual compound-substituent formatting keeps).
        ("CCCCOC(C)C", "1-(propan-2-yl)oxybutane"),
        # Longer R' side also branched: PubChem CID 28488 confirms the
        # structure ("2-methyl-2-propan-2-yloxypropane").
        ("CC(C)OC(C)(C)C", "2-(propan-2-yl)oxy-2-methylpropane"),
    ],
)
def test_ether(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Fully symmetric tie (PubChem CID-confirmed 'CC(C)OC(C)C' ->
        # '2-propan-2-yloxypropane'): either side works as parent, so any
        # tie-break gives the same answer.
        ("CC(C)OC(C)C", "2-(propan-2-yl)oxypropane"),
        # Tied total carbon count (4 each) and tied longest-chain length too
        # (both an isobutyl and a tert-butyl arm reduce to a
        # 2-methylpropane skeleton, chain length 3) -- the decider is the
        # locant set each side gives as parent: isobutyl-as-parent {1,2}
        # vs. tert-butyl-as-parent {2,2}, so isobutyl wins. PubChem
        # CID-confirmed: 'CC(C)COC(C)(C)C' ->
        # '2-methyl-1-[(2-methylpropan-2-yl)oxy]propane' (this project
        # retains 'tert-butyl' rather than PubChem's systematic
        # '2-methylpropan-2-yl', per `_substituents.py`'s existing
        # P-29.6.1 special case, and cites prefixes alphabetically --
        # 'butyloxy' before 'methyl' -- unaffected by this task).
        ("CC(C)COC(C)(C)C", "1-tert-butyloxy-2-methylpropane"),
        # Tied total carbon count (5 each) but different longest achievable
        # chain: a neopentyl arm's 5 carbons max out at chain length 3
        # around its quaternary carbon, while an isopentyl arm's 5 carbons
        # reach chain length 4 -- isopentyl wins outright, no locant-set
        # comparison needed. PubChem CID-confirmed: 'CC(C)(C)COCCC(C)C' ->
        # '1-(2,2-dimethylpropoxy)-3-methylbutane'.
        ("CC(C)(C)COCCC(C)C", "1-(2,2-dimethylpropyl)oxy-3-methylbutane"),
        # Tied total carbon count (4 each) but different longest achievable
        # chain: an isobutyl arm's own chain (used as parent) tops out at
        # length 3 (its 4th carbon has to be a methyl branch), while a
        # sec-butyl arm's own chain is a plain, unbranched 4-long butane --
        # sec-butyl wins outright, no locant-set comparison needed. PubChem
        # CID-confirmed: 'CC(C)COC(C)CC' -> '2-(2-methylpropoxy)butane'.
        ("CC(C)COC(C)CC", "2-(2-methylpropyl)oxybutane"),
    ],
)
def test_ether_both_sides_branched_and_tied(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem-verified (CID 57895829/54175646): the parent (longer)
        # chain carries the stereocenter, so the prefix mirrors
        # `_acetal.py`'s P-91.3 mechanism via `winning_chain_from_carbon_graph`.
        ("CC[C@H](C)OCC", "(2S)-2-ethoxybutane"),
        ("CC[C@@H](C)OCC", "(2R)-2-ethoxybutane"),
    ],
)
def test_stereocenter_on_parent_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unspecified_stereocenter_ignored():
    assert smiles_to_iupac("CCC(C)OCC") == "2-ethoxybutane"


def test_stereocenter_on_substituent_branch_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCO[C@H](C)CC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Direct ring-oxygen bond (P-63.2.2.1.1's own 'alkoxybenzene'
        # shape). Structure PubChem-confirmed: CID 7500 'c1ccccc1OCC' ->
        # 'ethoxybenzene', CID 7519 'c1ccccc1OC' -> 'methoxybenzene'.
        ("c1ccccc1OCC", "ethoxybenzene"),
        ("c1ccccc1OC", "methoxybenzene"),
        ("c1ccccc1OCCC", "propoxybenzene"),
        # Direct ring-oxygen bond with a branched R': the same
        # '(...)oxy' parenthesization as the plain two-chain path above
        # (P-63.2.2.1.1's worked example), even though PubChem's own
        # auto-generated name for the same structure omits the
        # parentheses ('c1ccccc1OC(C)C' -> 'propan-2-yloxybenzene').
        ("c1ccccc1OC(C)C", "(propan-2-yl)oxybenzene"),
        # Chain spacer between the ring and the ether oxygen -- the whole
        # branch is parenthesized when compound, matching this project's
        # existing benzene-ring-chain convention (`_nitro.py`/
        # `_azide.py`'s own PubChem-vs-PIN discrepancy, see
        # `test_nitro.py`) rather than PubChem's own un-parenthesized
        # auto-names ('ethoxymethylbenzene', '2-ethoxyethylbenzene').
        ("c1ccccc1COCC", "(ethoxymethyl)benzene"),
        ("c1ccccc1CCOCC", "(2-ethoxyethyl)benzene"),
        # A branched R' behind a chain spacer: the '(...)oxy' term
        # already carries round brackets, so the outer compound-branch
        # wrap escalates to square brackets instead of nesting round ones
        # (mirrors the plain two-chain path's own bracket-escalation
        # example in the module docstring). Structure PubChem-confirmed
        # ('c1ccccc1COC(C)C' -> 'propan-2-yloxymethylbenzene').
        ("c1ccccc1COC(C)C", "[(propan-2-yl)oxymethyl]benzene"),
    ],
)
def test_benzene_ring_parent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_benzene_ring_multiple_substituents_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COc1ccccc1OC")


def test_benzene_ring_stereocenter_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1O[C@H](C)CC")
