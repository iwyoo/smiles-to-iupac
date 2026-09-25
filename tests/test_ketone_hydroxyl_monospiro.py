import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 55263426 -- spiro[3.5]nonane, ketone and hydroxyl
        # each equidistant from both spiro neighbors in their own ring, so
        # the numbering direction is unambiguous and matches PubChem's own
        # computed name directly.
        ("OC1CCC2(CC1)CC(=O)C2", "7-hydroxyspiro[3.5]nonan-2-one"),
        # The remaining cases below are real registered PubChem structures
        # (CID 153873947, 165085875, 157976901, 145133375, 114821462) whose
        # own PubChem-computed "IUPACName" property does NOT take the
        # lowest locant available to the ketone suffix (or, when the
        # suffix locant is tied either direction, to the coexisting
        # hydroxyl) -- confirmed by hand for each below by walking both
        # spiro-ring traversal directions. Per P-31.1.4.3 (lowest locants
        # to the principal characteristic group suffix first, same rule
        # already proven for the von Baeyer branch, #1029), the name
        # asserted here is the one this project's `_candidate_key`
        # ordering also independently derives, not PubChem's own listed
        # string. E.g. CID 153873947's own PubChem name is
        # "8-hydroxyspiro[4.5]decan-4-one", but the ketone carbon sits
        # directly adjacent to the spiro atom in the smaller ring, so
        # numbering that ring in the other direction (equally valid under
        # P-24.2.1, which only fixes numbering to start in the smaller
        # ring, not which of its two spiro-adjacent atoms is C-1) gives
        # the ketone locant 1 instead of 4 -- strictly lower, and
        # therefore preferred.
        ("OC1CCC2(CC1)CCCC2=O", "8-hydroxyspiro[4.5]decan-1-one"),
        ("O=C1CCC2(CCC(O)CC2)C1", "8-hydroxyspiro[4.5]decan-2-one"),
        ("OC1CCC2(CC1)CCC2=O", "7-hydroxyspiro[3.5]nonan-1-one"),
        ("OC1CCCCC12CC(=O)C2", "5-hydroxyspiro[3.5]nonan-2-one"),
        # CID 114821462 -- ketone and hydroxyl on adjacent ring atoms; the
        # ketone's own locant (8) is tied between both directions (it sits
        # equidistant from each spiro neighbor), so the hydroxyl's own
        # lowest achievable locant (7, not PubChem's 9) breaks the tie.
        ("O=C1CCC2(CC1O)CCCC2", "7-hydroxyspiro[4.5]decan-8-one"),
    ],
)
def test_hydroxy_monospiro_ketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
