import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A *differing*-acid diester on a shared diol chain now uses the
        # P-65.6.3.3.3.2 method-1 PIN: the diyl backbone followed by each
        # differing acid name in alphanumerical order, no locants (the
        # backbone's own two-fold symmetry means either numbering
        # direction describes the identical molecule, mirroring the Blue
        # Book's own 'methylene acetate formate (PIN)'/'1,4-phenylene
        # acetate dichloroacetate (PIN)' examples, neither of which cites
        # an acid locant either). This project's systematic-name
        # convention (methanoate/ethanoate, not formate/acetate) carries
        # over here too.
        ("CC(=O)OCCOC(=O)CC", "ethane-1,2-diyl ethanoate propanoate"),
        ("CCC(=O)OCCOC(=O)C", "ethane-1,2-diyl ethanoate propanoate"),
        ("CCCC(=O)OCCOC(=O)C", "ethane-1,2-diyl butanoate ethanoate"),
        ("CC(=O)OCCOC(=O)CCC", "ethane-1,2-diyl butanoate ethanoate"),
        ("O=COCCOC(=O)C", "ethane-1,2-diyl ethanoate methanoate"),
        # An *identical*-acid diester instead uses the PIN's diyl+
        # multiplicative form (P-65.6.3.3.3.1, see module docstring for
        # the Blue Book's own 'ethane-1,2-diyl diacetate (PIN)' example).
        ("CC(=O)OCCOC(=O)C", "ethane-1,2-diyl diethanoate"),
        ("CC(=O)OCCCOC(=O)C", "propane-1,3-diyl diethanoate"),
        ("CCC(=O)OCCCCOC(=O)CC", "butane-1,4-diyl dipropanoate"),
        # An interior alcohol carbon with a plain terminal carbon beyond
        # it (e.g. propanediol/butanediol diacetate's own methyl end) is
        # still a genuinely unbranched carbon skeleton -- both esters plus
        # every plain chain end are found via `carbon_adjacency`'s own
        # connected-component/leaf search, not a naive single-direction
        # walk from one ester. PubChem CID 12198/66193.
        ("CC(OC(C)=O)COC(C)=O", "propane-1,2-diyl diethanoate"),
        ("CC(=O)OC(C)COC(=O)C", "propane-1,2-diyl diethanoate"),
        ("CC(OC(C)=O)C(C)OC(C)=O", "butane-2,3-diyl diethanoate"),
        # A genuine one-carbon alkyl branch (not reducible to a straight
        # chain) hanging off a non-terminal backbone carbon, for the
        # exactly-two-identical-esters case only.
        ("CC(COC(C)=O)COC(C)=O", "2-methylpropane-1,3-diyl diethanoate"),
    ],
)
def test_diester_acyloxy_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)OCC(Cl)OC(=O)C", "1-chloroethane-1,2-diyl diethanoate"),
        ("CC(=O)OCC(COC(C)=O)(COC(C)=O)COC(C)=O", "2,2-bis[(ethanoyloxy)methyl]propane-1,3-diyl diethanoate"),
        ("CC(COC(C)=O)COC(=O)CC", "2-methylpropane-1,3-diyl ethanoate propanoate"),
        ("CC(=O)OCC(O)COC(C)=O", "2-hydroxypropane-1,3-diyl diethanoate"),
        ("CC(=O)OCC=CCOC(C)=O", "but-2-ene-1,4-diyl diethanoate"),
        ("CC(=O)OCOC(C)=O", "methylene diethanoate"),
        ("CC(=O)OC(OC(C)=O)OC(C)=O", "methanetriyl triethanoate"),
    ],
)
def test_diester_acyloxy_generalized_backbones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)OCCOCCOC(C)=O",
        "CC(=O)OCC(=C)COC(C)=O",
    ],
)
def test_diester_acyloxy_out_of_scope(smiles):
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(smiles)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Three or more *identical*-acid esters on a single plain
        # unbranched chain generalize the same diyl+multiplicative
        # mechanism to triyl/tetrayl, per the Blue Book's own worked
        # example 'propane-1,2,3-triyl triacetate (PIN)'. Real PubChem
        # structures confirmed (glyceryl triacetate CID 5541, tripropionin
        # CID 8763, the four-ester case CID 539117) -- PubChem's own
        # auto-generated name uses the substitutive acyloxy-prefix style
        # instead (the same gap this module's 2-ester identical-acid path
        # already departs from), not the PIN form checked here.
        ("CC(=O)OCC(OC(C)=O)COC(C)=O", "propane-1,2,3-triyl triethanoate"),
        ("CCC(=O)OCC(OC(=O)CC)COC(=O)CC", "propane-1,2,3-triyl tripropanoate"),
        ("CC(=O)OCC(OC(C)=O)C(OC(C)=O)COC(C)=O", "butane-1,2,3,4-tetrayl tetraethanoate"),
        # A non-consecutive locant set (an interior position left plain)
        # is still a single unbranched chain, so it's in scope too.
        ("CC(=O)OCCC(OC(C)=O)COC(C)=O", "butane-1,2,4-triyl triethanoate"),
    ],
)
def test_diester_acyloxy_n_ary_identical_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Three or more esters, some identical/some differing, use the
        # same method-1 mechanism generalized to N>=3: each acid group
        # gets its own locant set (no longer omittable once a backbone
        # position is constitutionally distinct from the others, unlike
        # the 2-ester case above) -- matches the Blue Book's own
        # 'propane-1,2,3-triyl 1,2-diacetate 3-propanoate (PIN)' example
        # exactly (glyceryl diacetate monopropanoate).
        ("CC(=O)OCC(OC(C)=O)COC(=O)CC", "propane-1,2,3-triyl 1,2-diethanoate 3-propanoate"),
        # All three acids differ -- mirrors the Blue Book's own
        # 'propane-1,2,3-triyl 2-acetate 1-hexadecanoate
        # 3-[(9Z)-octadec-9-enoate] (PIN)' shape exactly (three singleton
        # groups, every one still gets its own locant).
        ("O=COCC(OC(C)=O)COC(=O)CC", "propane-1,2,3-triyl 2-ethanoate 1-methanoate 3-propanoate"),
    ],
)
def test_diester_acyloxy_n_ary_differing_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
