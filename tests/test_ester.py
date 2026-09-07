import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases (methyl formate, methyl/ethyl acetate, ethyl
        # propanoate), cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName).
        ("C(=O)OC", "methyl methanoate"),
        ("CC(=O)OC", "methyl ethanoate"),
        ("CC(=O)OCC", "ethyl ethanoate"),
        ("CCC(=O)OCC", "ethyl propanoate"),
        # A real positional choice for a substituent on the acyl chain: the
        # ester carbon still fixes C1 regardless.
        ("CC(C)C(=O)OC", "methyl 2-methylpropanoate"),
        # -oate combined with existing unsaturation support, cross-checked
        # against PubChem (methyl crotonate).
        ("CC=CC(=O)OC", "methyl but-2-enoate"),
        # -oate + halogen substituent prefix on the acyl chain, cross-checked
        # against PubChem.
        ("ClCC(=O)OC", "methyl 2-chloroethanoate"),
        # A longer, unbranched alcohol part.
        ("CC(=O)OCCC", "propyl ethanoate"),
        # A plain, unsubstituted benzene ring on the acyl chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/#271's
        # carboxylic-acid/ketone/alcohol chains): the ring is cited as a
        # "phenyl" substituent prefix. Cross-checked against PubChem CID
        # 7643 ("methyl 3-phenylpropanoate").
        ("c1ccccc1CCC(=O)OC", "methyl 3-phenylpropanoate"),
        # Two-carbon acyl chain: PubChem's own name for this SMILES
        # ('ethyl 2-phenylacetate', CID 7590) uses the retained 'acetate'
        # stem, but this module always uses the systematic 'ethanoate' stem
        # (see 'methyl ethanoate' above), so this is an accepted, reviewed
        # result rather than a PubChem-confirmed one -- same policy as the
        # carboxylic-acid module's '2-phenylethanoic acid' (PR #269).
        ("c1ccccc1CC(=O)OCC", "ethyl 2-phenylethanoate"),
        # P-31.1.1.1 'ene' elision before a multiplied 'yne' segment: the
        # elision must fire regardless of the yne count's own multiplying
        # prefix ('diyne' still counts as 'y'-initial), found via real-data
        # testing to be wrongly left unelided ('...dien-4,6-diynoate').
        ("BrC=CCCCCCCC=CCC#CC#CCCC(=O)OC", "methyl 18-bromooctadeca-9,17-dien-4,6-diynoate"),
    ],
)
def test_ester_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_diester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)CC(=O)OC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Structure PubChem-confirmed (CID 7915, 10908, 8038: 'propan-2-yl
        # acetate', 'tert-butyl acetate', '2-methylpropyl acetate'), but
        # this module always uses the systematic 'ethanoate' stem rather
        # than the retained 'acetate' one (see module docstring/
        # test_ester_names above), so the acyl word here is an accepted,
        # reviewed result rather than a PubChem-confirmed one.
        ("CC(=O)OC(C)C", "propan-2-yl ethanoate"),
        ("CC(=O)OC(C)(C)C", "tert-butyl ethanoate"),
        ("CC(=O)OCC(C)C", "2-methylpropyl ethanoate"),
    ],
)
def test_branched_alcohol_part(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_plain_cyclyl_alcohol_ester():
    # PubChem CID 12146 "cyclohexyl acetate"/CID 70273 "cyclopentyl
    # acetate"/CID 61375 "cyclohexyl propanoate" -- this project's own
    # systematic-stem convention gives 'ethanoate' not 'acetate' (see the
    # module's other tests).
    assert smiles_to_iupac("CC(=O)OC1CCCCC1") == "cyclohexyl ethanoate"
    assert smiles_to_iupac("CC(=O)OC1CCCC1") == "cyclopentyl ethanoate"
    assert smiles_to_iupac("CCC(=O)OC1CCCCC1") == "cyclohexyl propanoate"


def test_substituted_cyclyl_alcohol_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OC1CCC(C)CC1")


def test_amine_coexisting_now_supported_via_ester_amine():
    # A coexisting primary amine is now handled by `_ester_amine.py`
    # (P-41 Table 4.1: ester outranks amine) rather than rejected here --
    # see tests/test_ester_amine.py for that module's own coverage.
    assert smiles_to_iupac("NCC(=O)OC") == "methyl 2-aminoethanoate"


def test_benzoate_ester():
    # The acyl part hangs directly off the benzene ring, no intervening
    # chain carbon (P-65.6.3.2's retained-name axis) -- mirrors
    # `_carboxylic_acid.py`'s 'benzoic acid' with the 'benzoate' suffix
    # word. PubChem PUG REST-verified "methyl benzoate" (CID 7150).
    assert smiles_to_iupac("c1ccccc1C(=O)OC") == "methyl benzoate"


def test_benzoate_ester_substituted():
    # PubChem PUG REST-verified "methyl 2-methylbenzoate" (CID 33094).
    assert smiles_to_iupac("Cc1ccccc1C(=O)OC") == "methyl 2-methylbenzoate"


def test_ring_acyl_ester():
    # The acyl part hangs directly off a plain saturated ring, no
    # intervening chain carbon -- mirrors `_carboxylic_acid.py`'s
    # 'cyclohexanecarboxylic acid' with the 'carboxylate' suffix word.
    # PubChem PUG REST-verified "methyl cyclohexanecarboxylate"
    # (CID 20748).
    assert smiles_to_iupac("O=C(OC)C1CCCCC1") == "methyl cyclohexanecarboxylate"


def test_ring_acyl_ester_substituted():
    # PubChem PUG REST-verified "methyl 4-methylcyclohexane-1-
    # carboxylate" (CID 170993).
    assert (
        smiles_to_iupac("O=C(OC)C1CCC(C)CC1") == "methyl 4-methylcyclohexane-1-carboxylate"
    )


def test_ring_acyl_ester_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC)C1CCCC=C1")


def test_ring_acyl_ester_substituent_on_acyl_ring_atom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC)C1(C)CCCCC1")


def test_ring_acyl_chain_ester():
    # The acyl part (R-CO-) chain hangs off a single unbranched chain
    # attached to an otherwise-plain saturated ring -- the saturated
    # counterpart of `_name_phenyl_acyl_ester`. PubChem PUG REST-verified
    # "methyl 2-cyclohexylacetate" (CID 139743) -- this project's
    # systematic 'ethanoate' stem convention applies here too (see
    # `test_ester_names` above).
    assert smiles_to_iupac("O=C(OC)CC1CCCCC1") == "methyl 2-cyclohexylethanoate"


def test_ring_acyl_chain_ester_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC)CC1CCC(C)CC1")


def test_ring_acyl_chain_ester_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(OC)CC1CCCC=C1")


def test_ring_on_alcohol_chain_ester_still_raises():
    # A saturated ring reachable only from the alcohol side (not the acyl
    # side) is still out of scope -- this task only covers rings on the
    # acyl side (directly, or via a chain); a ring on the alcohol side via
    # a chain is a separate, unscoped gap.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)OCC1CCCCC1")


def test_phenyl_ester_oxygen_side():
    # A benzene ring attached directly via the ester oxygen (an aryl
    # ester, P-65.6.3.2) is now supported -- PubChem's own 'phenyl
    # acetate'/'phenyl formate'/'phenyl propanoate' (CID 31229/74626/
    # 12497), this project's systematic-name convention applied.
    assert smiles_to_iupac("CC(=O)Oc1ccccc1") == "phenyl ethanoate"
    assert smiles_to_iupac("O=COc1ccccc1") == "phenyl methanoate"
    assert smiles_to_iupac("CCC(=O)Oc1ccccc1") == "phenyl propanoate"


def test_phenyl_ester_oxygen_side_substituted_ring_raises():
    # A phenol ring with another substituent besides the ester oxygen is
    # still out of scope for this first pass.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)Oc1ccccc1C")


def test_phenyl_acyl_chain_ring_methyl():
    # PubChem 'methyl 2-(2-methylphenyl)acetate' -- a plain methyl ring
    # substituent, like halogens, no longer forces the benzene ring to
    # win the ring-vs-chain parent competition (see
    # test_carboxylic_acid.py's identical methylated-ring cases).
    assert smiles_to_iupac("Cc1ccccc1CC(=O)OC") == "methyl 2-(2-methylphenyl)ethanoate"


def test_phenyl_acyl_chain_ring_methyl_para():
    assert smiles_to_iupac("COC(=O)Cc1ccc(C)cc1") == "methyl 2-(4-methylphenyl)ethanoate"


def test_phenyl_acyl_chain_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CC(=O)OC") == "methyl 2-(4-ethylphenyl)ethanoate"


def test_phenyl_acyl_chain_ring_halogen_and_methyl():
    # PubChem "methyl 2-(2-chloro-5-methylphenyl)acetate" -- halogen and
    # methyl ring substituents mixed on the same ring.
    assert (
        smiles_to_iupac("ClC1=CC=C(C)C=C1CC(=O)OC")
        == "methyl 2-(2-chloro-5-methylphenyl)ethanoate"
    )


def test_phenyl_acyl_chain_ring_halogen():
    # PubChem CID 68493 "methyl 2-(4-chlorophenyl)acetate" -- this
    # project's own systematic-stem convention (see
    # test_carboxylic_acid.py's identical halogenated-ring cases) prefers
    # 'ethanoate' once substituted.
    assert smiles_to_iupac("COC(=O)Cc1ccc(Cl)cc1") == "methyl 2-(4-chlorophenyl)ethanoate"


def test_phenyl_acyl_chain_ring_dihalogen():
    assert smiles_to_iupac("COC(=O)Cc1ccc(Cl)c(Cl)c1") == "methyl 2-(3,4-dichlorophenyl)ethanoate"


def test_phenyl_acyl_chain_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)OC")


def test_acyl_carbon_off_longest_chain_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("COC(=O)C(CCC)CCCC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter on the acyl chain
        # (P-92): the stereodescriptor is cited immediately ahead of the
        # acyl part specifically, not the whole two-word ester name (CIP
        # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID 7156991.
        ("CC[C@@H](C)C(=O)OCC", "ethyl (2R)-2-methylbutanoate"),
        ("CC[C@H](C)C(=O)OCC", "ethyl (2S)-2-methylbutanoate"),
    ],
)
def test_ester_acyl_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ester_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(=O)OCC") == "ethyl 2-methylbutanoate"
