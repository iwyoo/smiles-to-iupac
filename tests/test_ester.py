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


def test_aryl_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)OC")


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


def test_phenyl_substituted_benzene_ring_ester_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)OC")


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
