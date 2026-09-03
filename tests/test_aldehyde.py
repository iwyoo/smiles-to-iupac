import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). Note the 'al'
        # locant is never cited (see module docstring): the -CHO carbon is
        # always the chain terminus, so 'propan-1-al' is never written.
        ("CC=O", "ethanal"),
        ("CCC=O", "propanal"),
        ("CC(C)C=O", "2-methylpropanal"),
        # Dialdehyde: multiplying prefix, still no locants at all (both ends
        # of the chain, unambiguous), cross-checked against PubChem.
        ("O=CCCCC=O", "pentanedial"),
        # -al combined with existing unsaturation support: the 'ene' locant
        # is still cited (counted from the -CHO end, suffix priority,
        # P-44.4.1.8), cross-checked against PubChem.
        ("C=CCCC=O", "pent-4-enal"),
        # -al combined with a halogen substituent prefix, cross-checked
        # against PubChem.
        ("ClCCC=O", "3-chloropropanal"),
        # A plain, unsubstituted benzene ring on the chain (P-2/P-3
        # aromatic-ring-substituent extension, mirroring PR #269/#270/
        # #271/#272's carboxylic-acid/ketone/alcohol/ester chains): the
        # ring is cited as a "phenyl" substituent prefix. Cross-checked
        # against PubChem CID 7707 ("3-phenylpropanal").
        ("c1ccccc1CCC=O", "3-phenylpropanal"),
        # Two-carbon chain: PubChem's own name for this SMILES
        # ('2-phenylacetaldehyde', CID 998) uses the retained
        # 'acetaldehyde' stem, but this module always uses the systematic
        # 'ethanal' stem (see 'ethanal' above), so this is an accepted,
        # reviewed result rather than a PubChem-confirmed one -- same
        # policy as the carboxylic-acid module's '2-phenylethanoic acid'
        # (PR #269).
        ("c1ccccc1CC=O", "2-phenylethanal"),
    ],
)
def test_aldehyde_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formaldehyde_raises():
    # A carbonyl carbon with zero carbon neighbors (formaldehyde) is out of
    # scope for this module (see module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=O")


def test_ketene_raises():
    # A ketene's carbonyl carbon (C=C=O) is itself doubly bonded to its
    # carbon neighbor (a cumulated double bond, sp-hybridized) -- it is not
    # an aldehyde's -CHO and must not be silently misnamed as one.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCCCCCCCCCCCC=C=O")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC=C=O")


def test_carboxylic_acid_not_misread_as_aldehyde():
    # A carboxylic acid is routed to the dedicated carboxylic-acid module
    # (see test_carboxylic_acid.py) instead of falling through here.
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_aryl_aldehyde_raises():
    # Benzaldehyde: an aromatic ring elsewhere in the molecule is out of
    # scope for this module (separate, in-progress aromatic-ring module's
    # territory).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=Cc1ccccc1")


def test_phenyl_substituted_benzene_ring_aldehyde_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC=O")


def test_phenyl_chain_aldehyde_with_hydroxyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCc1ccccc1CC=O")


def test_phenyl_chain_aldehyde_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC=O")


def test_alcohol_aldehyde_mix_names_hydroxy_prefix():
    # 'al' outranks 'ol' in Table 3.3, so a coexisting -OH is cited as the
    # 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OCC=O") == "2-hydroxyethanal"


def test_aldehyde_alcohol_mix_on_longer_chain():
    assert smiles_to_iupac("OCCCC=O") == "4-hydroxybutanal"


def test_aldehyde_enol_mix_raises():
    # A hydroxyl on a C=C carbon (an enol) is a tautomer of a more senior
    # carbonyl form and out of scope, same as `_alcohol.py`'s own enol check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC=O")


def test_ring_aldehyde_raises():
    # -CHO on a ring is the 'carbaldehyde' suffix (P-33.3.1.2), a different
    # naming pattern this module deliberately excludes (see module
    # docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CC1CCCCC1")


def test_bicyclic_carbon_skeleton_with_stray_aldehyde_raises():
    # A -CHO group whose carbon is not on any single longest chain of the
    # molecule (here, a branch off a longer chain) is out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCCC(C=O)CCCC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_carboxylic_acid.py` (CIP computed entirely by RDKit's
        # `rdCIPLabeler`, not reimplemented here). PubChem CID 76956407.
        ("C[C@@H](Cl)C=O", "(2R)-2-chloropropanal"),
        ("C[C@H](Cl)C=O", "(2S)-2-chloropropanal"),
        # Two specified stereocenters, ascending-locant group (P-91.3).
        # PubChem CID 92160220, name matches exactly.
        ("C[C@H](Cl)[C@H](Cl)C=O", "(2R,3S)-2,3-dichlorobutanal"),
    ],
)
def test_aldehyde_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_aldehyde_stereocenter_with_hydroxy_coexistence():
    # Two stereocenters coexisting with a standalone hydroxyl (already
    # supported by this module as the 'hydroxy' prefix). PubChem doesn't
    # have this exact stereoisomer registered (CID 0 for every @/@@
    # combination tried), so only the non-stereo parent structure is
    # verified (PubChem CID 24973902, "2-chloro-3-hydroxybutanal");
    # the R/S computation itself is RDKit's `rdCIPLabeler`, already
    # verified elsewhere (`_alcohol.py`/`_carboxylic_acid.py`).
    assert smiles_to_iupac("C[C@H](O)[C@H](Cl)C=O") == "(2S,3S)-2-chloro-3-hydroxybutanal"


def test_aldehyde_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CC(Cl)C=O") == "2-chloropropanal"


def test_aldehyde_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)C(Cl)C=O")
