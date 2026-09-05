import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The -COOH carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(=O)O", "methanoic acid"),
        ("CC(=O)O", "ethanoic acid"),
        ("CCC(=O)O", "propanoic acid"),
        # A real positional choice for a substituent, cross-checked against
        # PubChem: the -COOH carbon fixes C1 regardless.
        ("CC(C)CC(=O)O", "3-methylbutanoic acid"),
        # Dicarboxylic acid: both chain ends are -COOH carbons, multiplying
        # prefix + no locants cited, cross-checked against PubChem (adipic
        # acid's PIN).
        ("OC(=O)CCCCC(=O)O", "hexanedioic acid"),
        # -oic acid combined with existing unsaturation support, cross-checked
        # against PubChem (crotonic acid's PIN).
        ("CC=CC(=O)O", "but-2-enoic acid"),
        # Unsaturated dicarboxylic acid, cross-checked against PubChem
        # (fumaric/maleic acid's PIN).
        ("OC(=O)C=CC(=O)O", "but-2-enedioic acid"),
        # -oic acid + halogen substituent prefix, cross-checked against
        # PubChem.
        ("OC(=O)CCCl", "3-chloropropanoic acid"),
        # Diene acid: two multiplied 'ene' locants ahead of the acid suffix,
        # cross-checked against PubChem (sorbic acid's PIN).
        ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid"),
    ],
)
def test_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ether_raises():
    # An ether coexisting with a carboxylic acid (as opposed to a plain
    # ether on its own, now handled by the separate ether module) is still
    # out of scope for this module.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCOCC(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P-65.1.2.2.2: a -COOH directly attached to an otherwise
        # unsubstituted saturated monocyclic ring is named with the
        # 'carboxylic acid' suffix on the ring parent hydride, cross-checked
        # against PubChem (e.g. cyclohexanecarboxylic acid's PIN). The sole
        # substituent's ring locant is omitted (P-14.3.3), same as
        # _sulfonic_acid.py's analogous single-substituent ring case.
        ("C1CCC(CC1)C(=O)O", "cyclohexanecarboxylic acid"),
        ("C1CCCC1C(=O)O", "cyclopentanecarboxylic acid"),
        ("C1CCCCCC1C(=O)O", "cycloheptanecarboxylic acid"),
        ("C1CC1C(=O)O", "cyclopropanecarboxylic acid"),
        # P-65.1.2.2.2 + substituent(s) on a *different* ring atom than the
        # carboxylic acid: the ring is renumbered to give the -COOH carbon
        # locant 1 (never omitted once there's another substituent to
        # locate), then the other substituent(s) get the lowest remaining
        # locant set, cross-checked against PubChem PUG REST.
        ("CC1CCC(CC1)C(=O)O", "4-methylcyclohexane-1-carboxylic acid"),  # CID 20330
        ("ClC1CCC(CC1)C(=O)O", "4-chlorocyclohexane-1-carboxylic acid"),  # CID 12603304
        ("CC1CCC(C)C(C1)C(=O)O", "2,5-dimethylcyclohexane-1-carboxylic acid"),  # CID 14048268
    ],
)
def test_ring_carboxylic_acid_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_ring_carboxylic_acid_same_atom_substituent_raises():
    # A substituent sharing the same ring atom as the -COOH (a quaternary
    # ring carbon) is out of scope; a substituent on a *different* ring
    # atom is covered by test_ring_carboxylic_acid_names above.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(Cl)(CC1)C(=O)O")


def test_ring_carboxylic_acid_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCC1C(=O)O")


def test_ring_carboxylic_acid_standalone_hydroxyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC1CCCCC1C(=O)O")


def test_ring_carboxylic_acid_intervening_chain_carbon_raises():
    # A -COOH one chain carbon away from the ring (rather than directly
    # attached) is still out of scope for this first slice.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(CC1)CC(=O)O")


def test_ring_carboxylic_acid_two_carboxyls_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC(=O)C1CCC(C(=O)O)CC1")


def test_benzoic_acid():
    # -COOH directly on a benzene ring carbon: 'benzoic acid' is a fully
    # retained name (P-65.1.1), cross-checked against PubChem PUG REST.
    assert smiles_to_iupac("OC(=O)c1ccccc1") == "benzoic acid"  # CID 243


def test_substituted_benzoic_acid():
    # A substituent on a different ring atom than the -COOH: the retained
    # 'benzoic acid' name itself carries no locant, so only the other
    # substituent's position is cited.
    assert smiles_to_iupac("OC(=O)c1ccccc1C") == "2-methylbenzoic acid"  # CID 8373
    assert smiles_to_iupac("OC(=O)c1ccc(C)cc1") == "4-methylbenzoic acid"  # CID 7470


def test_phenyl_substituent_carboxylic_acid():
    # A plain, unsubstituted benzene ring hanging off a chain whose far end
    # carries the sole -COOH (P-65.1.1.2's territory doesn't apply here --
    # the -COOH is on the chain, not the ring itself). PubChem CID 999
    # structure-matches this to "phenylacetic acid", but this module's own
    # already-established convention uses the systematic stem once
    # substituted (see the halogen case just below) rather than the
    # retained 'acetic acid' name -- P-65.1.1 retains 'acetic acid' as PIN
    # only for unsubstituted CH3COOH itself.
    assert smiles_to_iupac("c1ccccc1CC(=O)O") == "2-phenylethanoic acid"


def test_phenyl_substituent_carboxylic_acid_matches_halogen_locant_pattern():
    # Same locant mechanism already verified for a halogen substituent
    # (test_carboxylic_acid_names' "3-chloropropanoic acid" case) applied
    # to a phenyl substituent instead -- one axis changed (substituent
    # identity), locant placement itself already proven correct.
    assert smiles_to_iupac("c1ccccc1CCC(=O)O") == "3-phenylpropanoic acid"


def test_phenyl_substituent_carboxylic_acid_directly_on_ring_is_benzoic_acid():
    # No intervening chain carbon between the ring and the -COOH carbon
    # is really the P-65.1.1 benzoic-acid case (test_benzoic_acid above),
    # reached here through the phenyl-chain module's own dispatch instead
    # of the ring module's.
    assert smiles_to_iupac("c1ccccc1C(=O)O") == "benzoic acid"


def test_phenyl_substituent_carboxylic_acid_substituted_ring_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-substituent.md) -- PubChem PUG
    # REST-verified "2-(2-methylphenyl)acetic acid".
    assert smiles_to_iupac("Cc1ccccc1CC(=O)O") == "2-(2-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_naphthalene_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2ccccc2c1CC(=O)O")


def test_phenyl_substituent_carboxylic_acid_ring_halogen():
    # The ring's other atoms may each carry a single halogen alongside the
    # chain attachment (PubChem CID 15880 IUPACName "2-(4-chlorophenyl)
    # acetic acid" -- this module's own retained-name convention prefers
    # the systematic 'ethanoic acid' stem once substituted, same as the
    # unsubstituted-ring case above).
    assert smiles_to_iupac("Clc1ccc(CC(=O)O)cc1") == "2-(4-chlorophenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_dihalogen():
    # PubChem CID 88209 "2-(2,4-dichlorophenyl)acetic acid".
    assert smiles_to_iupac("Clc1cc(Cl)ccc1CC(=O)O") == "2-(2,4-dichlorophenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_methyl():
    # PubChem PUG REST-verified "2-(4-methylphenyl)acetic acid" -- a
    # plain methyl ring substituent alongside the chain, same
    # non-halogen-specific mechanism
    # `ring_chain_attachment_with_halogens`/`halogenated_phenyl_substituent`
    # already used for halogens.
    assert smiles_to_iupac("CC1=CC=C(CC(=O)O)C=C1") == "2-(4-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_halogen_and_methyl():
    # PubChem "2-(2-chloro-5-methylphenyl)acetic acid" -- halogen and
    # methyl ring substituents mixed on the same ring.
    assert smiles_to_iupac("ClC1=CC=C(C)C=C1CC(=O)O") == "2-(2-chloro-5-methylphenyl)ethanoic acid"


def test_phenyl_substituent_carboxylic_acid_ring_alkyl_substituent_raises():
    # A non-halogen, non-methyl ring substituent alongside the chain is
    # still out of scope -- only plain halogens and a plain methyl are
    # handled so far (see tasks/aromatic-ring-methyl-substituent.md).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(C)Cc1ccc(C(C)C(=O)O)cc1")


def test_amine_coexisting_demotes_to_amino_prefix():
    # glycine (H2N-CH2-COOH): a coexisting primary amine is junior to -COOH
    # (Table 3.3) and demoted to the 'amino' prefix (see
    # _carboxylic_acid_amine.py) rather than rejected. This project's own
    # `_carboxylic_acid.py` always uses the systematic 'ethanoic acid' stem
    # rather than the retained 'acetic acid' PubChem uses for glycine's
    # actual PIN ('2-aminoacetic acid', CID 750) -- see
    # _carboxylic_acid_amine.py's module docstring for that pre-existing,
    # inherited divergence.
    assert smiles_to_iupac("NCC(=O)O") == "2-aminoethanoic acid"


def test_alcohol_mix_names_hydroxy_prefix():
    # '-oic acid' outranks 'ol' in Table 3.3, so a coexisting standalone -OH
    # is cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("OC(=O)CCO") == "3-hydroxypropanoic acid"


def test_carboxylic_acid_enol_mix_raises():
    # A hydroxyl on a C=C carbon (an enol) is a tautomer of a more senior
    # carbonyl form and out of scope, same as `_alcohol.py`'s own enol check.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(=O)O")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92), same pattern
        # as `_alcohol.py`'s `_name_acyclic_alcohol` (CIP computed
        # entirely by RDKit's `rdCIPLabeler`, not reimplemented here).
        # PubChem CID 2724540, name matches exactly (a rare case where
        # PubChem's auto-generated name equals the PIN).
        ("C[C@@H](Cl)C(=O)O", "(2R)-2-chloropropanoic acid"),
        ("C[C@H](Cl)C(=O)O", "(2S)-2-chloropropanoic acid"),
        # Two specified stereocenters, coexisting with a standalone
        # hydroxyl (P-91.3's own worked example, cited in
        # `_common.py`'s `specified_stereocenters` docstring). PubChem
        # CID 21586112, name matches exactly.
        ("C[C@H](O)[C@H](Cl)C(=O)O", "(2S,3S)-2-chloro-3-hydroxybutanoic acid"),
    ],
)
def test_carboxylic_acid_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_carboxylic_acid_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention (see `_common.py`'s `specified_stereocenters` docstring).
    assert smiles_to_iupac("CC(Cl)C(=O)O") == "2-chloropropanoic acid"


def test_carboxylic_acid_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[C@H](Cl)C(Cl)C(=O)O")
