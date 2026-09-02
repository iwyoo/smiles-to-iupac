import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName). The amide carbon is
        # always C1, and its own locant is never cited (P-14.3.3).
        ("C(N)=O", "methanamide"),
        ("CC(N)=O", "ethanamide"),
        ("CCC(N)=O", "propanamide"),
        # A real positional choice for a substituent: the amide carbon fixes
        # C1 regardless.
        ("CC(C)CC(N)=O", "3-methylbutanamide"),
        # -amide combined with existing unsaturation support, cross-checked
        # against PubChem (crotonamide's PIN).
        ("CC=CC(N)=O", "but-2-enamide"),
        # -amide + halogen substituent prefix, cross-checked against
        # PubChem.
        ("NC(=O)CCCl", "3-chloropropanamide"),
    ],
)
def test_amide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_methylethanamide():
    # PubChem structure match: "N-methylacetamide" (PubChem uses the
    # retained 'acetamide' stem; this module uses the systematic PIN stem
    # 'ethanamide' consistently with its own plain-amide tests above).
    assert smiles_to_iupac("CNC(C)=O") == "N-methylethanamide"


def test_n_n_dimethylethanamide():
    # PubChem structure match: "N,N-dimethylacetamide".
    assert smiles_to_iupac("CC(=O)N(C)C") == "N,N-dimethylethanamide"


def test_n_ethyl_n_methylethanamide():
    # PubChem structure match: "N-ethyl-N-methylacetamide".
    assert smiles_to_iupac("CC(=O)N(C)CC") == "N-ethyl-N-methylethanamide"


def test_n_substituted_amide_with_longer_acyl_chain():
    # The N-substituent must not be mistaken for the acyl chain even when
    # it happens to be longer than it.
    assert smiles_to_iupac("CC(=O)NCCCC") == "N-butylethanamide"


def test_n_n_disubstituted_amide_with_locant_leading_parent_name():
    # The N,N-prefix must be hyphen-separated from a parent name that
    # itself starts with a numeric locant, not just concatenated.
    assert (
        smiles_to_iupac("O=C(C(C(F)(F)F)C(F)(F)F)N(CCCl)CCCl")
        == "N,N-diethyl-3,3,3-trifluoro-2-(trifluoromethyl)propanamide"
    )


def test_n_substituted_amide_with_locant_leading_parent_name_raises():
    # This SMILES's N-substituent isn't actually a plain propyl group -- it
    # carries two hydroxyls of its own ("CC(O)C(O)N..."), which this
    # module's docstring already scopes out ("plain, unbranched,
    # unsubstituted" N-substituent only). Before `specified_stereocenters`
    # was wired in (tasks/amide-stereocenter-naming.md), this was silently
    # accepted and misnamed as "N-propyl-2,3,4,5,6,7-hexahydroxyheptanamide"
    # -- a genuine pre-existing bug on two fronts: the substituted
    # N-substituent was never rejected, and the acyl chain's five specified
    # stereocenters were silently dropped. It now correctly raises (either
    # via the new stereo safety net noticing the N-substituent's own two
    # stereocenters are left unspecified while the acyl chain's are
    # specified, or via the new explicit N-substituent-hydroxyl check).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(O)C(O)NC(=O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO")


def test_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC(C)C")


def test_unsaturated_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)NC=C")


def test_lactam_is_named_via_ketone_module():
    # A plain, unsubstituted lactam (ketone carbonyl directly bonded to the
    # ring's own N-H) is routed to `_ketone.py`'s hetero-ring ketone path
    # instead of being rejected here -- see
    # tasks/hetero-ring-ketone-lactam-routing.md and test_ketone.py's own
    # coverage. PubChem-verified: CID 12025.
    assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"


def test_n_substituted_lactam_raises():
    # An N-substituted lactam is still out of scope: the ring heteroatom
    # itself may only carry its own indicated hydrogen (see
    # `_ketone.py`'s hetero-ring ketone path, which is what actually
    # raises here -- this shape is routed there ahead of this module, same
    # as the plain lactam above).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C1CCCN1C")


def test_diamide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)CC(N)=O")


def test_ester_not_misnamed_as_amide():
    # -COO- (ester, `_ester.py`) is not amide-shaped: its carbonyl carbon's
    # other oxygen neighbor is carbon-bonded, not a nitrogen, so this must
    # not be routed here and misnamed.
    assert smiles_to_iupac("CC(=O)OC") == "methyl ethanoate"


def test_carboxylic_acid_not_misnamed_as_amide():
    assert smiles_to_iupac("CC(=O)O") == "ethanoic acid"


def test_aryl_amide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(=O)c1ccccc1")


def test_alcohol_mix_names_hydroxy_prefix():
    # 'amide' outranks 'ol' in Table 3.3, so a coexisting standalone -OH is
    # cited as the 'hydroxy' substituent prefix rather than rejected.
    assert smiles_to_iupac("NC(=O)CCO") == "3-hydroxypropanamide"


def test_amide_enol_mix_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OC=CC(N)=O")


def test_n_substituent_with_hydroxyl_raises():
    # A hydroxyl on the N-substituent itself (not the acyl chain) is
    # invisible to the carbon-only chain walk that measures N-substituent
    # length, so it must be checked separately -- only a plain,
    # unsubstituted alkyl N-substituent is in scope (module docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(O)NC(=O)C")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): an amide
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_carboxylic_acid.py`/`_aldehyde.py` (CIP computed entirely by
        # RDKit's `rdCIPLabeler`). PubChem CID 41097729.
        ("CC[C@@H](C)C(N)=O", "(2R)-2-methylbutanamide"),
        ("CC[C@H](C)C(N)=O", "(2S)-2-methylbutanamide"),
    ],
)
def test_amide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_amide_stereocenter():
    # The stereodescriptor sits outermost, ahead of the N-alkyl prefix
    # (P-91.3: always at the very front of the complete name). PubChem
    # CID 89121954.
    assert smiles_to_iupac("CC[C@@H](C)C(=O)NCC") == "(2R)-N-ethyl-2-methylbutanamide"


def test_amide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(N)=O") == "2-methylbutanamide"
