import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_unsubstituted_ammonium_name():
    # NH4+, PubChem CID 223 -- the Blue Book PIN and PubChem's own
    # auto-generated name agree exactly on "azanium" (P-73.1.1.2).
    assert smiles_to_iupac("[NH4+]") == "azanium"


def test_methanaminium_name():
    # CH3-NH3+, structure-verified via PubChem CID 644041
    # (ConnectivitySMILES "[NH3+]C"; PubChem's own auto-generated name
    # "methylazanium" isn't PIN format). Name derived from the Blue Book's
    # own P-73.1.1.2 rule ("final 'e' -> 'ium'") applied to this project's
    # existing, independently verified "methanamine" (P-33.1).
    assert smiles_to_iupac("C[NH3+]") == "methanaminium"


def test_propan_1_aminium_name():
    # CH3-CH2-CH2-NH3+, structure-verified via PubChem CID 3483736
    # ("propylazanium"). Derived from this project's existing
    # "propan-1-amine" the same way as above.
    assert smiles_to_iupac("CCC[NH3+]") == "propan-1-aminium"


def test_propan_2_aminium_name():
    # (CH3)2CH-NH3+, structure-verified via PubChem CID 3364502
    # ("propan-2-ylazanium"). Derived from this project's existing
    # "propan-2-amine" the same way as above.
    assert smiles_to_iupac("CC(C)[NH3+]") == "propan-2-aminium"


def test_chloromethanaminium_name():
    # Cl-CH2-NH3+, structure-verified via PubChem CID 20462735
    # ("chloromethylazanium"). Confirms halogen-prefix coexistence,
    # derived from this project's existing "chloromethanamine" the same
    # way as above.
    assert smiles_to_iupac("C(Cl)[NH3+]") == "chloromethanaminium"


def test_secondary_ammonium():
    assert smiles_to_iupac("C[NH2+]C") == "N-methylmethanaminium"


def test_secondary_ammonium_asymmetric():
    # Structure/base-name cross-check against _amine.py's own
    # "N-ethylethanamine" for the un-protonated diethylamine.
    assert smiles_to_iupac("CC[NH2+]CC") == "N-ethylethanaminium"


def test_tertiary_ammonium():
    assert smiles_to_iupac("CC[NH+](CC)CC") == "N,N-diethylethanaminium"


def test_quaternary_ammonium_symmetric():
    # Blue Book P-73.1.1.1 Table 7.3 worked example: (CH3)4N+ ->
    # "N,N,N-trimethylmethanaminium (PIN)" -- explicitly NOT the
    # alternative "tetramethylazanium" form (non-PIN, PubChem's own
    # auto-generated name for this structure).
    assert smiles_to_iupac("C[N+](C)(C)C") == "N,N,N-trimethylmethanaminium"


def test_quaternary_ammonium_asymmetric():
    assert smiles_to_iupac("CC[N+](C)(C)C") == "N,N,N-trimethylethanaminium"


def test_tertiary_ammonium_with_halogen_on_parent_chain():
    # Inherited from `_amine.py`'s N-prefix/halogen interleaving fix
    # (P-14.5.2) via the secondary/tertiary path's own name_amine call.
    assert smiles_to_iupac("ClCC[NH+](CC)CC") == "2-chloro-N,N-diethylethanaminium"


def test_quaternary_ammonium_with_halogen_on_parent_chain():
    # Same fix, reached via the quaternary path's direct call to
    # `_amine.py`'s `_name_acyclic_secondary_tertiary_amine`.
    assert smiles_to_iupac("ClCC[N+](CC)(CC)CC") == "2-chloro-N,N,N-triethylethanaminium"


def test_ring_ammonium_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C[NH2+]C1CCCCC1")


def test_aromatic_ammonium():
    # Anilinium (P-73, protonated aniline) -- now that `_amine.py` supports
    # aniline itself, the neutralize-and-'e'->'ium' path this module
    # already uses for chain amines applies here too. Blue Book worked
    # examples confirm "anilinium chloride (PIN)"/"N,N,N-trimethylanilinium
    # (PIN)" (`tmp/bluebook/P7.txt` lines 1795/3567).
    assert smiles_to_iupac("c1ccccc1[NH3+]") == "anilinium"


def test_doubly_charged_nitrogen_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[NH3++]")


def test_primary_ammonium_stereocenter():
    # Degree 1-3 ammonium delegates to `name_amine` on the neutralized
    # molecule, so it inherited full R/S support automatically once
    # `_amine.py` gained it (P-91.3/P-92).
    assert smiles_to_iupac("CC[C@@H](C)[NH3+]") == "(2R)-butan-2-aminium"


def test_quaternary_ammonium_stereocenter():
    # The quaternary (degree 4) path calls
    # `_name_acyclic_secondary_tertiary_amine` directly rather than
    # delegating to `name_amine`, so it needed its own `specified_stereocenters`
    # wiring -- CIP/locant cross-checked against PubChem's own (non-PIN
    # 'azanium'-style) name for the same structure: "[(2R)-butan-2-yl]-
    # trimethylazanium" (CID 102019526).
    assert smiles_to_iupac("CC[C@@H](C)[N+](C)(C)C") == "(2R)-N,N,N-trimethylbutan-2-aminium"
    assert smiles_to_iupac("CC[C@H](C)[N+](C)(C)C") == "(2S)-N,N,N-trimethylbutan-2-aminium"


def test_quaternary_ammonium_unspecified_stereocenter_unaffected():
    assert smiles_to_iupac("CCC(C)[N+](C)(C)C") == "N,N,N-trimethylbutan-2-aminium"


def test_nitrogen_ylide_not_misnamed_as_ammonium():
    # `has_ammonium_shape` only inspects the nitrogen's own local bonding,
    # which a genuine quaternary ammonium and a nitrogen ylide (a P-74.2
    # dipolar compound, e.g. trimethylammonium methylide) both match
    # identically -- without an explicit check this silently dropped the
    # carbanion and returned "N,N,N-trimethylmethanaminium" instead of the
    # real zwitterion name (`_ylide.py`, P-74.2.1.1.1).
    assert smiles_to_iupac("[CH2-][N+](C)(C)C") == "(N,N-dimethylmethanaminiumyl)methanide"
