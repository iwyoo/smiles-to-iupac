import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_paracyclophane():
    # [2.2]paracyclophane's PIN, confirmed directly against the Blue Book
    # text (P-26.3.2.1 and P-26.4.1.4, both marked "(PIN)"). Two equally
    # valid SMILES (this project's own, and Wikipedia's) both resolve to
    # the same canonical structure, C16H16.
    assert smiles_to_iupac("C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2") == "1,4(1,4)-dibenzenacyclohexaphane"
    assert smiles_to_iupac("C=1C=C2C=CC1CCC3=CC=C(C=C3)CC2") == "1,4(1,4)-dibenzenacyclohexaphane"


def test_metacyclophane():
    # [2.2]metacyclophane's PIN, confirmed directly against the Blue Book
    # text (P-26.4.1.4's own paired example alongside the para isomer
    # above, both marked "(PIN)"). PubChem CID 137543's own canonical
    # SMILES for "(2.2)Metacyclophane" is used here and independently
    # matches this project's own SMILES for the same structure, C16H16.
    assert smiles_to_iupac("C1CC2=CC(=CC=C2)CCC3=CC=CC1=C3") == "1,4(1,3)-dibenzenacyclohexaphane"
    assert smiles_to_iupac("c1cc2cc(c1)CCc1cccc(c1)CC2") == "1,4(1,3)-dibenzenacyclohexaphane"


def test_tetrabenzenacyclooctaphane():
    # [1.1.1.1]metacyclophane's PIN, confirmed directly against the Blue
    # Book text (P-26.4.1.4's third worked example, printed immediately
    # alongside the para/meta pair above, all marked "(PIN)"). This
    # project's own independently-built SMILES resolves (via InChIKey) to
    # PubChem CID 11740710, same formula and connectivity, C28H24.
    assert (
        smiles_to_iupac("C1c2cccc(c2)Cc2cccc(c2)Cc2cccc(c2)Cc2cccc1c2")
        == "1,3,5,7(1,3)-tetrabenzenacyclooctaphane"
    )


def test_substituted_tetrabenzenacyclooctaphane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cccc2c1CC1=CC=CC(=C1)CC1=CC=CC(=C1)CC1=CC=CC(=C1)C2")


def test_substituted_paracyclophane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC2=CC=C(CCC3=CC=C1C=C3)C=C2")


def test_substituted_metacyclophane_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CC2=CC(=CC=C2)CCC3=CC=CC1=C3")


def test_different_bridge_length_now_supported():
    # [3.3]paracyclophane: same shape as [2.2]paracyclophane but with
    # three-carbon bridges instead of two -- supported by the general
    # symmetric-phane algorithm, computed directly from the same formula verified
    # against the three Blue Book PIN worked examples, not a fresh guess.
    # Structure cross-checked against PubChem CID 137763 (its own computed
    # IUPACName is von-Baeyer-style, same limitation as the other cases).
    assert smiles_to_iupac("C1CCC2=CC=C(CCCC3=CC=C1C=C3)C=C2") == "1,5(1,4)-dibenzenacyclooctaphane"


def test_three_ring_symmetric_phane():
    # N=3 (a "triangular" assembly, not just the N=2/N=4 cases the Blue
    # Book's own worked examples happen to show): built from scratch via
    # RDKit's RWMol (three para-attached benzene rings, each connected to
    # the next by a two-carbon bridge), proving this is a genuine N-general
    # algorithm and not just a lookup covering the two N values that
    # happened to have worked examples. Structure cross-checked against
    # PubChem CID 4342704 (its own computed IUPACName is von-Baeyer-style,
    # same limitation as the other cases).
    assert (
        smiles_to_iupac("c1cc2ccc1CCc1ccc(cc1)CCc1ccc(cc1)CC2")
        == "1,4,7(1,4)-tribenzenacyclononaphane"
    )


def test_mixed_local_pattern_raises():
    # one ring para-attached, the other meta-attached (built from scratch
    # via RDKit's RWMol) -- a non-uniform pattern, out of scope for this
    # module (see docstring).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1cc2cc(c1)CCc1ccc(cc1)CC2")


def test_two_separate_paracyclophane_units_not_misread_as_one_n4_phane():
    # two independent [2.2]paracyclophane-like units (four rings total,
    # forming two separate 2-cycles, not one single 4-cycle spanning all
    # four rings) -- must not be misrecognized as a single N=4 symmetric
    # phane. In practice a multi-fragment input like this is already
    # rejected earlier by the general multi-fragment guard (P-13.6)
    # before this module's own dispatch is even reached, since two
    # complete, self-contained 2-cycles necessarily have no bond
    # connecting them to each other -- exercised here as a regression
    # guard on `_find_symmetric_phane`'s own single-cycle validation
    # regardless of which guard actually fires first.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac(
            "C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2.C1CC2=CC=C(CCC3=CC=C1C=C3)C=C2"
        )
