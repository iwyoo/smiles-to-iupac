import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Retained names, Blue Book P-72.2.2.2.2's own text.
        ("C[O-]", "methoxide"),
        ("CC[O-]", "ethoxide"),
        ("CCC[O-]", "propoxide"),
        ("CCCC[O-]", "butoxide"),
        # Non-terminal oxygen falls through to the systematic '-olate'
        # suffix, even for a chain length that has a retained name for its
        # terminal position: 'propan-2-olate (PIN)', not 'isopropoxide'
        # (Blue Book's own worked example; PubChem CID 3260420 structure
        # match, and PubChem's own generated name agrees exactly here).
        ("CC(C)[O-]", "propan-2-olate"),
        # Longer unbranched chain, no retained name at that length --
        # PubChem CID 13076479, name matches exactly.
        ("CCCCC[O-]", "pentan-1-olate"),
    ],
)
def test_alkoxide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_alkoxide_with_halogen_substituent():
    # PubChem CID 101079710, "fluoromethanolate" -- a halogen substituent
    # on the retained methoxide's plain shape falls through to the
    # systematic mononuclear-parent path (P-14.3.4.2(a): locant omitted).
    assert smiles_to_iupac("FC[O-]") == "fluoromethanolate"


def test_alkoxide_with_halogen_on_two_carbon_chain():
    # PubChem CID 17767529 (structure match; PubChem's own generated name
    # is "2-fluoroethanolate" -- this project's `_alcohol.py` already
    # established citing the '-1-' locant here for the neutral alcohol
    # ('2-fluoroethan-1-ol', PubChem CID divergence accepted project-wide),
    # so the anion mirrors that same convention with 'ate' appended.
    assert smiles_to_iupac("FCC[O-]") == "2-fluoroethan-1-olate"


def test_alkoxide_with_unsaturation():
    # PubChem CID 21252372, "prop-2-en-1-olate" -- exact match.
    assert smiles_to_iupac("C=CC[O-]") == "prop-2-en-1-olate"


def test_tert_butoxide():
    # Blue Book P-72.2.2.2.2's own text names 'tert-butoxide' as the
    # retained PIN for (CH3)3C-O(-) (structure confirmed via PubChem,
    # which itself returns the systematic '2-methylpropan-2-olate').
    assert smiles_to_iupac("CC(C)(C)[O-]") == "tert-butoxide"


def test_branched_alkoxide():
    # PubChem CID structure/name match: "2-methylpropan-1-olate".
    assert smiles_to_iupac("CC(C)C[O-]") == "2-methylpropan-1-olate"


def test_branched_alkoxide_longer_chain():
    # Cross-checked against _alcohol.py's own "3-methylbutan-1-ol" for the
    # same skeleton (isoamyl alcohol) -- same branch, 'ate' appended. Also
    # a regression check: this 5-carbon skeleton's longest chain is 4 atoms
    # (butoxide's own retained-name length), so the retained-name fast path
    # must not misfire here -- that path is reserved for a truly unbranched
    # chain using every carbon in the molecule.
    assert smiles_to_iupac("CC(C)CC[O-]") == "3-methylbutan-1-olate"


def test_branched_alkoxide_with_halogen():
    assert smiles_to_iupac("ClCC(C)C[O-]") == "3-chloro-2-methylpropan-1-olate"


def test_phenoxide():
    # -O(-) attached directly to a benzene ring carbon (P-63.8.1/
    # P-72.2.2.2.2). The retained name 'phenoxide' stands for the whole
    # ring+O(-) system (like 'phenol'), so the O(-)'s own ring locant is
    # never cited, only other substituents'. PubChem PUG REST:
    # "phenoxide" (CID 998, exact auto-generated-name match); the
    # substituted cases follow the Blue Book's own "substituted the same
    # way as the corresponding alcohols" text rather than PubChem's own
    # differently-styled "...phenolate" auto-generated names (module
    # docstring).
    assert smiles_to_iupac("[O-]c1ccccc1") == "phenoxide"
    assert smiles_to_iupac("[O-]c1ccc(C)cc1") == "4-methylphenoxide"
    assert smiles_to_iupac("[O-]c1ccccc1Cl") == "2-chlorophenoxide"


def test_phenoxide_non_alkyl_ring_substituent_raises():
    # A ring substituent containing another heteroatom (here, a
    # carboxylic acid) must not be silently walked as if it were a plain
    # carbon chain -- found via real-data testing: this exact SMILES (a
    # copper salt of a non-carboxylate anion, itself out of scope for
    # `_salt.py`) was misnamed '2-(propan-2-yl)phenoxide', treating the
    # -COOH branch's two oxygens as if they were carbons.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=C(O)c1ccccc1[O-].[Cu+]")


def test_two_alkoxide_groups_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CC[O-]")


def test_ether_oxygen_alongside_alkoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[O-]CCOC")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92) -- the
        # alkoxide oxygen itself is never a potential stereocenter (a
        # monovalent, negatively-charged terminal atom), confirmed via
        # RDKit `FindPotentialStereo`. Structure/CIP cross-checked against
        # the parent alcohol, PubChem CID 84682 "(2R)-butan-2-ol" (PubChem
        # doesn't resolve the bare anion SMILES to a CID).
        ("CC[C@@H](C)[O-]", "(2R)-butan-2-olate"),
        ("CC[C@H](C)[O-]", "(2S)-butan-2-olate"),
    ],
)
def test_acyclic_alkoxide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_alkoxide_stereocenter_with_coexisting_halogen():
    assert smiles_to_iupac("CC[C@@H](Cl)[O-]") == "(1R)-1-chloropropan-1-olate"


def test_alkoxide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)[O-]") == "butan-2-olate"


def test_phenyl_chain_alkoxide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_alcohol.py`'s
    # `test_phenyl_chain_alcohol`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "2-phenylethanolate"/
    # "3-phenylpropan-1-olate" (the two-carbon case keeps its locant here,
    # same known, out-of-scope-here limitation already established by this
    # module's own `test_fluoroethanolate_locant` for the halogen case).
    assert smiles_to_iupac("c1ccccc1CC[O-]") == "2-phenylethan-1-olate"
    assert smiles_to_iupac("c1ccccc1CCC[O-]") == "3-phenylpropan-1-olate"


def test_phenyl_directly_attached_alkoxide():
    # Same phenoxide-type case as test_phenoxide, reached via a different
    # ring-atom ordering in the SMILES.
    assert smiles_to_iupac("c1ccccc1[O-]") == "phenoxide"


def test_phenyl_substituted_benzene_ring_alkoxide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[O-]")


def test_phenyl_chain_alkoxide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[O-]")
