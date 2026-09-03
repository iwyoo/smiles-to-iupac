import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # Simplest cases, cross-checked against PubChem PUG REST
        # (compound/smiles/<smiles>/property/IUPACName) and matching the
        # Blue Book's own worked examples (P6a.pdf, P-66.3.1.1/.1.2.3):
        # 'butanehydrazide (PIN)', 'pentanehydrazide (PIN)'. The
        # hydrazide carbon is always C1, and its own locant is never
        # cited (P-14.3.3), same as `_amide.py`.
        ("CCCC(=O)NN", "butanehydrazide"),
        ("CCCCC(=O)NN", "pentanehydrazide"),
        ("CCC(=O)NN", "propanehydrazide"),
        # -hydrazide + halogen substituent prefix.
        ("ClCCC(=O)NN", "3-chloropropanehydrazide"),
        # -hydrazide combined with existing unsaturation support.
        ("CC=CC(=O)NN", "but-2-enehydrazide"),
        # 'hydrazide' outranks 'ol' in Table 3.3, so a coexisting
        # standalone -OH is cited as the 'hydroxy' prefix, same as
        # `_amide.py`.
        ("OCCC(=O)NN", "3-hydroxypropanehydrazide"),
    ],
)
def test_hydrazide_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_formohydrazide():
    # P-66.3.1.2.1: 'formohydrazide' (a retained name), not the
    # systematic 'methanehydrazide', is the actual PIN here. PubChem CID
    # 12229.
    assert smiles_to_iupac("C(=O)NN") == "formohydrazide"


def test_acetohydrazide():
    # 'acetohydrazide' is the PIN, not 'ethanehydrazide'. PubChem CID
    # 14039.
    assert smiles_to_iupac("CC(=O)NN") == "acetohydrazide"


def test_2_chloroacetohydrazide():
    # A substituent on acetohydrazide's terminal carbon is cited as an
    # ordinary prefix. PubChem CID 101883.
    assert smiles_to_iupac("ClCC(=O)NN") == "2-chloroacetohydrazide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # N-/N'-substituted hydrazide (P-66.3.1.1): the carbonyl-adjacent
        # nitrogen is "N", the terminal one "N'", mirroring `_amide.py`'s
        # "N-"/`_urea.py`'s "N-"/"N'-" convention. All cross-checked
        # against PubChem PUG REST (compound/smiles/<smiles>/property/
        # IUPACName).
        ("CC(=O)N(C)N", "N-methylacetohydrazide"),  # CID 19051
        ("CC(=O)NNC", "N'-methylacetohydrazide"),  # CID 122488
        ("CC(=O)N(C)NC", "N,N'-dimethylacetohydrazide"),  # CID 12596294
        ("CC(=O)NN(C)C", "N',N'-dimethylacetohydrazide"),  # CID 80385
        ("CC(=O)N(C)N(C)C", "N,N',N'-trimethylacetohydrazide"),  # CID 12362621
        ("CC(=O)N(CC)NC", "N-ethyl-N'-methylacetohydrazide"),  # CID 89037658
        ("CC(=O)N(C)NCC", "N'-ethyl-N-methylacetohydrazide"),  # CID 119096326
        # The systematic (non-retained-name) chain lengths behave the same
        # way as the retained 'aceto-' case above.
        ("CCCC(=O)N(C)N", "N-methylbutanehydrazide"),  # CID 53649946
        ("CCCC(=O)NNC", "N'-methylbutanehydrazide"),  # CID 88556637
    ],
)
def test_n_substituted_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_n_substituted_stereocenter_hydrazide():
    # Mechanical combination of two independently-verified axes: N-prefix
    # placement (see above) and stereodescriptor-outermost ordering
    # (`test_hydrazide_stereocenter` below) -- not itself registered on
    # PubChem, so treated as a reviewed/eyeballed result rather than an
    # independently-verified one (see the base '(2R)-2-methylbutanehydrazide'
    # case this extends, PubChem CID 30066157).
    assert smiles_to_iupac("CC[C@@H](C)C(=O)N(C)N") == "(2R)-N-methyl-2-methylbutanehydrazide"


def test_n_substituted_formohydrazide_raises():
    # Unlike the dinuclear (aceto-) and longer chain cases, the mononuclear
    # retained name itself doesn't survive N-substitution -- PubChem
    # switches to plain 'formamide' + amino-substituent naming instead
    # (`O=CN(C)N` -> "N-amino-N-methylformamide", CID 12962), a different
    # naming paradigm this module doesn't attempt.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("O=CN(C)N")


def test_branched_n_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC(=O)N(N)C(C)C")


def test_diacylhydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CCCC(=O)NNC(=O)C")


def test_amide_not_misnamed_as_hydrazide():
    # A plain primary amide (`_amide.py`) has only one nitrogen and must
    # not be routed here.
    assert smiles_to_iupac("CC(N)=O") == "ethanamide"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): a hydrazide
        # carbon is always chain-terminal (fixed C1), same pattern as
        # `_amide.py` (CIP computed entirely by RDKit's `rdCIPLabeler`).
        # PubChem CID 30066157.
        ("CC[C@@H](C)C(=O)NN", "(2R)-2-methylbutanehydrazide"),
        ("CC[C@H](C)C(=O)NN", "(2S)-2-methylbutanehydrazide"),
    ],
)
def test_hydrazide_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_hydrazide_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)C(=O)NN") == "2-methylbutanehydrazide"


def test_phenyl_chain_hydrazide():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#288's
    # carboxylic-acid/.../hydroperoxide chains): the ring is cited as a
    # "phenyl" substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCC(=O)NN") == "3-phenylpropanehydrazide"


def test_phenyl_chain_hydrazide_retained_name():
    # Dinuclear case: the retained 'acetohydrazide' name is still the PIN
    # even with the phenyl substituent, mirroring the module's own
    # '2-chloroacetohydrazide' pattern (module docstring).
    assert smiles_to_iupac("c1ccccc1CC(=O)NN") == "2-phenylacetohydrazide"


def test_phenyl_directly_attached_hydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C(=O)NN")


def test_phenyl_chain_hydrazide_n_alkyl_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1CC(=O)N(C)N")


def test_phenyl_substituted_benzene_ring_hydrazide_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC(=O)NN")


def test_phenyl_chain_hydrazide_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC(=O)NN")
