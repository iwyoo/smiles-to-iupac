import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methaneselenol():
    # PubChem PUG REST CID 440764, auto-generated name matches exactly.
    assert smiles_to_iupac("C[SeH]") == "methaneselenol"


def test_ethaneselenol():
    # PubChem PUG REST CID 5252527, auto-generated name matches exactly --
    # the same P-14.3.4.2(b) two-carbon locant omission as 'ethanethiol'.
    assert smiles_to_iupac("CC[SeH]") == "ethaneselenol"


def test_propane_1_selenol():
    # PubChem PUG REST CID 71373846, auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[SeH]") == "propane-1-selenol"


def test_2_methylpropane_1_selenol():
    # A branched chain. PubChem PUG REST CID 157368643, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC(C)C[SeH]") == "2-methylpropane-1-selenol"


def test_prop_2_ene_1_selenol():
    # An unsaturated chain. PubChem PUG REST CID 15821407, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("C=CC[SeH]") == "prop-2-ene-1-selenol"


def test_2_chloroethane_1_selenol():
    # Halogen coexistence (P-35.2.1) -- no PubChem-listed compound found
    # for this specific structure (every halogenated selenol SMILES tried
    # came back as CID 0), so this is a reviewed result, not an
    # independently verified one: the mechanism itself already has
    # independent confirmation via `_thiol.py`'s own identical
    # '2-chloroethane-1-thiol' case (same locant-citation rule once a
    # substituent is present on the two-carbon chain).
    assert smiles_to_iupac("ClCC[SeH]") == "2-chloroethane-1-selenol"


def test_ethane_1_2_diselenol():
    # PubChem PUG REST auto-generated name matches exactly.
    assert smiles_to_iupac("[SeH]CC[SeH]") == "ethane-1,2-diselenol"


def test_propane_1_3_diselenol():
    # PubChem PUG REST auto-generated name matches exactly.
    assert smiles_to_iupac("[SeH]CCC[SeH]") == "propane-1,3-diselenol"


def test_propane_1_2_3_triselenol():
    # `_thiol.py`'s own precedent generalizes group-count support beyond
    # two without a specific 3+-group PubChem worked example of its own
    # (the shared `_alcohol.py`-style locant/suffix machinery already
    # handles an arbitrary-length locant list) -- same reasoning applied
    # here.
    assert smiles_to_iupac("[SeH]CC([SeH])C[SeH]") == "propane-1,2,3-triselenol"


def test_butane_1_2_4_triselenol():
    assert smiles_to_iupac("[SeH]CC([SeH])CC[SeH]") == "butane-1,2,4-triselenol"


def test_cyclohexaneselenol():
    # PubChem structure match: "cyclohexaneselenol".
    assert smiles_to_iupac("C1CCCCC1[SeH]") == "cyclohexaneselenol"


def test_unsaturated_ring_selenol():
    # Monocyclic ring, single -SeH, single ring double bond (P-31.1.3):
    # same pattern already confirmed for _thiol.py (PR #350) -- the exact
    # ring structures aren't PubChem-registered (sparse selenol ring
    # coverage, same as the plain 'cyclohexaneselenol' case above), but
    # the acyclic ene+selenol combination and elision rule are already
    # PubChem-confirmed ('prop-2-ene-1-selenol', CID 15821407).
    assert smiles_to_iupac("[SeH]C1CCCC=C1") == "cyclohex-2-ene-1-selenol"
    assert smiles_to_iupac("[SeH]C1CC=CCC1") == "cyclohex-3-ene-1-selenol"


def test_unsaturated_ring_selenol_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC=C1C")


def test_unsaturated_ring_selenol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]C1CCCC#C1")


def test_cyclopentaneselenol():
    # PubChem structure match: "cyclopentaneselenol".
    assert smiles_to_iupac("C1CCCC1[SeH]") == "cyclopentaneselenol"


def test_substituted_ring_selenol_locant_cited():
    assert smiles_to_iupac("CC1CCCCC1[SeH]") == "2-methylcyclohexane-1-selenol"


def test_polycyclic_selenol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12(CCC(CC1)CC2)[SeH]")


def test_ring_substituent_chain_selenol():
    # A selenol entirely on a chain hanging off a plain saturated ring
    # (the ring itself bears no selenol) -- the ring is cited as a
    # "cyclo..." substituent prefix on the chain, mirroring
    # `_name_phenyl_chain_selenol`/`_thiol.py`'s `_name_ring_substituent_
    # chain_thiol`. No PubChem-registered structure for this exact
    # molecule -- a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`, same as `test_phenyl_chain_selenol_
    # internal_locant` above.
    assert smiles_to_iupac("C1CCCCC1C[SeH]") == "cyclohexylmethaneselenol"


def test_ring_substituent_chain_selenol_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCC(C)CC1C[SeH]")


def test_ring_substituent_chain_selenol_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CCCCC1C[SeH]")


def test_ring_with_selenol_chain_selenol_tie():
    # Ring and chain each carry exactly one selenol (P-44.1.1 tie,
    # P-44.1.2.2 resolves it in the ring's favor), mirroring `_thiol.py`'s
    # `test_ring_with_thiol_chain_thiol_tie`. No PubChem-registered
    # structure for this exact molecule (sparse selenol coverage), so
    # this is a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`.
    assert (
        smiles_to_iupac("[SeH]C1CCCCC1CC[SeH]") == "2-(2-selanylethyl)cyclohexane-1-selenol"
    )


def test_ring_with_selenol_chain_selenol_ring_wins_outright():
    # The ring carries two selenols against the chain's one -- P-44.1.1's
    # greater-count rule picks the ring outright, no tie-break needed.
    assert (
        smiles_to_iupac("[SeH]C1CCCCC1(C[SeH])[SeH]")
        == "1-(selanylmethyl)cyclohexane-1,2-diselenol"
    )


def test_selenol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[SeH]")


def test_acyclic_selenol_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): like a thiol's
    # -SH sulfur, a selenol's -SeH selenium is monovalent and can never
    # itself be a stereocenter, so this mirrors `_thiol.py` cleanly (CIP
    # computed entirely by RDKit's `rdCIPLabeler`). PubChem CID 175184280.
    assert smiles_to_iupac("C[C@@H](CC)C[SeH]") == "(2S)-2-methylbutane-1-selenol"


def test_cyclic_selenol_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_thiol.py`'s `_name_cyclic_thiol`. PubChem has no registered
    # stereoisomer for this exact ring (even the non-stereo parent isn't
    # registered), so this is a structural/regression check on the
    # already-proven mechanism ported verbatim from `_thiol.py`, not an
    # independent PubChem cross-check.
    assert smiles_to_iupac("[SeH][C@H]1CCCC[C@H]1C") == "(1S,2R)-2-methylcyclohexane-1-selenol"


def test_cyclic_selenol_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_thiol.py`'s
    # `_name_cyclic_thiol` -- the branch and the -SeH share the same ring
    # carbon (C1), same shape as the thiol precedent. PubChem has no
    # registered structure here either (not even the non-stereo parent),
    # so this is a structural/regression check on the already-proven
    # mechanism ported verbatim from `_thiol.py`.
    assert (
        smiles_to_iupac("[SeH]C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-selenol"
    )


def test_selenol_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)[SeH]") == "butane-2-selenol"


def test_selenol_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH][C@H]1CCCCC1Cl")


def test_phenyl_chain_selenol():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_thiol.py`'s
    # `test_phenyl_chain_thiol`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "3-phenylpropane-1-selenol".
    assert smiles_to_iupac("c1ccccc1CCC[SeH]") == "3-phenylpropane-1-selenol"


def test_phenyl_chain_selenol_internal_locant():
    # The -SeH locant is a genuine choice on the chain, same as the base
    # acyclic module. No PubChem-registered structure for this exact
    # molecule -- a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`.
    assert smiles_to_iupac("C(c1ccccc1)C(C)[SeH]") == "1-phenylpropane-2-selenol"


def test_benzeneselenol():
    # -SeH directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST (structure match; PubChem has no computed
    # IUPACName for the substituted case, so only the unsubstituted case
    # is independently name-verified here).
    assert smiles_to_iupac("c1ccccc1[SeH]") == "benzeneselenol"  # CID 69530


def test_substituted_benzeneselenol():
    # Mechanical extension of the verified unsubstituted case, mirroring
    # `_thiol.py`'s identical 'benzenethiol'/'2-methylbenzenethiol' proof.
    assert smiles_to_iupac("Cc1ccccc1[SeH]") == "2-methylbenzeneselenol"


def test_phenyl_substituted_benzene_ring_selenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[SeH]")


def test_phenyl_chain_selenol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[SeH]")


def test_phenyl_chain_diselenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C([SeH])CC[SeH]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[SeH]Cc1cccnc1", "(pyridin-3-yl)methaneselenol"),
        ("[SeH]Cc1cc[nH]c1", "(1H-pyrrol-3-yl)methaneselenol"),
    ],
)
def test_heteroaromatic_chain_selenol(smiles, expected):
    # Reviewed, not independently PubChem-verified (selenium
    # heteroaromatic-chain compounds are essentially unregistered there) --
    # reuses the exact same locant machinery already independently
    # verified for `_thiol.py`'s identical shape, just swapping the
    # chalcogen (see module docstring).
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_direct_attachment_selenol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]c1cccnc1")


def test_oxygen_or_sulfur_heteroaromatic_chain_selenol_not_yet_reached():
    # Furan/thiophene aren't reachable yet: `core.py`'s dispatch routes
    # these to `_ether.py`/`_sulfide.py` first (see module docstring) --
    # a separate, pre-existing dispatch gap, not specific to selenol.
    # Both still raise `UnsupportedStructure` rather than a wrong name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]Cc1cccs1")
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[SeH]Cc1ccco1")
