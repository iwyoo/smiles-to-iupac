import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanetellurol():
    # PubChem PUG REST CID 356643, auto-generated name matches exactly.
    assert smiles_to_iupac("C[TeH]") == "methanetellurol"


def test_ethanetellurol():
    # PubChem PUG REST CID 71407255, auto-generated name matches exactly --
    # the same P-14.3.4.2(b) two-carbon locant omission as 'ethaneselenol'/
    # 'ethanethiol'.
    assert smiles_to_iupac("CC[TeH]") == "ethanetellurol"


def test_propane_1_tellurol():
    # PubChem PUG REST CID 71405781, auto-generated name matches exactly.
    assert smiles_to_iupac("CCC[TeH]") == "propane-1-tellurol"


def test_2_methylpropane_1_tellurol():
    # A branched chain. PubChem PUG REST CID 101718199, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("CC(C)C[TeH]") == "2-methylpropane-1-tellurol"


def test_prop_2_ene_1_tellurol():
    # An unsaturated chain. PubChem PUG REST CID 101841305, auto-generated
    # name matches exactly.
    assert smiles_to_iupac("C=CC[TeH]") == "prop-2-ene-1-tellurol"


def test_2_chloroethane_1_tellurol():
    # Halogen coexistence (P-35.2.1) -- no PubChem-listed compound found
    # for this specific structure (CID 0), so this is a reviewed result,
    # not an independently verified one: the mechanism itself already has
    # independent confirmation via `_selenol.py`'s/`_thiol.py`'s own
    # identical '2-chloroethane-1-selenol'/'2-chloroethane-1-thiol' cases.
    assert smiles_to_iupac("ClCC[TeH]") == "2-chloroethanetellurol"


def test_ethane_1_2_ditellurol():
    # `_selenol.py`'s/`_thiol.py`'s own precedent generalizes group-count
    # support beyond one without a specific 2+-group PubChem worked
    # example of its own (the shared locant/suffix machinery already
    # handles an arbitrary-length locant list) -- same reasoning applied
    # here.
    assert smiles_to_iupac("[TeH]CC[TeH]") == "ethane-1,2-ditellurol"


def test_propane_1_3_ditellurol():
    assert smiles_to_iupac("[TeH]CCC[TeH]") == "propane-1,3-ditellurol"


def test_cyclohexanetellurol():
    # PubChem structure match: "cyclohexanetellurol".
    assert smiles_to_iupac("C1CCCCC1[TeH]") == "cyclohexanetellurol"


def test_unsaturated_ring_tellurol():
    # Monocyclic ring, single -TeH, single ring double bond (P-31.1.3):
    # same pattern already confirmed for _thiol.py/_selenol.py. The exact
    # ring structures aren't PubChem-registered (sparse tellurol ring
    # coverage, same as the plain 'cyclohexanetellurol' case above), but
    # the acyclic ene+tellurol combination is already PubChem-confirmed
    # ('prop-2-ene-1-tellurol', CID 101841305).
    assert smiles_to_iupac("[TeH]C1CCCC=C1") == "cyclohex-2-ene-1-tellurol"
    assert smiles_to_iupac("[TeH]C1CC=CCC1") == "cyclohex-3-ene-1-tellurol"


def test_unsaturated_ring_tellurol_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CCCC=C1C")


def test_unsaturated_ring_tellurol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]C1CCCC#C1")


def test_substituted_ring_tellurol_locant_cited():
    assert smiles_to_iupac("CC1CCCCC1[TeH]") == "2-methylcyclohexane-1-tellurol"


def test_von_baeyer_tellurol_now_supported():
    # A single tellurol on a von Baeyer bicyclic/polycyclic or monospiro
    # ring system is supported (see `test_von_baeyer_spiro_tellurol.py`).
    assert smiles_to_iupac("C12(CCC(CC1)CC2)[TeH]") == "bicyclo[2.2.2]octane-1-tellurol"


def test_ring_substituent_chain_tellurol():
    # A tellurol entirely on a chain hanging off a plain saturated ring
    # (the ring itself bears no tellurol) -- the ring is cited as a
    # "cyclo..." substituent prefix on the chain, mirroring
    # `_name_phenyl_chain_tellurol`/`_thiol.py`'s `_name_ring_substituent_
    # chain_thiol`. No PubChem-registered structure for this exact
    # molecule -- a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`.
    assert smiles_to_iupac("C1CCCCC1C[TeH]") == "cyclohexylmethanetellurol"


def test_ring_substituent_chain_tellurol_ring_with_substituent():
    assert smiles_to_iupac("C1CCC(C)CC1C[TeH]") == "(3-methylcyclohexyl)methanetellurol"


def test_ring_substituent_chain_tellurol_unsaturated_ring():
    assert smiles_to_iupac("C1=CCCCC1C[TeH]") == "(cyclohex-2-en-1-yl)methanetellurol"


def test_ring_with_tellurol_chain_tellurol_tie():
    # Ring and chain each carry exactly one tellurol (P-44.1.1 tie,
    # P-44.1.2.2 resolves it in the ring's favor), mirroring `_thiol.py`'s
    # `test_ring_with_thiol_chain_thiol_tie`. No PubChem-registered
    # structure for this exact molecule (sparse tellurol coverage), so
    # this is a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`.
    assert (
        smiles_to_iupac("[TeH]C1CCCCC1CC[TeH]") == "2-(2-tellanylethyl)cyclohexane-1-tellurol"
    )


def test_ring_with_tellurol_chain_tellurol_ring_wins_outright():
    # The ring carries two tellurols against the chain's one -- P-44.1.1's
    # greater-count rule picks the ring outright, no tie-break needed.
    assert (
        smiles_to_iupac("[TeH]C1CCCCC1(C[TeH])[TeH]")
        == "1-(tellanylmethyl)cyclohexane-1,2-ditellurol"
    )


def test_tellurol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("OCC[TeH]")


def test_acyclic_tellurol_stereocenter():
    # A single specified tetrahedral stereocenter (P-92): like -SH/-SeH, a
    # tellurol's -TeH tellurium is monovalent and can never itself be a
    # stereocenter, so this mirrors `_thiol.py`/`_selenol.py` cleanly (CIP
    # computed entirely by RDKit's `rdCIPLabeler`). PubChem has no
    # registered stereoisomer for this molecule (tellurium compounds are
    # sparse there), so this is a structural/regression check on the
    # already-proven mechanism, not an independent PubChem cross-check.
    assert smiles_to_iupac("C[C@@H](CC)C[TeH]") == "(2S)-2-methylbutane-1-tellurol"


def test_cyclic_tellurol_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_thiol.py`'s `_name_cyclic_thiol`. Same PubChem-sparsity caveat as
    # the acyclic case above.
    assert smiles_to_iupac("[TeH][C@H]1CCCC[C@H]1C") == "(1S,2R)-2-methylcyclohexane-1-tellurol"


def test_cyclic_tellurol_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_thiol.py`/`_selenol.py` -- the
    # branch and the -TeH share the same ring carbon (C1). Same
    # PubChem-sparsity caveat as the other tellurol stereocenter tests.
    assert (
        smiles_to_iupac("[TeH]C1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-tellurol"
    )


def test_tellurol_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention. PubChem CID 173009567.
    assert smiles_to_iupac("CCC(C)[TeH]") == "butane-2-tellurol"


def test_tellurol_partially_specified_stereocenters_cites_the_specified_elements():
    assert smiles_to_iupac("[TeH][C@H]1CCCCC1Cl") == '(1S)-2-chlorocyclohexane-1-tellurol'


def test_phenyl_chain_tellurol():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring `_thiol.py`'s
    # `test_phenyl_chain_thiol`): the ring is cited as a "phenyl"
    # substituent prefix. PubChem PUG REST: "3-phenylpropane-1-tellurol".
    assert smiles_to_iupac("c1ccccc1CCC[TeH]") == "3-phenylpropane-1-tellurol"


def test_phenyl_chain_tellurol_internal_locant():
    # The -TeH locant is a genuine choice on the chain, same as the base
    # acyclic module. No PubChem-registered structure for this exact
    # molecule -- a structural/regression check on the mechanism ported
    # verbatim from `_thiol.py`.
    assert smiles_to_iupac("C(c1ccccc1)C(C)[TeH]") == "1-phenylpropane-2-tellurol"


def test_benzenetellurol():
    # -TeH directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST (structure match; PubChem has no computed
    # IUPACName for the substituted case, so only the unsubstituted case
    # is independently name-verified here).
    assert smiles_to_iupac("c1ccccc1[TeH]") == "benzenetellurol"  # CID 5246059


def test_substituted_benzenetellurol():
    # Mechanical extension of the verified unsubstituted case, mirroring
    # `_thiol.py`'s/`_selenol.py`'s identical proofs.
    assert smiles_to_iupac("Cc1ccccc1[TeH]") == "2-methylbenzenetellurol"


def test_phenyl_substituted_benzene_ring_tellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccccc1CC[TeH]")


def test_phenyl_chain_tellurol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CC[TeH]")


def test_phenyl_chain_ditellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccccc1C([TeH])CC[TeH]")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[TeH]Cc1cccnc1", "(pyridin-3-yl)methanetellurol"),
        ("[TeH]Cc1cc[nH]c1", "(1H-pyrrol-3-yl)methanetellurol"),
    ],
)
def test_heteroaromatic_chain_tellurol(smiles, expected):
    # Reviewed, not independently PubChem-verified (tellurium
    # heteroaromatic-chain compounds are essentially unregistered there) --
    # reuses the exact same locant machinery already independently
    # verified for `_thiol.py`'s/`_selenol.py`'s identical shape, just
    # swapping the chalcogen (see module docstring).
    assert smiles_to_iupac(smiles) == expected


def test_heteroaromatic_direct_attachment_tellurol_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[TeH]c1cccnc1")


def test_oxygen_or_sulfur_heteroaromatic_chain_tellurol():
    assert smiles_to_iupac("[TeH]Cc1cccs1") == "(thiophen-2-yl)methanetellurol"
    assert smiles_to_iupac("[TeH]Cc1ccco1") == "(furan-2-yl)methanetellurol"
