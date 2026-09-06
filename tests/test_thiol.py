import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_methanethiol():
    # '-thiol' suffix construction mirrors '-ol' (_alcohol.py, see
    # test_alcohol.py, already cross-checked against PubChem) with the same
    # locant-citation rules; 'methanethiol'/'ethanethiol' are also the
    # well-known common names of these compounds.
    assert smiles_to_iupac("SC") == "methanethiol"


def test_ethanethiol():
    assert smiles_to_iupac("SCC") == "ethanethiol"


def test_propane_1_thiol():
    assert smiles_to_iupac("SCCC") == "propane-1-thiol"


def test_propane_2_thiol():
    assert smiles_to_iupac("CC(S)C") == "propane-2-thiol"


def test_butane_1_thiol_with_chloro_substituent():
    assert smiles_to_iupac("SCCCCCl") == "4-chlorobutane-1-thiol"


def test_pent_4_ene_1_thiol():
    assert smiles_to_iupac("SCCCC=C") == "pent-4-ene-1-thiol"


def test_ethane_1_2_dithiol():
    # PubChem CID 10902.
    assert smiles_to_iupac("SCCS") == "ethane-1,2-dithiol"


def test_propane_1_3_dithiol():
    # PubChem CID 8013.
    assert smiles_to_iupac("SCCCS") == "propane-1,3-dithiol"


def test_propane_1_2_dithiol():
    # PubChem CID 61217.
    assert smiles_to_iupac("CC(S)CS") == "propane-1,2-dithiol"


def test_butane_1_4_dithiol():
    # PubChem CID 79148.
    assert smiles_to_iupac("SCCCCS") == "butane-1,4-dithiol"


def test_cyclohexanethiol():
    # PubChem CID 15290.
    assert smiles_to_iupac("SC1CCCCC1") == "cyclohexanethiol"


def test_cyclohexane_1_2_dithiol():
    # PubChem CID 5242646.
    assert smiles_to_iupac("SC1CCCCC1S") == "cyclohexane-1,2-dithiol"


def test_unsaturated_ring_thiol():
    # Monocyclic ring, single -SH, single ring double bond (P-31.1.3): the
    # thiol always gets locant 1 (suffix priority), the ring double bond's
    # locant is minimized by choosing direction. Cross-checked against
    # PubChem (CID 21916868/53767965). Note 'thiol' starts with a
    # consonant so 'ene' doesn't elide, unlike the ketone/alcohol
    # '-en-1-one'/'-en-1-ol' pattern.
    assert smiles_to_iupac("SC1CCCC=C1") == "cyclohex-2-ene-1-thiol"
    assert smiles_to_iupac("SC1CC=CCC1") == "cyclohex-3-ene-1-thiol"


def test_unsaturated_ring_thiol_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC=C1C")


def test_unsaturated_ring_thiol_triple_bond_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CCCC#C1")


def test_2_methylcyclohexane_1_thiol():
    # PubChem CID 519947.
    assert smiles_to_iupac("CC1CCCCC1S") == "2-methylcyclohexane-1-thiol"


def test_cyclopentanethiol():
    # PubChem CID 15510.
    assert smiles_to_iupac("SC1CCCC1") == "cyclopentanethiol"


def test_polycyclic_thiol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SC1CC2CCC1CC2")


def test_ring_substituent_chain_thiol():
    # A thiol entirely on a chain hanging off a plain saturated ring (the
    # ring itself bears no thiol) -- the ring is cited as a "cyclo..."
    # substituent prefix on the chain, mirroring `_name_phenyl_chain_
    # thiol`/`_alcohol.py`'s `_name_ring_substituent_chain_alcohol`.
    # PubChem PUG REST-verified "cyclohexylmethanethiol" (CID 520209).
    assert smiles_to_iupac("SCC1CCCCC1") == "cyclohexylmethanethiol"


def test_ring_substituent_chain_thiol_internal_locant():
    # PubChem PUG REST: "1-cyclohexylethanethiol" (CID 18994009).
    assert smiles_to_iupac("SC(C)C1CCCCC1") == "1-cyclohexylethane-1-thiol"


def test_ring_substituent_chain_thiol_ring_with_substituent_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC1CCC(C)CC1")


def test_ring_substituent_chain_thiol_unsaturated_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCC1CCCC=C1")


def test_thiol_with_alcohol_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("SCCO")


def test_phenyl_chain_thiol():
    # A plain, unsubstituted benzene ring on the chain (P-2/P-3
    # aromatic-ring-substituent extension, mirroring PR #269-#277's
    # carboxylic-acid/ketone/alcohol/ester/aldehyde/amide/nitrile/acyl-
    # halide/sulfonic-acid chains): the ring is cited as a "phenyl"
    # substituent prefix.
    assert smiles_to_iupac("c1ccccc1CCCS") == "3-phenylpropane-1-thiol"


def test_phenyl_chain_thiol_internal_locant():
    # The -SH locant is a genuine choice on the chain, same as the base
    # acyclic module.
    assert smiles_to_iupac("C(c1ccccc1)C(C)S") == "1-phenylpropane-2-thiol"


def test_benzenethiol():
    # -SH directly on a benzene ring carbon, cross-checked against
    # PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1S") == "benzenethiol"  # CID 7969


def test_substituted_benzenethiol():
    # A substituent on a different ring atom than the -SH: the
    # mancude-ring numbering is free to start at the -SH carbon, so its
    # own locant is never cited, unlike the cycloalkane case.
    assert smiles_to_iupac("Cc1ccccc1S") == "2-methylbenzenethiol"  # CID 8712


def test_two_direct_ring_thiols_not_supported():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Sc1ccccc1S")


def test_phenyl_substituted_benzene_ring_thiol_ortho_methyl():
    # A ring methyl substituent is now supported (see
    # tasks/aromatic-ring-methyl-rollout-2.md).
    assert smiles_to_iupac("Cc1ccccc1CCS") == "2-(2-methylphenyl)ethane-1-thiol"


def test_phenyl_chain_thiol_ring_halogen():
    # PubChem PUG REST computes "3-(4-chlorophenyl)propane-1-thiol" for
    # this structure.
    assert smiles_to_iupac("Clc1ccc(CCCS)cc1") == "3-(4-chlorophenyl)propane-1-thiol"


def test_phenyl_chain_thiol_ring_methyl():
    # PubChem PUG REST-verified "3-(4-methylphenyl)propane-1-thiol".
    assert smiles_to_iupac("Cc1ccc(CCCS)cc1") == "3-(4-methylphenyl)propane-1-thiol"


def test_phenyl_chain_thiol_ring_ethyl():
    # PubChem PUG REST IUPACName match: any plain, fully saturated
    # acyclic alkyl ring substituent (not just methyl) is now supported,
    # reusing `name_branch` itself via `plain_alkyl_ring_substituents`.
    assert smiles_to_iupac("CCc1ccc(cc1)CS") == "(4-ethylphenyl)methanethiol"


def test_phenyl_chain_thiol_unsaturation_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C=Cc1ccccc1CCS")


def test_phenyl_chain_dithiol():
    # Two -SH groups on the chain hanging off the benzene ring, cross-
    # checked against PubChem PUG REST.
    assert smiles_to_iupac("c1ccccc1C(S)CCS") == "1-phenylpropane-1,3-dithiol"  # CID 57158210
    assert smiles_to_iupac("SCC(S)Cc1ccccc1") == "3-phenylpropane-1,2-dithiol"  # CID 154223538


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # A single specified tetrahedral stereocenter (P-92): unlike
        # `_sulfinic_acid.py`, a thiol's -SH sulfur is monovalent (bonded
        # only to carbon and H) and can never itself be a stereocenter, so
        # this mirrors `_carboxylic_acid.py`/`_aldehyde.py`/`_ketone.py`/
        # `_sulfonic_acid.py` cleanly (CIP computed entirely by RDKit's
        # `rdCIPLabeler`). PubChem CID 444090.
        ("CC[C@@H](C)S", "(2R)-butane-2-thiol"),
        ("CC[C@H](C)S", "(2S)-butane-2-thiol"),
    ],
)
def test_acyclic_thiol_stereocenter(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_cyclic_thiol_stereocenter():
    # Two stereocenters on the ring itself (P-92), same pattern as
    # `_alcohol.py`'s `_name_cyclic_alcohol`. PubChem CID 22211639 (name
    # carries a redundant 'cis-' relative descriptor this project drops
    # once full R/S is given, same policy as the existing ring-alcohol
    # stereocenter task).
    assert smiles_to_iupac("S[C@H]1CCCC[C@H]1C") == "(1S,2R)-2-methylcyclohexane-1-thiol"


def test_cyclic_thiol_branch_stereocenter():
    # A stereocenter on the ring's sole substituent branch rather than the
    # ring itself (P-92), same pattern as `_alcohol.py`'s
    # `_name_cyclic_alcohol`/`_ketone.py`'s `_name_cyclic_ketone` -- the
    # branch and the -SH share the same ring carbon (C1), same shape as
    # the alcohol precedent. PubChem has no cached record for this exact
    # structure (CID 0, sparse data gap), but the unstereo parent
    # ('1-ethylcyclohexane-1-thiol') matches PubChem exactly, and the
    # descriptor-construction code itself is identical to the
    # already-verified alcohol/ketone cases.
    assert (
        smiles_to_iupac("SC1(CCCCC1)[C@@H](C)CC")
        == "1-[(2S)-butan-2-yl]cyclohexane-1-thiol"
    )


def test_thiol_unspecified_stereocenter_unaffected():
    # A genuine stereocenter left unspecified (no @/@@) is named exactly
    # as before -- no stereo prefix, matching this project's long-standing
    # convention.
    assert smiles_to_iupac("CCC(C)S") == "butane-2-thiol"


def test_thiol_partially_specified_stereocenters_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("S[C@H]1CCCCC1Cl")
