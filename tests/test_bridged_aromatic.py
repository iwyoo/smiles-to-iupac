import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles",
    [
        # benzonorbornadiene: CAS 4453-90-1, PubChem CID 97391, NIST WebBook
        # gives the IUPAC name "1,4-Dihydro-1,4-methanonaphthalene".
        "C1=CC2CC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring
        # instead of the bridged ring -- the lowest-locants tie-break must
        # still normalize this to "1,4-", not some other numbering.
        "c1ccc2c(c1)C1C=CC2C1",
    ],
)
def test_bridged_naphthalene(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-methanonaphthalene"


@pytest.mark.parametrize(
    "smiles",
    [
        # 1,4-Epoxy-1,4-dihydronaphthalene ("7-oxabenzonorbornadiene"):
        # CAS 573-57-9, PubChem CID 97139 (structure/formula cross-checked
        # via ConnectivitySMILES "C1=CC=C2C3C=CC(C2=C1)O3" -- PubChem's own
        # computed IUPACName for this shape is a von Baeyer bridged-ring
        # name, not usable to verify the fusion+bridge name itself).
        "C1=CC2OC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2O1",
    ],
)
def test_bridged_naphthalene_epoxy(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-epoxynaphthalene"


def test_naphthalene_itself_is_unaffected():
    assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


def test_dihydronaphthalene_itself_is_unaffected():
    assert smiles_to_iupac("C1CC=Cc2ccccc12") == "1,2-dihydronaphthalene"


def test_fully_saturated_bridge_is_von_baeyer_not_fusion_name():
    # a fully saturated bridgehead pair (no "ene" double bond in the
    # reduced ring) is a different molecule from `test_bridged_naphthalene`
    # above -- P-25.4's fusion+bridge name doesn't apply, but the aromatic
    # ring is still a valid von Baeyer Kekule cyclohexatriene (P-31.1.4.2).
    # PubChem CID 138272 confirms this exact name.
    assert smiles_to_iupac("C1CC2CC1c1ccccc12") == "tricyclo[6.2.1.0^2,7]undeca-2,4,6-triene"


def test_substituted_bridge_atom_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)C1c1ccccc12")


def test_substituted_aromatic_ring_raises():
    # a methyl (not a halogen) substituent remains out of scope.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2C1")


def test_halogen_on_bridged_naphthalene_aromatic_ring():
    # chlorobenzonorbornadiene: built from scratch from the existing
    # unsubstituted test fixture plus one Cl on the intact aromatic ring.
    # Structure cross-checked against PubChem CID 12473502 (PubChem's own
    # computed IUPACName is von-Baeyer-style, same limitation as the
    # unsubstituted case -- see module docstring).
    assert smiles_to_iupac("Clc1ccc2c(c1)C1C=CC2C1") == "6-chloro-1,4-dihydro-1,4-methanonaphthalene"
    # different ring position -> different, still lowest-locant, result.
    assert smiles_to_iupac("c1cc(Cl)c2c(c1)C1C=CC2C1") == "5-chloro-1,4-dihydro-1,4-methanonaphthalene"


def test_halogen_on_bridged_naphthalene_epoxy_aromatic_ring():
    # the epoxy analogue of the chlorinated case above -- same mechanism,
    # not bridge-prefix-specific. Structure cross-checked against PubChem
    # CID 14208771 (same von-Baeyer-name limitation).
    assert smiles_to_iupac("Clc1ccc2c(c1)C1C=CC2O1") == "6-chloro-1,4-dihydro-1,4-epoxynaphthalene"


def test_ether_itself_is_unaffected():
    assert smiles_to_iupac("COC") == "methoxymethane"


def test_substituted_epoxy_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2O1")


@pytest.mark.parametrize(
    "smiles",
    [
        # 9,10-Dihydro-9,10-methanoanthracene: literature name confirmed in
        # J. Org. Chem. ("9,10-Dihydro-9,10-methanoanthracene and Its
        # Perhydro Derivatives"); structure/formula (C15H12) cross-checked
        # against PubChem CID 12651785's ConnectivitySMILES.
        "C12c3ccccc3C(c3ccccc31)C2",
        # Same molecule, atom order starting from one of the intact
        # aromatic rings instead of a bridgehead -- the lowest-locants
        # tie-break must still normalize this to "9,10-", not some other
        # numbering (anthracene's meso positions are always 9,10 by
        # definition, but this still exercises the candidate search).
        "c12ccccc1C1c3ccccc3C2C1",
    ],
)
def test_bridged_anthracene(smiles):
    assert smiles_to_iupac(smiles) == "9,10-dihydro-9,10-methanoanthracene"


def test_anthracene_itself_is_unaffected():
    assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1") == "anthracene"


def test_bridged_anthracene_epoxy_raises():
    # Unlike the naphthalene case, an 'epoxy' bridge on anthracene's 9,10
    # positions is deliberately out of scope here -- no independently
    # verifiable name was found while scoping this (see
    # `find_bridged_anthracene_core`'s docstring), only the structure
    # (PubChem CID cross-checked via ConnectivitySMILES), which isn't
    # enough under this project's test-writing policy to assert a name.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12c3ccccc3C(c3ccccc31)O2")


def test_substituted_bridged_anthracene_bridge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12c3ccccc3C(c3ccccc31)C2(C)")


def test_substituted_bridged_anthracene_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1c3ccccc3C2C1")


def test_bridged_anthracene_on_terminal_ring():
    # A 1,4-type bridge on one of anthracene's own *terminal* rings, not
    # the already-tested 9,10 meso positions -- structurally the same
    # "one fusion neighbor + one ene neighbor per bridgehead" shape
    # `test_bridged_naphthalene` exercises, just on a 3-ring parent, so
    # this is the case `find_bridged_aromatic_core`'s generalization now
    # reaches via `_anthracene_candidates` (previously only reachable via
    # `find_bridged_naphthalene_core`'s hardcoded 11-atom, 2-ring-only
    # match). PubChem CID 606064 confirms both the structure and this
    # exact name as a registered synonym ("1,4-Dihydro-1,4-methanoanthracene",
    # alongside the alternate CAS index name "Naphtho[2,3-b]norbornadiene").
    assert smiles_to_iupac("C1=CC2CC1c1cc3ccccc3cc12") == "1,4-dihydro-1,4-methanoanthracene"


def test_halogen_on_bridged_anthracene_terminal_ring():
    assert (
        smiles_to_iupac("Clc1ccc2cc3c(cc2c1)C1C=CC3C1") == "6-chloro-1,4-dihydro-1,4-methanoanthracene"
    )


def test_bridged_phenanthrene():
    # A 1,4-type bridge on one terminal ring of phenanthrene (a *bent*
    # 3-ring chain, unlike anthracene's straight one) -- exercises
    # `_phenanthrene_candidates` via the same generalized core-finder.
    # PubChem CID 606062 confirms the structure (its own registered
    # synonym uses the alternate CAS norbornadiene-base convention,
    # "Naphtho[1,2-b]norbornadiene", rather than this fusion+bridge
    # style -- P-25.4.1.2 requires the bridge to attach to an
    # already-named *fused ring system*, so the mancude phenanthrene
    # system, not the non-mancude bridged bicyclic, is the correct P-25.4
    # base component; the terminal-ring-anthracene case above
    # independently confirms this project's mechanism produces a real,
    # registered name for the identical shape on a different parent).
    assert smiles_to_iupac("C1=CC2CC1c1ccc3ccccc3c12") == "1,4-dihydro-1,4-methanophenanthrene"


def test_bridged_tetracene():
    # A 1,4-type bridge on one terminal ring of a straight 4-ring chain
    # (tetracene) -- exercises the polyacene branch of
    # `_retained_chain_name`/`_straight_chain_candidates` via the same
    # generalized core-finder. PubChem CID 13082583 confirms the
    # structure (same CAS-convention-vs-P-25.4 caveat as the phenanthrene
    # case above).
    assert smiles_to_iupac("C1=CC2CC1c1cc3cc4ccccc4cc3cc12") == "1,4-dihydro-1,4-methanotetracene"


@pytest.mark.parametrize(
    "smiles",
    [
        # PubChem CID 68694281, structure/formula (C10H8S) cross-checked
        # via ConnectivitySMILES "C1=CC=C2C3C=CC(C2=C1)S3" -- same
        # von-Baeyer-name limitation as the epoxy case above.
        "C1=CC2SC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2S1",
    ],
)
def test_bridged_naphthalene_sulfano(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-sulfanonaphthalene"


@pytest.mark.parametrize(
    "smiles",
    [
        # PubChem CID 138429, structure/formula (C10H9N) cross-checked via
        # ConnectivitySMILES "C1=CC=C2C3C=CC(C2=C1)N3" -- same
        # von-Baeyer-name limitation as the epoxy/sulfano cases above.
        "C1=CC2NC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2N1",
    ],
)
def test_bridged_naphthalene_azano(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-azanonaphthalene"


def test_halogen_on_bridged_naphthalene_sulfano_aromatic_ring():
    # the sulfano analogue of the chlorinated methano/epoxy cases above --
    # same mechanism, not bridge-prefix-specific.
    assert smiles_to_iupac("Clc1ccc2c(c1)C1C=CC2S1") == "6-chloro-1,4-dihydro-1,4-sulfanonaphthalene"


def test_substituted_sulfano_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2S1")


def test_bridged_anthracene_terminal_ring_sulfano():
    # A 1,4-type sulfano bridge on one terminal ring of anthracene (the
    # same terminal-ring shape as `test_bridged_naphthalene_sulfano`, on a
    # 3-ring chain instead of a 2-ring one) -- exercises the same
    # generalized core-finder's anthracene-parent branch. Not to be
    # confused with `find_bridged_anthracene_core`'s deliberately
    # methano-only 9,10-meso bridge (see that function's own docstring).
    assert smiles_to_iupac("C1=CC2SC1c1cc3ccccc3cc12") == "1,4-dihydro-1,4-sulfanoanthracene"


@pytest.mark.parametrize(
    "smiles",
    [
        # Primary-source worked example (`tmp/bluebook/P2.txt` ~5892):
        # "1,4-dihydro-1,4-ethanoanthracene (PIN) (not
        # 1,2,3,4-tetrahydro-1,4-ethenoanthracene)" -- confirms both the
        # bridge prefix and the required 'dihydro' added-hydrogen prefix
        # for the naphthalene-shaped case here (same terminal-ring shape
        # as the methano/epoxy/sulfano/azano tests above, a 2-atom bridge
        # instead of 1).
        "C1=CC2CCC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2CC1",
    ],
)
def test_bridged_naphthalene_ethano(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-ethanonaphthalene"


def test_halogen_on_bridged_naphthalene_ethano_aromatic_ring():
    # the ethano analogue of the chlorinated methano/epoxy/sulfano cases
    # above -- same mechanism, not bridge-length-specific.
    assert smiles_to_iupac("Clc1ccc2c(c1)C1C=CC2CC1") == "6-chloro-1,4-dihydro-1,4-ethanonaphthalene"


def test_substituted_ethano_bridge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)CC1c1ccccc12")


def test_substituted_ethano_aromatic_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccc2c(c1)C1C=CC2CC1")


def test_bridged_anthracene_terminal_ring_ethano():
    # Directly confirmed by the primary source's own worked example (see
    # module docstring): "1,4-dihydro-1,4-ethanoanthracene (PIN)".
    assert smiles_to_iupac("C1=CC2CCC1c1cc3ccccc3cc12") == "1,4-dihydro-1,4-ethanoanthracene"


@pytest.mark.parametrize(
    "smiles",
    [
        # P-25.4.2.1.1's own consecutive worked-example table gives
        # "methano"/"ethano"/"propano" together (`tmp/bluebook/P2.txt`
        # ~5306) -- the same 'dihydro' structural-necessity reasoning as
        # the methano/ethano cases (#899) applies unchanged, this being
        # one more length of the same acyclic carbon-chain bridge, now
        # generalized (#905) instead of hardcoded per length. Structure
        # existence confirmed via PubChem (CID 14689353's connectivity).
        "C1=CC2CCCC1c1ccccc12",
        # Same molecule, atom order starting from the intact aromatic ring.
        "c1ccc2c(c1)C1C=CC2CCC1",
    ],
)
def test_bridged_naphthalene_propano(smiles):
    assert smiles_to_iupac(smiles) == "1,4-dihydro-1,4-propanonaphthalene"


def test_substituted_propano_bridge_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC2C(C)CCC1c1ccccc12")


def test_bridged_anthracene_terminal_ring_propano():
    assert smiles_to_iupac("C1=CC2CCCC1c1cc3ccccc3cc12") == "1,4-dihydro-1,4-propanoanthracene"
