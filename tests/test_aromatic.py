import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # benzene, the retained PIN (P-25.1.2.1); RDKit sanitizes both
        # lowercase-aromatic and Kekulized SMILES to the same structure.
        ("c1ccccc1", "benzene"),
        ("C1=CC=CC=C1", "benzene"),
        # naphthalene: the only possible 2-ring ortho-fused shape, numbered
        # 1,2,3,4,4a,5,6,7,8,8a by the general algorithm (P-25.3.3.1.1) --
        # matches the well-known numbering, a good sanity check even though
        # naphthalene isn't a P-25.3.3 exception.
        ("c1ccc2ccccc2c1", "naphthalene"),
        # anthracene: straight 3-ring chain, but P-25.3.3 retains its
        # traditional (non-walk) numbering rather than the general
        # algorithm's -- see _aromatic.py's module docstring for the
        # structural proof (derived from, and cross-checked against,
        # PubChem's own round-tripped IUPAC names below).
        ("c1ccc2cc3ccccc3cc2c1", "anthracene"),
        # phenanthrene: bent/angular 3-ring chain (verified via RDKit that
        # its ring-adjacency graph is angular, not straight, by checking the
        # two fusion bonds of the middle ring are NOT opposite hexagon
        # edges) -- also a P-25.3.3 exception, traditional numbering.
        ("c1ccc2c(c1)ccc1ccccc12", "phenanthrene"),
        # tetracene: straight 4-ring chain (P-25.1.2.1 polyacene), a real
        # test of the general algorithm producing the standard numbering
        # (meso-like positions 5,6,11,12) since tetracene is not a P-25.3.3
        # exception. Cross-checked against PubChem CID 7080.
        ("C1=CC=C2C=C3C=C4C=CC=CC4=CC3=CC2=C1", "tetracene"),
        # pentacene: straight 5-ring chain, same general algorithm one ring
        # further. Cross-checked against PubChem CID 8671.
        ("C1=CC=C2C=C3C=C4C=C5C=CC=CC5=CC4=CC3=CC2=C1", "pentacene"),
        # Halogen substituents (P-35.2.1), reusing the existing
        # halogen_substituents/name_branch machinery, cross-checked against
        # PubChem: chlorobenzene (CID 7964), 2-chloronaphthalene (CID 7056).
        ("C1=CC=C(C=C1)Cl", "chlorobenzene"),
        ("C1=CC=C2C=C(C=CC2=C1)Cl", "2-chloronaphthalene"),
        # 2-chloroanthracene (CID 28308): exercises the traditional
        # numbering's non-fusion "1-4, 5-8" terminal-ring locants under
        # substitution, not just the unsubstituted parent name.
        ("C1=CC=C2C=C3C=C(C=CC3=CC2=C1)Cl", "2-chloroanthracene"),
        # 1-bromoanthracene (CID 12529827) and 1,8-dichloroanthracene
        # (CID 618890): pin down which physical atom is locant 1 in the
        # traditional numbering (adjacent to a terminal-ring fusion carbon).
        ("C1=CC=C2C=C3C(=CC2=C1)C=CC=C3Br", "1-bromoanthracene"),
        ("C1=CC2=CC3=C(C=C2C(=C1)Cl)C(=CC=C3)Cl", "1,8-dichloroanthracene"),
        # 9,10-dibromoanthracene (CID 68226): both traditional "meso"
        # carbons substituted at once -- these are the two ring-2 nonfusion
        # atoms, which are graph-theoretically nonadjacent (verified via
        # RDKit), so this also regression-tests that the anthracene-specific
        # locant lookup doesn't accidentally wander onto a substituent atom
        # instead of the real ring neighbor when both meso positions bear a
        # substituent.
        ("C1=CC=C2C(=C1)C(=C3C=CC=CC3=C2Br)Br", "9,10-dibromoanthracene"),
        # 3-chlorophenanthrene (CID 3013936) and 9-bromophenanthrene
        # (CID 11309): phenanthrene's traditional numbering under
        # substitution, including the K-region (9,10) position.
        ("C1=CC=C2C(=C1)C=CC3=C2C=C(C=C3)Cl", "3-chlorophenanthrene"),
        ("C1=CC=C2C(=C1)C=C(C3=CC=CC=C23)Br", "9-bromophenanthrene"),
        # 2-chlorotetracene (CID 12336463): substituent on a general-
        # algorithm (non-exception) 4-ring straight chain.
        ("C1=CC=C2C=C3C=C4C=C(C=CC4=CC3=CC2=C1)Cl", "2-chlorotetracene"),
    ],
)
def test_smiles_to_iupac_aromatic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_peri_fused_raises():
    # pyrene: verified via RDKit that some atom is shared by three rings
    # (peri-fusion), explicitly out of scope (P-25.3.1.3, ortho-fusion only).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1cc2ccc3cccc4ccc(c1)c2c34")


def test_branched_fusion_raises():
    # triphenylene: verified via RDKit that its ring-adjacency graph has a
    # degree-3 node (three rings all ortho-fused to one central ring), not
    # a simple chain -- out of scope (P-25.3.1.3).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2c(c1)c1ccccc1c1ccccc21")


def test_heteroaromatic_raises():
    # pyridine: not even detected as in scope (find_aromatic_fused_core
    # requires an all-carbon ring), falls through to the existing
    # heteroatom rejection in validate_atoms_and_bonds.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccncc1")


def test_heteroaromatic_fused_raises():
    # quinoline: a fused heteroaromatic, same reasoning as pyridine above.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1ccc2ncccc2c1")


def test_angular_four_ring_raises():
    # chrysene (CID 9171): a bent/angular 4-ring chain (polyphene, n=4) --
    # explicitly out of scope for this PR (only the straight 4-ring
    # tetracene is supported; an angular n>=4 chain has no single retained
    # name and would need genuine benzo[x,y-z]fusion[...] construction).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CC4=CC=CC=C43")


def test_saturated_rings_still_resolve_unaffected():
    # Regression check: the new aromatic dispatch branch in core.py must
    # not shadow or interfere with existing saturated-ring detection.
    assert smiles_to_iupac("C1CCCCC1") == "cyclohexane"
    assert smiles_to_iupac("C1CC2CCC1C2") == "bicyclo[2.2.1]heptane"
    assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "tricyclo[3.3.1.1^3,7]decane"
    assert smiles_to_iupac("C1C2C3C2C4C1C34") == "tetracyclo[3.2.0.0^2,7.0^4,6]heptane"
