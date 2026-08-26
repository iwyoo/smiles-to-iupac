import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # adamantane (flagship case): four bridgeheads, each pair connected by
        # a one-carbon bridge (tetrahedral K4 arrangement). Cross-checked
        # against the NIST WebBook / PubChem name for CAS 281-23-2, and
        # independently re-derived by hand from P-23.2.2/P-23.2.3/P-23.2.5.1/
        # P-23.2.5.2 (see _tricyclic.py's module docstring): with main
        # bridgeheads numbered 1 and 5, the other two bridgeheads land at the
        # midpoints of the two 3-atom main-ring segments, positions 3 and 7.
        ("C1C2CC3CC1CC(C2)C3", "tricyclo[3.3.1.1^3,7]decane"),
        # twistane: built from scratch (not copied from a SMILES database) to
        # match the descriptor tricyclo[4.4.0.0^3,8]decane reported by
        # Wikipedia/ChemSpider (CAS 253-14-5) -- main bridgeheads directly
        # bonded (c=0) and the secondary bridgeheads (at positions 3 and 8,
        # one on each 4-atom main-ring segment) also directly bonded (d=0).
        ("C1CC2CC3CCC2CC13", "tricyclo[4.4.0.0^3,8]decane"),
        # 1-methyladamantane: methyl on a main bridgehead. All four
        # bridgeheads are equivalent in unsubstituted adamantane, so the
        # lowest-locant tie-break (P-14.4/P-45.2) picks the numbering that
        # puts the substituted one at position 1.
        ("CC12CC3CC(CC(C3)C1)C2", "1-methyltricyclo[3.3.1.1^3,7]decane"),
        # Same, with a halogen substituent -- proves `halogens` threads
        # through find_tricyclic_core/name_tricycloalkane the same way it
        # does through _bicyclic.py (see this repo's history for the exact
        # bug this guards against: halogens dropped by a bicyclic/tricyclic
        # merge that didn't thread the map through).
        ("ClC12CC3CC(CC(C3)C1)C2", "1-chlorotricyclo[3.3.1.1^3,7]decane"),
        # Built directly from the Blue Book's own P-23.2.5.2 worked example:
        # both secondary-bridge attachment points (locants 2 and 5) fall on
        # the *same* 4-atom main-ring segment, so the branch-atom multigraph
        # has a doubled edge between the two main bridgeheads (the two
        # length-2 bridges b, c) and another doubled edge between the two
        # secondary bridgeheads (the interior of the length-4 segment a, and
        # the independent secondary bridge d itself) -- not the six distinct
        # single bridges (K4) adamantane/twistane have. This is the exact
        # case `test_secondary_bridgeheads_on_same_main_bridge_raises` used
        # to assert was out of scope; the general composite-bridge search in
        # `name_tricycloalkane` now covers it too.
        ("C1CC2CCC1C1CCC2CC1", "tricyclo[4.2.2.2^2,5]dodecane"),
        # perhydroanthracene: three fused cyclohexanes in a row (ortho-fused
        # "chain"), C14H24. Abstractly the same doubled-main-bridgeheads
        # branch-atom multigraph as the case above, just with different
        # bridge lengths -- both main bridgeheads and both secondary
        # bridgeheads directly bonded (c=0, d=0), giving an 8+4-atom main
        # ring using every skeletal atom. Not independently cross-checked
        # against a von Baeyer name database (PubChem's computed IUPAC name
        # for this compound, CID 93034, uses hydro-prefixed fusion
        # nomenclature -- "tetradecahydroanthracene" -- not a tricyclo[...]
        # von Baeyer name, so there was nothing to compare against); verified
        # instead by hand-deriving the branch-atom multigraph and by the
        # a+b+c+d+2 atom count matching the C14H24 formula.
        ("C1CCC2CC3CCCCC3CC2C1", "tricyclo[8.4.0.0^3,8]tetradecane"),
        # perhydrophenanthrene: same three-fused-cyclohexane formula, angular
        # instead of linear fusion -- same abstract branch-atom multigraph
        # again, but the secondary bridgeheads land at different locants (2
        # and 7 instead of 3 and 8), confirming linear vs. angular fusion is
        # distinguished correctly.
        ("C1CCC2C(C1)CCC1CCCCC21", "tricyclo[8.4.0.0^2,7]tetradecane"),
    ],
)
def test_smiles_to_iupac_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_bicyclic_is_not_tricyclic():
    # bicyclo[2.2.2]octane: genuinely bicyclic (cyclomatic number 2), still
    # resolved by _bicyclic.py unaffected by tricyclic support.
    assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"


def test_two_separate_rings_still_raises():
    # two cyclohexane rings joined by a single bond: share no atom at all,
    # not a polycyclic ring system by any of these modules' definitions.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCC1C1CCCCC1")


def test_pentagonal_prism_is_not_tricyclic():
    # pentaprismane: a genuinely hexacyclic system (cyclomatic number 6, ten
    # branch atoms of degree 3) -- not mistaken for tricyclic, and
    # correctly resolved via _polycyclic.py's hexacyclic support instead
    # (see tests/test_hexacyclic.py).
    assert smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12") == "hexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # [1.1.1]propellane: two branch atoms, both of degree 4, directly
        # bonded to each other in addition to three one-carbon bridges. The
        # direct bond is the independent secondary bridge (P-23.2.5.1),
        # length 0, attached at the two main bridgeheads themselves (1 and
        # 3, since each of the two main-ring bridges contributes one atom:
        # 1 -> 2 -> 3). Verified against Wikipedia ("1,1,1-Propellane"),
        # ChemSpider (CID 125285), and the ACS "Molecule of the Week"
        # writeup, all of which give tricyclo[1.1.1.0^1,3]pentane; also
        # cross-checked C5H6 via RDKit's computed molecular formula.
        ("C1C23CC12C3", "tricyclo[1.1.1.0^1,3]pentane"),
        # [2.2.2]propellane: same topology with three two-carbon bridges
        # instead of one-carbon ones, so the main bridgeheads land at 1 and
        # 4 instead of 1 and 3. This SMILES was built from scratch with
        # RDKit's RWMol (two bridgehead atoms bonded directly, plus three
        # explicit two-atom C-C bridges between them) rather than copied
        # from a database, then confirmed to reduce to that exact
        # propellane branch-atom graph. Verified against Wikipedia
        # ("2,2,2-Propellane") and Wikidata (Q4596979), which give
        # tricyclo[2.2.2.0^1,4]octane; also cross-checked C8H12 via RDKit's
        # computed molecular formula.
        ("C1CC23CCC12CC3", "tricyclo[2.2.2.0^1,4]octane"),
        # Secondary bridge with one atom of its own instead of a direct
        # bond: the two main bridgeheads (both degree 4) are joined by three
        # one-carbon main-ring/main-bridge segments *and* a fourth,
        # one-carbon secondary bridge (rather than being bonded directly).
        # Verified against PubChem CID 59850615, which gives exactly
        # tricyclo[1.1.1.1^1,3]hexane for this SMILES; C6H8 cross-checked
        # via RDKit's computed molecular formula.
        ("C1C23CC1(C2)C3", "tricyclo[1.1.1.1^1,3]hexane"),
    ],
)
def test_smiles_to_iupac_propellane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
