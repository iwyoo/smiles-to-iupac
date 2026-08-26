import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # cubane (flagship case): C8H8, eight branch atoms all of degree 3
        # (verified via RDKit -- CalcMolFormula and per-atom GetDegree()),
        # the cube graph. Cross-checked externally against PubChem CID
        # 136090's computed IUPAC name via the PUG REST API
        # (https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/136090/
        # property/CanonicalSMILES,IUPACName/JSON), not taken from memory.
        ("C12C3C4C1C1C2C3C41", "pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane"),
        # 1-methylcubane: all eight vertices of unsubstituted cubane are
        # equivalent, so the lowest-locant tie-break (P-14.4/P-45.2) puts the
        # substituted one at position 1, exactly as it does for adamantane's
        # equivalent bridgeheads (see tests/test_tricyclic.py).
        ("CC12C3C4C1C1C2C3C41", "1-methylpentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane"),
        # Same, with a halogen substituent -- proves `halogens` threads
        # through find_polycyclic_core/name_polycycloalkane the same way it
        # does through the tricyclic/tetracyclic cases this engine also
        # covers (see this repo's history for the exact bug this class of
        # test guards against).
        ("ClC12C3C4C1C1C2C3C41", "1-chloropentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane"),
    ],
)
def test_smiles_to_iupac_pentacyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_tricyclic_is_not_pentacyclic():
    # adamantane: genuinely tricyclic (cyclomatic number 3, four branch
    # atoms), still resolved as such unaffected by pentacyclic support.
    assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "tricyclo[3.3.1.1^3,7]decane"


def test_tetracyclic_is_not_pentacyclic():
    # quadricyclane: genuinely tetracyclic (cyclomatic number 4, six branch
    # atoms), still resolved as such unaffected by pentacyclic support.
    assert smiles_to_iupac("C1C2C3C2C4C1C34") == "tetracyclo[3.2.0.0^2,7.0^4,6]heptane"


def test_pentacyclic_propellane_like_degree_four_branch_atom_raises():
    # Cubane with one extra direct bond added between two face-diagonal
    # branch atoms: still eight skeletal atoms, but now cyclomatic number 6
    # (hexacyclic) with two branch atoms pushed to degree 4 by the added
    # bond -- a propellane-like topology, out of scope regardless of ring
    # count (mirrors _tricyclic.py's own degree-4 exclusion for the
    # non-propellane-shaped case). Built from scratch via RDKit's RWMol
    # (cubane's cube graph plus one bond between vertices 0 and 2), verified
    # to give degree sequence two 4s / six 3s, cyclomatic number 6.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C12C3C4C1C15C2C31C45")
