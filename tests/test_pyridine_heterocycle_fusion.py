import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._pyridine_heterocycle_fusion import (
    has_pyridine_heterocycle_fusion_name,
    name_pyridine_heterocycle_fusion,
)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 12210217
        ("C1=CC2=C(C=CO2)N=C1", "furo[3,2-b]pyridine"),
        # PubChem CID 12421098
        ("C1=CC2=C(N=C1)OC=C2", "furo[2,3-b]pyridine"),
        # PubChem CID 289928
        ("C1=CC2=C(N=C1)SC=C2", "thieno[2,3-b]pyridine"),
        # PubChem CID 12210218
        ("C1=CC2=C(C=CS2)N=C1", "thieno[3,2-b]pyridine"),
        # PubChem CID 12826108 -- confirms the base letter isn't fixed at 'b'
        ("C1=CN=CC2=C1C=CO2", "furo[2,3-c]pyridine"),
    ],
)
def test_pyridine_heterocycle_fusion_matches_pubchem(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    assert has_pyridine_heterocycle_fusion_name(mol)
    assert name_pyridine_heterocycle_fusion(mol) == expected


def test_pyrrole_out_of_scope_needs_indicated_hydrogen():
    # 1H-pyrrolo[2,3-b]pyridine (PubChem CID 9222) -- excluded because its
    # correct name needs an indicated-hydrogen citation ("1H-") this
    # module doesn't implement (see module docstring); rejecting is
    # correct, a silently-incomplete name would not be.
    mol = Chem.MolFromSmiles("C1=CC2=C(NC=C2)N=C1")
    assert not has_pyridine_heterocycle_fusion_name(mol)
    with pytest.raises(UnsupportedStructure):
        name_pyridine_heterocycle_fusion(mol)


def test_bridgehead_nitrogen_unsupported():
    # Indolizine (PubChem CID 9230): the shared fusion atom itself is
    # pyridine's own nitrogen, a bridgehead-N shape this module rejects
    # (the attached ring has no O/S heteroatom of its own once the shared
    # N is accounted for as pyridine's).
    mol = Chem.MolFromSmiles("C1=CC2=CC=CN2C=C1")
    assert not has_pyridine_heterocycle_fusion_name(mol)


def test_fusion_not_reachable_from_attached_heteroatom_unsupported():
    # furo[3,4-b]pyridine (PubChem CID 18442745): furan fused at its own
    # 3,4-bond, neither atom adjacent to furan's own oxygen. Matches the
    # coarse ring-shape check (same as the sibling
    # `_two_component_heterocycle_fusion.py` module's `has_...`), but
    # `name_...` rejects it once it looks at which atoms are fused.
    mol = Chem.MolFromSmiles("C1=CC2=COC=C2N=C1")
    assert has_pyridine_heterocycle_fusion_name(mol)
    with pytest.raises(UnsupportedStructure):
        name_pyridine_heterocycle_fusion(mol)


def test_wrong_ring_sizes_unsupported():
    # Naphthalene: two 6-membered rings, not a 5+6 pyridine/heterocycle pair.
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    assert not has_pyridine_heterocycle_fusion_name(mol)


def test_substituent_out_of_scope():
    mol = Chem.MolFromSmiles("Cc1ccc2occc2n1")
    assert not has_pyridine_heterocycle_fusion_name(mol)
