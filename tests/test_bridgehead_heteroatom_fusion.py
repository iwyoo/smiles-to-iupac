import pytest
from rdkit import Chem

from smiles_to_iupac._common import UnsupportedStructure
from smiles_to_iupac._bridgehead_heteroatom_fusion import (
    has_bridgehead_heteroatom_fusion_name,
    name_bridgehead_heteroatom_fusion,
)


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # PubChem CID 817024 -- the Blue Book's own P-25.3.2.5.1 worked
        # example: imidazole (N, N) vs thiazole (N, S) tie through
        # P-25.3.2.4(a)-(d); only (e) "greater variety of heteroatoms"
        # decides it (thiazole's N+S beats imidazole's N+N).
        ("C1=CN2C=CSC2=N1", "imidazo[2,1-b]thiazole"),
        # PubChem CID 19832717 -- pyrazole (N, N) vs thiazole (N, S), same
        # (e)-decided seniority as above.
        ("C1=CSC2=CC=NN21", "pyrazolo[5,1-b]thiazole"),
        # PubChem CID 78960 -- imidazole (5-ring) vs pyridine (6-ring):
        # decided by (c) "the larger ring at the first point of
        # difference" before (d)/(e) are even reached.
        ("C1=CC2=NC=CN2C=C1", "imidazo[1,2-a]pyridine"),
        # PubChem CID 577018 -- imidazole vs pyrimidine, (c) decides again.
        ("C1=CN2C=CN=C2N=C1", "imidazo[1,2-a]pyrimidine"),
        # PubChem CID 67507 -- pyrazole vs pyridine, (c) decides.
        ("C1=CC2=CC=NN2C=C1", "pyrazolo[1,5-a]pyridine"),
        # PubChem CID 136599 -- imidazole vs pyridazine; the fusion bond
        # lands on pyridazine's *wraparound* bond (locants {1, n}), which
        # is the earliest letter by *position* even though its own locant
        # set isn't numerically the lowest - this is the case that caught
        # an earlier version of this module's letter-selection logic
        # picking by raw locant-tuple order instead of bond position.
        ("C1=CC2=NC=CN2N=C1", "imidazo[1,2-b]pyridazine"),
        # PubChem CID 2771670 -- imidazole vs pyrazine.
        ("C1=CN2C=CN=C2C=N1", "imidazo[1,2-a]pyrazine"),
        # PubChem CID 11636795 -- pyrazole vs pyrimidine.
        ("C1=CN2C(=CC=N2)N=C1", "pyrazolo[1,5-a]pyrimidine"),
        # PubChem CID 577096 -- imidazole vs pyrimidine, a different
        # fusion bond than CID 577018 above (tests the attached
        # component's own numbering-direction choice independently of the
        # base's).
        ("C1=CN2C=NC=C2N=C1", "imidazo[1,5-a]pyrimidine"),
        # PubChem CID 9548860 -- imidazole vs pyrimidine, letter 'c'.
        ("C1=CN=CN2C1=NC=C2", "imidazo[1,2-c]pyrimidine"),
    ],
)
def test_bridgehead_heteroatom_fusion_matches_pubchem(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    assert has_bridgehead_heteroatom_fusion_name(mol)
    assert name_bridgehead_heteroatom_fusion(mol) == expected


def test_substituent_out_of_scope():
    # 2-methylimidazo[1,2-a]pyridine (PubChem CID 136742).
    mol = Chem.MolFromSmiles("CC1=CN2C=CC=CC2=N1")
    assert not has_bridgehead_heteroatom_fusion_name(mol)


def test_three_heteroatom_base_out_of_scope():
    # imidazo[2,1-b][1,3,4]thiadiazole (PubChem CID 18625232): the base
    # ring (1,3,4-thiadiazole) has 3 heteroatoms (S, N, N), not the 2 this
    # module's `_RING_TYPES` patterns support - correctly rejected rather
    # than guessed at.
    mol = Chem.MolFromSmiles("C1=CN2C(=N1)SC=N2")
    assert not has_bridgehead_heteroatom_fusion_name(mol)
    with pytest.raises(UnsupportedStructure):
        name_bridgehead_heteroatom_fusion(mol)


def test_no_shared_heteroatom_out_of_scope():
    # Naphthalene: shared fusion atoms are both plain carbon, not this
    # module's bridgehead-heteroatom shape at all.
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    assert not has_bridgehead_heteroatom_fusion_name(mol)


def test_furo_pyridine_not_misrouted_here():
    # furo[3,2-b]pyridine (from #595): the shared fusion atoms are both
    # carbon (furan's own oxygen isn't at the fusion bond) - this module
    # must not claim it.
    mol = Chem.MolFromSmiles("C1=CC2=C(C=CO2)N=C1")
    assert not has_bridgehead_heteroatom_fusion_name(mol)
