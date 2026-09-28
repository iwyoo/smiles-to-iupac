from rdkit import Chem

from smiles_to_iupac._quinoline_bicyclic_numbering import peripheral_numbering


def _locant_of(mol, atomic_num):
    numbering = peripheral_numbering(mol)
    (atom,) = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == atomic_num]
    return numbering, numbering[atom]


def test_benzo_g_quinoline_matches_bluebook_diagram():
    # tmp/bluebook/P2.pdf p.106: 10,5-[2,3]furanobenzo[g]quinoline diagram
    # gives these exact meso locants and N=1.
    mol = Chem.MolFromSmiles("C1=CC=C2C=C3C(=CC2=C1)C=CC=N3")
    numbering, n_locant = _locant_of(mol, 7)
    assert n_locant == "1"
    assert set(numbering.values()) == {
        "1", "2", "3", "4", "4a", "5", "5a", "6", "7", "8", "9", "9a", "10", "10a",
    }
    meso = {atom for atom, loc in numbering.items() if loc in ("5", "10")}
    assert all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in meso)


def test_benzo_g_isoquinoline_algorithm_consistent():
    # Same algorithm as the diagram-verified quinoline case above; no
    # second Blue Book diagram exists to cross-check isoquinoline's N
    # locant specifically, so this only pins the derived output.
    mol = Chem.MolFromSmiles("C1=CC=C2C=C3C=NC=CC3=CC2=C1")
    numbering, n_locant = _locant_of(mol, 7)
    assert n_locant == "2"
    assert set(numbering.values()) == {
        "1", "2", "3", "4", "4a", "5", "5a", "6", "7", "8", "9", "9a", "10", "10a",
    }


def test_angular_benzo_f_and_h_quinoline_out_of_scope():
    f = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2C=CC=N3")
    h = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=CC3=C2N=CC=C3")
    assert peripheral_numbering(f) is None
    assert peripheral_numbering(h) is None


def test_bare_quinoline_not_a_g_fusion():
    assert peripheral_numbering(Chem.MolFromSmiles("c1ccc2ncccc2c1")) is None
