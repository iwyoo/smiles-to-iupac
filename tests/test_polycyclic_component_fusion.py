import pytest
from rdkit import Chem

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


def test_benzo_g_indole():
    # PubChem CID 98617, structure-verified (formula C12H9N, properly
    # aromatic) -- an earlier, different CID (170319311, bare
    # "benzo[g]indole") is a malformed database entry and must not be
    # used (see module docstring).
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2NC=C3") == "1H-benzo[g]indole"


def test_benzo_e_1_benzofuran():
    # PubChem CID 9192, structure-verified.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2C=CO3") == "benzo[e][1]benzofuran"


def test_benzo_g_1_benzofuran():
    # PubChem CID 67475, structure-verified -- same fusion-bond letter
    # ('g') as the indole case above, an independent cross-check that the
    # periphery-walk algorithm is base-heteroatom-agnostic.
    assert smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2OC=C3") == "benzo[g][1]benzofuran"


def _fuse_new_ring(base_smiles, atom_a, atom_b):
    mol = Chem.RWMol(Chem.MolFromSmiles(base_smiles))
    new_idxs = []
    for _ in range(4):
        atom = Chem.Atom(6)
        atom.SetIsAromatic(True)
        new_idxs.append(mol.AddAtom(atom))
    mol.AddBond(atom_a, new_idxs[0], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[0], new_idxs[1], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[1], new_idxs[2], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[2], new_idxs[3], Chem.BondType.AROMATIC)
    mol.AddBond(new_idxs[3], atom_b, Chem.BondType.AROMATIC)
    result = mol.GetMol()
    Chem.SanitizeMol(result)
    return Chem.MolToSmiles(result)


def test_benzo_f_indole_reviewed_extension():
    # The middle bond of the same periphery stretch (C5-C6, atoms 0 and 1
    # in the indole reference SMILES "c1ccc2[nH]ccc2c1") -- no PubChem-
    # listed compound of its own, a reviewed extension by the same
    # mechanism already confirmed twice over for 'e' and 'g' (see module
    # docstring).
    smi = _fuse_new_ring("c1ccc2[nH]ccc2c1", 0, 1)
    assert smiles_to_iupac(smi) == "1H-benzo[f]indole"


def test_bare_indole_and_benzofuran_still_work():
    # Sanity check: bare indole/benzofuran must not be misrouted here.
    assert smiles_to_iupac("c1ccc2[nH]ccc2c1") == "1H-indole"
    assert smiles_to_iupac("c1ccc2occc2c1") == "1-benzofuran"


def test_substituted_variant_raises():
    # A methyl substituent on the new ring is out of scope for this
    # module (any substituent at all is rejected).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1=CC=C2C(=C1)C=CC3=C2NC(C)=C3")
