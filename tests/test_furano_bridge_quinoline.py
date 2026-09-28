from rdkit import Chem
from rdkit.Chem import RWMol, BondType

from smiles_to_iupac import smiles_to_iupac


def _furanobenzoquinoline_smiles():
    # Constructed from the Blue Book's own diagram (tmp/bluebook/P2.pdf
    # p.106): benzo[g]quinoline with furan's C2/C3 bonded to its meso
    # locants 10 and 5. No PubChem CID exists for this worked example.
    base = Chem.MolFromSmiles("C1=CC=C2C=C3C(=CC2=C1)C=CC=N3")
    furan = Chem.MolFromSmiles("c1ccoc1")
    combo = Chem.CombineMols(base, furan)
    rw = RWMol(combo)
    n_base = base.GetNumAtoms()
    rw.AddBond(4, n_base + 4, BondType.SINGLE)
    rw.AddBond(7, n_base + 0, BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def test_furano_bridge_benzo_g_quinoline():
    assert smiles_to_iupac(_furanobenzoquinoline_smiles()) == "10,5-[2,3]furanobenzo[g]quinoline"
