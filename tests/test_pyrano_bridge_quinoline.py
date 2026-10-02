from rdkit import Chem
from rdkit.Chem import RWMol, Atom, BondType

from smiles_to_iupac import smiles_to_iupac


def _epipyranobenzoquinoline_smiles():
    # Constructed from the Blue Book's own diagram (tmp/bluebook/P2.pdf
    # p.112): benzo[g]quinoline with pyran's para C2/C5 bonded to its
    # meso locants 5 and 10. No PubChem CID exists for this worked example.
    base = Chem.MolFromSmiles("c1ccc2cc3ncccc3cc2c1")
    rw = RWMol(base)
    o1 = rw.AddAtom(Atom(8))
    c2 = rw.AddAtom(Atom(6))
    c3 = rw.AddAtom(Atom(6))
    c4 = rw.AddAtom(Atom(6))
    c5 = rw.AddAtom(Atom(6))
    c6 = rw.AddAtom(Atom(6))
    rw.AddBond(o1, c2, BondType.SINGLE)
    rw.AddBond(c2, c3, BondType.DOUBLE)
    rw.AddBond(c3, c4, BondType.SINGLE)
    rw.AddBond(c4, c5, BondType.DOUBLE)
    rw.AddBond(c5, c6, BondType.SINGLE)
    rw.AddBond(c6, o1, BondType.SINGLE)
    rw.AddBond(c2, 11, BondType.SINGLE)
    rw.AddBond(c5, 4, BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def test_epipyrano_bridge_benzo_g_quinoline():
    assert smiles_to_iupac(_epipyranobenzoquinoline_smiles()) == "12H-5,10-[2,5]epipyranobenzo[g]quinoline"
