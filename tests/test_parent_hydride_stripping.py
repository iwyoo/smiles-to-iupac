from rdkit import Chem
from rdkit.Chem import RWMol

from smiles_to_iupac._parent_hydride_stripping import strip_substituents
from smiles_to_iupac._steroid_parent_hydrides import _PLAIN_CANONICAL_TO_NAME, _locant_map


def _androstane_with_two_hydroxyls():
    """A synthetic androstane-3,17-diol-shaped molecule: two -OH oxygens
    added at the C3 and C17 ring carbons (one -CH2- position each, valence
    reduced by one to make room), built programmatically off the already-
    verified plain androstane parent so its own ring-atom indices are known
    -- proves `strip_substituents` handles N=2 simultaneous removals, not
    just #1144's original single-substituent case."""
    mol = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
    rw = RWMol(mol)
    # Atom 4 and atom 15 are both ring -CH2- carbons in this SMILES's own
    # atom order (independently confirmed via RDKit's atom map, not a
    # locant -- `_locant_map` below recovers the real steroid locants).
    c1, c2 = 4, 15
    for c in (c1, c2):
        atom = rw.GetAtomWithIdx(c)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
    o1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c1, o1, Chem.BondType.SINGLE)
    o2 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c2, o2, Chem.BondType.SINGLE)
    diol = rw.GetMol()
    Chem.SanitizeMol(diol)
    return diol, c1, c2, o1, o2


def test_strip_substituents_removes_two_atoms_and_reindexes():
    diol, c1, c2, o1, o2 = _androstane_with_two_hydroxyls()
    stripped, old_to_new = strip_substituents(diol, {o1, o2})
    assert stripped is not None
    assert stripped.GetNumAtoms() == diol.GetNumAtoms() - 2
    # Both original ring carbons survive, correctly reindexed.
    assert c1 in old_to_new and c2 in old_to_new
    assert old_to_new[c1] != old_to_new[c2]


def test_strip_substituents_two_atom_removal_matches_androstane_skeleton():
    diol, c1, c2, o1, o2 = _androstane_with_two_hydroxyls()
    stripped, old_to_new = strip_substituents(diol, {o1, o2})
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    assert name == "androstane"

    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    locant_c1 = locant_of_atom[old_to_new[c1]]
    locant_c2 = locant_of_atom[old_to_new[c2]]
    # Both stripped substituents' own ring-attachment locants are
    # correctly recovered simultaneously, not just one at a time.
    assert {locant_c1, locant_c2} <= set(range(1, 18))
    assert locant_c1 != locant_c2


def test_strip_substituents_single_atom_still_works():
    # Zero-regression check on the original N=1 case #1144 generalized from
    # (mirrors `steroid_suffix_name`'s own real usage: a ketone oxygen
    # double-bonded to a ring carbon).
    mol = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
    rw = RWMol(mol)
    ketone_c = 4
    atom = rw.GetAtomWithIdx(ketone_c)
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(atom.GetTotalNumHs() - 2)
    ketone_o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(ketone_c, ketone_o, Chem.BondType.DOUBLE)
    ketone_mol = rw.GetMol()
    Chem.SanitizeMol(ketone_mol)

    stripped, old_to_new = strip_substituents(ketone_mol, {ketone_o})
    assert stripped is not None
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    assert name == "androstane"
    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    assert locant_of_atom[old_to_new[ketone_c]] in range(1, 18)
