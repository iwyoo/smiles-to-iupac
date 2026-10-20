"""Names of polyatomic inorganic molecules cited as adduct components (P-14.8.2, ref. 12): a homonuclear molecule takes
the multiplied element name, a binary molecule is named stoichiometrically with the element cited later in the
electronegativity sequence as the '-ide', and the acids with retained names keep them."""

from collections import Counter

from rdkit import Chem

from ._numerals import numerical_term

_SEQUENCE = (86, 54, 36, 18, 10, 2, 5, 14, 6, 51, 33, 15, 7, 1, 52, 34, 16, 85, 53, 35, 17, 8, 9)
_ELEMENT = {
    5: "boron", 6: "carbon", 7: "nitrogen", 8: "oxygen", 9: "fluorine", 14: "silicon", 15: "phosphorus", 16: "sulfur",
    17: "chlorine", 33: "arsenic", 34: "selenium", 35: "bromine", 51: "antimony", 52: "tellurium", 53: "iodine",
    85: "astatine",
}
_IDE = {
    5: "boride", 6: "carbide", 7: "nitride", 8: "oxide", 9: "fluoride", 14: "silicide", 15: "phosphide",
    16: "sulfide", 17: "chloride", 33: "arsenide", 34: "selenide", 35: "bromide", 51: "antimonide", 52: "telluride",
    53: "iodide", 85: "astatide",
}
_CHALCOGENS = {8, 16, 34, 52}
_RETAINED_ACIDS = {
    "N=[N+]=[N-]": "hydrazoic acid",
    "N=C=O": "isocyanic acid",
    "OC#N": "cyanic acid",
    "C#[N+][O-]": "fulminic acid",
    "N=C=S": "isothiocyanic acid",
    "SC#N": "thiocyanic acid",
}


def _skeleton(frag):
    """Connectivity and hydrogen counts of a molecule, whatever its bond orders and charge separation."""
    graph = Chem.RWMol()
    for atom in frag.GetAtoms():
        copy = Chem.Atom(atom.GetAtomicNum())
        copy.SetNumExplicitHs(atom.GetTotalNumHs())
        copy.SetNoImplicit(True)
        graph.AddAtom(copy)
    for bond in frag.GetBonds():
        graph.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)
    graph.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(graph)


_RETAINED_SKELETONS = {_skeleton(Chem.MolFromSmiles(smiles)): name for smiles, name in _RETAINED_ACIDS.items()}


def molecule_name(frag):
    """Name of a neutral polyatomic molecule of nonmetals (P-14.8.2), else None."""
    atoms = list(frag.GetAtoms())
    if len(atoms) < 2 or any(a.GetIsotope() for a in atoms) or sum(a.GetFormalCharge() for a in atoms):
        return None
    if any(a.GetAtomicNum() not in _ELEMENT for a in atoms):
        return None
    if any(a.GetTotalNumHs() for a in atoms):
        return _RETAINED_SKELETONS.get(_skeleton(frag))
    counts = Counter(a.GetAtomicNum() for a in atoms)
    if len(counts) == 1:
        (z, n), = counts.items()
        return numerical_term(n) + _ELEMENT[z]
    if len(counts) != 2:
        return None
    first, second = sorted(counts, key=_SEQUENCE.index)
    if first == 6 and (counts[6] != 1 or second not in _CHALCOGENS):
        return None
    head = _ELEMENT[first] if counts[first] == 1 else numerical_term(counts[first]) + _ELEMENT[first]
    tail = numerical_term(counts[second]) + _IDE[second]
    return f"{head} {tail.replace('monooxide', 'monoxide')}"


def is_inorganic_molecule(frag):
    return not any(a.GetAtomicNum() == 6 for a in frag.GetAtoms()) or molecule_name(frag) is not None
