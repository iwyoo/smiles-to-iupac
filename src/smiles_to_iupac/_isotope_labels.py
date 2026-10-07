"""Isotopic modification of a parent hydride named by the chain engine (P-82.2.1, P-82.3, P-82.5.2): the molecule is
named without its nuclide labels, then the descriptor '(locants-nuclide count)' is inserted before the parent hydride
name, with the locants of the parent's own numbering."""

from rdkit import Chem

from ._common import UnsupportedStructure

_HYDROGEN_MASS = {2: "2H", 3: "3H"}
HYDROGEN_ISOTOPES = {nuclide: mass for mass, nuclide in _HYDROGEN_MASS.items()}


def split_isotopes(mol):
    """(unlabelled mol, labels) where labels maps an atom index of the new mol to {"skeleton": nuclide or None,
    "H": {nuclide: count}}; None when `mol` has no isotope label."""
    if not any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    editable = Chem.RWMol(mol)
    labels = {}
    removed = []
    for atom in mol.GetAtoms():
        if not atom.GetIsotope():
            continue
        if atom.GetAtomicNum() == 1:
            if atom.GetDegree() != 1 or atom.GetIsotope() not in _HYDROGEN_MASS:
                raise UnsupportedStructure("this isotopically labelled hydrogen is not a terminal hydrogen atom")
            neighbor = atom.GetNeighbors()[0].GetIdx()
            counts = labels.setdefault(neighbor, {"skeleton": None, "H": {}})["H"]
            nuclide = _HYDROGEN_MASS[atom.GetIsotope()]
            counts[nuclide] = counts.get(nuclide, 0) + 1
            removed.append(atom.GetIdx())
        else:
            labels.setdefault(atom.GetIdx(), {"skeleton": None, "H": {}})["skeleton"] = f"{atom.GetIsotope()}{atom.GetSymbol()}"
            editable.GetAtomWithIdx(atom.GetIdx()).SetIsotope(0)
    new_index = {}
    kept = 0
    for idx in range(mol.GetNumAtoms()):
        if idx not in removed:
            new_index[idx] = kept
            kept += 1
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    clean = editable.GetMol()
    Chem.SanitizeMol(clean)
    return clean, {new_index[a]: entry for a, entry in labels.items()}


def _nuclide_sort_key(nuclide):
    mass = int("".join(ch for ch in nuclide if ch.isdigit()))
    symbol = "".join(ch for ch in nuclide if ch.isalpha())
    return symbol, mass


def descriptor(labels, position_of, single_position):
    """The '(…)' isotopic descriptor for the labelled parent atoms; `single_position` drops the locants for a
    one-atom parent (P-82.6.1.1)."""
    if any(a not in position_of for a in labels):
        raise UnsupportedStructure("an isotopically modified atom outside the parent hydride is not supported yet")
    located = {}
    for atom, entry in labels.items():
        if entry["skeleton"]:
            located.setdefault(entry["skeleton"], []).append(position_of[atom])
        for nuclide, count in entry["H"].items():
            located.setdefault(nuclide, []).extend([position_of[atom]] * count)
    groups = []
    for nuclide in sorted(located, key=_nuclide_sort_key):
        locants = sorted(located[nuclide])
        count = len(locants)
        subscript = str(count) if count > 1 or nuclide[-1] == "H" else ""
        symbol = f"{nuclide}{subscript}"
        groups.append(symbol if single_position else f"{','.join(str(x) for x in locants)}-{symbol}")
    return "(" + ",".join(groups) + ")"
