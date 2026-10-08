"""Isotopic modification named by the chain engine (P-82.2.1, P-82.3, P-82.5.2, P-82.6): the molecule is named without
its nuclide labels, then the descriptor '(locants-nuclide count)' is inserted before the parent hydride name with the
locants of the parent's own numbering; a nuclide on a characteristic-group atom takes a letter locant (P-82.6.2)."""

from rdkit import Chem

from ._common import UnsupportedStructure

_HYDROGEN_MASS = {2: "2H", 3: "3H"}
HYDROGEN_ISOTOPES = {nuclide: mass for mass, nuclide in _HYDROGEN_MASS.items()}


def split_isotopes(mol):
    """(unlabelled mol, labels, index map) where labels maps an atom index of the new mol to {"skeleton": nuclide or None,
    "H": {nuclide: count}}; None when `mol` has no isotope label."""
    if not any(a.GetIsotope() for a in mol.GetAtoms()):
        return None
    editable = Chem.RWMol(mol)
    labels = {}
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
        else:
            labels.setdefault(atom.GetIdx(), {"skeleton": None, "H": {}})["skeleton"] = f"{atom.GetIsotope()}{atom.GetSymbol()}"
            editable.GetAtomWithIdx(atom.GetIdx()).SetIsotope(0)
    parameters = Chem.RemoveHsParameters()
    parameters.removeIsotopes = True
    parameters.removeDefiningBondStereo = True
    clean = Chem.RemoveHs(editable.GetMol(), parameters)
    new_index = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            new_index[atom.GetIdx()] = len(new_index)
    if clean.GetNumAtoms() != len(new_index):
        raise UnsupportedStructure("this molecule keeps hydrogen atoms that cannot be dropped to name its isotopic labels")
    return clean, {new_index[a]: entry for a, entry in labels.items()}, new_index


def _nuclide_sort_key(nuclide):
    mass = int("".join(ch for ch in nuclide if ch.isdigit()))
    symbol = "".join(ch for ch in nuclide if ch.isalpha())
    return symbol, mass


def modification_key(labels, position_of):
    """Sort key of a candidate parent or chain by its isotopic modification (P-44.4.1.11.1-6, P-46.2): more modified
    atoms, then more nuclides of higher atomic number and of higher mass number, then the lowest locants for the
    modified atoms and for those nuclides; the smaller key is senior."""
    table = Chem.GetPeriodicTable()

    def nuclides(entry):
        found = [entry["skeleton"]] if entry["skeleton"] else []
        for nuclide, count in entry["H"].items():
            found.extend([nuclide] * count)
        return found

    order = sorted(
        {n for entry in labels.values() for n in nuclides(entry)},
        key=lambda n: (-table.GetAtomicNumber(_nuclide_sort_key(n)[0]), -_nuclide_sort_key(n)[1]),
    )
    modified = [(position_of[atom], n) for atom, entry in labels.items() if atom in position_of for n in nuclides(entry)]
    return (
        -len(modified),
        tuple(-sum(1 for _, m in modified if m == n) for n in order),
        tuple(sorted(locant for locant, _ in modified)),
        tuple(tuple(sorted(locant for locant, m in modified if m == n)) for n in order),
    )


def descriptor(labels, position_of, single_position, extra=(), capacity=None, sole=frozenset()):
    """The '(…)' isotopic descriptor for the labelled parent atoms; `single_position` drops the locants for a
    one-atom parent (P-82.6.1.1). `extra`: (nuclide, locant text or None, count, repeatable) for atoms outside the
    parent, whose letter locant (P-82.6.2) is cited once however many atoms carry the nuclide; the count is a
    subscript only where several atoms can be substituted at that position (P-82.2.1). `capacity`: {atom: number of
    hydrogens the atom carries unmodified}, so a position with one hydrogen takes no subscript for a single atom."""
    if any(a not in position_of for a in labels):
        raise UnsupportedStructure("an isotopically modified atom outside the parent hydride is not supported yet")
    located = {}
    several = set()
    for atom, entry in labels.items():
        if entry["skeleton"]:
            located.setdefault(entry["skeleton"], []).append(position_of[atom])
        for nuclide, count in entry["H"].items():
            located.setdefault(nuclide, []).extend([position_of[atom]] * count)
            if capacity is None or capacity.get(atom, 2) > 1:
                several.add(nuclide)
    letters = {}
    for nuclide, locant, count, repeatable in extra:
        entry = letters.setdefault(nuclide, {"locants": [], "count": 0, "repeatable": False})
        entry["repeatable"] |= repeatable
        if locant is not None and locant not in entry["locants"]:
            entry["locants"].append(locant)
        entry["count"] += count
    groups = []
    for nuclide in sorted(set(located) | set(letters), key=_nuclide_sort_key):
        numbers = sorted(located.get(nuclide, []))
        outside = letters.get(nuclide, {"locants": [], "count": 0, "repeatable": False})
        count = len(numbers) + outside["count"]
        counted = nuclide[-1] == "H" and (nuclide in several or outside["repeatable"])
        subscript = str(count) if count > 1 or counted else ""
        symbol = f"{nuclide}{subscript}"
        locants = [str(x) for x in numbers] + sorted(outside["locants"])
        unlocated = (single_position or nuclide in sole) and not outside["locants"]
        groups.append(symbol if unlocated or not locants else f"{','.join(locants)}-{symbol}")
    return "(" + ",".join(groups) + ")"
