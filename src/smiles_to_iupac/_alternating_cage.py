"""Von Baeyer and spiro systems of one skeletal heteroatom (P-52.1.6.1: bicyclo[4.2.1]nonasilane) or of two alternating
ones (P-23.5, P-24.2.4.3, P-52.1.6.2): the preselected name cites the hydrocarbon descriptor, the number of atoms of the
element cited first in the 'a' term, that term, and the parent hydride of the other element: bicyclo[3.3.1]tetrasiloxane,
spiro[5.5]pentasiloxane, 1N-tricyclo[3.3.1.1^2,4]pentasilazane. The skeleton is numbered as the hydrocarbon; locant 1
is named by '1Si-'/'1N-' where either could be there."""

import re

from rdkit import Chem

from ._numerals import numerical_term

_SENIORITY = ["O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl"]
_A_TERM = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "aluma",
    "Ga": "galla", "In": "inda", "Tl": "thalla",
}
_PARENT = {
    "O": "oxane", "S": "thiane", "Se": "selenane", "Te": "telluranane", "N": "azane", "P": "phosphane", "As": "arsane",
    "Sb": "stibane", "Bi": "bismuthane", "Si": "silane", "Ge": "germane", "Sn": "stannane", "Pb": "plumbane",
    "B": "borane",
}
_HYDRIDE = {
    "O": "oxidane", "S": "sulfane", "Se": "selane", "Te": "tellane", "N": "azane", "P": "phosphane", "As": "arsane",
    "Sb": "stibane", "Bi": "bismuthane", "Si": "silane", "Ge": "germane", "Sn": "stannane", "Pb": "plumbane", "B": "borane",
}
_VALENCE = {"O": 2, "S": 2, "Se": 2, "Te": 2, "N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4, "B": 3}
_RETAINED_HEAD = {"cubane": "pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]"}
_HEAD = re.compile(r"^((?:bi|tri|tetra|penta|hexa|hepta|octa|nona)cyclo\[[^\]]+\]|(?:di|tri|tetra|penta)?spiro\[[^\]]+\])")


def _shape(mol):
    """(first-cited element, second element, count of the first, bridgehead elements) or None."""
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings() < 2:
        return None
    atoms = list(mol.GetAtoms())
    elements = {a.GetSymbol() for a in atoms}
    if len(elements) == 1 and elements <= set(_HYDRIDE):
        (element,) = elements
        if any(
            not a.IsInRing() or a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons()
            or a.GetDegree() + a.GetTotalNumHs() != _VALENCE[element]
            for a in atoms
        ) or any(b.GetBondTypeAsDouble() not in (1.0, 1.5) for b in mol.GetBonds()):
            return None
        return element, None, len(atoms), set()
    if len(elements) != 2 or not elements <= set(_PARENT):
        return None
    if any(not a.IsInRing() or a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in atoms):
        return None
    if any(
        bond.GetBondTypeAsDouble() != 1.0 or bond.GetBeginAtom().GetSymbol() == bond.GetEndAtom().GetSymbol()
        for bond in mol.GetBonds()
    ):
        return None
    if any(a.GetDegree() + a.GetTotalNumHs() != _VALENCE[a.GetSymbol()] for a in atoms):
        return None
    junior, senior = sorted(elements, key=_SENIORITY.index, reverse=True)
    heads = {a.GetSymbol() for a in atoms if a.GetDegree() >= 3}
    return junior, senior, sum(a.GetSymbol() == junior for a in atoms), heads


def _carbon_skeleton_name(mol):
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetAtomicNum(6)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
        atom.SetIsAromatic(False)
    for bond in editable.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    hydrocarbon = editable.GetMol()
    Chem.SanitizeMol(hydrocarbon)
    return smiles_to_iupac(Chem.MolToSmiles(hydrocarbon))


def _retained_head(carbon_name, extra=None):
    """The von Baeyer or spiro prefix and descriptor of a carbon skeleton name, or None for any other name."""
    head = {**_RETAINED_HEAD, **(extra or {})}.get(carbon_name)
    if head is None:
        match = _HEAD.match(carbon_name)
        head = match.group(1) if match else None
    return head


def _locant_one_elements(mol):
    """The elements that can occupy locant 1 over every von Baeyer numbering of the cage the name's descriptor allows
    (P-23.2), empty when the cage is not a bicyclic or polycyclic von Baeyer system."""
    from ._bicyclic import find_bicyclic_core, iter_bicyclic_numberings
    from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates

    rings = mol.GetNumBonds() - mol.GetNumAtoms() + 1
    if rings == 2:
        core = find_bicyclic_core(mol)
        orders = list(iter_bicyclic_numberings(core)) if core is not None else []
    else:
        core = find_polycyclic_core(mol, rings) if rings in (3, 4, 5, 6) else None
        candidates = list(iter_polycyclic_candidates(core, rings)) if core is not None else []
        best = min((outer for _, _, outer in candidates), default=None)
        orders = [order for order, _, outer in candidates if outer == best]
    return {mol.GetAtomWithIdx(order[0]).GetSymbol() for order in orders}


def _exchange_is_another_compound(mol, first, second):
    """Whether exchanging the two elements gives a compound different from `mol` (P-23.5.2: only then must the atom at
    locant 1 be specified)."""
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in (first, second):
            atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(second if symbol == first else first))
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    exchanged = editable.GetMol()
    try:
        Chem.SanitizeMol(exchanged)
    except Chem.MolSanitizeException:
        return False
    return Chem.MolToSmiles(exchanged) != Chem.MolToSmiles(mol)


def has_alternating_cage_shape(mol) -> bool:
    return _shape(mol) is not None


def name_alternating_cage(mol):
    found = _shape(mol)
    junior, senior, count, heads = found
    carbon_name = _carbon_skeleton_name(mol)
    if senior is None:
        head = _retained_head(carbon_name, {"adamantane": "tricyclo[3.3.1.1^3,7]"})
        return None if head is None else f"{head}{numerical_term(count)}{_HYDRIDE[junior]}"
    if carbon_name == "adamantane":
        # P-68.2.1.2 prints the replacement name of an adamantane skeleton as the preselected one, against P-52.1.6.2
        return None
    head = _retained_head(carbon_name)
    if head is None:
        return None
    parent = _PARENT[senior]
    a_term = _A_TERM[junior]
    if parent[0] == "a":
        joined = a_term + parent[1:]
    elif parent[0] in "eiou":
        joined = a_term[:-1] + parent
    else:
        joined = a_term + parent
    prefix = ""
    if len(heads) == 1:
        (bridgehead,) = heads
        other = senior if bridgehead == junior else junior
        if _VALENCE[other] >= max(a.GetDegree() for a in mol.GetAtoms()):
            prefix = f"1{bridgehead}-"
    elif len(heads) == 2 and _exchange_is_another_compound(mol, junior, senior):
        elements = _locant_one_elements(mol)
        if len(elements) == 1:
            prefix = f"1{next(iter(elements))}-"
    multiplier = numerical_term(count) if count > 1 else ""
    return f"{prefix}{head}{multiplier}{joined}"
