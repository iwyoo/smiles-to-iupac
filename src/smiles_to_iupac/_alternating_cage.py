"""Von Baeyer and spiro systems made of two alternating skeletal heteroatoms (P-23.5, P-24.2.4.3, P-52.1.6.2): the
preselected name cites the hydrocarbon descriptor, then the number of atoms of the element cited first in the 'a'
term, that term, and the mononuclear parent hydride of the other element: bicyclo[3.3.1]tetrasiloxane,
tricyclo[3.3.1.1^3,7]tetrasiloxane, 1N-tricyclo[3.3.1.1^2,4]pentasilazane, spiro[5.5]pentasiloxane. The skeleton is
numbered as the corresponding hydrocarbon; the atom at locant 1 is named by '1Si-'/'1N-' where either could be there."""

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
_VALENCE = {"O": 2, "S": 2, "Se": 2, "Te": 2, "N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4, "B": 3}
_RETAINED_HEAD = {"adamantane": "tricyclo[3.3.1.1^3,7]", "cubane": "pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]"}
_HEAD = re.compile(r"^((?:bi|tri|tetra|penta|hexa|hepta|octa|nona)cyclo\[[^\]]+\]|(?:di|tri|tetra|penta)?spiro\[[^\]]+\])")


def _shape(mol):
    """(first-cited element, second element, count of the first, bridgehead elements) or None."""
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetRingInfo().NumRings() < 2:
        return None
    atoms = list(mol.GetAtoms())
    elements = {a.GetSymbol() for a in atoms}
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


def has_alternating_cage_shape(mol) -> bool:
    return _shape(mol) is not None


def name_alternating_cage(mol):
    from .core import smiles_to_iupac

    found = _shape(mol)
    junior, senior, count, heads = found
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetAtomicNum(6)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    hydrocarbon = editable.GetMol()
    Chem.SanitizeMol(hydrocarbon)
    carbon_name = smiles_to_iupac(Chem.MolToSmiles(hydrocarbon))
    head = _RETAINED_HEAD.get(carbon_name)
    if head is None:
        match = _HEAD.match(carbon_name)
        if match is None:
            return None
        head = match.group(1)
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
    multiplier = numerical_term(count) if count > 1 else ""
    return f"{prefix}{head}{multiplier}{joined}"
