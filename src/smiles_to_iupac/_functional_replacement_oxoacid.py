"""Mononuclear oxoacids of P, As, Sb and B modified by functional replacement (P-67.1.2.1, P-67.1.2.3, P-67.1.4.1.1.4):
the acid stem follows the number of chalcogen ligands on the centre (one 'phosphin-', two 'phosphon-', three 'phosphor-'),
'ic' for the =X form of the pentavalent centre and 'ous' for the trivalent one, carbon groups replace its hydrogens as
prefixes, the replacing S, Se, Te are infixes ('phosphinothious acid'), and the hydrogen-bearing elements are cited
when the ligands differ ('phosphorothioic O,O,O-acid')."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._numerals import numerical_term
from ._substituents import cited_stereo_around, format_mononuclear_prefixes, name_branch

_STEMS = {
    15: ("phosphin", "phosphon", "phosphor"),
    33: ("arsin", "arson", "arsor"),
    51: ("stibin", "stibon", "stibor"),
    5: ("borin", "boron", "bor"),
}
_PENTAVALENT = {15, 33, 51}
_CHALCOGENS = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_ORDER = {8: 0, 16: 1, 34: 2, 52: 3}
_HALIDE_INFIX = {9: "fluorid", 17: "chlorid", 35: "bromid", 53: "iodid"}
# pseudohalide groups bonded to the acid centre by their first atom, with the infix that replaces an OH (P-67.1.2.4.1.3)
_PSEUDOHALIDES = (
    ("[O;D2]-C#N", "cyanatid"),
    ("[N;D2]=C=O", "isocyanatid"),
    ("[S;D2]-C#N", "thiocyanatid"),
    ("[N;D2]=C=S", "isothiocyanatid"),
    ("[Se;D2]-C#N", "selenocyanatid"),
    ("[C;D2]#N", "cyanid"),
    ("[N;D2]#[C;D1]", "isocyanid"),
    ("[N;D2]=[N+]=[N-]", "azid"),
)
_PSEUDOHALIDE_PATTERNS = [(Chem.MolFromSmarts(smarts), infix) for smarts, infix in _PSEUDOHALIDES]


def _replacer(mol, center, neighbor):
    """(infix, atoms) of a halogen or pseudohalogen group bonded to the acid centre, else None."""
    if neighbor.GetAtomicNum() in _HALIDE_INFIX and neighbor.GetDegree() == 1:
        return _HALIDE_INFIX[neighbor.GetAtomicNum()], {neighbor.GetIdx()}
    for pattern, infix in _PSEUDOHALIDE_PATTERNS:
        for match in mol.GetSubstructMatches(pattern):
            if match[0] == neighbor.GetIdx():
                return infix, set(match)
    return None


def _ligands(mol, center):
    """([=X atoms], [XH atoms], [carbon neighbours]) of the acid centre, or None if it has any other neighbour."""
    oxo, hydroxy, carbon, replacers = [], [], [], []
    for neighbor in center.GetNeighbors():
        z = neighbor.GetAtomicNum()
        order = mol.GetBondBetweenAtoms(center.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble()
        if z in _CHALCOGENS and neighbor.GetDegree() == 1 and not neighbor.GetFormalCharge() and not neighbor.GetIsotope():
            if order == 2.0 and neighbor.GetTotalNumHs() == 0:
                oxo.append(neighbor)
            elif order == 1.0 and neighbor.GetTotalNumHs() == 1:
                hydroxy.append(neighbor)
            else:
                return None
        elif z == 6 and order == 1.0:
            carbon.append(neighbor)
        elif order == 1.0 and _replacer(mol, center, neighbor) is not None:
            replacers.append(neighbor)
        else:
            return None
    return oxo, hydroxy, carbon, replacers


def _acid_parts(mol):
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS]
    if len(centers) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    center = centers[0]
    if center.GetFormalCharge() or center.GetIsotope() or center.IsInRing():
        return None
    ligands = _ligands(mol, center)
    if ligands is None:
        return None
    oxo, hydroxy, carbon, replacers = ligands
    z = center.GetAtomicNum()
    pentavalent = z in _PENTAVALENT and len(oxo) == 1
    if len(oxo) > 1 or (len(oxo) == 1 and not pentavalent) or not 1 <= len(hydroxy) <= 3:
        return None
    if len(hydroxy) + len(replacers) + len(carbon) + center.GetTotalNumHs() != 3:
        return None
    chalcogens = [*oxo, *hydroxy]
    if all(a.GetAtomicNum() == 8 for a in chalcogens) and pentavalent and not replacers:
        return None
    ligand_atoms = {x.GetIdx() for x in chalcogens} | {center.GetIdx()}
    for neighbor in replacers:
        ligand_atoms |= _replacer(mol, center, neighbor)[1]
    if any(a.GetAtomicNum() not in (6, *HALOGEN_PREFIXES) for a in mol.GetAtoms() if a.GetIdx() not in ligand_atoms):
        return None
    return center, oxo, hydroxy, carbon, pentavalent, replacers


def has_functional_replacement_oxoacid_shape(mol) -> bool:
    return _acid_parts(mol) is not None


def name_functional_replacement_oxoacid(mol) -> str:
    parts = _acid_parts(mol)
    if parts is None:
        raise UnsupportedStructure("this is not a mononuclear oxoacid modified by functional replacement")
    center, oxo, hydroxy, carbon, pentavalent, replacers = parts
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    with cited_stereo_around(mol, center.GetIdx()):
        entries = [name_branch(graph, c.GetIdx(), center.GetIdx(), halogens, mol=mol) for c in carbon]

    chalcogens = [*oxo, *hydroxy]
    replaced = sorted((a.GetAtomicNum() for a in chalcogens if a.GetAtomicNum() != 8), key=lambda z: _INFIX[z])
    pieces = []
    for z in sorted(set(replaced), key=lambda z: _INFIX[z]):
        count = replaced.count(z)
        pieces.append(((numerical_term(count) if count > 1 else "") + _INFIX[z], _INFIX[z]))
    infixes = [_replacer(mol, center, n)[0] for n in replacers]
    for text in sorted(set(infixes)):
        count = infixes.count(text)
        pieces.append(((numerical_term(count) if count > 1 else "") + text, text))
    pieces.sort(key=lambda piece: piece[1])
    infix = "o".join(text for text, _ in pieces)
    z_center = center.GetAtomicNum()
    stem = _STEMS[z_center][len(hydroxy) + len(replacers) - 1]
    if z_center == 5 or pentavalent:
        word = stem + ("o" + infix if infix else "") + "ic acid"
    else:
        word = stem + ("o" + infix[:-1] if infix else "") + "ous acid"
    if len(chalcogens) >= 3 and len({a.GetAtomicNum() for a in chalcogens}) > 1:
        # the replacing atoms that carry hydrogen, else the oxygens that do (P-67.1.2.4.1)
        cited = [a.GetAtomicNum() for a in hydroxy if a.GetAtomicNum() != 8] or [a.GetAtomicNum() for a in hydroxy]
        word = word[: -len("acid")] + ",".join(_CHALCOGENS[z] for z in sorted(cited, key=lambda z: _ORDER[z])) + "-acid"
    return format_mononuclear_prefixes(entries) + word
