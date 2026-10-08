"""Acyl groups of mononuclear P, As and Sb oxoacids modified by functional replacement (P-67.1.4.1.1.4, P-67.1.4.1.1.5):
the infix name of the acid with the ending 'ic' changed to 'oyl', 'phosphorodichloridoyl', 'phosphorothioyl',
'N,N-dimethylphosphoramidoyl', with hydroxy, alkoxy and sulfanyl groups cited as prefixes of the acyl group."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_CENTRES = {15, 33, 51}
_OXO = {8, 16, 34, 52}
_HALOGENS = {9, 17, 35, 53}
_PREFIX_ELEMENTS = {8, 16, 34, 52}
_MONO = {15: "phosphono", 33: "arsono", 51: "stibono"}
_ACID_TAIL = re.compile(r" (?:[A-Za-z]+(?:,[A-Za-z]+)*-)?acid$")


_CHALCOGEN_RANK = {8: 0, 16: 1, 34: 2, 52: 3}
_HALOGEN_RANK = {9: 10, 17: 11, 35: 12, 53: 13}
_UNLISTED = 99


def _ligand_rank(mol, centre, n):
    atom = mol.GetAtomWithIdx(n)
    z = atom.GetAtomicNum()
    if z in _CHALCOGEN_RANK:
        return _CHALCOGEN_RANK[z]
    if z in _HALOGEN_RANK:
        return _HALOGEN_RANK[z]
    if z == 7:
        return 30 if all(m.GetAtomicNum() in (1, 6) for m in atom.GetNeighbors() if m.GetIdx() != centre) else 20
    if z == 6 and atom.GetDegree() == 2 and any(m.GetAtomicNum() == 7 for m in atom.GetNeighbors()):
        return 21
    return _UNLISTED


def centre_seniority_key(mol, atom):
    """P-67.1.5.2: the ligands of an acid centre in the order O, S, Se, Te, then F, Cl, Br, I, then the pseudohalides,
    then amides and hydrazides; the centre whose sorted list is the smallest is the senior one."""
    ranks = sorted(_ligand_rank(mol, atom.GetIdx(), n.GetIdx()) for n in atom.GetNeighbors())
    return tuple(ranks + [_UNLISTED] * (4 - len(ranks)))


def is_senior_centre(mol, atom):
    """Whether `atom` outranks every other P, As or Sb acid centre of `mol`."""
    key = centre_seniority_key(mol, atom)
    return all(
        centre_seniority_key(mol, other) > key
        for other in mol.GetAtoms()
        if other.GetAtomicNum() in _CENTRES and other.GetIdx() != atom.GetIdx()
    )


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _side(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _is_oxo(mol, centre, n):
    atom = mol.GetAtomWithIdx(n)
    if _bond(mol, centre, n) != 2.0 or atom.GetDegree() != 1 or atom.GetFormalCharge():
        return False
    return atom.GetAtomicNum() in _OXO or (atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() == 1)


def _ligand_role(mol, graph, centre, n):
    """'infix' for a halide, pseudohalide or amide-type group, 'organyl' for carbon, 'prefix' for an oxygen-type group."""
    atom = mol.GetAtomWithIdx(n)
    z = atom.GetAtomicNum()
    if atom.GetFormalCharge() or _bond(mol, centre, n) != 1.0:
        return None
    if z in _HALOGENS and atom.GetDegree() == 1:
        return "infix"
    if z == 6:
        if atom.GetDegree() == 2 and any(
            m.GetAtomicNum() == 7 and _bond(mol, n, m.GetIdx()) == 3.0 for m in atom.GetNeighbors()
        ):
            return "infix"
        return "organyl"
    if z == 7:
        others = [m for m in atom.GetNeighbors() if m.GetIdx() != centre]
        if all(m.GetAtomicNum() == 6 and _bond(mol, n, m.GetIdx()) == 1.0 for m in others):
            return "infix"
        if len(others) == 1 and others[0].GetAtomicNum() == 7 and others[0].GetFormalCharge() == 0:
            return "infix"
        if len(others) == 1 and others[0].GetAtomicNum() == 6 and _bond(mol, n, others[0].GetIdx()) == 2.0:
            return "infix"
        return None
    if z in _PREFIX_ELEMENTS:
        return "prefix"
    return None


def _retained_thio_group(mol, centre, oxo, ligands):
    """'thiophosphono' for -P(S)(OH)2 or -P(O)(OH)(SH) and 'trithiophosphono' for -P(S)(SH)2 (P-67.1.4.1.1.1)."""
    symbols = [mol.GetAtomWithIdx(oxo).GetSymbol()]
    for n in ligands:
        atom = mol.GetAtomWithIdx(n)
        if atom.GetAtomicNum() not in (8, 16) or atom.GetDegree() != 1 or atom.GetTotalNumHs() != 1:
            return None
        symbols.append(atom.GetSymbol())
    if set(symbols) - {"O", "S"}:
        return None
    thio = symbols.count("S")
    mono = _MONO[mol.GetAtomWithIdx(centre).GetAtomicNum()]
    return {1: "thio", 3: "trithio"}.get(thio, "") + mono if thio in (1, 3) else None


def infix_acyl_name(mol, graph, centre, attach, halogens, aromatic_atoms=None):
    """(name, compound) of the acyl group of the P, As or Sb `centre` attached through `attach` when functional replacement
    modifies the acid, else None."""
    from ._substituents import format_mononuclear_prefixes, name_branch
    from .core import smiles_to_iupac

    atom = mol.GetAtomWithIdx(centre)
    if atom.GetAtomicNum() not in _CENTRES or atom.GetFormalCharge() or atom.IsInRing() or atom.GetDegree() != 4:
        return None
    oxo, ligands = None, []
    for n in graph[centre]:
        if n == attach:
            continue
        if _is_oxo(mol, centre, n) and oxo is None:
            oxo = n
        else:
            ligands.append(n)
    if oxo is None or len(ligands) != 2:
        return None
    roles = {n: _ligand_role(mol, graph, centre, n) for n in ligands}
    if None in roles.values():
        return None
    retained = _retained_thio_group(mol, centre, oxo, ligands)
    if retained is not None:
        return retained, False
    replaced_oxo = mol.GetAtomWithIdx(oxo).GetAtomicNum() != 8
    if not replaced_oxo and "infix" not in roles.values():
        return None
    acid_ligands = [n for n in ligands if roles[n] in ("infix", "organyl")]
    keep = {centre, oxo}
    for n in acid_ligands:
        keep |= _side(graph, n, centre)
    if keep & _side(graph, attach, centre):
        return None
    editable = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in keep), reverse=True):
        editable.RemoveAtom(idx)
    new_centre = sum(1 for a in mol.GetAtoms() if a.GetIdx() < centre and a.GetIdx() in keep)
    for _ in range(3 - len(acid_ligands)):
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(new_centre, oxygen, Chem.BondType.SINGLE)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
        name = smiles_to_iupac(Chem.MolToSmiles(acid))
    except (UnsupportedStructure, ValueError, Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    stem = _ACID_TAIL.sub("", name)
    if stem == name or " " in stem:
        return None
    if stem.endswith("oic"):
        stem = stem[:-3] + "oyl"
    elif stem.endswith("ic"):
        stem = stem[:-2] + "oyl"
    else:
        return None
    entries = [
        name_branch(graph, n, centre, halogens, aromatic_atoms, mol=mol) for n in ligands if roles[n] == "prefix"
    ]
    return (format_mononuclear_prefixes(entries) if entries else "") + stem, bool(entries) or not stem.startswith(("phosph", "ars", "stib"))
