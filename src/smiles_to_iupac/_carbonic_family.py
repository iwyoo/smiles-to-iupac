"""Carbonic, carbamic and cyanic acids and their functional-replacement analogues (P-65.2.1, P-65.2.2): an =X position
and two -Y-H positions, or one amino ('carbam'), halide or pseudohalide position; replaced oxygen is cited by infixes
and the hydrogen-bearing atoms by italic letters (P-65.2.1.2). Two non-hydroxy positions make an acyl halide (P-65.5.3).
"""

import re

from rdkit import Chem

from ._acid_lexicon import _INFIX, _MULTIPLIER, _peroxo_word
from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._retained_acids import single_site_prefixes
from ._substituents import format_substituent_prefixes, name_branch

_SYMBOL = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_HALIDE_INFIX = {9: "fluorid", 17: "chlorid", 35: "bromid", 53: "iodid"}
_X_INFIX = {"S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido", "NNH2": "hydrazono"}
PSEUDO_INFIX = {
    "N3": "azid",
    "CN": "cyanid",
    "NC": "isocyanid",
    "NCO": "isocyanatid",
    "NCS": "isothiocyanatid",
    "NCSe": "isoselenocyanatid",
    "NCTe": "isotellurocyanatid",
}
_CUMULENE_END = {8: "NCO", 16: "NCS", 34: "NCSe", 52: "NCTe"}


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def pseudohalide_at(mol, root, parent):
    """'N3', 'CN', 'NC', 'NCO', 'NCS', 'NCSe', 'NCTe' for the pseudohalogen group rooted at `root`, else None."""
    atom = mol.GetAtomWithIdx(root)
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() != parent]
    z = atom.GetAtomicNum()
    if len(rest) != 1:
        return None
    other = rest[0]
    order = _bond(mol, root, other.GetIdx())
    if z == 6 and other.GetAtomicNum() == 7 and order == 3.0 and other.GetDegree() == 1 and not atom.GetFormalCharge():
        return "CN"
    if z != 7:
        return None
    if other.GetAtomicNum() == 6 and order == 3.0 and other.GetFormalCharge() == -1 and other.GetDegree() == 1 and atom.GetFormalCharge() == 1:
        return "NC"
    if other.GetAtomicNum() == 7 and not atom.GetFormalCharge() and other.GetFormalCharge() == 1:
        end = [n for n in other.GetNeighbors() if n.GetIdx() != root]
        if len(end) == 1 and end[0].GetAtomicNum() == 7 and end[0].GetDegree() == 1 and end[0].GetFormalCharge() == -1:
            return "N3"
    if other.GetAtomicNum() == 6 and order == 2.0 and not atom.GetFormalCharge():
        end = [n for n in other.GetNeighbors() if n.GetIdx() != root]
        if len(end) == 1 and end[0].GetDegree() == 1 and _bond(mol, other.GetIdx(), end[0].GetIdx()) == 2.0:
            return _CUMULENE_END.get(end[0].GetAtomicNum())
    return None


def _chain(mol, center, atom, ion=False):
    """(symbols from the centre outward, atoms) of a -Y-H or -Y-Y-H position (an anionic -Y(-) or -Y-Y(-)
    position when `ion`), or None."""
    z = atom.GetAtomicNum()
    if z not in _SYMBOL or _bond(mol, center, atom.GetIdx()) != 1.0:
        return None
    charge, hydrogens = (-1, 0) if ion else (0, 1)
    if atom.GetFormalCharge() == charge and atom.GetDegree() == 1 and atom.GetTotalNumHs() == hydrogens:
        return [_SYMBOL[z]], {atom.GetIdx()}
    onward = [n for n in atom.GetNeighbors() if n.GetIdx() != center]
    if (
        not atom.GetFormalCharge()
        and len(onward) == 1
        and onward[0].GetAtomicNum() in _SYMBOL
        and onward[0].GetDegree() == 1
        and onward[0].GetFormalCharge() == charge
        and onward[0].GetTotalNumHs() == hydrogens
        and not atom.GetTotalNumHs()
    ):
        return [_SYMBOL[z], _SYMBOL[onward[0].GetAtomicNum()]], {atom.GetIdx(), onward[0].GetIdx()}
    return None


def _group_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _infix_text(terms):
    """Infixes from (word, count, compound) triples in alphabetical order; the last word carries 'ic'."""
    terms = sorted(terms)
    out = []
    for i, (word, count, compound) in enumerate(terms):
        last = i == len(terms) - 1
        if compound and count > 1:
            text = f"{'bis' if count == 2 else 'tris'}({word}{'ic' if last else ''})"
        else:
            text = _MULTIPLIER[count] + word
            if last:
                text = text[:-1] + "ic" if word in ("imido", "hydrazono") else text + "ic"
            if compound:
                text = f"({text})"
        out.append(text)
    return "".join(out)


def _letters(x, chains):
    """Italic letters of the atoms bearing hydrogen when several chalcogens are present: the chains that contain
    a replaced atom, or all of them when only the =X atom is replaced (P-65.2.1.2)."""
    elements = {a for chain in chains for a in chain} | ({x} if x in _SYMBOL.values() else set())
    if len(elements) < 2:
        return ""
    marked = ["".join(c) for c in chains if any(a != "O" for a in c)] or ["".join(c) for c in chains]
    return ",".join(sorted(marked))


def carbonic_acid_name(x, ligands, amino=False):
    """Name of a carbonic acid analogue: `x` the =X atom ('O', 'S', ..., 'NH'), `ligands` the other positions as
    ('Y', symbols) or ('hal', infix) tuples; `amino` when one position is an amino group."""
    simple, compound, chains = {}, {}, []
    if x != "O":
        simple[_X_INFIX[x]] = simple.get(_X_INFIX[x], 0) + 1
    for kind, value in ligands:
        if kind == "Y":
            chains.append(value)
            if len(value) == 2:
                word = _peroxo_word(value)
                compound[word] = compound.get(word, 0) + 1
            elif value[0] != "O":
                simple[_INFIX[value[0]]] = simple.get(_INFIX[value[0]], 0) + 1
        else:
            simple[value] = simple.get(value, 0) + 1
    terms = [(w, c, False) for w, c in simple.items()] + [(w, c, w != "peroxo") for w, c in compound.items()]
    stem = "carbam" if amino else "carbon"
    if not terms:
        return stem + "ic acid"
    text = _infix_text(terms)
    letters = _letters(x, chains)
    return stem + ("" if text[0] in "aeiouy" else "o") + text + (f" {letters}-acid" if letters else " acid")


def cyanic_acid_name(chain):
    """P-65.2.2: cyanic acid and the chalcogen analogues named with prefixes."""
    if chain == ["O"]:
        return "cyanic acid"
    if len(chain) == 1:
        return _INFIX[chain[0]] + "cyanic acid"
    word = _peroxo_word(chain)[: -len("o")]
    letters = _letters("O", [chain])
    return word + "cyanic" + (f" {letters}-acid" if letters else " acid")


_ANION_ENDING = re.compile(r"ic(?: [A-Za-z,]+-acid| acid)?(\)?)$")


def _anionic(name):
    return _ANION_ENDING.sub(lambda m: "ate" + (m.group(1) or ""), name)


def carbonic_anion_center(mol, idx):
    """Whether atom `idx` is the centre of a carbonic or cyanic acid group (for the anion test of salts)."""
    return _find(mol) == idx


def _find(mol):
    centers, cyano = [], []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetTotalNumHs():
            continue
        neighbors = list(atom.GetNeighbors())
        triples = [n for n in neighbors if _bond(mol, atom.GetIdx(), n.GetIdx()) == 3.0 and n.GetAtomicNum() == 7 and n.GetDegree() == 1]
        doubles = [n for n in neighbors if _bond(mol, atom.GetIdx(), n.GetIdx()) == 2.0]
        if triples and len(neighbors) == 2:
            cyano.append(atom.GetIdx())
        elif len(neighbors) == 3 and len(doubles) == 1:
            centers.append(atom.GetIdx())
    lone = [
        c
        for c in cyano
        if not any(n.GetIdx() in centers for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 6)
    ]
    found = centers + lone
    return found[0] if len(found) == 1 else None


def _imine_slot(mol, center, n):
    """('X', hydrazine tail atom or None, N-substituted imine nitrogen or None) for a double-bonded neighbor."""
    i, z = n.GetIdx(), n.GetAtomicNum()
    if z in _SYMBOL and n.GetDegree() == 1 and not n.GetFormalCharge():
        return _SYMBOL[z], None, None
    if z != 7 or n.GetFormalCharge():
        raise UnsupportedStructure("unsupported =X atom on a carbonic acid")
    far = [m for m in n.GetNeighbors() if m.GetIdx() != center]
    if not far and n.GetTotalNumHs() == 1:
        return "NH", None, None
    if len(far) == 1 and far[0].GetAtomicNum() == 7 and far[0].GetDegree() == 1 and far[0].GetTotalNumHs() == 2 and not n.GetTotalNumHs():
        return "NNH2", far[0].GetIdx(), None
    if far and not n.IsInRing() and all(_bond(mol, i, m.GetIdx()) == 1.0 for m in far) and not n.GetTotalNumHs():
        return "NH", None, i
    raise UnsupportedStructure("this imine nitrogen is not supported on a carbonic acid")


def name_carbonic_family(mol):
    """Name of a molecule that is one carbonic, carbamic or cyanic acid analogue with organyl groups only on nitrogen."""
    center = _find(mol)
    if center is None or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("no single carbonic or cyanic acid group")
    graph = adjacency(mol)
    atom = mol.GetAtomWithIdx(center)
    if any(_bond(mol, center, n.GetIdx()) == 3.0 for n in atom.GetNeighbors()):
        (y,) = [n for n in atom.GetNeighbors() if _bond(mol, center, n.GetIdx()) == 1.0]
        chain = _chain(mol, center, y)
        ionic = _chain(mol, center, y, ion=True) if chain is None else None
        found = chain or ionic
        if found is None or len(found[1]) + 2 != mol.GetNumAtoms():
            raise UnsupportedStructure("this cyanic acid derivative is not an acid")
        name = cyanic_acid_name(found[0])
        return _anionic(name) if ionic is not None else name
    x, ligands, amino, imine_n = None, [], None, None
    owned = {center}
    anions = hydrogens = 0
    for n in atom.GetNeighbors():
        i, z = n.GetIdx(), n.GetAtomicNum()
        if _bond(mol, center, i) == 2.0:
            x, tail, imine_n = _imine_slot(mol, center, n)
            owned |= {i, *([tail] if tail is not None else [])}
            continue
        chain = _chain(mol, center, n)
        ionic = _chain(mol, center, n, ion=True) if chain is None else None
        pseudo = pseudohalide_at(mol, i, center)
        if ionic is not None:
            ligands.append(("Y", ionic[0]))
            owned |= ionic[1]
            anions += 1
        elif chain is not None:
            ligands.append(("Y", chain[0]))
            owned |= chain[1]
            hydrogens += 1
        elif z in _HALIDE_INFIX and n.GetDegree() == 1:
            ligands.append(("hal", _HALIDE_INFIX[z]))
            owned.add(i)
        elif pseudo is not None:
            ligands.append(("hal", PSEUDO_INFIX[pseudo]))
            owned |= _group_atoms(graph, i, center)
        elif z == 7 and not n.IsInRing() and not n.GetFormalCharge() and amino is None and all(_bond(mol, i, m.GetIdx()) == 1.0 for m in n.GetNeighbors()):
            amino = i
            owned.add(i)
        else:
            raise UnsupportedStructure("this ligand of the carbonic acid group is not supported")
    if x is None:
        raise UnsupportedStructure("no =X atom")
    if any(kind != "Y" for kind, _ in ligands) and (len([1 for kind, _ in ligands if kind != "Y"]) + (amino is not None)) > 1:
        raise UnsupportedStructure("two non-hydroxy positions make an acyl halide, not an acid")
    if amino is not None and any(mol.GetAtomWithIdx(m).GetAtomicNum() == 7 for m in graph[amino] if m != center):
        raise UnsupportedStructure("an amino group bearing another nitrogen is a hydrazine parent")
    branch = set()
    for root in (amino, imine_n):
        if root is not None:
            for m in graph[root]:
                if m != center:
                    branch |= _group_atoms(graph, m, root)
    if set(range(mol.GetNumAtoms())) - owned != branch:
        raise UnsupportedStructure("atoms outside the carbonic acid group and its nitrogen substituents")
    name = carbonic_acid_name(x, ligands, amino=amino is not None)
    if anions:
        name = _anionic(name)
        if hydrogens:
            name = f"{'hydrogen' if hydrogens == 1 else 'dihydrogen'} {name}"
    return _nitrogen_prefixes(mol, graph, center, amino, imine_n, x) + name


def _nitrogen_prefixes(mol, graph, center, amino, imine_n, x):
    halogens = halogen_substituents(mol)
    entries = {}
    for root, locant in ((amino, "N"), (imine_n, "N'" if amino is not None else "N")):
        if root is None:
            continue
        for m in graph[root]:
            if m != center:
                name, compound = name_branch(graph, m, root, halogens, mol=mol)
                entries.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    if not entries:
        return ""
    if amino is not None and imine_n is None and x == "O":
        return single_site_prefixes(entries)
    return format_substituent_prefixes(entries)
