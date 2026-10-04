"""Prefix names for heteroatom groups and acyl groups hanging off a skeleton (P-63, P-65.1.7), computed
bottom-up so `name_branch`/`_name_acyl_part` can cite them as plain leaves: hydroxy, oxo, amino,
(alkyl)amino, nitro, cyano, sulfanyl, alkoxy, alkylsulfanyl and acyl/acyloxy (acid name '-ic acid' to '-yl').
"""

import re

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, alpha_sort_key
from ._hetero_prefixes import (
    ANIONIC_PREFIXES,
    CHALCOGEN_PREFIXES,
    _functional_carbon,
    _has_senior_principal_group,
    phosphoryl_name,
    require_plain_chalcogen_kids,
    require_senior_group,
)
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._substituents import name_branch

_NATIVE_ROOTS = frozenset({6, 7, 8, 9, 16, 17, 34, 35, 52, 53})
_RETAINED_ALKYL_END = re.compile(r"(meth|eth|prop|but)yl$")


_enclose = enclose


def _is_compound(name):
    return any(ch.isdigit() for ch in name) or "(" in name


def _bond_order(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def is_nitro_nitrogen(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 7 or atom.GetDegree() != 3:
        return False
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    return len(oxygens) == 2 and atom.GetFormalCharge() in (0, 1)


def nitro_atoms(mol):
    found = set()
    for atom in mol.GetAtoms():
        if is_nitro_nitrogen(mol, atom.GetIdx()):
            found.add(atom.GetIdx())
            found.update(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)
    return found


def _acyl_prefix(mol, subtree, root, parent):
    """'-yl' name of the acyl group rooted at `root`, from its parent acid's name."""
    from .core import smiles_to_iupac

    rw = Chem.RWMol(mol)
    new_o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(root, new_o, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - subtree, reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    Chem.SanitizeMol(sub)
    acid = smiles_to_iupac(Chem.MolToSmiles(sub))
    if not acid.endswith("ic acid"):
        raise UnsupportedStructure("an acyl substituent whose parent acid has no 'ic acid' name is not supported yet")
    return acid[: -len("ic acid")] + "yl"


def _multiplied_amino(children, tail):
    parts = sorted(children, key=lambda c: alpha_sort_key(c[0]))
    names = [(_enclose(n) if c else n) for n, c in parts]
    if len(set(parts)) == 1 and len(parts) > 1:
        n, c = parts[0]
        return multiplying_prefix(len(parts), compound=c) + (_enclose(n) if c else n) + tail
    first, rest = names[0], names[1:]
    return first + "".join(_enclose(r) if not r.startswith(("(", "[")) else r for r in rest) + tail


def _phosphoryl_parts(mol, node, parent, kids, named, bond_order):
    """Substituent names of a P(=O)(X)(Y) group bonded through oxygen, from the already-named oxygen children."""
    atom = mol.GetAtomWithIdx(node)
    if atom.GetFormalCharge() != 0 or atom.GetDegree() != 4:
        raise UnsupportedStructure("a charged or non-tetracoordinate phosphorus substituent is not supported yet")
    oxo = [k for k in kids if bond_order(mol, node, k) == 2.0 and mol.GetAtomWithIdx(k).GetDegree() == 1]
    rest = [k for k in kids if k not in oxo]
    if len(oxo) != 1 or len(rest) != 2 or any(bond_order(mol, node, k) != 1.0 for k in rest):
        raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
    parts = []
    for k in rest:
        if k not in named or mol.GetAtomWithIdx(k).GetAtomicNum() != 8:
            raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
        if "phospho" in named[k][0]:
            raise UnsupportedStructure("a polyphosphate chain substituent is not supported yet")
        parts.append(named[k])
    return parts


def functional_names(mol, graph, seeds, blocked, halogens, aromatic_atoms=frozenset(), chain_seeds=False):
    """({atom: (name, is_compound)}, {atom: display string}, covered atoms) for every functional group rooted at an
    atom reachable outward from `seeds` ((parent, root) pairs) without entering `blocked`."""
    parent_of = {}
    order = []
    for parent, root in seeds:
        if root in parent_of or root in blocked:
            continue
        parent_of[root] = parent
        queue = [root]
        while queue:
            node = queue.pop(0)
            order.append(node)
            for n in graph[node]:
                if n not in parent_of and n not in blocked and n != parent_of[node]:
                    parent_of[n] = node
                    queue.append(n)

    named = {}
    shown = dict(halogens)
    children_of = {}
    for node in order:
        children_of.setdefault(parent_of[node], []).append(node)

    def subtree(node):
        found = {node}
        stack = [node]
        while stack:
            for child in children_of.get(stack.pop(), []):
                found.add(child)
                stack.append(child)
        return found

    nitro = nitro_atoms(mol)
    root_set = set() if chain_seeds else {root for _, root in seeds}

    from ._diester_ring_diyl import _system_of

    ring_info = mol.GetRingInfo()
    ring_entries = set()
    for node in order:
        if ring_info.NumAtomRings(node) and node not in nitro:
            if parent_of[node] not in _system_of(mol, node)[1]:
                ring_entries.add(node)
    skip = set()
    for node in ring_entries:
        skip |= subtree(node) - {node}
    delegated = {}
    for _, root in seeds:
        if root in parent_of and root not in skip and mol.GetAtomWithIdx(root).GetAtomicNum() not in _NATIVE_ROOTS:
            try:
                delegated[root] = name_branch(graph, root, parent_of[root], shown, aromatic_atoms, mol=mol)
            except UnsupportedStructure:
                continue
            skip |= subtree(root) - {root}

    def child_name(child, via):
        if child in named:
            return named[child]
        if mol.GetAtomWithIdx(child).GetAtomicNum() == 6:
            return name_branch(graph, child, via, shown, aromatic_atoms, mol=mol)
        raise UnsupportedStructure("this substituent group is not supported yet")

    def record(node, name, compound):
        named[node] = (name, compound)
        shown[node] = _enclose(name) if compound else name

    phosphoryl_nodes = set()
    for node in reversed(order):
        if node in skip:
            continue
        if node in ring_entries:
            name, compound = name_branch(graph, node, parent_of[node], shown, aromatic_atoms, mol=mol)
            record(node, name, compound)
            continue
        if node in delegated:
            record(node, *delegated[node])
            continue
        atom = mol.GetAtomWithIdx(node)
        z = atom.GetAtomicNum()
        parent = parent_of[node]
        kids = [n for n in graph[node] if n != parent]
        if z in HALOGEN_PREFIXES:
            continue
        if z == 15 and mol.GetAtomWithIdx(parent).GetAtomicNum() == 8:
            group = phosphoryl_name(_phosphoryl_parts(mol, node, parent, kids, named, _bond_order))
            phosphoryl_nodes.add(node)
            record(node, group, group != "phosphono")
            continue
        if node in nitro:
            if is_nitro_nitrogen(mol, node):
                record(node, "nitro", False)
            continue
        if atom.HasProp("_anion") and not kids and _bond_order(mol, node, parent) == 1.0:
            record(node, ANIONIC_PREFIXES[z], False)
            continue
        if z == 8 and node not in named:
            if _bond_order(mol, node, parent) == 2.0 and not kids:
                record(node, "oxo", False)
            elif not kids:
                record(node, "hydroxy", False)
            elif len(kids) == 1 and kids[0] in phosphoryl_nodes:
                rname, rcompound = named[kids[0]]
                record(node, "phosphonooxy" if rname == "phosphono" else _enclose(rname) + "oxy", True)
            elif len(kids) == 1 and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() in (6, 9, 17, 35, 53):
                rname, rcompound = child_name(kids[0], node)
                if rname == "phenyl":
                    name = "phenoxy"
                elif _RETAINED_ALKYL_END.search(rname):
                    name = rname[:-2] + "oxy"
                elif rname[0].isdigit() or rname[0] == "(":
                    name = _enclose(rname) + "oxy"
                else:
                    name = rname + "oxy"
                record(node, name, _is_compound(name) or name != "phenoxy" and not _RETAINED_ALKYL_END.search(rname))
            else:
                raise UnsupportedStructure("this oxygen-bearing substituent is not supported yet")
        elif z in CHALCOGEN_PREFIXES:
            word = CHALCOGEN_PREFIXES[z]
            if not kids:
                require_senior_group(mol, z)
            require_plain_chalcogen_kids(mol, z, kids)
            if not kids and _bond_order(mol, node, parent) == 2.0:
                record(node, word[:-2] + "ylidene", False)
            elif not kids:
                record(node, word, False)
            elif len(kids) == 1 and atom.GetDegree() == 2 and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() == 6:
                rname, rcompound = child_name(kids[0], node)
                name = (_enclose(rname) if rcompound else rname) + word
                record(node, name, True)
            else:
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
        elif z == 7:
            bond = _bond_order(mol, node, parent)
            if bond == 3.0 and not kids:
                continue
            if bond == 2.0 and len(kids) == 1 and named.get(kids[0]) == ("amino", False) and _bond_order(mol, node, kids[0]) == 1.0:
                record(node, "hydrazinylidene", False)
                continue
            if (
                bond == 2.0
                and len(kids) <= 1
                and not atom.GetFormalCharge()
                and all(_bond_order(mol, node, k) == 1.0 for k in kids)
                and _has_senior_principal_group(mol)
            ):
                if not kids:
                    record(node, "imino", False)
                else:
                    rname, rcompound = child_name(kids[0], node)
                    record(node, (_enclose(rname) if rcompound else rname) + "imino", True)
                continue
            if bond != 1.0 or any(_bond_order(mol, node, k) != 1.0 for k in kids):
                raise UnsupportedStructure("an imine/azo/nitroso-type substituent is not supported yet")
            if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 for k in kids):
                raise UnsupportedStructure("a hetero-substituted nitrogen substituent is not supported yet")
            if not kids:
                record(node, "amino", False)
            else:
                children = [child_name(k, node) for k in kids]
                name = _multiplied_amino(children, "amino")
                record(node, name, True)
        elif z == 6 and node not in named:
            carbonyl = [k for k in kids if mol.GetAtomWithIdx(k).GetAtomicNum() == 8 and _bond_order(mol, node, k) == 2.0]
            triple_n = [k for k in kids if mol.GetAtomWithIdx(k).GetAtomicNum() == 7 and _bond_order(mol, node, k) == 3.0]
            if triple_n and len(kids) == 1:
                record(node, "cyano", False)
            elif carbonyl and (node in root_set or mol.GetAtomWithIdx(parent).GetAtomicNum() != 6):
                others = [k for k in kids if k not in carbonyl]
                if len(carbonyl) != 1 or len(others) > 1:
                    raise UnsupportedStructure("this carbonyl substituent is not supported yet")
                if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 for k in others):
                    record(node, *_functional_carbon(graph, node, parent, shown, aromatic_atoms, mol))
                    continue
                name = "formyl" if not others else _acyl_prefix(mol, subtree(node), node, parent)
                record(node, name, _is_compound(name))
            elif carbonyl:
                others = [k for k in kids if k not in carbonyl]
                if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 for k in others):
                    record(node, *_functional_carbon(graph, node, parent, shown, aromatic_atoms, mol))
    covered = set()
    for node in named:
        covered |= subtree(node)
    return named, shown, covered
