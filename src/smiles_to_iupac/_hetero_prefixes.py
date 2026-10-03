"""Substituent prefixes joined through a heteroatom or a functional carbon
(P-29.3.3, P-29.4, P-35, P-63.2, P-66): hydroxy, oxo, alkoxy/aryloxy, sulfanyl,
amino, nitro, cyano, formyl, carboxy, carbamoyl, alkoxycarbonyl and acyl.
`name_branch` calls `hetero_branch_name` before it walks a carbon chain.
"""

from ._common import UnsupportedStructure, alpha_sort_key
from ._numerals import alkane_name

_ALKOXY_STEMS = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy", "phenyl": "phenoxy"}
_SIMPLE_NAMES = {
    "methoxy", "ethoxy", "propoxy", "butoxy", "tert-butoxy", "phenoxy", "amino", "anilino", "hydroxy", "oxo",
    "nitro", "nitroso", "cyano", "sulfanyl", "formyl", "carboxy", "carbamoyl",
}
_MULTIPLE_TARGETS = {7, 8, 16}


def _carbonyl_oxygen(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return next(
        (
            b.GetOtherAtom(atom).GetIdx()
            for b in atom.GetBonds()
            if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8
        ),
        None,
    )


def is_functional_carbon(mol, idx):
    """A carbon that is a substituent group of its own (-C#N, -COOH, -COOR,
    -CONR2, ...) rather than a chain member. Ketone and aldehyde carbons stay
    in the chain and are cited as 'oxo' (P-64.3)."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    for b in atom.GetBonds():
        other = b.GetOtherAtom(atom)
        if b.GetBondTypeAsDouble() == 3.0 and other.GetAtomicNum() == 7:
            return True
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() in _MULTIPLE_TARGETS:
            return any(
                n.GetIdx() != other.GetIdx() and n.GetAtomicNum() in _MULTIPLE_TARGETS for n in atom.GetNeighbors()
            )
    return False


def _enclose(name, compound):
    return f"({name})" if compound else name


def _alkoxy(rname):
    if rname == "tert-butyl":
        return "tert-butoxy"
    for stem, short in _ALKOXY_STEMS.items():
        if rname.endswith(stem) and not rname.endswith("cyclo" + stem):
            return rname[: -len(stem)] + short
    return rname + "oxy"


def _compound(name):
    return name not in _SIMPLE_NAMES


def _acyl_name(mol, graph, carbon, from_atom):
    from ._substituents import name_branch

    others = [n for n in graph[carbon] if n != from_atom and mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() == 1.0]
    (alkyl,) = others
    if mol.GetAtomWithIdx(alkyl).GetIsAromatic():
        name, _ = name_branch(graph, alkyl, carbon, mol=mol)
        if name == "phenyl":
            return "benzoyl"
        raise UnsupportedStructure("this aroyl group is not supported yet")
    count = 1
    node, parent = alkyl, carbon
    while True:
        atom = mol.GetAtomWithIdx(node)
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or is_functional_carbon(mol, node):
            raise UnsupportedStructure("this acyl group is not supported yet")
        onward = [n for n in graph[node] if n != parent]
        if any(mol.GetBondBetweenAtoms(node, n).GetBondTypeAsDouble() != 1.0 for n in onward):
            raise UnsupportedStructure("an unsaturated acyl group is not supported yet")
        count += 1
        if not onward:
            break
        if len(onward) > 1:
            raise UnsupportedStructure("a branched acyl group is not supported yet")
        node, parent = onward[0], node
    return alkane_name(count)[:-1] + "oyl"


def _group_names(graph, mol, atoms, parent, halogens, aromatic_atoms):
    from ._substituents import name_branch

    out = []
    for a in atoms:
        if is_functional_carbon(mol, a) and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and any(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() == 2.0
            for n in graph[a]
        ):
            out.append((_acyl_name(mol, graph, a, parent), True))
        else:
            out.append(name_branch(graph, a, parent, halogens, aromatic_atoms, mol=mol))
    return out


def _amino(names):
    if not names:
        return "amino"
    if names == [("phenyl", False)]:
        return "anilino"
    ordered = sorted(names, key=lambda item: alpha_sort_key(item[0]))
    if len(ordered) == 2 and ordered[0] == ordered[1]:
        name, compound = ordered[0]
        return (f"bis({name})" if compound else "di" + name) + "amino"
    parts = [_enclose(ordered[0][0], ordered[0][1])]
    parts += [f"({n})" for n, _ in ordered[1:]]
    return "".join(parts) + "amino"


def hetero_branch_name(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """(name, is_compound) of a heteroatom- or functional-carbon-rooted
    substituent, or None when `root` is an ordinary carbon."""
    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    if z == 6:
        if is_functional_carbon(mol, root) or _carbonyl_oxygen(mol, root) is not None:
            return _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol)
        return None
    if atom.GetFormalCharge() and z != 7:
        raise UnsupportedStructure("a charged atom in a substituent group is not supported yet")
    if atom.IsInRing():
        return None
    others = [n for n in graph[root] if n != coming_from]
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    if z == 8:
        if order == 2.0:
            return "oxo", False
        if not others:
            return "hydroxy", False
        (other,) = others
        if mol.GetAtomWithIdx(other).GetAtomicNum() != 6:
            raise UnsupportedStructure("this oxygen-linked group is not supported yet")
        if is_functional_carbon(mol, other):
            (acyl,) = _group_names(graph, mol, [other], root, halogens, aromatic_atoms)
            return acyl[0] + "oxy", True
        from ._substituents import name_branch

        rname, _ = name_branch(graph, other, root, halogens, aromatic_atoms, mol=mol)
        name = _alkoxy(rname)
        return name, _compound(name)
    if z == 16:
        if order != 1.0 or len(others) > 1 or atom.GetDegree() > 2:
            raise UnsupportedStructure("this sulfur-linked group is not supported yet")
        if not others:
            return "sulfanyl", False
        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() != 6:
            raise UnsupportedStructure("a chalcogen chain (disulfanyl, ...) is not supported yet")
        from ._substituents import name_branch

        rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
        return _enclose(rname, rcomp) + "sulfanyl", True
    if z == 7:
        oxygens = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
        if len(oxygens) == len(others) and others:
            if len(oxygens) == 2 and atom.GetTotalDegree() == 3:
                return "nitro", False
            if len(oxygens) == 1 and mol.GetBondBetweenAtoms(root, oxygens[0]).GetBondTypeAsDouble() == 2.0:
                return "nitroso", False
        if order != 1.0 or atom.GetFormalCharge() or any(
            mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others
        ):
            raise UnsupportedStructure("this nitrogen-linked group is not supported yet")
        if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in others):
            raise UnsupportedStructure("this nitrogen-linked group is not supported yet")
        name = _amino(_group_names(graph, mol, others, root, halogens, aromatic_atoms))
        return name, _compound(name)
    raise UnsupportedStructure("this heteroatom-linked substituent is not supported yet")


def _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol):
    from ._substituents import name_branch

    atom = mol.GetAtomWithIdx(root)
    others = [n for n in graph[root] if n != coming_from]
    triple_n = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 3.0]
    if triple_n and len(others) == 1:
        return "cyano", False
    carbonyl = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0]
    if len(carbonyl) != 1 or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")
    rest = [n for n in others if n != carbonyl[0]]
    if not rest:
        return "formyl", False
    (x,) = rest
    z = mol.GetAtomWithIdx(x).GetAtomicNum()
    if z == 8:
        tail = [n for n in graph[x] if n != root]
        if not tail:
            return "carboxy", False
        rname, _ = name_branch(graph, tail[0], x, halogens, aromatic_atoms, mol=mol)
        return _alkoxy(rname) + "carbonyl", True
    if z == 7:
        subs = [n for n in graph[x] if n != root]
        if not subs:
            return "carbamoyl", False
        name = _amino(_group_names(graph, mol, subs, x, halogens, aromatic_atoms))
        return name[: -len("amino")] + "carbamoyl", True
    if z == 6:
        return _acyl_name(mol, graph, root, coming_from), True
    raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")
