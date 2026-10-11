"""Divalent acyl groups that bridge two identical units in a multiplicative name (P-65.1.7.4.1-3): the group is named
from the diacid it comes from, 'butanedioyl', 'butanebis(thioyl)', 'cyclohexane-1,2-dicarbonyl'."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._multiplicative_linker import Part

_TERMINAL_X = {7, 8, 16, 34, 52}
_PLAIN_ACYL = re.compile(
    r"^oxalyl$|^(?:cyclo)?(?:meth|eth|prop|but|pent|hex|hept|oct|non|dec|undec|dodec)an(?:e|-[\d,]+-)?(?:di|tri|tetra)?"
    r"(?:oyl|carbonyl|carbothioyl|carboselenoyl|carbotelluroyl|carboximidoyl|bis\([a-z]+\))$"
)


def _terminals(mol, atom):
    return [n.GetIdx() for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetDegree() == 1]


def _acyl_neighbours(mol, carbon):
    """The two single-bonded neighbours of an acyl carbon C(=X), else None."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
        return None
    terminals, others = [], []
    for n in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble()
        if bond == 2.0 and n.GetDegree() == 1 and n.GetAtomicNum() in _TERMINAL_X and not n.GetFormalCharge():
            terminals.append(n.GetIdx())
        elif bond == 1.0:
            others.append(n.GetIdx())
    return others if len(terminals) == 1 and len(others) == 2 else None


def diacyl_linker(mol, graph, central_atoms, edges, first_atoms, protected):
    """(acyl name, inner arm count, whether it has substituents) of the divalent acyl group around the central component `central_atoms`, whose
    outward bonds are `edges` [(inside, outside, order)] and whose first outward components are `first_atoms`
    [atom set per side]; None when the centre is not a diacyl group. `protected` holds the atoms of the units."""
    from ._acid_derivatives import acyl_name
    from ._acid_groups import acid_group_at
    from .core import smiles_to_iupac

    if len(edges) != 2 or any(order != 1 for _, _, order in edges):
        return None
    ends, drops = [], set()
    for (inside, outside, _), atoms in zip(edges, first_atoms):
        if inside in central_atoms and outside in (_acyl_neighbours(mol, inside) or ()):
            ends.append((inside, outside))
            drops.add(0)
        elif atoms <= {outside, *_terminals(mol, outside)} and inside in (others := _acyl_neighbours(mol, outside) or ()):
            ends.append((outside, next(o for o in others if o != inside)))
            drops.add(1)
        else:
            return None
    if len(drops) != 1:
        return None
    if len({mol.GetAtomWithIdx(t).GetAtomicNum() for c, _ in ends for t in _terminals(mol, c) if mol.GetBondBetweenAtoms(c, t).GetBondTypeAsDouble() == 2.0}) != 1:
        return None
    cuts = {frozenset(pair) for pair in ends}
    region = set(central_atoms) | {c for c, _ in ends}
    stack = list(region)
    while stack:
        a = stack.pop()
        for n in graph[a]:
            if n not in region and frozenset((a, n)) not in cuts:
                region.add(n)
                stack.append(n)
    if region & protected or any(o in region for _, o in ends):
        return None
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in (mol.GetAtomWithIdx(i) for i in region)):
        return None
    editable = Chem.RWMol(mol)
    for carbon, _ in ends:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - region, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
    except Exception:
        return None
    if sum(1 for a in acid.GetAtoms() if acid_group_at(acid, a.GetIdx()) is not None) != 2:
        return None
    try:
        name = acyl_name(smiles_to_iupac(Chem.MolToSmiles(acid)))
    except UnsupportedStructure:
        return None
    if " " in name or re.search(r"acid|ate\b", name):
        return None
    return name, drops.pop(), _PLAIN_ACYL.match(name) is None


def acyl_part(name, substituted):
    return Part(name, substituted or bool(re.match(r"[\d(\[{]", name)), bool(re.search(r"[\d(]", name)))
