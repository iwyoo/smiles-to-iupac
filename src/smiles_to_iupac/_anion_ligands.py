"""Anionic heteroatom-donor ligands (P-69.2.2, `tmp/bluebook/P6a.txt` lines
8619-8680): the '-ide/-ate/-ite -> -ido/-ato/-ito' rule applied to the
ligand's parent acid or hydride, e.g. 'acetato', 'benzenethiolato',
'dimethylazanido', 'nitrato', 'azido', 'thiocyanato-kappaS', and the chelates
'pentane-2,4-dionato-kappa2O,O\'', 'oxalato-kappa2O,O\'', 'carbonato-kappa2O,O\''.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._substituents import format_mononuclear_prefixes, name_branch

_K = "κ"


def _neutral_name(mol, atoms, edit=None):
    keep = sorted(atoms)
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        rw.RemoveAtom(idx)
    index = {old: new for new, old in enumerate(keep)}
    if edit:
        edit(rw, index)
    frag = rw.GetMol()
    for a in frag.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNoImplicit(False)
        a.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(frag)
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(frag))


def _acid_to_ato(name: str) -> str:
    if name.endswith("ic acid"):
        stem = name[: -len("ic acid")] + "ato"
        return {"ethanoato": "acetato", "methanoato": "formato"}.get(stem, stem)
    raise UnsupportedStructure("this carboxylate ligand is not supported yet")


def _bond(mol, a, b):
    bond = mol.GetBondBetweenAtoms(a, b)
    return None if bond is None else bond.GetBondTypeAsDouble()


def _prefixed(mol, graph, donor, metal_idx, suffix, atoms):
    entries = []
    for n in graph[donor]:
        if n == metal_idx:
            continue
        if mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or mol.GetAtomWithIdx(n).GetIsAromatic():
            raise UnsupportedStructure("only alkyl substituents are supported on this anionic donor")
        entries.append(name_branch(graph, n, donor, {}, mol=mol))
    return (format_mononuclear_prefixes(entries) if entries else "") + suffix


def monodentate_anion(mol, metal, donor, atoms):
    graph = adjacency(mol)
    z, d, m = donor.GetAtomicNum(), donor.GetIdx(), metal.GetIdx()
    own = [n.GetIdx() for n in donor.GetNeighbors() if n.GetIdx() != m]
    hydrogens = donor.GetNumExplicitHs()
    sym = lambda i: mol.GetAtomWithIdx(i).GetAtomicNum()
    if z == 8 and len(own) == 1 and hydrogens == 0:
        x = own[0]
        if sym(x) == 6:
            terminal_o = [n for n in graph[x] if n != d and sym(n) == 8 and _bond(mol, x, n) == 2.0]
            if terminal_o and len(graph[x]) == 2 and mol.GetAtomWithIdx(x).GetTotalNumHs() == 1:
                return "formato"
            if terminal_o and len(graph[x]) == 3:
                return _acid_to_ato(_neutral_name(mol, atoms))
        if sym(x) == 7:
            others = [n for n in graph[x] if n != d]
            if len(others) == 2 and all(sym(n) == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1 for n in others):
                return "nitrato"
            if len(others) == 1 and sym(others[0]) == 8 and _bond(mol, x, others[0]) == 2.0:
                return f"nitrito-{_K}O"
    if z == 7:
        chain = [n for n in own if sym(n) == 7]
        if len(own) == 1 and len(chain) == 1:
            far = [n for n in graph[chain[0]] if n != d]
            if len(far) == 1 and sym(far[0]) == 7 and len(atoms) == 3:
                return "azido"
        if len(own) == 1 and sym(own[0]) == 6 and len(atoms) == 3:
            tail = [n for n in graph[own[0]] if n != d]
            if len(tail) == 1 and sym(tail[0]) == 16:
                return f"thiocyanato-{_K}N"
        if len(own) == 1 and sym(own[0]) == 6 and len(atoms) == 2 and _bond(mol, d, own[0]) == 3.0:
            return f"cyanido-{_K}N"
        if hydrogens + sum(_bond(mol, d, n) for n in own) == 2 and donor.GetFormalCharge() in (0, -1) and own:
            return _prefixed(mol, graph, d, m, "azanido", atoms)
    if z == 16 and len(own) == 1 and hydrogens == 0:
        x = own[0]
        if sym(x) == 6:
            tail = [n for n in graph[x] if n != d]
            if len(tail) == 1 and sym(tail[0]) == 7 and _bond(mol, x, tail[0]) == 3.0 and len(atoms) == 3:
                return f"thiocyanato-{_K}S"
            thiol = _neutral_name(mol, atoms)
            if thiol.endswith("thiol"):
                return thiol + "ato"
    if z == 15 and own and hydrogens + sum(_bond(mol, d, n) for n in own) == 2:
        return _prefixed(mol, graph, d, m, "phosphanido", atoms)
    return None


def _kappa(symbols):
    seen, cited = {}, []
    for symbol in sorted(symbols):
        cited.append(symbol + "'" * seen.get(symbol, 0))
        seen[symbol] = seen.get(symbol, 0) + 1
    return f"{_K}{len(symbols)}{','.join(cited)}"


def bidentate_anion(mol, metal, donors_in, atoms):
    graph = adjacency(mol)
    m = metal.GetIdx()
    ids = [d.GetIdx() for d in donors_in]
    if len(ids) != 2 or any(d.GetAtomicNum() != 8 for d in donors_in):
        return None
    sym = lambda i: mol.GetAtomWithIdx(i).GetAtomicNum()
    kappa = _kappa(["O", "O"])
    carbons = {d: [n for n in graph[d] if n != m] for d in ids}
    if any(len(v) != 1 or sym(v[0]) != 6 for v in carbons.values()):
        return None
    c1, c2 = carbons[ids[0]][0], carbons[ids[1]][0]
    if c1 == c2:
        others = [n for n in graph[c1] if n not in ids]
        if len(others) == 1 and sym(others[0]) == 8 and mol.GetAtomWithIdx(others[0]).GetDegree() == 1:
            return f"carbonato-{kappa}"
        if len(others) == 1 and sym(others[0]) == 6 and len(atoms) >= 4:
            return f"{_acid_to_ato(_neutral_name(mol, atoms))}-{kappa}"
        return None
    if mol.GetBondBetweenAtoms(c1, c2) is not None and len(atoms) == 6:
        return f"oxalato-{kappa}"
    middle = set(graph[c1]) & set(graph[c2])
    if len(middle) == 1:
        (c3,) = middle
        if sym(c3) == 6 and mol.GetAtomWithIdx(c3).GetDegree() == 2 and mol.GetAtomWithIdx(c3).GetTotalNumHs() <= 1:
            o1_double = _bond(mol, c1, ids[0])
            o2_double = _bond(mol, c2, ids[1])

            def edit(rw, index):
                rw.GetBondBetweenAtoms(index[c1], index[ids[0]]).SetBondType(Chem.BondType.DOUBLE)
                rw.GetBondBetweenAtoms(index[c2], index[ids[1]]).SetBondType(Chem.BondType.DOUBLE)
                rw.GetBondBetweenAtoms(index[c1], index[c3]).SetBondType(Chem.BondType.SINGLE)
                rw.GetBondBetweenAtoms(index[c2], index[c3]).SetBondType(Chem.BondType.SINGLE)
                rw.GetAtomWithIdx(index[c3]).SetNumExplicitHs(2)
                rw.GetAtomWithIdx(index[c3]).SetNoImplicit(True)

            diketone = _neutral_name(mol, atoms, edit)
            if diketone.endswith("dione"):
                return f"{diketone[:-1]}ato-{kappa}"
    return None
