"""Polyesters of Group 13 or 14 hydrides (pseudoesters, P-65.6.3.4.1, P-68.1.5.2.3, P-68.2.6.2): the atom, or the chain of
such atoms with the heteroatoms between them, is the multivalent group cited before the multiplied acid anions,
'alumanetriyl tri(octadecanoate)', '1,1,3,3-tetramethyldistannoxane-1,3-diyl dipropanoate'."""

from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diester_anions import acid_anions, cite_anions
from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
from ._substituents import format_mononuclear_prefixes, name_branch

_HOSTS = (13, 14, 31, 32, 49, 50, 81, 82)
_BRIDGES = (8, 16, 34, 52)
_VALENCE_WORD = {2: "diyl", 3: "triyl", 4: "tetrayl"}


def _host_matches(mol):
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        carbonyls = [
            o
            for o in atom.GetNeighbors()
            if o.GetAtomicNum() == 8
            and o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        for o in atom.GetNeighbors():
            hosts = [n for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() in _HOSTS]
            if (
                len(carbonyls) == 1
                and o.GetAtomicNum() == 8
                and o.GetDegree() == 2
                and not o.GetFormalCharge()
                and len(hosts) == 1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            ):
                matches.append((atom, carbonyls[0], o, hosts[0]))
    return matches


def _host(mol):
    matches = _host_matches(mol)
    hosts = {m[3].GetIdx() for m in matches}
    if len(matches) < 2 or len(hosts) != 1:
        return None
    host = mol.GetAtomWithIdx(next(iter(hosts)))
    if host.IsInRing() or host.GetFormalCharge() or host.GetIsotope() or len(Chem.GetMolFrags(mol)) > 1:
        return None
    return host, matches


def _is_bridge(atom):
    return (
        atom.GetAtomicNum() in _BRIDGES
        and atom.GetDegree() == 2
        and not atom.IsInRing()
        and not atom.GetFormalCharge()
        and all(n.GetAtomicNum() in _HOSTS for n in atom.GetNeighbors())
    )


def _chain_matches(mol):
    """The acyloxy groups of several hosts joined into one acyclic chain, possibly through chalcogen bridges, or None."""
    matches = _host_matches(mol)
    hosts = {m[3].GetIdx() for m in matches}
    if len(hosts) < 2 or len(Chem.GetMolFrags(mol)) > 1:
        return None
    if any(a.IsInRing() or a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms() if a.GetIdx() in hosts):
        return None
    reached, stack = set(), [min(hosts)]
    while stack:
        node = stack.pop()
        if node in reached:
            continue
        reached.add(node)
        stack.extend(
            n.GetIdx()
            for n in mol.GetAtomWithIdx(node).GetNeighbors()
            if n.GetIdx() not in reached and (n.GetAtomicNum() in _HOSTS or _is_bridge(n))
        )
    return matches if hosts <= reached else None


def _chain_group(mol, matches):
    """The name of the multivalent group left by removing the acyloxy groups of a chain, or None."""
    editable = Chem.RWMol(mol)
    for _, _, ester, host in matches:
        atom = editable.GetAtomWithIdx(host.GetIdx())
        hydrogens = atom.GetTotalNumHs()
        editable.RemoveBond(host.GetIdx(), ester.GetIdx())
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(hydrogens)
        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + 1)
    cut = editable.GetMol()
    fragments = Chem.GetMolFrags(cut)
    if len(fragments) != len(matches) + 1:
        return None
    hosts = {m[3].GetIdx() for m in matches}
    skeleton = next(f for f in fragments if hosts & set(f))
    return _group_name(Chem.MolFragmentToSmiles(cut, atomsToUse=list(skeleton)), len(matches))


@lru_cache(maxsize=256)
def _group_name(smiles, valence):
    from .core import smiles_to_iupac

    try:
        name = smiles_to_iupac(smiles)
    except (UnsupportedStructure, ValueError):
        return None
    return name if name.endswith(_VALENCE_WORD.get(valence, "?")) else None


def has_chain_polyester_shape(mol) -> bool:
    chain = _chain_matches(mol)
    return chain is not None and _chain_group(mol, chain) is not None


def has_hydride_polyester_shape(mol) -> bool:
    chain = _chain_matches(mol)
    if chain is not None:
        return _chain_group(mol, chain) is not None
    found = _host(mol)
    if found is None:
        return False
    host, matches = found
    if any(a.GetAtomicNum() in _HOSTS and a.GetIdx() != host.GetIdx() for a in mol.GetAtoms()):
        return False
    return len(matches) <= MONONUCLEAR_HYDRIDES[host.GetAtomicNum()][2]


def name_hydride_polyester(mol) -> str:
    chain = _chain_matches(mol)
    if chain is not None:
        group = _chain_group(mol, chain)
        if group is None:
            raise UnsupportedStructure("this chain of Group 13 and 14 atoms is not named as a multivalent group")
        return f"{group} {cite_anions(acid_anions(mol, chain), [1] * len(chain), False)}"
    found = _host(mol)
    if found is None:
        raise UnsupportedStructure("a single hydride atom carrying several acyloxy groups is required")
    host, matches = found
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    acyloxy = {m[2].GetIdx() for m in matches}
    entries = [
        name_branch(graph, n, host.GetIdx(), halogens, mol=mol, unsaturated=True)
        for n in graph[host.GetIdx()]
        if n not in acyloxy
    ]
    stem = MONONUCLEAR_HYDRIDES[host.GetAtomicNum()][0]
    group = (format_mononuclear_prefixes(entries) if entries else "") + stem + _VALENCE_WORD[len(matches)]
    anions = cite_anions(acid_anions(mol, matches), [1] * len(matches), False)
    return f"{group} {anions}"
