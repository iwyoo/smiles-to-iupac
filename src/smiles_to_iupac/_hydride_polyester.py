"""Polyesters of one Group 13 or 14 hydride atom (pseudoesters, P-65.6.3.4.1, P-68.1.5.2.3): the atom with its organyl
groups is the multivalent group cited before the multiplied acid anions, 'alumanetriyl tri(octadecanoate)'."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diester_anions import acid_anions, cite_anions
from ._hetero_prefixes import MONONUCLEAR_HYDRIDES
from ._substituents import format_mononuclear_prefixes, name_branch

_HOSTS = (13, 14, 31, 32, 49, 50, 81, 82)
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


def has_hydride_polyester_shape(mol) -> bool:
    found = _host(mol)
    if found is None:
        return False
    host, matches = found
    if any(a.GetAtomicNum() in _HOSTS and a.GetIdx() != host.GetIdx() for a in mol.GetAtoms()):
        return False
    return len(matches) <= MONONUCLEAR_HYDRIDES[host.GetAtomicNum()][2]


def name_hydride_polyester(mol) -> str:
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
