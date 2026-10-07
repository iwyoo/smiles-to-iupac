"""Thiourea, selenourea and tellurourea (H2N-C(=E)-NH2, E = S/Se/Te) with any number of N-substituents
(P-66.1.4.1); unlike `_urea.py` no semicarbazide (amino-substituted nitrogen) handling. The N/N' citation helpers
are shared with `_urea.py`."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
)
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch


def _n_substituent_roots(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx)


def _substituent_atoms(graph, roots, core_atoms):
    """Every atom reachable from `roots` without passing through the core, or
    None when a substituent loops back into it (a ring-fused urea)."""
    reached = set()
    stack = list(roots)
    while stack:
        idx = stack.pop()
        if idx in reached:
            continue
        reached.add(idx)
        stack.extend(n for n in graph[idx] if n not in core_atoms)
    for idx in reached:
        if sum(1 for n in graph[idx] if n in core_atoms) > (1 if idx in roots else 0):
            return None
    return reached


def n_substituent_names(mol, graph, core_atoms, n1_idx, n2_idx, carbon_idx):
    """([names of n1's substituents], [names of n2's]) for every atom outside `core_atoms`
    (P-66.1.4.1); raises when a ring fuses into the core or a non-halogen
    acyclic heteroatom is present."""
    n1_roots = tuple(r for r in _n_substituent_roots(mol, n1_idx, carbon_idx) if r not in core_atoms)
    n2_roots = tuple(r for r in _n_substituent_roots(mol, n2_idx, carbon_idx) if r not in core_atoms)
    outside = _substituent_atoms(graph, n1_roots + n2_roots, core_atoms)
    if outside is None or len(outside) + len(core_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a ring-fused urea (e.g. hydantoin) or a characteristic group outside the urea core "
            "and its N-substituents is not supported yet"
        )
    halogens = halogen_substituents(mol)
    for idx in outside:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 and idx not in halogens and not atom.IsInRing():
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the urea core and its N-substituents "
                "is not supported yet"
            )
    return (
        [name_branch(graph, c, n1_idx, halogens, mol=mol) for c in n1_roots],
        [name_branch(graph, c, n2_idx, halogens, mol=mol) for c in n2_roots],
    )


def n_prefix(n1_names, n2_names):
    """The 'N'/'N'' substituent prefix block: the nitrogen with more
    substituents takes the unprimed locant, then the one whose substituent
    comes first alphanumerically (P-14.3.5, P-14.5.2)."""
    if len(n1_names) != len(n2_names):
        unprimed, primed = (n1_names, n2_names) if len(n1_names) > len(n2_names) else (n2_names, n1_names)
    else:
        key = lambda names: min((alpha_sort_key(name) for name, _ in names), default="")
        unprimed, primed = sorted((n1_names, n2_names), key=key)
    positions = {"N": unprimed, "N'": primed}
    return format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))


@dataclass(frozen=True)
class Chalcogenourea:
    atomic_num: int
    symbol: str
    word: str

    def _core(self, mol):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
                continue
            if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
                continue
            neighbors = atom.GetNeighbors()
            chalcogens = [n for n in neighbors if n.GetAtomicNum() == self.atomic_num]
            nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
            if len(chalcogens) != 1 or len(nitrogens) != 2:
                continue
            (chalcogen,) = chalcogens
            if (
                chalcogen.GetDegree() != 1
                or mol.GetBondBetweenAtoms(atom.GetIdx(), chalcogen.GetIdx()).GetBondTypeAsDouble() != 2.0
            ):
                continue
            if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
                continue
            if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
                continue
            if any(
                nn.GetAtomicNum() != 6 for n in nitrogens for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx()
            ):
                continue
            return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
        return None

    def has_shape(self, mol) -> bool:
        return self._core(mol) is not None

    def name(self, mol) -> str:
        word = self.word
        core = self._core(mol)
        if core is None:
            raise UnsupportedStructure(
                f"no {word} (H2N-C(={self.symbol})-NH2 or an N-substituted derivative) "
                f"shape found; this module only handles {word} and simple "
                f"N-substituted {word}s"
            )
        carbon_idx, (n1_idx, n2_idx) = core

        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        (chalcogen_idx,) = (
            n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == self.atomic_num
        )
        graph = adjacency(mol)
        n1_names, n2_names = n_substituent_names(
            mol, graph, {carbon_idx, chalcogen_idx, n1_idx, n2_idx}, n1_idx, n2_idx, carbon_idx
        )
        return f"{n_prefix(n1_names, n2_names)}{word}"
