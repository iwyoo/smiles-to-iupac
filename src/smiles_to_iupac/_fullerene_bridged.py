"""Bridges between atoms of a fullerene, a modified fullerene or the rings fused to it, cited as bridge prefixes with the
locants of the bridged atoms in front of the parent (P-27.6.2). The bridge prefixes and their order are those of P-25.4,
with the fullerene and the fused rings as the fused parent whose numbering is given."""

from itertools import product

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_bridged import PHANE_NODES, _assemble, _decompositions, _kekule, _metrics

_MAX_COMBINATIONS = 200


class Bridges:
    def __init__(self, mol, graph, fused, clusters):
        aromatic = Chem.Mol(mol)
        plain = _kekule(aromatic)
        singles = [(b, f) for cluster in clusters for b in cluster for f in graph[b] if f in fused]
        self.mol = _kekule(aromatic, singles) or plain
        self.fused = frozenset(fused)
        self.atoms = frozenset(graph)
        cache = {"aromatic_rings": [r for r in aromatic.GetRingInfo().AtomRings() if all(aromatic.GetAtomWithIdx(a).GetIsAromatic() for a in r)]}
        rings = [set(r) for r in aromatic.GetRingInfo().AtomRings() if set(r) <= self.fused]
        token = PHANE_NODES.set(10**6)
        try:
            readings = [_decompositions(self.mol, graph, self.fused, set(cluster), cache, rings) for cluster in clusters]
        finally:
            PHANE_NODES.reset(token)
        combinations = [sum(combo, []) for combo in list(product(*readings))[:_MAX_COMBINATIONS]]
        best = min(_metrics(parts) for parts in combinations)
        self.parts = [parts for parts in combinations if _metrics(parts) == best]

    def assembled(self, position):
        """The bridge prefixes for the atoms of the fused parent labelled by `position`, as (key, text) of every reading,
        the senior first."""
        found = []
        for parts in self.parts:
            for parent in _assemble(self.mol, self.atoms, self.fused, sorted(self.fused), "", position, parts):
                found.append((parent.key, parent.stem))
        if not found:
            raise UnsupportedStructure("no bridge name")
        return sorted(found)
