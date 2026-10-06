"""Selenides and tellurides R-E-R' (E = Se/Te, P-63.2.1) named as 'R-selanyl'/'R-tellanyl' substituted parents: the
larger carbon component (or a lone benzene ring) is the parent, the other side the chalcogenyl prefix. Saturated
acyclic chains and one plain benzene ring only; no stereocenters on the ring path."""

from dataclasses import dataclass

from rdkit import Chem

from ._acyclic import longest_chain_length, name_from_carbon_graph, winning_chain_with_key
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    component_subgraph,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._multiplicative_text import enclose
from ._substituents import name_branch


@dataclass(frozen=True)
class Chalcogenide:
    atomic_num: int
    symbol: str
    element: str
    word: str
    prefix: str

    def _prefix(self, name: str) -> str:
        return name + self.prefix

    def has_shape(self, mol) -> bool:
        atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == self.atomic_num]
        if len(atoms) != 1:
            return False
        (atom,) = atoms
        return atom.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())

    def _validate_atoms(self, mol, ring_atoms=frozenset()):
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if atom.GetAtomicNum() not in (6, self.atomic_num):
                raise UnsupportedStructure(
                    f"heteroatoms other than the {self.word} {self.element} are not "
                    f"supported yet (P-63.2.1 is restricted to a plain -{self.symbol}- "
                    f"{self.word} here)"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and idx not in ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                )
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in non_single_bonds(mol)):
            raise UnsupportedStructure(
                "unsaturation is not supported by this module (see P-31 for "
                f"alkenes/alkynes; not yet combined with a {self.word} here)"
            )

    def _name_benzene_ring_chain(self, mol, ring_atoms) -> str:
        self._validate_atoms(mol, ring_atoms)
        if specified_stereocenters(mol) is not None:
            raise UnsupportedStructure(
                f"a specified stereocenter is not supported yet for a benzene-ring-substituent {self.word}"
            )

        graph = adjacency(mol)
        (hetero_idx,) = (idx for idx in graph if mol.GetAtomWithIdx(idx).GetAtomicNum() == self.atomic_num)

        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            raise UnsupportedStructure(
                "a benzene ring with more than one exocyclic substituent "
                f"alongside a {self.word} is not supported yet"
            )
        ring_atom, chain_root = attachment

        if chain_root == hetero_idx:
            (r_prime,) = [n for n in graph[hetero_idx] if n != ring_atom]
            sub_name, sub_compound = name_branch(graph, r_prime, hetero_idx, {}, mol=mol)
            if sub_compound:
                sub_name = enclose(sub_name)
            return f"{self._prefix(sub_name)}benzene"

        blocked_graph = {node: [n for n in neighbors if n != hetero_idx] for node, neighbors in graph.items()}
        del blocked_graph[hetero_idx]
        reached, _ = bfs(blocked_graph, ring_atom)
        (r_prime,) = [n for n in graph[hetero_idx] if n not in reached]

        sub_name, sub_compound = name_branch(graph, r_prime, hetero_idx, {}, mol=mol)
        if sub_compound:
            sub_name = enclose(sub_name)
        prefix_term = self._prefix(sub_name)
        branch_name, is_compound = name_branch(graph, chain_root, ring_atom, {hetero_idx: prefix_term}, mol=mol)
        if not is_compound:
            return f"{branch_name}benzene"
        if "(" in branch_name:
            return f"[{branch_name}]benzene"
        return f"({branch_name})benzene"

    def name(self, mol) -> str:
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure(
                "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
            )
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_benzene_ring_chain(mol, ring_atoms)
            raise UnsupportedStructure("rings are not supported by this module yet")
        if ring_info.NumRings() > 0:
            raise UnsupportedStructure("rings are not supported by this module yet")

        self._validate_atoms(mol)

        (hetero,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == self.atomic_num)
        hetero_idx = hetero.GetIdx()
        n1, n2 = (n.GetIdx() for n in hetero.GetNeighbors())

        full_graph = adjacency(mol)
        carbon_graph = carbon_adjacency(mol)
        size1 = len(component_subgraph(carbon_graph, n1))
        size2 = len(component_subgraph(carbon_graph, n2))

        if size1 == size2:
            graph_a = component_subgraph(carbon_graph, n1)
            graph_b = component_subgraph(carbon_graph, n2)
            len_a, len_b = longest_chain_length(graph_a), longest_chain_length(graph_b)
            if len_a > len_b:
                parent_root, sub_root = n1, n2
            elif len_b > len_a:
                parent_root, sub_root = n2, n1
            else:
                name_a, compound_a = name_branch(full_graph, n1, hetero_idx, {}, mol=mol)
                name_b, compound_b = name_branch(full_graph, n2, hetero_idx, {}, mol=mol)
                sub_from_a = enclose(name_a) if compound_a else name_a
                sub_from_b = enclose(name_b) if compound_b else name_b
                key_a, _, _ = winning_chain_with_key(
                    full_graph, graph_a, {hetero_idx: self._prefix(sub_from_b)}, mol=mol
                )
                key_b, _, _ = winning_chain_with_key(
                    full_graph, graph_b, {hetero_idx: self._prefix(sub_from_a)}, mol=mol
                )
                parent_root, sub_root = (n1, n2) if key_a <= key_b else (n2, n1)
        elif size1 > size2:
            parent_root, sub_root = n1, n2
        else:
            parent_root, sub_root = n2, n1

        sub_name, sub_compound = name_branch(full_graph, sub_root, hetero_idx, {}, mol=mol)
        if sub_compound:
            sub_name = enclose(sub_name)

        parent_carbon_graph = component_subgraph(carbon_graph, parent_root)
        terminals = {hetero_idx: self._prefix(sub_name)}
        return name_from_carbon_graph(full_graph, parent_carbon_graph, terminals, mol=mol)
