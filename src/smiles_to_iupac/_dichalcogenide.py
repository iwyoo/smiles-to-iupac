"""Polysulfides, polyselenides and polytellurides R-E(n)-R' (E = S/Se/Te, n >= 2) named as '(R-disulfanyl)' /
'(R-triselanyl)' ... substituted parents: the larger carbon component (or a lone benzene ring) is the parent, the other
side the multiplied chalcogenyl prefix. Only saturated acyclic chains and one plain benzene ring bonded directly to E;
no suffix form exists (P-63.2.1.2). A terminal -EH on one end is allowed ('disulfanylmethane')."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    component_subgraph,
    group_substituents,
    is_plain_benzene_ring,
    longest_chains,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._multiplicative_text import enclose
from ._numerals import alkane_name, multiplying_prefix
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain_forced_compound_terminals


@dataclass(frozen=True)
class Dichalcogenide:
    atomic_num: int
    symbol: str
    element: str
    word: str
    base: str
    stem: str
    perol: str

    def _prefix_word(self, n: int) -> str:
        return multiplying_prefix(n) + self.stem

    def _chain(self, mol):
        atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == self.atomic_num]
        if len(atoms) < 2:
            return None
        idxs = {atom.GetIdx() for atom in atoms}

        chain_graph = {}
        for atom in atoms:
            neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in idxs]
            if len(neighbors) not in (1, 2):
                return None
            if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n).GetBondTypeAsDouble() != 1.0 for n in neighbors):
                return None
            chain_graph[atom.GetIdx()] = neighbors

        termini = [idx for idx, neighbors in chain_graph.items() if len(neighbors) == 1]
        if len(termini) != 2:
            return None
        order = [termini[0]]
        previous, current = None, termini[0]
        while len(order) < len(atoms):
            next_atoms = [n for n in chain_graph[current] if n != previous]
            if not next_atoms:
                return None
            previous, current = current, next_atoms[0]
            order.append(current)
        if len(order) != len(atoms) or order[-1] != termini[1]:
            return None

        for idx in order[1:-1]:
            if mol.GetAtomWithIdx(idx).GetDegree() != 2:
                return None
        for idx in (order[0], order[-1]):
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetDegree() == 2:
                (other,) = [n for n in atom.GetNeighbors() if n.GetIdx() not in idxs]
                if other.GetAtomicNum() != 6:
                    return None
            elif atom.GetDegree() == 1:
                if atom.GetTotalNumHs() != 1:
                    return None
            else:
                return None
        return order

    def has_shape(self, mol) -> bool:
        order = self._chain(mol)
        if order is None:
            return False
        e1, e2 = mol.GetAtomWithIdx(order[0]), mol.GetAtomWithIdx(order[-1])
        return not (e1.GetDegree() == 1 and e2.GetDegree() == 1)

    def _validate_and_find(self, mol, aromatic_ring_atoms=frozenset()):
        if not self.has_shape(mol):
            raise UnsupportedStructure(
                f"no plain poly{self.base} (R-{self.symbol}(n)-R') skeleton found; this module "
                f"only handles poly{self.base}s"
            )
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in (6, self.atomic_num):
                raise UnsupportedStructure(
                    f"heteroatoms other than the {self.word}'s own two {self.element}s are not supported yet"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                )
        if any(a not in aromatic_ring_atoms or b not in aromatic_ring_atoms for a, b, _ in non_single_bonds(mol)):
            raise UnsupportedStructure(
                "unsaturation is not supported by this module (see P-31 for "
                f"alkenes/alkynes; not yet combined with a {self.word} here)"
            )
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0 and not (
            ring_info.NumRings() == 1 and set(ring_info.AtomRings()[0]) == set(aromatic_ring_atoms)
        ):
            raise UnsupportedStructure(
                "rings are not supported yet, other than the separate benzene-ring-substituent path"
            )
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        order = self._chain(mol)
        e1_idx, e2_idx = order[0], order[-1]
        inner_neighbor_1, inner_neighbor_2 = order[1], order[-2]
        e1, e2 = mol.GetAtomWithIdx(e1_idx), mol.GetAtomWithIdx(e2_idx)
        c1 = None
        if e1.GetDegree() == 2:
            (c1,) = [n.GetIdx() for n in e1.GetNeighbors() if n.GetIdx() != inner_neighbor_1]
        c2 = None
        if e2.GetDegree() == 2:
            (c2,) = [n.GetIdx() for n in e2.GetNeighbors() if n.GetIdx() != inner_neighbor_2]
        return e1_idx, e2_idx, c1, c2, len(order)

    @staticmethod
    def _name_from_substituents(chain_length, grouped, bare_name):
        if chain_length == 1 and grouped:
            (name,) = grouped
            if name == bare_name:
                return name + alkane_name(chain_length)
            return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
        total_subs = sum(len(info["locants"]) for info in grouped.values())
        if chain_length == 2 and total_subs == 1:
            (name,) = grouped
            if name == bare_name:
                display_name = name
            else:
                display_name = enclose(name) if grouped[name]["compound"] else name
            return display_name + alkane_name(chain_length)
        prefix = format_substituent_prefixes(grouped)
        return prefix + alkane_name(chain_length)

    def _candidate_key(self, chain_length, substituents, bare_name):
        grouped = group_substituents(substituents)
        locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
        name = self._name_from_substituents(chain_length, grouped, bare_name)
        return (-total_count, locant_set, citation_locants, name), name

    def _name_parent_chain(self, full_graph, carbon_graph, terminals, bare_name, mol=None):
        chains = longest_chains(carbon_graph)
        chain_length = len(chains[0])

        best_key = None
        best_name = None
        best_chain = None
        for chain in chains:
            for candidate in (chain, list(reversed(chain))):
                substituents = substituents_for_chain_forced_compound_terminals(
                    full_graph, candidate, terminals, mol=mol
                )
                key, name = self._candidate_key(chain_length, substituents, bare_name)
                if best_key is None or key < best_key:
                    best_key, best_name, best_chain = key, name, candidate
        return best_chain, best_name

    def _name_benzene_ring_chain(self, mol, ring_atoms) -> str:
        e1_idx, e2_idx, c1, c2, n = self._validate_and_find(mol, aromatic_ring_atoms=ring_atoms)
        if c1 is None or c2 is None:
            raise UnsupportedStructure(
                f"a {self.perol} (-{self.symbol}(n)-{self.symbol}H) combined with a benzene-ring-substituent "
                f"poly{self.base} is out of scope for this module"
            )
        if specified_stereocenters(mol) is not None:
            raise UnsupportedStructure(
                f"a specified stereocenter is not supported yet for a benzene-ring-substituent {self.word}"
            )

        full_graph = adjacency(mol)
        attachment = ring_chain_attachment(full_graph, ring_atoms, set())
        if attachment is None:
            raise UnsupportedStructure(
                "a benzene ring with more than one exocyclic substituent "
                f"alongside a {self.word} chain is not supported yet"
            )
        ring_atom, root = attachment

        if root == e1_idx:
            other_e, other_root = e2_idx, c2
        elif root == e2_idx:
            other_e, other_root = e1_idx, c1
        else:
            raise UnsupportedStructure(
                f"a chain spacer between the benzene ring and the {self.word}'s "
                f"near {self.element} is not supported yet (only a direct ring-{self.element} "
                "bond is, see module docstring)"
            )

        sub_name, sub_compound = name_branch(full_graph, other_root, other_e, {}, mol=mol)
        if sub_compound:
            raise UnsupportedStructure(f"a branched alkylpoly{self.stem} substituent is not supported yet")
        return f"({sub_name}{self._prefix_word(n)})benzene"

    def name(self, mol) -> str:
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_benzene_ring_chain(mol, ring_atoms)
        e1_idx, e2_idx, c1, c2, n = self._validate_and_find(mol)
        bare_name = self._prefix_word(n)
        full_graph = adjacency(mol)
        carbon_graph = carbon_adjacency(mol)

        if c1 is None or c2 is None:
            parent_root, parent_e = (c2, e2_idx) if c1 is None else (c1, e1_idx)
            terminals = {parent_e: bare_name}
        else:
            size1 = len(component_subgraph(carbon_graph, c1))
            size2 = len(component_subgraph(carbon_graph, c2))

            if size1 == size2:
                _, compound_a = name_branch(full_graph, c1, e1_idx, {}, mol=mol)
                _, compound_b = name_branch(full_graph, c2, e2_idx, {}, mol=mol)
                if compound_a and compound_b:
                    raise UnsupportedStructure(
                        f"a {self.word} tied in skeletal-atom count with both sides "
                        "branched is not supported yet"
                    )
                if compound_a:
                    parent_root, parent_e, sub_root, sub_e = c2, e2_idx, c1, e1_idx
                else:
                    parent_root, parent_e, sub_root, sub_e = c1, e1_idx, c2, e2_idx
            elif size1 > size2:
                parent_root, parent_e, sub_root, sub_e = c1, e1_idx, c2, e2_idx
            else:
                parent_root, parent_e, sub_root, sub_e = c2, e2_idx, c1, e1_idx

            sub_name, sub_compound = name_branch(full_graph, sub_root, sub_e, {}, mol=mol)
            if sub_compound:
                raise UnsupportedStructure(f"a branched alkylpoly{self.stem} substituent is not supported yet")
            terminals = {parent_e: sub_name + bare_name}

        parent_carbon_graph = component_subgraph(carbon_graph, parent_root)
        chain, name = self._name_parent_chain(full_graph, parent_carbon_graph, terminals, bare_name, mol=mol)

        stereo = specified_stereocenters(mol)
        if stereo is None:
            return name

        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if any(atom not in position_of for atom, _ in stereo):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        labels = sorted((position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{name}"
