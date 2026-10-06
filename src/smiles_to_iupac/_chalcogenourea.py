"""Thiourea, selenourea and tellurourea (H2N-C(=E)-NH2, E = S/Se/Te) with plain N-alkyl or N-phenyl substituents
(P-66.1.4.1); unlike `_urea.py` no semicarbazide (amino-substituted nitrogen) handling."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    plain_phenyl_substituent_atoms,
    reject_unsaturated_substituents,
)
from ._multiplicative_text import enclose
from ._substituents import alpha_sort_key, name_branch


def _n_substituent_carbons(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx)


def _substituent_names(full_graph, nitrogen_idx, substituent_carbons, aromatic_atoms=frozenset(), mol=None):
    return [name_branch(full_graph, c, nitrogen_idx, {}, aromatic_atoms, mol=mol) for c in substituent_carbons]


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        reached, _ = bfs(carbon_graph, root)
        atoms.update(reached)
    return atoms


def _di_name(name, is_compound):
    return enclose(name) if is_compound else name


def _n_letter_entry(letter, name, is_compound):
    return f"{letter}-({name})" if is_compound else f"{letter}-{name}"


def _n_prefix(letter, entries):
    if not entries:
        return ""
    if len(entries) == 1:
        ((name, is_compound),) = entries
        return _n_letter_entry(letter, name, is_compound)
    (name_a, compound_a), (name_b, compound_b) = entries
    if name_a == name_b:
        return f"{letter},{letter}-di{_di_name(name_a, compound_a)}"
    (a, ca), (b, cb) = sorted(entries, key=lambda e: alpha_sort_key(e[0]))
    return f"{_n_letter_entry(letter, a, ca)}-{_n_letter_entry(letter, b, cb)}"


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
        n1_carbons = _n_substituent_carbons(mol, n1_idx, carbon_idx)
        n2_carbons = _n_substituent_carbons(mol, n2_idx, carbon_idx)

        full_graph = adjacency(mol)
        phenyl_atoms = plain_phenyl_substituent_atoms(mol, full_graph, n1_carbons + n2_carbons)
        if phenyl_atoms and (
            (any(c in phenyl_atoms for c in n1_carbons) and len(n1_carbons) > 1)
            or (any(c in phenyl_atoms for c in n2_carbons) and len(n2_carbons) > 1)
        ):
            raise UnsupportedStructure(
                "a phenyl N-substituent alongside another substituent on the same nitrogen is not supported yet"
            )

        if mol.GetRingInfo().NumRings() > 0:
            all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
            if all_ring_atoms - phenyl_atoms:
                raise UnsupportedStructure(
                    f"a ring-fused {word} (e.g. hydantoin) or a ring "
                    "N-substituent other than a plain, unsubstituted benzene "
                    "ring is out of scope for this module"
                )

        carbon_graph = carbon_adjacency(mol)
        n1_chain_atoms = _substituent_chain_atoms(carbon_graph, n1_carbons)
        n2_chain_atoms = _substituent_chain_atoms(carbon_graph, n2_carbons)
        known_atoms = {carbon_idx, chalcogen_idx, n1_idx, n2_idx} | n1_chain_atoms | n2_chain_atoms
        for atom in mol.GetAtoms():
            if atom.GetIdx() not in known_atoms:
                raise UnsupportedStructure(
                    "a heteroatom or other characteristic group outside the "
                    f"{word} core and its plain N-alkyl substituents is not "
                    "supported yet"
                )

        reject_unsaturated_substituents(mol, n1_chain_atoms - phenyl_atoms)
        reject_unsaturated_substituents(mol, n2_chain_atoms - phenyl_atoms)

        aromatic_atoms = frozenset(phenyl_atoms)
        n1_names = _substituent_names(full_graph, n1_idx, n1_carbons, aromatic_atoms, mol=mol)
        n2_names = _substituent_names(full_graph, n2_idx, n2_carbons, aromatic_atoms, mol=mol)

        if not n1_names and not n2_names:
            return word

        if n1_names and n2_names:
            if len(n1_names) != 1 or len(n2_names) != 1:
                raise UnsupportedStructure(
                    f"a different substituent count on each of {word}'s two "
                    "nitrogens is not supported yet (no confirmed worked "
                    "example settles the locant tie-break for that case)"
                )
            (name_a, compound_a), (name_b, _) = n1_names[0], n2_names[0]
            if name_a == name_b:
                return f"N,N'-di{_di_name(name_a, compound_a)}{word}"
            (first, first_compound), (second, second_compound) = sorted(
                (n1_names[0], n2_names[0]), key=lambda e: alpha_sort_key(e[0])
            )
            first_entry = _n_letter_entry("N", first, first_compound)
            second_entry = _n_letter_entry("N'", second, second_compound)
            return f"{first_entry}-{second_entry}{word}"

        names = n1_names or n2_names
        return f"{_n_prefix('N', names)}{word}"
