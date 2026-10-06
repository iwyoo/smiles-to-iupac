"""Cyanate esters 'R cyanate', 'R thiocyanate', 'R selenocyanate', 'R tellurocyanate' (R-X-C#N, Chapter P-6). R is a
plain saturated acyclic alkyl group or one plain benzene ring bonded directly to X, named with `name_branch`."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, is_plain_benzene_ring, non_single_bonds
from ._substituents import name_branch

_YNE_ORDER = 3.0


@dataclass(frozen=True)
class Cyanate:
    atomic_num: int
    symbol: str
    element: str
    word: str

    def _cores(self, mol):
        cores = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != self.atomic_num or atom.GetDegree() != 2:
                continue
            neighbors = atom.GetNeighbors()
            if any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in neighbors
            ):
                continue
            carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
            if len(carbons) != 2:
                continue
            for nitrile_c, alkyl_c in ((carbons[0], carbons[1]), (carbons[1], carbons[0])):
                if nitrile_c.GetDegree() != 2:
                    continue
                other_neighbors = [n for n in nitrile_c.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                if len(other_neighbors) != 1:
                    continue
                (nitrogen,) = other_neighbors
                if nitrogen.GetAtomicNum() != 7 or nitrogen.GetDegree() != 1:
                    continue
                if (
                    mol.GetBondBetweenAtoms(nitrile_c.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble()
                    != _YNE_ORDER
                ):
                    continue
                if nitrogen.GetFormalCharge() != 0 or nitrogen.GetIsotope() != 0:
                    continue
                cores.append((atom.GetIdx(), nitrile_c.GetIdx(), alkyl_c.GetIdx()))
        return cores

    def has_shape(self, mol) -> bool:
        return bool(self._cores(mol))

    def name(self, mol) -> str:
        cores = self._cores(mol)
        if len(cores) != 1:
            raise UnsupportedStructure(
                f"exactly one {self.word} (-{self.symbol}-C#N) group is required; zero or "
                "multiple such groups are not supported yet"
            )
        heteroatom_idx, nitrile_c_idx, alkyl_c_idx = cores[0]

        ring_info = mol.GetRingInfo()
        ring_atoms = set()
        if ring_info.NumRings() == 1:
            candidate_ring_atoms = set(ring_info.AtomRings()[0])
            if not is_plain_benzene_ring(mol, candidate_ring_atoms):
                raise UnsupportedStructure("a non-benzene ring is out of scope for this module")
            if alkyl_c_idx not in candidate_ring_atoms:
                raise UnsupportedStructure(
                    "a benzene ring reached through a chain spacer (rather "
                    f"than bonded directly to the {self.word} {self.element}) is not "
                    "supported yet"
                )
            ring_atoms = candidate_ring_atoms
        elif ring_info.NumRings() > 1:
            raise UnsupportedStructure("more than one ring is out of scope for this module")
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        (nitrogen_idx,) = (
            n.GetIdx() for n in mol.GetAtomWithIdx(nitrile_c_idx).GetNeighbors() if n.GetIdx() != heteroatom_idx
        )
        excluded = {heteroatom_idx, nitrile_c_idx, nitrogen_idx}
        for atom in mol.GetAtoms():
            if atom.GetIdx() in excluded:
                continue
            if atom.GetAtomicNum() != 6:
                raise UnsupportedStructure(
                    f"heteroatoms other than this single {self.word} group are not supported yet"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atom.GetIsAromatic() and atom.GetIdx() not in ring_atoms:
                raise UnsupportedStructure("aromatic rings are out of scope for this module")

        non_single = [
            b
            for b in non_single_bonds(mol)
            if nitrile_c_idx not in (b[0], b[1]) and (b[0] not in ring_atoms or b[1] not in ring_atoms)
        ]
        if non_single:
            raise UnsupportedStructure("unsaturation in the R group is not supported yet")

        r_name, _ = name_branch(adjacency(mol), alkyl_c_idx, heteroatom_idx, {}, ring_atoms, mol=mol)
        return f"{r_name} {self.word}"
