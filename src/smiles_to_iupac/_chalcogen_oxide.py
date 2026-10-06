"""Sulfones/selenones/tellurones R-E(=O)(=O)-R' and sulfoxides/selenoxides/telluroxides R-E(=O)-R' (P-63.6), named as
'1-(alkane-onyl)alkane' / '1-(alkane-inyl)alkane' from two unbranched, saturated, acyclic chains."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import UnsupportedStructure, non_single_bonds, unbranched_chain_length
from ._numerals import alkane_name


@dataclass(frozen=True)
class ChalcogenOxide:
    atomic_num: int
    element: str
    word: str
    prefix: str
    oxygens: int

    def _atoms(self, mol):
        matches = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != self.atomic_num or atom.GetDegree() != 2 + self.oxygens:
                continue
            neighbors = atom.GetNeighbors()
            carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
            oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
            if len(carbons) != 2 or len(oxygens) != self.oxygens:
                continue
            if any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 2.0 or o.GetDegree() != 1
                for o in oxygens
            ):
                continue
            matches.append(atom)
        return matches

    def has_shape(self, mol) -> bool:
        return bool(self._atoms(mol))

    def name(self, mol) -> str:
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure(
                "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
            )
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in (6, 8, self.atomic_num):
                oxygen_word = "oxygens" if self.oxygens > 1 else "oxygen"
                raise UnsupportedStructure(
                    f"heteroatoms other than the {self.word}'s own {self.element} and "
                    f"{oxygen_word} are not supported yet (P-63.6 is restricted to a "
                    f"plain acyclic {self.word} here)"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                )
        centers = self._atoms(mol)
        if len(centers) != 1:
            raise UnsupportedStructure(f"more than one {self.word} group is out of scope for this module")
        (center,) = centers
        group_idxs = {center.GetIdx()} | {n.GetIdx() for n in center.GetNeighbors() if n.GetAtomicNum() == 8}
        if any(a not in group_idxs and b not in group_idxs for a, b, _ in non_single_bonds(mol)):
            raise UnsupportedStructure(
                "unsaturation is not supported by this module (P-63.6's scope "
                "here is limited to two saturated chains)"
            )
        if mol.GetRingInfo().NumRings() > 0:
            raise UnsupportedStructure("rings are not supported by this module yet")

        center_idx = center.GetIdx()
        c1_idx, c2_idx = (n.GetIdx() for n in center.GetNeighbors() if n.GetAtomicNum() == 6)

        len1 = unbranched_chain_length(mol, c1_idx, center_idx)
        len2 = unbranched_chain_length(mol, c2_idx, center_idx)
        if len1 is None or len2 is None:
            raise UnsupportedStructure("a branched R or R' group is out of scope for this module (see module docstring)")

        parent_len, acyl_len = (len1, len2) if len1 >= len2 else (len2, len1)

        acyl_prefix = f"({alkane_name(acyl_len)}{self.prefix})"
        parent = alkane_name(parent_len)
        if parent_len <= 2:
            return acyl_prefix + parent
        return f"1-{acyl_prefix}{parent}"
