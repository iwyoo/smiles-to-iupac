"""A carboxylic acid combined with sulfonic/sulfinic/seleninic acid groups on one acyclic chain (P-65.1.1/P-65.3.1):
the carboxylic acid is the principal characteristic group, the other acid cited as 'sulfo'/'sulfino'/'selenino'."""

from dataclasses import dataclass


from ._coexisting_groups import name_via_senior_acyclic, name_via_senior_phenyl_chain
from ._oxo_acid import OxoAcid
from ._common import (
    UnsupportedStructure,
    is_plain_benzene_ring,
    non_single_bonds,
    validate_allowed_atoms,
)
from ._carboxylic_acid import _name_acyclic_carboxylic_acid, _name_phenyl_chain_carboxylic_acid


@dataclass(frozen=True)
class CarboxylicAcidCombo:
    atomic_num: int
    word: str
    prefix: str
    element: str
    acid: OxoAcid

    def _carboxyl_carbons(self, mol):
        matches = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
            if len(oxygens) != 2:
                continue
            carbonyls = [
                o
                for o in oxygens
                if o.GetDegree() == 1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            ]
            hydroxyls = [
                o
                for o in oxygens
                if o.GetDegree() == 1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
                and o.GetTotalNumHs() == 1
            ]
            if len(carbonyls) != 1 or len(hydroxyls) != 1:
                continue
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) > 1:
                continue
            matches.append(atom)
        return matches

    def has_shape(self, mol) -> bool:
        return bool(self._carboxyl_carbons(mol)) and bool(self.acid._hetero_atoms(mol))

    def _validate_and_collect(self, mol, aromatic_ring_atoms=frozenset()):
        carboxyl_carbons = self._carboxyl_carbons(mol)
        if len(carboxyl_carbons) != 1:
            raise UnsupportedStructure(
                "exactly one carboxylic acid is required; this module only "
                "handles a single carboxylic acid combined with one or more "
                f"{self.word} acids"
            )
        (carboxyl_carbon,) = carboxyl_carbons
        carboxyl_oxygens = {n.GetIdx() for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 8}

        acids = self.acid._hetero_atoms(mol)
        if not acids:
            raise UnsupportedStructure(
                f"no {self.word} acid found; this module only handles a "
                f"carboxylic acid combined with at least one {self.word} acid "
                "(see _carboxylic_acid.py for a plain carboxylic acid)"
            )
        acid_idxs = {hetero.GetIdx() for hetero in acids}
        acid_oxygens = {o.GetIdx() for hetero in acids for o in hetero.GetNeighbors() if o.GetAtomicNum() == 8}

        accounted_oxygen_idxs = carboxyl_oxygens | acid_oxygens

        validate_allowed_atoms(
            mol,
            f"heteroatoms other than the carboxylic/{self.word} acid groups' own "
            "oxygens (P-65.1.1/P-65.3.1) and halogen substituents (P-35.2.1) "
            "are not supported yet",
            [
                (
                    8,
                    accounted_oxygen_idxs,
                    f"an oxygen that isn't part of the carboxylic/{self.word} acid "
                    "groups is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)",
                ),
                (
                    self.atomic_num,
                    acid_idxs,
                    f"a {self.element} atom not shaped like a {self.word} acid group is out of scope for this module",
                ),
            ],
            aromatic_ring_atoms=aromatic_ring_atoms,
        )

        return carboxyl_carbon.GetIdx(), carboxyl_oxygens, acid_idxs

    def _name_phenyl_chain(self, mol, ring_atoms):
        _carboxyl_carbon, _carboxyl_oxygens, acid_idxs = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)

        acid_carbons = {
            n.GetIdx()
            for hetero in acid_idxs
            for n in mol.GetAtomWithIdx(hetero).GetNeighbors()
            if n.GetAtomicNum() == 6
        }
        acid_oxygens = {
            n.GetIdx()
            for hetero in acid_idxs
            for n in mol.GetAtomWithIdx(hetero).GetNeighbors()
            if n.GetAtomicNum() == 8
        }
        return name_via_senior_phenyl_chain(
            _name_phenyl_chain_carboxylic_acid,
            "carboxylic_acid",
            f"{self.word}_acid",
            mol,
            ring_atoms,
            {hetero: f"{self.prefix}" for hetero in acid_idxs},
            required_atoms=acid_carbons,
            extra_accounted_atoms=acid_idxs | acid_oxygens,
        )

    def name(self, mol) -> str:
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_phenyl_chain(mol, ring_atoms)

        carboxyl_carbon, carboxyl_oxygens, acid_idxs = self._validate_and_collect(mol)

        if mol.GetRingInfo().NumRings() != 0:
            raise UnsupportedStructure(
                f"a carboxylic acid/{self.word} acid combination on a ring is out of scope for this acyclic-only module"
            )
        all_non_single = non_single_bonds(mol)
        excluded_from_unsaturation_check = {carboxyl_carbon} | acid_idxs
        if any(
            a not in excluded_from_unsaturation_check and b not in excluded_from_unsaturation_check
            for a, b, _ in all_non_single
        ):
            raise UnsupportedStructure(
                "chain unsaturation (ene/yne) alongside a carboxylic acid/"
                f"{self.word} acid combination is out of scope for this module"
            )

        acid_carbons = {
            n.GetIdx()
            for hetero in acid_idxs
            for n in mol.GetAtomWithIdx(hetero).GetNeighbors()
            if n.GetAtomicNum() == 6
        }
        return name_via_senior_acyclic(
            _name_acyclic_carboxylic_acid,
            "carboxylic_acid",
            f"{self.word}_acid",
            (mol, {carboxyl_carbon}, carboxyl_oxygens, set(), []),
            {hetero: f"{self.prefix}" for hetero in acid_idxs},
            required_atoms=acid_carbons,
        )
