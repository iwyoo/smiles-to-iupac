"""Thioic, selenoic and telluroic acids 'alkane-ethioic O-acid' / 'alkaneselenoic Se-acid' / ... (P-65.1.5): one
saturated chain terminating in -C(=O)EH / -C(=E)OH (E = S/Se/Te), or a branched chain hanging off one benzene ring
bearing only halogen/alkyl substituents. A group directly on the ring stays out of scope."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    non_single_bonds,
    ring_branch_attachments,
    separate_aromatic_monocycles,
    unbranched_chain_length,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch


@dataclass(frozen=True)
class ChalcogenoicAcid:
    atomic_num: int
    symbol: str
    word: str

    def _carbons(self, mol):
        matches = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, self.atomic_num)]
            if len(chalcogens) != 2:
                continue
            double_bonded = []
            single_bonded = []
            for chalcogen in chalcogens:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), chalcogen.GetIdx())
                if chalcogen.GetDegree() != 1:
                    double_bonded = single_bonded = None
                    break
                if bond.GetBondTypeAsDouble() == 2.0:
                    double_bonded.append(chalcogen)
                elif bond.GetBondTypeAsDouble() == 1.0 and chalcogen.GetTotalNumHs() == 1:
                    single_bonded.append(chalcogen)
                else:
                    double_bonded = single_bonded = None
                    break
            if not double_bonded or len(double_bonded) != 1 or len(single_bonded) != 1:
                continue
            (double_atom,) = double_bonded
            (single_atom,) = single_bonded
            if {double_atom.GetAtomicNum(), single_atom.GetAtomicNum()} != {8, self.atomic_num}:
                continue
            label = "O" if single_atom.GetAtomicNum() == 8 else self.symbol
            matches.append((atom, label, double_atom, single_atom))
        return matches

    def has_shape(self, mol) -> bool:
        return bool(self._carbons(mol))

    def _validate_and_collect(self, mol, aromatic_ring_atoms=frozenset()):
        allowed = {6, 8, self.atomic_num, *HALOGEN_PREFIXES}
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in allowed:
                raise UnsupportedStructure(
                    f"heteroatoms other than the {self.word} acid's own chalcogens "
                    "and halogen substituents (P-35.2.1) are not supported yet "
                    f"(P-65.1.5 is restricted to a plain acyclic {self.word} acid "
                    "here)"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atomic_num in HALOGEN_PREFIXES and atom.GetDegree() != 1:
                raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
            if atomic_num == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                )
        matches = self._carbons(mol)
        if len(matches) != 1:
            raise UnsupportedStructure(f"more than one {self.word} acid group is out of scope for this module")
        (acid_carbon, label, double_atom, single_atom) = matches[0]
        acid_atom_idxs = {acid_carbon.GetIdx(), double_atom.GetIdx(), single_atom.GetIdx()}
        if any(
            a not in acid_atom_idxs
            and b not in acid_atom_idxs
            and not (a in aromatic_ring_atoms and b in aromatic_ring_atoms)
            for a, b, _ in non_single_bonds(mol)
        ):
            raise UnsupportedStructure(
                "unsaturation is not supported by this module (P-65.1.5's scope "
                "here is limited to a single saturated chain)"
            )
        return acid_carbon, label, acid_atom_idxs

    def _name_phenyl_chain(self, mol, ring_atoms):
        acid_carbon, label, acid_atom_idxs = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
        attachment = ring_branch_attachments(mol, graph, rings)
        if not attachment:
            raise UnsupportedStructure(
                "a benzene ring with more than one non-halogen, non-alkyl "
                f"exocyclic substituent alongside a chain {self.word} acid is not "
                "supported yet"
            )
        acid_carbon_idx = acid_carbon.GetIdx()
        hetero_idxs = acid_atom_idxs - {acid_carbon_idx}
        chain, branches = longest_branched_chain(
            graph, acid_carbon_idx, ring_atoms, hetero_idxs, halogens=halogen_substituents(mol)
        )
        if len(chain) < 2:
            raise UnsupportedStructure(
                f"a {self.word} acid directly attached to the benzene ring uses a "
                "separate naming construction, out of scope for this "
                "acyclic-chain-parent module"
            )

        chain_length = len(chain)
        substituents = {
            position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
            for position, roots in branches.items()
        }
        grouped = group_substituents(substituents)
        prefix = format_substituent_prefixes(grouped)
        return f"{prefix}{alkane_name(chain_length)}{self.word} {label}-acid"

    def name(self, mol) -> str:
        aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
        if aromatic_rings is not None:
            return self._name_phenyl_chain(mol, set().union(*aromatic_rings))
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure(
                "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
            )
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_phenyl_chain(mol, ring_atoms)
        acid_carbon, label, acid_atom_idxs = self._validate_and_collect(mol)
        if mol.GetRingInfo().NumRings() > 0:
            raise UnsupportedStructure("rings are not supported by this module yet")

        chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(chain_neighbors) > 1:
            raise UnsupportedStructure(f"a {self.word}-acid carbon with more than one carbon neighbor is not valid")

        if chain_neighbors:
            (chain_root,) = chain_neighbors
            chain_length = unbranched_chain_length(mol, chain_root.GetIdx(), acid_carbon.GetIdx())
            if chain_length is None:
                raise UnsupportedStructure("a branched R group is out of scope for this module (see module docstring)")
            chain_length += 1
        else:
            chain_length = 1

        return f"{alkane_name(chain_length)}{self.word} {label}-acid"
