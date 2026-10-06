"""Thioate and selenoate anions R-CO-E(-)/R-C(=E)O(-) (E = S/Se, P-72.2.2.2.1.1) named 'alkanethioate'/'alkaneselenoate':
acyclic chains and chains hanging off a benzene ring bearing only halogen/alkyl substituents."""

from dataclasses import dataclass


from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    all_chains,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_branch_attachments,
    separate_aromatic_monocycles,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    substituents_for_chain,
)


@dataclass(frozen=True)
class Ate:
    atomic_num: int
    symbol: str
    word: str

    def _matches(self, mol):
        matches = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, self.atomic_num)]
            if len(chalcogens) != 2:
                continue
            carbonyls = [
                o
                for o in chalcogens
                if o.GetDegree() == 1
                and o.GetFormalCharge() == 0
                and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            ]
            anions = [
                o
                for o in chalcogens
                if o.GetDegree() == 1
                and o.GetFormalCharge() == -1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            ]
            if len(carbonyls) != 1 or len(anions) != 1:
                continue
            if {carbonyls[0].GetAtomicNum(), anions[0].GetAtomicNum()} != {8, self.atomic_num}:
                continue
            matches.append((atom, carbonyls[0], anions[0]))
        return matches

    def has_shape(self, mol) -> bool:
        return bool(self._matches(mol))

    def _find_group(self, mol):
        matches = self._matches(mol)
        if len(matches) != 1:
            raise UnsupportedStructure(
                f"exactly one {self.word} (-CO{self.symbol}-/-C{self.symbol}O-) group is required; zero or "
                "multiple such groups are not supported yet (P-72.2.2.2.1.1)"
            )
        group_carbon, carbonyl_atom, anion_atom = matches[0]
        total_chalcogens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (8, self.atomic_num))
        if total_chalcogens != 2:
            raise UnsupportedStructure(
                f"a chalcogen outside the single {self.word} group's carbonyl/anion "
                "pair is out of scope for this module"
            )
        carbon_neighbors = [n for n in group_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                f"a {self.word} carbon with more than one carbon neighbor is not a valid {self.word} group"
            )
        return group_carbon, carbonyl_atom, anion_atom

    def _name_from_substituents(self, chain_length, ene_locants, yne_locants, grouped):
        return format_substituent_prefixes(grouped) + name_from_substituents(
            chain_length, ene_locants, yne_locants, f"{self.word}"
        )

    def _candidate_key(self, chain_length, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
        return (
            (
                combined_locant_set,
                ene_locant_set,
                -total_count,
                locant_set,
                citation_locants,
                name,
            ),
            name,
        )

    def _name_acyclic(self, mol, group_carbon_idx, excluded_atoms, bonds, stereo=None):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        chains = all_chains(carbon_adjacency(mol))
        stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

        eligible = []
        for chain in chains:
            if group_carbon_idx not in chain:
                continue
            chain_set = set(chain)
            if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
                continue
            eligible.append(chain)
        if not eligible:
            if stereo is not None and any(group_carbon_idx in c for c in chains):
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the "
                    "principal chain is not supported yet (see P-92)"
                )
            raise UnsupportedStructure(
                f"the {self.word} carbon (and/or multiple bonds) does not lie on a single longest carbon chain"
            )

        best_key = None
        best_name = None
        best_position_of = None
        chain_length = max(len(c) for c in eligible)
        eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
        for chain in eligible:
            for candidate in (chain, list(reversed(chain))):
                if candidate[0] != group_carbon_idx:
                    continue
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
                substituents = substituents_for_chain(graph, candidate, halogens, excluded_atoms, mol=mol)
                key, name = self._candidate_key(chain_length, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, name, position_of

        if stereo is not None:
            labels = sorted((best_position_of[atom], code) for atom, code in stereo)
            prefix = ",".join(f"{locant}{code}" for locant, code in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def _validate_and_prepare(self, mol, aromatic_ring_atoms=frozenset()):
        group_carbon, carbonyl_atom, anion_atom = self._find_group(mol)
        excluded_atoms = {carbonyl_atom.GetIdx(), anion_atom.GetIdx()}
        stereo = specified_stereocenters(mol)

        has_carbon = False
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in {6, 8, self.atomic_num, *HALOGEN_PREFIXES}:
                raise UnsupportedStructure(
                    f"heteroatoms other than the {self.word}'s own chalcogens "
                    "(P-72.2.2.2.1.1) and halogen substituents (P-35.2.1) are "
                    "not supported yet"
                )
            if atom.GetIdx() in excluded_atoms:
                continue
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure(
                    "a charged or isotopically modified atom other than the "
                    f"single {self.word} anion chalcogen is not supported yet"
                )
            if atomic_num == 6:
                has_carbon = True
                if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                    raise UnsupportedStructure(
                        "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                    )
            elif atomic_num not in (8, self.atomic_num) and atom.GetDegree() != 1:
                raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
        if not has_carbon:
            raise UnsupportedStructure(
                "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
            )
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        all_non_single = [
            b
            for b in non_single_bonds(mol)
            if b[0] not in excluded_atoms
            and b[1] not in excluded_atoms
            and (b[0] not in aromatic_ring_atoms or b[1] not in aromatic_ring_atoms)
        ]
        bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
        if len(bonds) != len(all_non_single):
            raise UnsupportedStructure(
                "a bond order other than single, double, or triple is not supported (see P-31.1.1.1)"
            )
        return group_carbon, excluded_atoms, bonds, stereo

    def _name_phenyl_chain(self, mol, ring_atoms):
        group_carbon, excluded_atoms, bonds, stereo = self._validate_and_prepare(mol, aromatic_ring_atoms=ring_atoms)
        if stereo:
            raise UnsupportedStructure(
                f"a specified stereocenter alongside a benzene-ring-substituent {self.word} chain is not supported yet"
            )
        non_ring_unsaturation = [b for b in bonds if b[0] not in ring_atoms or b[1] not in ring_atoms]
        if non_ring_unsaturation:
            raise UnsupportedStructure(
                f"chain unsaturation alongside a benzene-ring-substituent {self.word} chain is not supported yet"
            )

        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
        attachment = ring_branch_attachments(mol, graph, rings)
        if not attachment:
            raise UnsupportedStructure(
                "a benzene ring with more than one non-halogen, non-alkyl "
                f"exocyclic substituent alongside a chain {self.word} is not "
                "supported yet"
            )
        chain, branches = longest_branched_chain(
            graph, group_carbon.GetIdx(), ring_atoms, excluded_atoms, halogens=halogen_substituents(mol)
        )
        if len(chain) < 2:
            raise UnsupportedStructure(
                f"a {self.word} group directly attached to the benzene ring uses a "
                "separate construction, out of scope for this "
                "acyclic-chain-parent module"
            )

        chain_length = len(chain)
        substituents = {
            position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
            for position, roots in branches.items()
        }
        grouped = group_substituents(substituents)
        return self._name_from_substituents(chain_length, [], [], grouped)

    def name(self, mol) -> str:
        aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
        if aromatic_rings is not None:
            return self._name_phenyl_chain(mol, set().union(*aromatic_rings))
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_phenyl_chain(mol, ring_atoms)
        if mol.GetRingInfo().NumRings() > 0:
            raise UnsupportedStructure(f"a {self.word} group on/in a ring is out of scope for this acyclic-only module")

        group_carbon, excluded_atoms, bonds, stereo = self._validate_and_prepare(mol)
        return self._name_acyclic(mol, group_carbon.GetIdx(), excluded_atoms, bonds, stereo)
