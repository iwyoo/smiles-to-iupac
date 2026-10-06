"""Chalcogen analogues of ketones (P-64.6.1): '-thione' (C=S), '-selone' (C=Se), '-tellone' (C=Te).
One implementation mirrors `_ketone.py` with a non-eliding consonant suffix (P-16.3.3); per-element words live in
`Chalcogenone`. Scope: acyclic chains, saturated monocycles (one ring C=C), benzene/cycloalkyl-substituted chains and a
single C=X on a von Baeyer/monospiro skeleton, with no other heteroatom but halogens."""

from dataclasses import dataclass

from rdkit import Chem

from ._bicyclic import find_bicyclic_core
from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    YNE_BOND_ORDER,
    UnsupportedStructure,
    adjacency,
    all_chains,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    most_multiple_bonds,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ordered_chain,
    ring_bond_locants,
    ring_chain_attachment,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    ring_branch_stereo_display,
    substituents_for_chain,
    substituents_for_ring,
)


@dataclass(frozen=True)
class Chalcogenone:
    atomic_num: int
    symbol: str
    element: str
    suffix: str
    ylidene: str
    ether: str
    hydride: str
    aldehyde: str

    def has_shape(self, mol) -> bool:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != self.atomic_num:
                continue
            if atom.GetDegree() == 1:
                (bond,) = atom.GetBonds()
                if bond.GetBondTypeAsDouble() == 2.0:
                    return True
        return False

    def _validate_and_collect(self, mol, aromatic_ring_atoms=frozenset()):
        word = self.suffix
        allowed = {6, self.atomic_num, *HALOGEN_PREFIXES}
        found = set()
        has_carbon = False
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in allowed:
                raise UnsupportedStructure(
                    f"heteroatoms other than a {word} {self.element} (P-64.6.1) and "
                    "halogen substituents (P-35.2.1) are not supported yet -- "
                    "in particular, a coexisting ketone C=O or hydroxyl -OH "
                    "needs Table 3.3 seniority-coexistence handling not yet "
                    f"implemented for {word}s"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atomic_num == 6:
                has_carbon = True
                if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                    raise UnsupportedStructure(
                        "aromatic rings are out of scope for this module (see the separate aromatic-ring module)"
                    )
            elif atomic_num == self.atomic_num:
                if atom.GetDegree() != 1:
                    raise UnsupportedStructure(
                        f"a {self.element} bonded to more than one heavy atom (e.g. a "
                        f"{self.ether}) is out of scope; only an isolated "
                        f"{word} is supported (P-64.6.1)"
                    )
                (bond,) = atom.GetBonds()
                (carbon,) = atom.GetNeighbors()
                if carbon.GetAtomicNum() != 6:
                    raise UnsupportedStructure(f"a {word} {self.element} must be attached to a carbon atom")
                if bond.GetBondTypeAsDouble() != 2.0:
                    raise UnsupportedStructure(
                        f"a {self.element} that isn't a {word} (C={self.symbol}) is out of scope for "
                        f"this module (e.g. a {self.hydride})"
                    )
                if carbon.GetIsAromatic():
                    raise UnsupportedStructure(f"a {word} on an aromatic ring is out of scope for this module")
                carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
                if len(carbon_neighbors) != 2:
                    raise UnsupportedStructure(
                        f"a {word} carbon with fewer than two carbon neighbors "
                        f"(a {self.aldehyde}) is a different suffix, which "
                        "this module does not attempt to disambiguate"
                    )
                found.add(atom.GetIdx())
            else:
                if atom.GetDegree() != 1:
                    raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
        if not has_carbon:
            raise UnsupportedStructure(
                "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
            )
        if not found:
            raise UnsupportedStructure(
                f"no {word} (C={self.symbol}) group found; this module only handles {word}s"
            )
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")
        return found

    def _name_from_substituents(self, chain_length, locants, ene_locants, yne_locants, grouped):
        return format_substituent_prefixes(grouped) + name_from_substituents(
            chain_length,
            ene_locants,
            yne_locants,
            multiplied_word(len(locants), self.suffix),
            locants,
            substituted=bool(grouped),
        )

    def _candidate_key(self, chain_length, locants, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
        group_locant_set = lowest_locant_set(locants)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._name_from_substituents(chain_length, locants, ene_locants, yne_locants, grouped)
        return (
            (
                group_locant_set,
                combined_locant_set,
                ene_locant_set,
                -total_count,
                locant_set,
                citation_locants,
                name,
            ),
            name,
        )

    @staticmethod
    def _locants(position_of, groups, graph):
        locants = []
        for s in groups:
            (carbon,) = graph[s]
            if carbon not in position_of:
                return None
            locants.append(position_of[carbon])
        return locants

    def _name_acyclic(self, mol, groups, bonds, stereo=None):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        chains = all_chains(carbon_adjacency(mol))
        stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

        eligible = []
        for chain in chains:
            position_of = {atom: i + 1 for i, atom in enumerate(chain)}
            if self._locants(position_of, groups, graph) is None:
                continue
            chain_set = set(chain)
            if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
                continue
            eligible.append(chain)
        if not eligible:
            if stereo is not None and any(
                self._locants({a: i + 1 for i, a in enumerate(c)}, groups, graph) is not None for c in chains
            ):
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the "
                    "principal chain is not supported yet (see P-92)"
                )
            raise UnsupportedStructure(
                f"not every {self.suffix}-bearing carbon (and/or multiple bond) lies "
                "on a single longest carbon chain; a shorter principal chain "
                f"capturing more C={self.symbol} groups is not supported yet"
            )

        best_key = None
        best_name = None
        best_position_of = None
        chain_length = max(len(c) for c in eligible)
        eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
        for chain in eligible:
            for candidate in (chain, list(reversed(chain))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                locants = self._locants(position_of, groups, graph)
                ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
                substituents = substituents_for_chain(graph, candidate, halogens, groups, mol=mol)
                key, name = self._candidate_key(chain_length, locants, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, name, position_of

        if stereo is not None:
            labels = sorted((best_position_of[atom], code) for atom, code in stereo)
            prefix = ",".join(f"{locant}{code}" for locant, code in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def _ring_name_from_substituents(self, ring_size, locants, ene_locants, yne_locants, grouped):
        total_subs = sum(len(info["locants"]) for info in grouped.values())
        prefix = format_substituent_prefixes(grouped)
        return ring_name_from_substituents(
            ring_size,
            ene_locants,
            yne_locants,
            prefix,
            total_subs,
            multiplied_word(len(locants), self.suffix),
            locants,
        )

    def _ring_candidate_key(self, ring_size, locants, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
        group_locant_set = lowest_locant_set(locants)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._ring_name_from_substituents(ring_size, locants, ene_locants, yne_locants, grouped)
        return group_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name

    def _name_cyclic(self, mol, groups, stereo=None, bonds=()):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_info = mol.GetRingInfo()
        ring_atoms = list(ring_info.AtomRings()[0])
        ring_order = ring_cycle(graph, ring_atoms)
        ring_size = len(ring_order)
        branch_stereo = None
        if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
            branch_stereo = ring_branch_stereo_display(graph, ring_order, groups, stereo, halogens, mol=mol)
            if branch_stereo is None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the ring "
                    "itself is not supported yet (see P-92)"
                )
        if bonds and any(substituents_for_ring(graph, ring_order, halogens, groups, mol=mol).values()):
            raise UnsupportedStructure(
                "a substituent alongside both a ring double/triple bond and a "
                f"{self.suffix} is not supported yet (see module docstring)"
            )

        best_key = None
        best_name = None
        best_position_of = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                locants = self._locants(position_of, groups, graph)
                if locants is None:
                    raise UnsupportedStructure(
                        f"a {self.suffix} not on the ring itself (e.g. on a substituent branch) is not supported yet"
                    )
                substituents = substituents_for_ring(graph, candidate, halogens, groups, mol=mol)
                if branch_stereo is not None:
                    branch_ring_atom, display = branch_stereo
                    substituents[position_of[branch_ring_atom]] = [(display, False)]
                ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
                key = self._ring_candidate_key(ring_size, locants, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, key[-1], position_of

        if stereo is not None and branch_stereo is None:
            labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
            prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def _name_phenyl_chain(self, mol, ring_atoms):
        groups = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
        if len(groups) != 1:
            raise UnsupportedStructure(
                f"more than one {self.suffix} group alongside a benzene-ring substituent is not supported yet"
            )
        if specified_stereocenters(mol):
            raise UnsupportedStructure(
                "a specified stereocenter alongside a benzene-ring-substituent "
                f"{self.suffix} chain is not supported yet"
            )
        non_ring_unsaturation = [
            b
            for b in non_single_bonds(mol)
            if b[0] not in groups and b[1] not in groups and b[0] not in ring_atoms and b[1] not in ring_atoms
        ]
        if non_ring_unsaturation:
            raise UnsupportedStructure(
                f"chain unsaturation alongside a benzene-ring-substituent {self.suffix} chain is not supported yet"
            )

        graph = adjacency(mol)
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            raise UnsupportedStructure(
                "a benzene ring with more than one exocyclic substituent "
                f"alongside a chain {self.suffix} is not supported yet"
            )
        ring_atom, chain_root = attachment
        (heteroatom,) = groups
        (carbonyl_carbon,) = graph[heteroatom]
        if carbonyl_carbon == chain_root:
            raise UnsupportedStructure(
                f"a {self.suffix} carbon directly attached to the benzene ring (an "
                f"aryl {self.suffix}) is out of scope for this module (see the "
                "separate aromatic-ring module)"
            )

        chain, branches = longest_branched_chain_through(
            graph, carbonyl_carbon, ring_atoms, groups, halogens=halogen_substituents(mol)
        )
        branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

        chain_length = len(chain)
        halogens = halogen_substituents(mol)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locants = self._locants(position_of, groups, graph)
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            key, name = self._candidate_key(chain_length, locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_ring_substituent_chain(self, mol, groups):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

        attachment = ring_chain_attachment(graph, ring_atoms, groups)
        if attachment is None:
            raise UnsupportedStructure("a ring with more than one exocyclic branch is not supported yet")
        ring_atom, chain_root = attachment
        (heteroatom,) = groups
        (carbonyl_carbon,) = graph[heteroatom]

        chain, branches = longest_branched_chain_through(
            graph, carbonyl_carbon, ring_atoms, groups, halogens=halogen_substituents(mol)
        )
        branches_by_atom = {
            chain[position - 1]: [r for r in roots if r != ring_atom] for position, roots in branches.items()
        }
        branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}

        ring_name = "cyclo" + alkyl_name(len(ring_atoms))
        chain_length = len(chain)

        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locants = self._locants(position_of, groups, graph)
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
            key, name = self._candidate_key(chain_length, locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_ring_with_chain(self, mol, groups):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

        attachment = ring_chain_attachment(graph, ring_atoms, groups)
        if attachment is None:
            raise UnsupportedStructure("a ring with more than one exocyclic branch is not supported yet")
        ring_atom, chain_root = attachment
        chain = ordered_chain(graph, chain_root, ring_atom, groups)
        if chain is None:
            raise UnsupportedStructure("a branched substituent chain hanging off the ring is not supported yet")

        chain_set = set(chain)
        chain_groups = {o for o in groups if next(iter(graph[o])) in chain_set}
        ring_groups = groups - chain_groups
        if len(ring_groups) < len(chain_groups):
            ring_name, ring_is_compound = name_branch(
                graph, ring_atom, chain_root, {**halogens, **{o: self.ylidene for o in ring_groups}}, mol=mol
            )
            chain_length = len(chain)
            best_key = None
            best_name = None
            for candidate in (chain, list(reversed(chain))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                locants = self._locants(position_of, chain_groups, graph)
                substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
                key, name = self._candidate_key(chain_length, locants, [], [], substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, name
            return best_name

        chain_name, chain_is_compound = name_branch(
            graph, chain_root, ring_atom, {**halogens, **{o: self.ylidene for o in chain_groups}}, mol=mol
        )

        ring_order = ring_cycle(graph, list(ring_atoms))
        ring_size = len(ring_order)
        best_key = None
        best_name = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                locants = self._locants(position_of, ring_groups, graph)
                substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
                key = self._ring_candidate_key(ring_size, locants, [], [], substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]
        return best_name

    def name(self, mol) -> str:
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            if is_plain_benzene_ring(mol, ring_atoms):
                return self._name_phenyl_chain(mol, ring_atoms)
        groups = self._validate_and_collect(mol)
        stereo = specified_stereocenters(mol)
        graph = adjacency(mol)
        all_non_single = [b for b in non_single_bonds(mol) if b[0] not in groups and b[1] not in groups]
        bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
        if len(bonds) != len(all_non_single):
            raise UnsupportedStructure(
                "a bond order other than single, double, or triple is not supported (see P-31.1.1.1)"
            )

        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings == 0:
            return self._name_acyclic(mol, groups, bonds, stereo)
        if num_rings == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
            if any(order == YNE_BOND_ORDER for _, _, order in ring_bonds):
                raise UnsupportedStructure(
                    f"a ring triple bond (cycloalkyne) alongside a {self.suffix} is "
                    "not supported yet -- only a ring double bond is in scope "
                    "for this first pass (see P-31.1.3)"
                )
            ring_groups = {s for s in groups if next(iter(graph[s])) in ring_atoms}
            if not ring_groups:
                return self._name_acyclic(mol, groups, bonds, stereo)
            if not bonds and ring_groups and ring_groups != groups:
                if stereo is not None:
                    raise UnsupportedStructure(
                        f"a stereocenter alongside a ring-vs-chain {self.suffix} comparison is not supported yet (see P-92)"
                    )
                return self._name_ring_with_chain(mol, groups)
            return self._name_cyclic(mol, groups, stereo, ring_bonds)
        return self._name_von_baeyer_or_spiro(mol, groups, stereo, bonds)

    def _name_von_baeyer_or_spiro(self, mol, groups, stereo, bonds):
        if len(groups) != 1:
            raise UnsupportedStructure(
                f"more than one {self.suffix} on a von Baeyer bicyclic/polycyclic or "
                "monospiro ring system is not supported yet"
            )

        (heteroatom,) = groups
        graph = adjacency(mol)
        (carbonyl_carbon,) = graph[heteroatom]

        bicyclic_core = find_bicyclic_core(mol)
        polycyclic_core = None
        von_baeyer_ring_count = None
        if bicyclic_core is None:
            for candidate_ring_count in (3, 4, 5, 6):
                polycyclic_core = find_polycyclic_core(mol, candidate_ring_count)
                if polycyclic_core is not None:
                    von_baeyer_ring_count = candidate_ring_count
                    break
        if bicyclic_core is not None or polycyclic_core is not None:
            return name_von_baeyer_suffix(
                mol,
                carbonyl_carbon,
                groups,
                self.suffix,
                self.suffix,
                bicyclic_core,
                polycyclic_core,
                von_baeyer_ring_count,
                elide_e=False,
                stereo=stereo,
                bonds=bonds,
            )

        if bonds:
            raise UnsupportedStructure("an unsaturated monospiro ring system is not supported yet (see P-31.1.5)")

        spiro_atom = find_monospiro_atom(mol)
        if spiro_atom is not None:
            return name_monospiro_suffix(
                mol, carbonyl_carbon, groups, self.suffix, self.suffix, spiro_atom, elide_e=False, stereo=stereo
            )

        raise UnsupportedStructure(
            f"polycyclic and fused-ring {self.suffix}s are not supported yet (P-23/"
            "P-25 numbering integration with a suffix group is future work)"
        )
