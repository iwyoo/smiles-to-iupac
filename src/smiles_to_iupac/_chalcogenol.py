"""Thiols, selenols and tellurols R-SH / R-SeH / R-TeH (P-63.1.1), the chalcogen analogues of alcohols: '-thiol'/'-selenol'/
'-tellurol' suffix on acyclic chains, monocycles, benzene (benzenethiol), benzene/cycloalkyl-substituted chains and von
Baeyer/spiro skeletons, with no other heteroatom but halogens."""

from dataclasses import dataclass


from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    all_chains,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    heteroaromatic_monocycle_name,
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
    ring_branch_attachments,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    two_separate_rings_with_plain_aromatic_substituent,
)
from ._bicyclic import find_bicyclic_core
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

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


@dataclass(frozen=True)
class Chalcogenol:
    atomic_num: int
    symbol: str
    element: str
    word: str
    phenol: str
    base: str
    stem: str
    article: str
    coexisting: str

    def _validate_and_collect(self, mol, aromatic_ring_atoms=frozenset()):
        hetero_atoms = set()
        has_carbon = False
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atom.GetIdx() in aromatic_ring_atoms:
                if atomic_num == 6:
                    has_carbon = True
                continue
            if atomic_num not in {6, self.atomic_num, *HALOGEN_PREFIXES}:
                raise UnsupportedStructure(
                    f"heteroatoms other than a {self.word} {self.element} (P-63.1.1) and "
                    "halogen substituents (P-35.2.1) are not supported yet -- "
                    f"in particular, a coexisting {self.coexisting} nitrogen needs "
                    "Table 3.3 seniority-coexistence handling not yet "
                    f"implemented for {self.word}s"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
            if atomic_num == 6:
                has_carbon = True
                if atom.GetIsAromatic():
                    raise UnsupportedStructure("aromatic rings are out of scope for this module")
            elif atomic_num == self.atomic_num:
                if atom.GetDegree() != 1:
                    raise UnsupportedStructure(
                        f"a {self.element} bonded to more than one heavy atom (e.g. a "
                        f"{self.base}) is out of scope; only an isolated {self.word} "
                        f"(-{self.symbol}H) is supported (Table 3.3, P-63.1.1)"
                    )
                (bond,) = atom.GetBonds()
                if bond.GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure(f"a {self.element} double-bonded to carbon is not a {self.word}")
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        f"{self.article} -{self.symbol}- atom that isn't a simple {self.word} (-{self.symbol}H) is out of "
                        "scope for this module"
                    )
                (neighbor,) = atom.GetNeighbors()
                if neighbor.GetAtomicNum() != 6:
                    raise UnsupportedStructure(f"a {self.word} must be attached to a carbon atom")
                hetero_atoms.add(atom.GetIdx())
            else:
                if atom.GetDegree() != 1:
                    raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
        if not has_carbon:
            raise UnsupportedStructure(
                "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
            )
        if not hetero_atoms:
            raise UnsupportedStructure(
                f"no {self.word} (-{self.symbol}H) group found; this module only handles {self.word}s"
            )
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")
        return hetero_atoms

    def _reject_ene_carbons(self, graph, hetero_atoms, bonds):
        unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
        for hetero_idx in hetero_atoms:
            (carbon,) = graph[hetero_idx]
            if carbon in unsaturated_atoms:
                raise UnsupportedStructure(
                    f"a {self.word} on a carbon that is also part of a C=C/C#C bond is out of scope for this module"
                )

    def _name_from_substituents(self, chain_length, hetero_locants, ene_locants, yne_locants, grouped):
        prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
        return prefix + name_from_substituents(
            chain_length,
            ene_locants,
            yne_locants,
            multiplied_word(len(hetero_locants), f"{self.word}"),
            hetero_locants,
            substituted=bool(grouped),
        )

    def _candidate_key(self, chain_length, hetero_locants, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
        hetero_locant_set = lowest_locant_set(hetero_locants)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._name_from_substituents(chain_length, hetero_locants, ene_locants, yne_locants, grouped)
        return (
            (
                hetero_locant_set,
                combined_locant_set,
                ene_locant_set,
                -total_count,
                locant_set,
                citation_locants,
                name,
            ),
            name,
        )

    def _hetero_locants(self, position_of, hetero_atoms, graph):
        locants = []
        for hetero in hetero_atoms:
            (carbon,) = graph[hetero]
            if carbon not in position_of:
                return None
            locants.append(position_of[carbon])
        return locants

    def _ring_name_from_substituents(self, ring_size, hetero_locants, ene_locants, yne_locants, grouped):
        total_subs = sum(len(info["locants"]) for info in grouped.values())
        prefix = format_substituent_prefixes(grouped)
        return ring_name_from_substituents(
            ring_size,
            ene_locants,
            yne_locants,
            prefix,
            total_subs,
            multiplied_word(len(hetero_locants), f"{self.word}"),
            hetero_locants,
        )

    def _ring_candidate_key(self, ring_size, hetero_locants, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
        hetero_locant_set = lowest_locant_set(hetero_locants)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._ring_name_from_substituents(ring_size, hetero_locants, ene_locants, yne_locants, grouped)
        return hetero_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name

    def _ring_branch_stereo_display(
        self, graph, ring_order, hetero_atoms, stereo, halogens, mol=None, aromatic_atoms=frozenset()
    ):
        return ring_branch_stereo_display(
            graph, ring_order, hetero_atoms, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms
        )

    def _name_cyclic(self, mol, hetero_atoms, stereo=None, bonds=(), ring_atoms=None, aromatic_atoms=frozenset()):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        if ring_atoms is None:
            ring_atoms = list(mol.GetRingInfo().AtomRings()[0])
        else:
            ring_atoms = list(ring_atoms)
        ring_order = ring_cycle(graph, ring_atoms)
        ring_size = len(ring_order)
        branch_stereo = None
        if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
            branch_stereo = self._ring_branch_stereo_display(
                graph, ring_order, hetero_atoms, stereo, halogens, mol=mol, aromatic_atoms=aromatic_atoms
            )
            if branch_stereo is None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the ring itself is not supported yet (see P-92)"
                )
        if bonds and any(
            substituents_for_ring(
                graph, ring_order, halogens, hetero_atoms, mol=mol, aromatic_atoms=aromatic_atoms
            ).values()
        ):
            raise UnsupportedStructure(
                "a substituent alongside both a ring double/triple bond and a "
                f"{self.word} is not supported yet (see module docstring)"
            )

        best_key = None
        best_name = None
        best_position_of = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                hetero_locants = self._hetero_locants(position_of, hetero_atoms, graph)
                substituents = substituents_for_ring(
                    graph, candidate, halogens, hetero_atoms, mol=mol, aromatic_atoms=aromatic_atoms
                )
                if branch_stereo is not None:
                    branch_ring_atom, display = branch_stereo
                    substituents[position_of[branch_ring_atom]] = [(display, False)]
                ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
                key = self._ring_candidate_key(ring_size, hetero_locants, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, key[-1], position_of

        if stereo is not None and branch_stereo is None:
            labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
            prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def _benzene_name_from_substituents(self, hetero_locants, grouped):
        group_word = multiplied_word(len(hetero_locants), f"{self.word}")
        if not grouped:
            return "benzene" + group_word
        return f"{format_substituent_prefixes(grouped)}benzene{group_word}"

    def _benzene_candidate_key(self, hetero_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
        hetero_locant_set = lowest_locant_set(hetero_locants)
        name = self._benzene_name_from_substituents(hetero_locants, grouped)
        return hetero_locant_set, locant_set, citation_locants, name

    def _name_benzene(self, mol, ring_atoms, exempt_atoms=None):
        hetero_atoms = self._validate_and_collect(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
        if len(hetero_atoms) != 1:
            raise UnsupportedStructure(f"more than one {self.word} directly on the benzene ring is not supported yet")
        if specified_stereocenters(mol):
            raise UnsupportedStructure(f"a specified stereocenter alongside benzene{self.word} is not supported yet")

        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_order = ring_cycle(graph, list(ring_atoms))
        ring_size = len(ring_order)

        best_key = None
        best_name = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                hetero_locants = self._hetero_locants(position_of, hetero_atoms, graph)
                substituents = substituents_for_ring(graph, candidate, halogens, hetero_atoms, mol=mol)
                key = self._benzene_candidate_key(hetero_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]
        return best_name

    def _name_phenyl_chain(self, mol, ring_atoms):
        hetero_atoms = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
        if specified_stereocenters(mol):
            raise UnsupportedStructure(
                f"a specified stereocenter alongside a benzene-ring-substituent {self.word} chain is not supported yet"
            )
        non_ring_unsaturation = [b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms]
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
        if any(chain_root in hetero_atoms for _, chain_root in attachment):
            raise UnsupportedStructure(
                f"a {self.word} directly on the benzene ring ({self.phenol}-type) uses "
                "a separate construction, out of scope for this chain-parent "
                "module"
            )
        anchor_hetero = next(iter(hetero_atoms))
        (anchor_carbon,) = graph[anchor_hetero]
        chain, branches = longest_branched_chain_through(
            graph, anchor_carbon, ring_atoms, hetero_atoms, halogens=halogen_substituents(mol)
        )
        chain_set = set(chain)
        for hetero in hetero_atoms:
            (carbon,) = graph[hetero]
            if carbon not in chain_set:
                raise UnsupportedStructure(
                    f"a {self.word} outside the single unbranched chain hanging off "
                    "the benzene ring is not supported yet"
                )
        branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            hetero_locants = self._hetero_locants(position_of, hetero_atoms, graph)
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            key, name = self._candidate_key(chain_length, hetero_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_ring_substituent_chain(self, mol, hetero_atoms):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

        attachment = ring_chain_attachment(graph, ring_atoms, hetero_atoms)
        if attachment is None:
            raise UnsupportedStructure("a ring with more than one exocyclic branch is not supported yet")
        ring_atom, chain_root = attachment
        anchor_hetero = next(iter(hetero_atoms))
        (anchor_carbon,) = graph[anchor_hetero]
        chain, branches = longest_branched_chain_through(
            graph, anchor_carbon, ring_atoms, hetero_atoms, halogens=halogen_substituents(mol)
        )
        chain_set = set(chain)
        for hetero in hetero_atoms:
            (carbon,) = graph[hetero]
            if carbon not in chain_set:
                raise UnsupportedStructure(
                    f"a {self.word} outside the single branched chain hanging off the ring is not supported yet"
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
            hetero_locants = self._hetero_locants(position_of, hetero_atoms, graph)
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
            key, name = self._candidate_key(chain_length, hetero_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_ring_with_chain(self, mol, hetero_atoms):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

        attachment = ring_chain_attachment(graph, ring_atoms, hetero_atoms)
        if attachment is None:
            raise UnsupportedStructure("a ring with more than one exocyclic branch is not supported yet")
        ring_atom, chain_root = attachment
        chain = ordered_chain(graph, chain_root, ring_atom, hetero_atoms)
        if chain is None:
            raise UnsupportedStructure("a branched substituent chain hanging off the ring is not supported yet")

        chain_set = set(chain)
        chain_hetero_atoms = {hetero for hetero in hetero_atoms if next(iter(graph[hetero])) in chain_set}
        ring_hetero_atoms = hetero_atoms - chain_hetero_atoms
        if len(ring_hetero_atoms) < len(chain_hetero_atoms):
            ring_name, ring_is_compound = name_branch(
                graph,
                ring_atom,
                chain_root,
                {**halogens, **{hetero: f"{self.stem}" for hetero in ring_hetero_atoms}},
                mol=mol,
            )
            chain_length = len(chain)
            best_key = None
            best_name = None
            for candidate in (chain, list(reversed(chain))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                hetero_locants = self._hetero_locants(position_of, chain_hetero_atoms, graph)
                substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
                key, name = self._candidate_key(chain_length, hetero_locants, [], [], substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, name
            return best_name

        chain_name, chain_is_compound = name_branch(
            graph,
            chain_root,
            ring_atom,
            {**halogens, **{hetero: f"{self.stem}" for hetero in chain_hetero_atoms}},
            mol=mol,
        )

        ring_order = ring_cycle(graph, list(ring_atoms))
        ring_size = len(ring_order)
        best_key = None
        best_name = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                hetero_locants = self._hetero_locants(position_of, ring_hetero_atoms, graph)
                substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
                key = self._ring_candidate_key(ring_size, hetero_locants, [], [], substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]
        return best_name

    def _name_von_baeyer_or_spiro(self, mol, hetero_atoms, stereo, bonds):
        if len(hetero_atoms) != 1:
            raise UnsupportedStructure(
                f"more than one {self.word} on a von Baeyer bicyclic/polycyclic or "
                "monospiro ring system is not supported yet"
            )
        (hetero_atom,) = hetero_atoms
        graph = adjacency(mol)
        (group_carbon,) = graph[hetero_atom]

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
                group_carbon,
                hetero_atoms,
                f"{self.word}",
                f"{self.word}",
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
                mol,
                group_carbon,
                hetero_atoms,
                f"{self.word}",
                f"{self.word}",
                spiro_atom,
                elide_e=False,
                stereo=stereo,
            )

        raise UnsupportedStructure(
            f"polycyclic and fused-ring {self.word}s are not supported yet (P-23/"
            "P-25 numbering integration with a suffix group is future work)"
        )

    def name(self, mol) -> str:
        aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
        if aromatic_rings is not None:
            union = set().union(*aromatic_rings)
            anchors = list(self._validate_and_collect(mol, aromatic_ring_atoms=union))
            host = ring_hosting_anchors(mol, adjacency(mol), aromatic_rings, anchors)
            if host is not None:
                return self._name_benzene(mol, host, union)
            return self._name_phenyl_chain(mol, union)
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            is_benzene = is_plain_benzene_ring(mol, ring_atoms)
            is_heteroaromatic = not is_benzene and (
                heteroaromatic_monocycle_name(mol, ring_cycle(adjacency(mol), list(ring_atoms))) is not None
            )
            if is_benzene or is_heteroaromatic:
                if is_benzene:
                    ring_hetero_atoms = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
                    if len(ring_hetero_atoms) == 1:
                        (only_hetero,) = ring_hetero_atoms
                        (only_hetero_carbon,) = adjacency(mol)[only_hetero]
                        if only_hetero_carbon in ring_atoms:
                            return self._name_benzene(mol, ring_atoms)
                return self._name_phenyl_chain(mol, ring_atoms)

        aromatic_shape = None
        if ring_info.NumRings() == 2:
            aromatic_shape = two_separate_rings_with_plain_aromatic_substituent(mol, adjacency(mol))
        aromatic_atoms = aromatic_shape[1] if aromatic_shape is not None else frozenset()

        hetero_atoms = self._validate_and_collect(mol, aromatic_ring_atoms=aromatic_atoms)
        stereo = specified_stereocenters(mol)
        graph = adjacency(mol)
        all_non_single = [b for b in non_single_bonds(mol) if not (b[0] in aromatic_atoms and b[1] in aromatic_atoms)]
        bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
        if len(bonds) != len(all_non_single):
            raise UnsupportedStructure(
                "a bond order other than single, double, or triple is not supported (see P-31.1.1.1)"
            )
        self._reject_ene_carbons(graph, hetero_atoms, bonds)

        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings == 2 and aromatic_shape is not None:
            ring_atoms, _, _, _ = aromatic_shape
            chain_hetero_atoms = hetero_atoms - {
                hetero for hetero in hetero_atoms if next(iter(graph[hetero])) in ring_atoms
            }
            if chain_hetero_atoms:
                raise UnsupportedStructure(
                    f"a {self.word} on a chain hanging off the ring, with the ring "
                    f"itself bearing no {self.word} of its own, alongside this "
                    "two-ring aromatic-substituent shape is not supported yet"
                )
            ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
            return self._name_cyclic(
                mol, hetero_atoms, stereo, ring_bonds, ring_atoms=ring_atoms, aromatic_atoms=aromatic_atoms
            )
        if num_rings > 1:
            return self._name_von_baeyer_or_spiro(mol, hetero_atoms, stereo, bonds)
        if num_rings == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
            if any(order == _YNE_ORDER for _, _, order in ring_bonds):
                raise UnsupportedStructure(
                    f"a ring triple bond (cycloalkyne) alongside a {self.word} is not "
                    "supported yet -- only a ring double bond is in scope for "
                    "this first pass (see P-31.1.3)"
                )
            ring_hetero_atoms = {hetero for hetero in hetero_atoms if next(iter(graph[hetero])) in ring_atoms}
            if not ring_hetero_atoms:
                return self._name_acyclic(mol, hetero_atoms, bonds, stereo)
            if not bonds and ring_hetero_atoms and ring_hetero_atoms != hetero_atoms:
                if stereo is not None:
                    raise UnsupportedStructure(
                        f"a stereocenter alongside a ring-vs-chain {self.word} "
                        "comparison is not supported yet (see P-92)"
                    )
                return self._name_ring_with_chain(mol, hetero_atoms)
            if ring_hetero_atoms != hetero_atoms:
                raise UnsupportedStructure(
                    f"a {self.word} on a substituent branch chain rather than the ring itself is not supported yet"
                )
            return self._name_cyclic(mol, hetero_atoms, stereo, ring_bonds)

        return self._name_acyclic(mol, hetero_atoms, bonds, stereo)

    def _name_acyclic(
        self, mol, hetero_atoms, bonds, stereo=None, extra_names=None, required_atoms=frozenset(), carbon_graph=None
    ):
        graph = adjacency(mol)
        halogens = {**halogen_substituents(mol), **(extra_names or {})}
        chains = all_chains(carbon_graph if carbon_graph is not None else carbon_adjacency(mol))
        stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

        eligible = []
        for chain in chains:
            position_of = {atom: i + 1 for i, atom in enumerate(chain)}
            if self._hetero_locants(position_of, hetero_atoms, graph) is None:
                continue
            chain_set = set(chain)
            if not required_atoms <= chain_set:
                continue
            if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
                continue
            eligible.append(chain)
        if not eligible:
            if stereo is not None and any(
                self._hetero_locants({a: i + 1 for i, a in enumerate(c)}, hetero_atoms, graph) is not None
                and required_atoms <= set(c)
                for c in chains
            ):
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the "
                    "principal chain is not supported yet (see P-92)"
                )
            raise UnsupportedStructure(
                f"not every {self.word}-bearing carbon (and/or multiple bond) lies on "
                "a single longest carbon chain; a shorter principal chain "
                f"capturing more -{self.symbol}H groups, or {self.article} -{self.symbol}H expressed as a "
                f"'{self.stem}' substituent prefix, is not supported yet"
            )

        best_key = None
        best_name = None
        best_position_of = None
        chain_length = max(len(c) for c in eligible)
        eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
        for chain in eligible:
            for candidate in (chain, list(reversed(chain))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                hetero_locants = self._hetero_locants(position_of, hetero_atoms, graph)
                ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
                substituents = substituents_for_chain(graph, candidate, halogens, hetero_atoms, mol=mol)
                key, name = self._candidate_key(chain_length, hetero_locants, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, name, position_of

        if stereo is not None:
            labels = sorted((best_position_of[atom], code) for atom, code in stereo)
            prefix = ",".join(f"{locant}{code}" for locant, code in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def has_shape(self, mol) -> bool:
        return any(atom.GetAtomicNum() == self.atomic_num and not atom.GetIsAromatic() for atom in mol.GetAtoms())
