"""Sulfonic/selenonic/telluronic acids R-E(=O)(=O)OH and sulfinic/seleninic/tellurinic acids R-E(=O)OH (P-65.3.1) as
'alkane-sulfonic acid' style names: acyclic chains, monocycles, benzene and heteroaromatic rings bearing the group,
ring-substituent chains and von Baeyer/spiro skeletons. `oxygens` is 3 or 2; a sole stereogenic chalcogen of an -inic
acid is cited as '(R)-'/'(S)-' (P-93.3.4.1)."""

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
    heteroatom_stereo_prefix,
    heteroaromatic_monocycle_name,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_bond_locants,
    ring_chain_attachment,
    ring_branch_attachments,
    ring_hosting_anchors,
    separate_aromatic_monocycles,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._bicyclic import find_bicyclic_core
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    substituents_for_ring,
    format_substituent_prefixes,
    name_branch,
    ring_branch_stereo_display,
    substituents_for_chain,
)

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


@dataclass(frozen=True)
class OxoAcid:
    atomic_num: int
    word: str
    element: str
    formula: str
    oxygens: int

    def _hetero_atoms(self, mol):
        matches = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != self.atomic_num or atom.GetDegree() != 1 + self.oxygens:
                continue
            neighbors = atom.GetNeighbors()
            carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
            oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
            if len(carbons) != 1 or len(oxygens) != self.oxygens:
                continue
            double_os = [
                o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            ]
            hydroxyl_os = [
                o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            ]
            if len(double_os) != self.oxygens - 1 or len(hydroxyl_os) != 1:
                continue
            if any(o.GetDegree() != 1 for o in double_os):
                continue
            (hydroxyl_o,) = hydroxyl_os
            if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
                continue
            matches.append(atom)
        return matches

    def has_shape(self, mol) -> bool:
        return bool(self._hetero_atoms(mol))

    def _validate_and_collect(self, mol, aromatic_ring_atoms=frozenset()):
        hetero_atoms = self._hetero_atoms(mol)
        if not hetero_atoms:
            raise UnsupportedStructure(
                f"no {self.word} acid (-{self.formula}) group found; this module only handles {self.word} acids"
            )
        if len(hetero_atoms) > 1:
            raise UnsupportedStructure(f"more than one {self.word} acid group is out of scope for this module")
        group_atom_idxs = set()
        for hetero in hetero_atoms:
            group_atom_idxs.add(hetero.GetIdx())
            group_atom_idxs.update(n.GetIdx() for n in hetero.GetNeighbors() if n.GetAtomicNum() == 8)

        has_carbon = False
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 6:
                has_carbon = True
                if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                    raise UnsupportedStructure("aromatic rings are out of scope for this module")
            elif atomic_num in HALOGEN_PREFIXES:
                if atom.GetDegree() != 1:
                    raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
            elif atom.GetIdx() not in group_atom_idxs and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    f"heteroatoms other than a {self.word} acid group (P-65.3.1), "
                    "a heteroaromatic ring's own heteroatom (P-29.3.4.1), and "
                    "halogen substituents (P-35.2.1) are not supported "
                    "yet -- in particular a coexisting carboxylic acid or "
                    "other characteristic group needs acid-vs-acid Table 3.3 "
                    "seniority handling not yet implemented here"
                )
            if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if not has_carbon:
            raise UnsupportedStructure(
                "a structure with no carbon atom has no hydrocarbon parent hydride to substitute"
            )
        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        (hetero_atom,) = hetero_atoms
        (carbon,) = (n for n in hetero_atom.GetNeighbors() if n.GetAtomicNum() == 6)
        return hetero_atom.GetIdx(), carbon.GetIdx()

    def _reject_ene_carbon(self, graph, acid_carbon, bonds):
        unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
        if acid_carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                f"a {self.word} acid on a carbon that is also part of a C=C/C#C bond is out of scope for this module"
            )

    def _name_from_substituents(self, chain_length, acid_locant, ene_locants, yne_locants, grouped):
        prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
        return prefix + name_from_substituents(
            chain_length, ene_locants, yne_locants, f"{self.word} acid", [acid_locant], substituted=bool(grouped)
        )

    def _candidate_key(self, chain_length, acid_locant, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._name_from_substituents(chain_length, acid_locant, ene_locants, yne_locants, grouped)
        return (
            (
                acid_locant,
                combined_locant_set,
                ene_locant_set,
                -total_count,
                locant_set,
                citation_locants,
                name,
            ),
            name,
        )

    def _ring_name_from_substituents(self, ring_size, acid_locant, ene_locants, yne_locants, grouped):
        total_subs = sum(len(info["locants"]) for info in grouped.values())
        prefix = format_substituent_prefixes(grouped)
        return ring_name_from_substituents(
            ring_size, ene_locants, yne_locants, prefix, total_subs, f"{self.word} acid", [acid_locant]
        )

    def _ring_candidate_key(self, ring_size, acid_locant, ene_locants, yne_locants, substituents):
        grouped = group_substituents(substituents)
        locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
        combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
        ene_locant_set = lowest_locant_set(ene_locants)
        name = self._ring_name_from_substituents(ring_size, acid_locant, ene_locants, yne_locants, grouped)
        return acid_locant, combined_locant_set, ene_locant_set, locant_set, citation_locants, name

    def _ring_branch_stereo_display(self, graph, ring_order, excluded, stereo, halogens, mol=None):
        return ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)

    def _name_cyclic(self, mol, hetero_idx, acid_carbon, stereo=None, bonds=()):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        excluded = {hetero_idx}
        ring_info = mol.GetRingInfo()
        ring_atoms = list(ring_info.AtomRings()[0])
        ring_order = ring_cycle(graph, ring_atoms)
        ring_size = len(ring_order)
        branch_stereo = None
        if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
            branch_stereo = self._ring_branch_stereo_display(graph, ring_order, excluded, stereo, halogens, mol=mol)
            if branch_stereo is None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the ring itself is not supported yet (see P-92)"
                )
        if bonds and any(substituents_for_ring(graph, ring_order, halogens, excluded, mol=mol).values()):
            raise UnsupportedStructure(
                "a substituent alongside both a ring double/triple bond and a "
                f"{self.word} acid is not supported yet (see module docstring)"
            )

        best_key = None
        best_name = None
        best_position_of = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                acid_locant = position_of[acid_carbon]
                substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
                if branch_stereo is not None:
                    branch_ring_atom, display = branch_stereo
                    substituents[position_of[branch_ring_atom]] = [(display, False)]
                ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
                key = self._ring_candidate_key(ring_size, acid_locant, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, key[-1], position_of

        if stereo is not None and branch_stereo is None:
            labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
            prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
            return f"({prefix})-{best_name}"
        return best_name

    def _benzene_name_from_substituents(self, grouped):
        if not grouped:
            return f"benzene{self.word} acid"
        return f"{format_substituent_prefixes(grouped)}benzene-1-{self.word} acid"

    def _benzene_candidate_key(self, acid_locant, substituents):
        grouped = group_substituents(substituents)
        locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
        name = self._benzene_name_from_substituents(grouped)
        return acid_locant, locant_set, citation_locants, name

    def _name_benzene(self, mol, ring_atoms, exempt_atoms=None):
        hetero_idx, acid_carbon = self._validate_and_collect(mol, aromatic_ring_atoms=exempt_atoms or ring_atoms)
        if specified_stereocenters(mol):
            raise UnsupportedStructure(
                f"a specified stereocenter alongside benzene{self.word} acid is not supported yet"
            )

        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        excluded = {hetero_idx}
        ring_order = ring_cycle(graph, list(ring_atoms))
        ring_size = len(ring_order)

        best_key = None
        best_name = None
        for start in range(ring_size):
            rotated = ring_order[start:] + ring_order[:start]
            for candidate in (rotated, list(reversed(rotated))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                acid_locant = position_of[acid_carbon]
                substituents = substituents_for_ring(graph, candidate, halogens, excluded, mol=mol)
                key = self._benzene_candidate_key(acid_locant, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]
        return best_name

    def _name_phenyl_chain(self, mol, ring_atoms):
        hetero_idx, acid_carbon = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
        if specified_stereocenters(mol):
            raise UnsupportedStructure(
                "a specified stereocenter alongside a benzene-ring-substituent "
                f"{self.word} acid chain is not supported yet"
            )
        excluded = {hetero_idx}
        non_ring_unsaturation = [
            b
            for b in non_single_bonds(mol)
            if hetero_idx not in (b[0], b[1]) and b[0] not in ring_atoms and b[1] not in ring_atoms
        ]
        if non_ring_unsaturation:
            raise UnsupportedStructure(
                f"chain unsaturation alongside a benzene-ring-substituent {self.word} acid chain is not supported yet"
            )

        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
        attachment = ring_branch_attachments(mol, graph, rings)
        if not attachment:
            raise UnsupportedStructure(
                "a benzene ring with more than one non-halogen, non-alkyl "
                f"exocyclic substituent alongside a chain {self.word} acid is "
                "not supported yet"
            )
        chain, branches = longest_branched_chain_through(
            graph, acid_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol)
        )
        branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            acid_locant = position_of[acid_carbon]
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            key, name = self._candidate_key(chain_length, acid_locant, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_ring_substituent_chain(self, mol, hetero_idx, acid_carbon):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        ring_atoms = set(mol.GetRingInfo().AtomRings()[0])
        excluded = {hetero_idx}

        attachment = ring_chain_attachment(graph, ring_atoms, excluded)
        if attachment is None:
            raise UnsupportedStructure("a ring with more than one exocyclic branch is not supported yet")
        ring_atom, chain_root = attachment

        chain, branches = longest_branched_chain_through(
            graph, acid_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol)
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
            acid_locant = position_of[acid_carbon]
            substituents = {
                position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
                for atom, roots in branches_by_atom.items()
            }
            substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
            key, name = self._candidate_key(chain_length, acid_locant, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    def _name_von_baeyer_or_spiro(self, mol, hetero_idx, acid_carbon, bonds, stereo):
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
                acid_carbon,
                {hetero_idx},
                f"{self.word} acid",
                f"{self.word} acid",
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
                acid_carbon,
                {hetero_idx},
                f"{self.word} acid",
                f"{self.word} acid",
                spiro_atom,
                elide_e=False,
                stereo=stereo,
            )

        raise UnsupportedStructure(
            f"polycyclic and fused-ring {self.word} acids are not supported yet "
            "(P-23/P-25 numbering integration with a suffix group is future "
            "work)"
        )

    def name(self, mol) -> str:
        if self.oxygens == 3:
            return self._name_core(mol)
        centers = self._hetero_atoms(mol)
        if len(centers) != 1:
            return self._name_core(mol)
        hetero_idx = centers[0].GetIdx()
        stereo = specified_stereocenters(mol)
        if stereo is None or all(atom != hetero_idx for atom, _ in stereo):
            return self._name_core(mol)
        prefix = heteroatom_stereo_prefix(mol, hetero_idx)
        achiral = Chem.Mol(mol)
        achiral.GetAtomWithIdx(hetero_idx).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        return prefix + self._name_core(achiral, hetero_prefix=prefix)

    def _name_core(self, mol, hetero_prefix="") -> str:
        aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
        if aromatic_rings is not None:
            union = set().union(*aromatic_rings)
            hetero_idx, _ = self._validate_and_collect(mol, aromatic_ring_atoms=union)
            anchors = [hetero_idx]
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
                hetero_idx, acid_carbon = self._validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
                if acid_carbon in ring_atoms:
                    if is_heteroaromatic:
                        raise UnsupportedStructure(
                            f"a {self.word} acid directly attached to a heteroaromatic ring is not supported yet"
                        )
                    return self._name_benzene(mol, ring_atoms)
                return self._name_phenyl_chain(mol, ring_atoms)
        hetero_idx, acid_carbon = self._validate_and_collect(mol)
        stereo = specified_stereocenters(mol)
        graph = adjacency(mol)
        all_non_single = non_single_bonds(mol)
        bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and hetero_idx not in (b[0], b[1])]
        if len(bonds) != len(all_non_single) - (self.oxygens - 1):
            raise UnsupportedStructure(
                "a bond order other than single, double, or triple is not supported (see P-31.1.1.1)"
            )
        self._reject_ene_carbon(graph, acid_carbon, bonds)

        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings > 1:
            if hetero_prefix:
                raise UnsupportedStructure(
                    "a specified stereocenter alongside a von Baeyer bicyclic/"
                    f"polycyclic or monospiro {self.word} acid is not supported yet "
                    "(see P-92)"
                )
            return self._name_von_baeyer_or_spiro(mol, hetero_idx, acid_carbon, bonds, stereo)
        if num_rings == 1:
            ring_atoms = set(ring_info.AtomRings()[0])
            ring_bonds = [b for b in bonds if b[0] in ring_atoms and b[1] in ring_atoms]
            if any(order == _YNE_ORDER for _, _, order in ring_bonds):
                raise UnsupportedStructure(
                    f"a ring triple bond (cycloalkyne) alongside a {self.word} "
                    "acid is not supported yet -- only a ring double bond is "
                    "in scope for this first pass (see P-31.1.3)"
                )
            if acid_carbon not in ring_atoms:
                if not bonds:
                    if stereo is not None:
                        raise UnsupportedStructure(
                            "a stereocenter on a substituent branch rather "
                            "than the ring itself is not supported yet (see "
                            "P-92)"
                        )
                    return self._name_ring_substituent_chain(mol, hetero_idx, acid_carbon)
                raise UnsupportedStructure(
                    f"a {self.word} acid on a substituent branch chain rather than the ring itself is not supported yet"
                )
            return self._name_cyclic(mol, hetero_idx, acid_carbon, stereo, ring_bonds)

        return self._name_acyclic(mol, hetero_idx, acid_carbon, bonds, stereo)

    def _name_acyclic(
        self, mol, hetero_idx, acid_carbon, bonds, stereo=None, extra_names=None, required_atoms=frozenset()
    ):
        graph = adjacency(mol)
        halogens = {**halogen_substituents(mol), **(extra_names or {})}
        excluded = {hetero_idx}
        chains = all_chains(carbon_adjacency(mol))
        stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

        eligible = []
        for chain in chains:
            chain_set = set(chain)
            if acid_carbon not in chain_set:
                continue
            if not required_atoms <= chain_set:
                continue
            if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
                continue
            eligible.append(chain)
        if not eligible:
            if stereo is not None and any(acid_carbon in c and required_atoms <= set(c) for c in chains):
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than the "
                    "principal chain is not supported yet (see P-92)"
                )
            raise UnsupportedStructure(
                f"the {self.word}-acid-bearing carbon (and/or a multiple bond) "
                "does not lie on a single longest carbon chain; a shorter "
                "principal chain is not supported yet"
            )

        best_key = None
        best_name = None
        best_position_of = None
        chain_length = max(len(c) for c in eligible)
        eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
        for chain in eligible:
            for candidate in (chain, list(reversed(chain))):
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                acid_locant = position_of[acid_carbon]
                ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
                substituents = substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
                key, name = self._candidate_key(chain_length, acid_locant, ene_locants, yne_locants, substituents)
                if best_key is None or key < best_key:
                    best_key, best_name, best_position_of = key, name, position_of

        if stereo is not None:
            labels = sorted((best_position_of[atom], code) for atom, code in stereo)
            prefix = ",".join(f"{locant}{code}" for locant, code in labels)
            return f"({prefix})-{best_name}"
        return best_name
