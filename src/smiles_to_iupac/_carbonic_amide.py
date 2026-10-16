"""Amides of cyanic acid and of di- and polycarbonic acids (P-66.1.6.1.4, P-66.1.6.2, P-66.1.6.3):
(propan-2-yl)cyanamide, 2-imidodicarbonic diamide, N1-methyl-2-imido-1-thiodicarbonic diamide, 1,2,3-trithiodicarbonic
diamide."""

from rdkit import Chem

from ._chalcogenourea import has_group_senior_to_urea, is_urea_substituent_atom, n_substituent_names
from ._cited_group import cited_group
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents, is_nitro_nitrogen
from ._numerals import numerical_term
from ._substituents import cited_whole_stereo, format_substituent_prefixes, name_branch

_CHALCOGEN = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}
_CARBON_AND_HALOGEN = {6, 9, 17, 35, 53}
_MAXIMUM_CARBONIC_UNITS = 4
_IMIDIC_DERIVATIVE = Chem.MolFromSmarts("[CX3;!R](=[NX2;!R])[NX3,OX2,SX2]")


def _has_stereo(mol):
    return any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )


def _plain(mol, stereo=False):
    return (
        len(Chem.GetMolFrags(mol)) == 1
        and not any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        and (stereo or not _has_stereo(mol))
    )


def plain_arm_atom(atom):
    """An atom that adds no characteristic group to an N-substituent: carbon, halogen, a divalent ether oxygen or sulfide
    sulfur between carbons, or a ring nitrogen of a heterocycle."""
    z = atom.GetAtomicNum()
    if z in _CARBON_AND_HALOGEN:
        return True
    if z == 7:
        return atom.IsInRing()
    return z in (8, 16) and atom.GetDegree() == 2 and all(
        n.GetAtomicNum() == 6 and (bond.GetBondTypeAsDouble() == 1.0 or bond.GetIsAromatic()) for bond in atom.GetBonds() for n in [bond.GetOtherAtom(atom)]
    )


def _arm_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n != blocked:
                seen.add(n)
                stack.append(n)
    return seen


_JUNIOR_ROOTS = {8, 16, 34, 52, 9, 17, 35, 53}


def _junior_arm_atom(atom):
    return is_urea_substituent_atom(atom) or atom.IsInRing()


def _hydrocarbon_arms(mol, graph, nitrogen, skip, junior=False):
    """(name, compound) of each carbon group on `nitrogen` outside `skip` and the atoms they cover, or None when a
    group carries a characteristic group; with `junior` the groups junior to the amide are prefixes (P-66.1.6.1.1.2)."""
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    names, covered = [], set()
    for root in graph[nitrogen]:
        if root in skip:
            continue
        root_atom = mol.GetAtomWithIdx(root)
        if root_atom.GetAtomicNum() != 6 and not (junior and root_atom.GetAtomicNum() in _JUNIOR_ROOTS):
            return None
        arm = _arm_atoms(graph, root, nitrogen)
        if nitrogen in arm or covered & arm or not all(
            plain_arm_atom(mol.GetAtomWithIdx(i)) or (junior and _junior_arm_atom(mol.GetAtomWithIdx(i))) for i in arm
        ):
            return None
        covered |= arm
        if _has_stereo(mol):
            names.append(cited_group(mol, graph, root, nitrogen))
        else:
            names.append(name_branch(graph, root, nitrogen, halogens, aromatic, mol=mol, unsaturated=True))
    return names, covered


def cyanamide_name(mol):
    if not _plain(mol, stereo=True):
        return None
    graph = adjacency(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 2 or atom.IsInRing():
            continue
        triple = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
        if len(triple) != 1:
            continue
        terminal = triple[0].GetOtherAtom(atom)
        amino = next(n for n in atom.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
        if terminal.GetAtomicNum() != 7 or terminal.GetDegree() != 1 or amino.GetAtomicNum() != 7:
            continue
        if amino.IsInRing() or amino.GetIsAromatic() or any(
            b.GetBondTypeAsDouble() != 1.0
            and not (b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(amino).GetAtomicNum() in (6, 16, 34, 52))
            for b in amino.GetBonds()
        ):
            continue
        betas = [n for n in amino.GetNeighbors() if n.GetAtomicNum() == 7 and not n.IsInRing()]
        if len(betas) > 1:
            continue
        if not betas:
            if mol.HasSubstructMatch(_IMIDIC_DERIVATIVE):
                continue
            try:
                (names,) = n_substituent_names(
                    mol, graph, {atom.GetIdx(), terminal.GetIdx(), amino.GetIdx()}, (amino.GetIdx(),), atom.GetIdx(), junior_groups=True
                )
            except UnsupportedStructure:
                continue
            grouped = {}
            for name, compound in names:
                grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(1)
            return (format_substituent_prefixes(grouped, omit_locants=True) if grouped else "") + "cyanamide"
        found = _hydrocarbon_arms(mol, graph, amino.GetIdx(), {atom.GetIdx(), *(b.GetIdx() for b in betas)})
        if found is None:
            continue
        names, covered = found
        entries = [("N" if betas else 1, name, compound) for name, compound in names]
        extra = 3
        if betas:
            beta = betas[0]
            if beta.IsInRing() or beta.GetIsAromatic() or any(b.GetBondTypeAsDouble() != 1.0 for b in beta.GetBonds()):
                continue
            tail = _hydrocarbon_arms(mol, graph, beta.GetIdx(), {amino.GetIdx()})
            if tail is None:
                continue
            entries += [("N'", name, compound) for name, compound in tail[0]]
            covered = covered | tail[1]
            extra = 4
        if len(covered) + extra != mol.GetNumAtoms():
            continue
        grouped = {}
        for locant, name, compound in entries:
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
        prefix = format_substituent_prefixes(grouped, omit_locants=not betas) if grouped else ""
        return prefix + ("cyanohydrazide" if betas else "cyanamide")
    return None


def carbonohydrazonic_diamide_name(mol):
    """C(=N-NH2)(NH2)2 with carbon groups on its nitrogens (P-66.4.2.2): N and N' on the two amide nitrogens, N'' on
    the terminal hydrazone nitrogen, as the locants of an amidrazone follow the nitrogens in the order of the group."""
    if not _plain(mol):
        return None
    graph = adjacency(mol)
    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6 or carbon.IsInRing() or carbon.GetDegree() != 3:
            continue
        bonds = {n.GetIdx(): mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in carbon.GetNeighbors()}
        if any(n.GetAtomicNum() != 7 or n.IsInRing() or n.GetIsAromatic() for n in carbon.GetNeighbors()):
            continue
        imines = [i for i, order in bonds.items() if order == 2.0]
        amides = [i for i, order in bonds.items() if order == 1.0]
        if len(imines) != 1 or len(amides) != 2:
            continue
        imine = mol.GetAtomWithIdx(imines[0])
        terminals = [n for n in imine.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
        if len(terminals) != 1 or terminals[0].GetAtomicNum() != 7 or terminals[0].IsInRing() or terminals[0].GetIsAromatic():
            continue
        terminal = terminals[0]
        nitrogens = [mol.GetAtomWithIdx(i) for i in amides] + [terminal]
        if any(b.GetBondTypeAsDouble() != 1.0 for n in nitrogens for b in n.GetBonds() if b.GetOtherAtomIdx(n.GetIdx()) != carbon.GetIdx()):
            continue
        if any(n.GetAtomicNum() == 7 for amide in nitrogens[:2] for n in amide.GetNeighbors() if n.GetIdx() != carbon.GetIdx()):
            continue
        found = []
        used = {carbon.GetIdx(), imine.GetIdx(), terminal.GetIdx(), *amides}
        for nitrogen, skip in ((nitrogens[0], {carbon.GetIdx()}), (nitrogens[1], {carbon.GetIdx()}), (terminal, {imine.GetIdx()})):
            arms = _hydrocarbon_arms(mol, graph, nitrogen.GetIdx(), skip)
            if arms is None:
                break
            found.append(arms[0])
            used |= arms[1]
        else:
            if len(used) != mol.GetNumAtoms():
                continue
            first, second = sorted(found[:2], key=lambda names: (-len(names), sorted(name for name, _ in names)))
            grouped = {}
            for locant, names in (("N", first), ("N'", second), ("N''", found[2])):
                for name, compound in names:
                    grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
            head = format_substituent_prefixes(grouped) if grouped else ""
            return f"{head}carbonohydrazonic diamide"
    return None


def _carbonyl_chalcogen(atom):
    for b in atom.GetBonds():
        other = b.GetOtherAtom(atom)
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() in _CHALCOGEN and other.GetDegree() == 1:
            return other
    return None


def _carbonic_carbon(atom):
    """A carbon bearing a doubly bonded chalcogen and otherwise only nitrogen, chalcogen or peroxy links: a unit of a
    condensed urea, not an acyl group."""
    return (
        atom.GetAtomicNum() == 6
        and atom.GetDegree() == 3
        and _carbonyl_chalcogen(atom) is not None
        and all(n.GetAtomicNum() != 6 for n in atom.GetNeighbors())
    )


def _linker(mol, carbon, neighbor):
    """(linker atoms, next carbon) when `neighbor` bridges `carbon` to another carbonyl carbon, else None."""
    if neighbor.IsInRing() or neighbor.GetIsAromatic():
        return None
    z = neighbor.GetAtomicNum()
    if z == 7:
        if neighbor.GetDegree() == 2 and neighbor.GetTotalNumHs() == 1:
            return ([neighbor], next(n for n in neighbor.GetNeighbors() if n.GetIdx() != carbon.GetIdx()))
        if neighbor.GetDegree() != 3 or neighbor.GetTotalNumHs():
            return None
        onward = [n for n in neighbor.GetNeighbors() if n.GetIdx() != carbon.GetIdx() and _carbonic_carbon(n)]
        return ([neighbor], onward[0]) if len(onward) == 1 else None
    if z not in _CHALCOGEN or neighbor.GetDegree() != 2:
        return None
    forward = [n for n in neighbor.GetNeighbors() if n.GetIdx() != carbon.GetIdx()][0]
    if forward.GetAtomicNum() == z and z != 8 and forward.GetDegree() == 2:
        beyond = [n for n in forward.GetNeighbors() if n.GetIdx() != neighbor.GetIdx()]
        return ([neighbor, forward], beyond[0])
    return ([neighbor], forward)


def _chain(mol, start):
    """[carbon, linker atoms, carbon, ...] from the terminal carbonyl carbon `start` along the carbonic chain, else
    None."""
    chain, previous, current = [start], None, start
    while True:
        linkers = []
        for n in current.GetNeighbors():
            if previous is not None and n.GetIdx() in previous:
                continue
            found = _linker(mol, current, n)
            if found is not None and _carbonyl_chalcogen(found[1]) is not None and found[1].GetAtomicNum() == 6:
                linkers.append(found)
        if len(linkers) != 1:
            return chain if len(chain) > 1 else None
        atoms, nxt = linkers[0]
        chain.append(atoms)
        chain.append(nxt)
        previous = {a.GetIdx() for a in atoms} | {current.GetIdx()}
        current = nxt


def _terminal_nitrogen(mol, carbon, chalcogen, used):
    ends = [
        n
        for n in carbon.GetNeighbors()
        if n.GetIdx() != chalcogen.GetIdx()
        and n.GetIdx() not in used
        and n.GetAtomicNum() == 7
        and not n.IsInRing()
        and not n.GetIsAromatic()
    ]
    return ends[0] if len(ends) == 1 else None


def polycarbonic_amide_name(mol):
    if not _plain(mol):
        return None
    graph = adjacency(mol)
    for start in mol.GetAtoms():
        if start.GetAtomicNum() != 6 or start.IsInRing() or start.GetDegree() != 3 or _carbonyl_chalcogen(start) is None:
            continue
        chain = _chain(mol, start)
        if chain is None or len(chain) // 2 + 1 > _MAXIMUM_CARBONIC_UNITS:
            continue
        if chain[-1].GetDegree() != 3 or start.GetIdx() > chain[-1].GetIdx():
            continue
        name = _assemble(mol, graph, chain)
        if name is not None:
            return name
    return None


def _assemble(mol, graph, chain):
    carbons = chain[0::2]
    linkers = chain[1::2]
    used = {c.GetIdx() for c in carbons} | {a.GetIdx() for group in linkers for a in group}
    chalcogens = [_carbonyl_chalcogen(c) for c in carbons]
    if any(c.GetDegree() != 3 for c in carbons):
        return None
    used |= {x.GetIdx() for x in chalcogens}
    ends = []
    for carbon, chalcogen in ((carbons[0], chalcogens[0]), (carbons[-1], chalcogens[-1])):
        nitrogen = _terminal_nitrogen(mol, carbon, chalcogen, used)
        if nitrogen is None:
            return None
        ends.append(nitrogen)
    for middle, chalcogen in zip(carbons[1:-1], chalcogens[1:-1]):
        if sum(1 for n in middle.GetNeighbors() if n.GetIdx() != chalcogen.GetIdx()) != 2:
            return None
    arms, hydrazide = [], []
    carbon_set = {c.GetIdx() for c in carbons}
    if has_group_senior_to_urea(mol, used | {n.GetIdx() for n in ends}):
        return None
    linker_arms = {}
    for index, group in enumerate(linkers):
        if len(group) == 1 and group[0].GetAtomicNum() == 7 and group[0].GetDegree() == 3:
            found = _hydrocarbon_arms(mol, graph, group[0].GetIdx(), carbon_set, junior=True)
            if found is None:
                return None
            linker_arms[index] = found[0]
            used |= found[1]
    for nitrogen in ends:
        if any(b.GetBondTypeAsDouble() != 1.0 for b in nitrogen.GetBonds()):
            return None
        betas = [n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7]
        hydrazide.append(bool(betas))
        if len(betas) > 1:
            return None
        found = _hydrocarbon_arms(mol, graph, nitrogen.GetIdx(), carbon_set | {b.GetIdx() for b in betas}, junior=True)
        if found is None:
            return None
        names = {"": found[0]}
        used |= {nitrogen.GetIdx()} | found[1]
        if betas:
            beta = betas[0]
            if beta.IsInRing() or beta.GetIsAromatic() or beta.GetDegree() > 3 or any(b.GetBondTypeAsDouble() != 1.0 for b in beta.GetBonds()):
                return None
            tail = _hydrocarbon_arms(mol, graph, beta.GetIdx(), {nitrogen.GetIdx()})
            if tail is None:
                return None
            names["'"] = tail[0]
            used |= {beta.GetIdx()} | tail[1]
        arms.append(names)
    if hydrazide[0] != hydrazide[1] or len(used) != mol.GetNumAtoms():
        return None

    best = None
    for forward in (True, False):
        order = list(range(len(carbons)))
        if not forward:
            order.reverse()
        positions, replacements, position = [], [], 1
        linker_positions = {}
        for rank, index in enumerate(order):
            if rank:
                position += 1
                joined = min(index, order[rank - 1])
                replacements.extend(_linker_prefix(linkers[joined], position))
                linker_positions[joined] = position
                position += 1
            positions.append(position)
            symbol = _CHALCOGEN[chalcogens[index].GetAtomicNum()]
            if symbol:
                replacements.append((position, symbol, True))
        substituents = [
            (f"N{linker_positions[joined]}", name, compound) for joined, names in linker_arms.items() for name, compound in names
        ]
        for rank, index in enumerate(order):
            if index in (0, len(carbons) - 1):
                for prime, names in arms[0 if index == 0 else 1].items():
                    substituents.extend((f"N{prime}{positions[rank]}", name, compound) for name, compound in names)
        key = (
            tuple(sorted(p for p, _, _ in replacements)),
            tuple(sorted(_locant_number(loc) for loc, _, _ in substituents)),
            tuple(sorted((name, loc) for loc, name, _ in substituents)),
        )
        if best is None or key < best[0]:
            best = (key, replacements, substituents)
    _, replacements, substituents = best
    prefixes = {}
    for locant, name, compound in substituents:
        prefixes.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    head = "-".join(
        part for part in (format_substituent_prefixes(prefixes) if prefixes else "", _replacement_text(replacements)) if part
    )
    return f"{head}{numerical_term(len(carbons))}carbonic {'dihydrazide' if hydrazide[0] else 'diamide'}"


def _locant_number(locant):
    return int(locant.lstrip("N'"))


def _linker_prefix(group, position):
    z = group[0].GetAtomicNum()
    if len(group) == 2:
        return [(position, numerical_term(2) + _CHALCOGEN[z] + "peroxy", False)]
    if z == 7:
        return [(position, "imido", True)]
    symbol = _CHALCOGEN[z]
    return [(position, symbol, True)] if symbol else []


def _replacement_text(replacements):
    grouped = {}
    for locant, name, multipliable in replacements:
        grouped.setdefault((name, multipliable), []).append(locant)
    pieces = []
    for (name, multipliable), locants in sorted(grouped.items(), key=lambda item: item[0][0]):
        locants.sort()
        word = (numerical_term(len(locants)) if multipliable and len(locants) > 1 else "") + name
        pieces.append(",".join(str(x) for x in locants) + "-" + word)
    return "-".join(pieces)


_DICARBOXAMIDE = {8: "carboxamide", 16: "carbothioamide", 34: "carboselenoamide", 52: "carbotelluroamide"}


def _nitro_atoms(mol):
    atoms = set()
    for atom in mol.GetAtoms():
        if is_nitro_nitrogen(mol, atom.GetIdx()):
            atoms |= {atom.GetIdx(), *(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)}
    return atoms


def _carbamoyl(mol, nitrogen, other):
    """(carbonyl carbon, amide nitrogen) of the carbamoyl group C(=O)N or its chalcogen analogue on `nitrogen`, else None."""
    for carbon in nitrogen.GetNeighbors():
        if carbon.GetIdx() == other.GetIdx() or carbon.GetAtomicNum() != 6 or carbon.IsInRing() or carbon.GetDegree() != 3:
            continue
        chalcogen = _carbonyl_chalcogen(carbon)
        amides = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != nitrogen.GetIdx()]
        if chalcogen is not None and len(amides) == 1 and not amides[0].IsInRing():
            return carbon, amides[0]
    return None


def hydrazine_dicarboxamide_name(mol):
    """Hydrazine-1,2-dicarboxamide with substituents on the hydrazine and amide nitrogens (P-66.1.7); two chalcogen
    analogues of one kind give the dicarbothioamide, of two kinds the more senior is the suffix and the other a prefix
    (P-65.1.5.1)."""
    nitro = _nitro_atoms(mol)
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        (a.GetFormalCharge() and a.GetIdx() not in nitro) or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()
    ):
        return None
    graph = adjacency(mol)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() != 7 or b.GetAtomicNum() != 7 or bond.IsInRing() or bond.GetBondTypeAsDouble() != 1.0:
            continue
        if a.GetFormalCharge() or b.GetFormalCharge() or a.GetIdx() in nitro or b.GetIdx() in nitro:
            continue
        found = []
        for nitrogen, other in ((a, b), (b, a)):
            carbamoyl = _carbamoyl(mol, nitrogen, other)
            if carbamoyl is None or nitrogen.GetDegree() > 3 or any(x.GetBondTypeAsDouble() != 1.0 for x in nitrogen.GetBonds()):
                break
            found.append((nitrogen, carbamoyl))
        else:
            try:
                with cited_whole_stereo(mol):
                    name = _hydrazine_dicarboxamide(mol, graph, found, nitro)
            except UnsupportedStructure:
                name = None
            if name is not None:
                return name
    return None


def _hydrazine_dicarboxamide(mol, graph, found, nitro):
    halogens = halogen_substituents(mol)
    aromatic = frozenset(x.GetIdx() for x in mol.GetAtoms() if x.GetIsAromatic())
    core = set()
    for nitrogen, (carbon, amide) in found:
        core |= {nitrogen.GetIdx(), carbon.GetIdx(), amide.GetIdx(), _carbonyl_chalcogen(carbon).GetIdx()}
    elements = [_carbonyl_chalcogen(carbon).GetAtomicNum() for _, (carbon, _) in found]
    senior = min(range(2), key=lambda i: elements[i])
    mixed = elements[0] != elements[1]
    sides = []
    for nitrogen, (carbon, amide) in found:
        branches = []
        covered = set()
        for atom, skip in ((nitrogen, {x.GetIdx() for x, _ in found} | {carbon.GetIdx()}), (amide, {carbon.GetIdx()})):
            names = []
            for root in graph[atom.GetIdx()]:
                if root in skip:
                    continue
                if mol.GetAtomWithIdx(root).GetAtomicNum() != 6:
                    return None
                arm = _arm_atoms(graph, root, atom.GetIdx())
                if atom.GetIdx() in arm or covered & arm or core & arm:
                    return None
                if any(mol.GetAtomWithIdx(i).GetAtomicNum() not in _CARBON_AND_HALOGEN and i not in nitro for i in arm):
                    return None
                covered |= arm
                names.append(name_branch(graph, root, atom.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True))
            branches.append(names)
        sides.append((branches, covered))
    if len(core) + sum(len(c) for _, c in sides) != mol.GetNumAtoms():
        return None
    best = None
    for first, second in (((senior, 1 - senior),) if mixed else ((0, 1), (1, 0))):
        entries = []
        for position, index in ((1, first), (2, second)):
            hydrazine_names, amide_names = sides[index][0]
            entries.extend((position, name, compound) for name, compound in hydrazine_names)
            if mixed and index != senior:
                junior_nitrogen, (junior_carbon, _) = found[index]
                name, compound = name_branch(graph, junior_carbon.GetIdx(), junior_nitrogen.GetIdx(), halogens, aromatic, mol=mol)
                entries.append((position, name, compound))
            else:
                entries.extend((("N" if mixed else f"N{position}"), name, compound) for name, compound in amide_names)
        key = (
            tuple(sorted(loc for loc, _, _ in entries if isinstance(loc, int))),
            tuple(sorted(int(loc[1:] or 0) for loc, _, _ in entries if isinstance(loc, str))),
            tuple(sorted((alpha_sort_key(name), str(loc)) for loc, name, _ in entries)),
        )
        if best is None or key < best[0]:
            best = (key, entries)
    prefixes = {}
    for locant, name, compound in best[1]:
        prefixes.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    head = format_substituent_prefixes(prefixes) if prefixes else ""
    if mixed:
        return f"{head}hydrazine-1-{_DICARBOXAMIDE[elements[senior]]}"
    return f"{head}hydrazine-1,2-di{_DICARBOXAMIDE[elements[0]]}"


_HYDRAZIDE_INFIX = {8: "carbo", 16: "carbothio", 34: "carboseleno", 52: "carbotelluro"}


def diazenecarbohydrazide_name(mol):
    """Diazenecarbohydrazide and its chalcogen analogues with N, N' and 2 substituents (P-68.3.1.3.4)."""
    if not _plain(mol):
        return None
    graph = adjacency(mol)
    for carbon in mol.GetAtoms():
        chalcogen = _carbonyl_chalcogen(carbon) if carbon.GetAtomicNum() == 6 and carbon.GetDegree() == 3 else None
        if chalcogen is None or carbon.IsInRing():
            continue
        nitrogens = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7]
        if len(nitrogens) != 2:
            continue
        for azo, acyl in (nitrogens, nitrogens[::-1]):
            far = [n for n in azo.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != carbon.GetIdx()]
            beta = [n for n in acyl.GetNeighbors() if n.GetAtomicNum() == 7]
            if len(far) != 1 or len(beta) != 1 or mol.GetBondBetweenAtoms(azo.GetIdx(), far[0].GetIdx()).GetBondTypeAsDouble() != 2.0:
                continue
            if azo.IsInRing() or acyl.IsInRing() or beta[0].IsInRing() or far[0].IsInRing() or far[0].GetDegree() > 2:
                continue
            if any(b.GetBondTypeAsDouble() != 1.0 for atom in (acyl, beta[0]) for b in atom.GetBonds()):
                continue
            entries, used = [], {carbon.GetIdx(), chalcogen.GetIdx(), azo.GetIdx(), acyl.GetIdx(), far[0].GetIdx(), beta[0].GetIdx()}
            valid = True
            for atom, skip, locant in (
                (far[0], {azo.GetIdx()}, 2),
                (acyl, {carbon.GetIdx(), beta[0].GetIdx()}, "N"),
                (beta[0], {acyl.GetIdx()}, "N'"),
            ):
                found = _hydrocarbon_arms(mol, graph, atom.GetIdx(), skip)
                if found is None:
                    valid = False
                    break
                entries += [(locant, name, compound) for name, compound in found[0]]
                used |= found[1]
            if not valid or len(used) != mol.GetNumAtoms():
                continue
            grouped = {}
            for locant, name, compound in entries:
                grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
            head = format_substituent_prefixes(grouped) if grouped else ""
            return f"{head}diazene{_HYDRAZIDE_INFIX[chalcogen.GetAtomicNum()]}hydrazide"
    return None
