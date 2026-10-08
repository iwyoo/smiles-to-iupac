"""Amides of cyanic acid and of di- and polycarbonic acids (P-66.1.6.1.4, P-66.1.6.2, P-66.1.6.3):
(propan-2-yl)cyanamide, 2-imidodicarbonic diamide, N1-methyl-2-imido-1-thiodicarbonic diamide, 1,2,3-trithiodicarbonic
diamide."""

from rdkit import Chem

from ._common import adjacency, halogen_substituents, is_nitro_nitrogen
from ._numerals import numerical_term
from ._substituents import format_substituent_prefixes, name_branch

_CHALCOGEN = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}
_CARBON_AND_HALOGEN = {6, 9, 17, 35, 53}
_MAXIMUM_CARBONIC_UNITS = 4


def _plain(mol):
    return (
        len(Chem.GetMolFrags(mol)) == 1
        and not any(a.GetFormalCharge() or a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        and not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())
        and not any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
    )


def _arm_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n not in seen and n != blocked:
                seen.add(n)
                stack.append(n)
    return seen


def _hydrocarbon_arms(mol, graph, nitrogen, skip):
    """(name, compound) of each carbon group on `nitrogen` outside `skip` and the atoms they cover, or None when a
    group is not a plain hydrocarbon (or halogenated) substituent."""
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    names, covered = [], set()
    for root in graph[nitrogen]:
        if root in skip:
            continue
        if mol.GetAtomWithIdx(root).GetAtomicNum() != 6:
            return None
        arm = _arm_atoms(graph, root, nitrogen)
        if nitrogen in arm or covered & arm or any(mol.GetAtomWithIdx(i).GetAtomicNum() not in _CARBON_AND_HALOGEN for i in arm):
            return None
        covered |= arm
        names.append(name_branch(graph, root, nitrogen, halogens, aromatic, mol=mol, unsaturated=True))
    return names, covered


def cyanamide_name(mol):
    if not _plain(mol):
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
        if amino.IsInRing() or amino.GetIsAromatic() or any(b.GetBondTypeAsDouble() != 1.0 for b in amino.GetBonds()):
            continue
        betas = [n for n in amino.GetNeighbors() if n.GetAtomicNum() == 7]
        if len(betas) > 1:
            continue
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


def _carbonyl_chalcogen(atom):
    for b in atom.GetBonds():
        other = b.GetOtherAtom(atom)
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() in _CHALCOGEN and other.GetDegree() == 1:
            return other
    return None


def _linker(mol, carbon, neighbor):
    """(linker atoms, next carbon) when `neighbor` bridges `carbon` to another carbonyl carbon, else None."""
    if neighbor.IsInRing() or neighbor.GetIsAromatic():
        return None
    z = neighbor.GetAtomicNum()
    if z == 7:
        if neighbor.GetDegree() != 2 or neighbor.GetTotalNumHs() != 1:
            return None
        forward = [n for n in neighbor.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
        return ([neighbor], forward[0])
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
    for nitrogen in ends:
        if any(b.GetBondTypeAsDouble() != 1.0 for b in nitrogen.GetBonds()):
            return None
        betas = [n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 7]
        hydrazide.append(bool(betas))
        if len(betas) > 1:
            return None
        found = _hydrocarbon_arms(mol, graph, nitrogen.GetIdx(), carbon_set | {b.GetIdx() for b in betas})
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
        for rank, index in enumerate(order):
            if rank:
                position += 1
                replacements.extend(_linker_prefix(linkers[min(index, order[rank - 1])], position))
                position += 1
            positions.append(position)
            symbol = _CHALCOGEN[chalcogens[index].GetAtomicNum()]
            if symbol:
                replacements.append((position, symbol, True))
        substituents = []
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


def _nitro_atoms(mol):
    atoms = set()
    for atom in mol.GetAtoms():
        if is_nitro_nitrogen(mol, atom.GetIdx()):
            atoms |= {atom.GetIdx(), *(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)}
    return atoms


def _carbamoyl(mol, nitrogen, other):
    """(carbonyl carbon, amide nitrogen) of the carbamoyl group C(=O)N on `nitrogen`, else None."""
    for carbon in nitrogen.GetNeighbors():
        if carbon.GetIdx() == other.GetIdx() or carbon.GetAtomicNum() != 6 or carbon.IsInRing() or carbon.GetDegree() != 3:
            continue
        oxygen = _carbonyl_chalcogen(carbon)
        amides = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != nitrogen.GetIdx()]
        if oxygen is not None and oxygen.GetAtomicNum() == 8 and len(amides) == 1 and not amides[0].IsInRing():
            return carbon, amides[0]
    return None


def hydrazine_dicarboxamide_name(mol):
    """Hydrazine-1,2-dicarboxamide with substituents on the hydrazine and amide nitrogens (P-66.1.7)."""
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
            name = _hydrazine_dicarboxamide(mol, graph, found, nitro)
            if name is not None:
                return name
    return None


def _hydrazine_dicarboxamide(mol, graph, found, nitro):
    halogens = halogen_substituents(mol)
    aromatic = frozenset(x.GetIdx() for x in mol.GetAtoms() if x.GetIsAromatic())
    core = set()
    for nitrogen, (carbon, amide) in found:
        core |= {nitrogen.GetIdx(), carbon.GetIdx(), amide.GetIdx(), _carbonyl_chalcogen(carbon).GetIdx()}
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
    for first, second in ((0, 1), (1, 0)):
        entries = []
        for position, index in ((1, first), (2, second)):
            hydrazine_names, amide_names = sides[index][0]
            entries.extend((position, name, compound) for name, compound in hydrazine_names)
            entries.extend((f"N{position}", name, compound) for name, compound in amide_names)
        key = (
            tuple(sorted(loc for loc, _, _ in entries if isinstance(loc, int))),
            tuple(sorted(int(loc[1:]) for loc, _, _ in entries if isinstance(loc, str))),
            tuple(sorted((name, str(loc)) for loc, name, _ in entries)),
        )
        if best is None or key < best[0]:
            best = (key, entries)
    prefixes = {}
    for locant, name, compound in best[1]:
        prefixes.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    head = format_substituent_prefixes(prefixes) if prefixes else ""
    return f"{head}hydrazine-1,2-dicarboxamide"
