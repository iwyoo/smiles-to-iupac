"""Ring assemblies of two identical fused ring systems joined by one single bond (P-28.2.1, P-29.3.5):
[1,1'-binaphthalene]-2,2'-diol, 2'H-1,2'-biindole-type names, [1,2'-binaphthalen]-4'-yl."""

import re


from ._common import UnsupportedStructure, assembly_join
from ._multiplicative import _bare_key
from ._multiplicative_ring import _SUFFIX_WORDS
from ._common import multiplied_word
from ._multiplicative_text import PrimedLocant
from ._ring_diyl_numbering import SUFFIX_ATOMS, is_hydro_fusion_system, system_numberings
from ._substituents import format_substituent_prefixes, name_branch


def _order(locant):
    return float(locant[1]), locant[0]


def _primed(position, primes):
    """P-16.9.4: primes follow the number, not the letter of a fusion locant (4'a, 3'a)."""
    return re.sub(r"^(\d+)", lambda m: m.group(1) + chr(39) * primes, str(position))


def _cite(locant):
    return _primed(locant[1], locant[0])


def _hetero_mancude_ring(mol, atoms):
    """A monocycle with a heteroatom and at least one ring double bond (or aromatic atoms): pyridine, azepine, pyran."""
    ring_bonds = [b for b in mol.GetBonds() if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms]
    return (
        len(atoms) == len(ring_bonds)
        and any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in atoms)
        and any(b.GetBondTypeAsDouble() == 2.0 or b.GetIsAromatic() for b in ring_bonds)
        and not any(b.GetBondTypeAsDouble() == 3.0 for b in ring_bonds)
    )


def _monocycle_numberings(mol, graph, atoms):
    from ._common import ring_cycle
    from ._ring_diyl_numbering import monocycle_numberings

    return monocycle_numberings(mol, ring_cycle(graph, list(atoms)), set(), 1)


def _systems(mol):
    from ._diester_ring_diyl import _system_of

    found = []
    for ring in mol.GetRingInfo().AtomRings():
        rings, atoms = _system_of(mol, ring[0])
        if not any(atoms == other[1] for other in found):
            found.append((rings, atoms))
    return found


_HYDRO_HEAD = re.compile(r"^(?:[\d,a-z]+-)?(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca)?hydro-?")


def _stem(numbering):
    if numbering.added:
        return None
    if numbering.stem:
        return re.sub(r"^\d+[a-z]*H(?:,\d+[a-z]*H)*-", "", numbering.stem)
    text = numbering.text((), 0, frozenset(), "")
    text = _HYDRO_HEAD.sub("", text)
    return re.sub(r"^\d+[a-z]*H(?:,\d+[a-z]*H)*-", "", text)


def _hydro(numbering):
    if numbering.hydro:
        return tuple(numbering.hydro)
    return tuple(numbering.unsat_key[1]) if len(numbering.unsat_key) == 2 else ()


def _fully_hydro(numbering):
    return numbering.fully_hydro or bool(re.match(r"^(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca)hydro", numbering.text((), 0, frozenset(), "")))


def _lit_hydrogens(mol, numbering):
    """Indicated-hydrogen positions that still carry a hydrogen: a junction atom whose only hydrogen was
    replaced by the bond between the components needs none (P-28.2.3)."""
    atom_of = {p: a for a, p in numbering.position_of.items()}
    return [p for p in numbering.ih_positions if mol.GetAtomWithIdx(atom_of[p]).GetTotalNumHs() > 0]


def _lit(mol, numbering, positions, monocycles):
    """Indicated-hydrogen positions of a component. A fused system keeps those its numbering names, minus junction atoms
    that lost their hydrogen (P-28.2.3). A monocycle takes every saturated ring carbon or NH that still has a hydrogen or
    carries =O (a pseudoketone, P-64.3.1): the junction atoms and the positions of double bonds need none."""
    if not monocycles:
        return _lit_hydrogens(mol, numbering) + list(positions)
    found = []
    for atom_idx, position in numbering.position_of.items():
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() not in (6, 7):
            continue
        oxo = False
        double = False
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0 or bond.GetIsAromatic():
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() in (8, 16, 34, 52) and other.GetDegree() == 1:
                    oxo = True
                else:
                    double = True
        if double:
            continue
        if atom.GetTotalNumHs() > 0 or (oxo and atom.GetAtomicNum() == 6):
            found.append(position)
    return sorted(found)


def _skeleton_key(mol, atoms):
    """Canonical SMILES of the ring system alone, every bond single and every atom saturated, so that the same parent
    hydride in different hydrogenation states compares equal."""
    from rdkit import Chem

    editable = Chem.RWMol(mol)
    for bond in editable.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in editable.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    for index in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(index)
    return Chem.MolToSmiles(editable.GetMol())


def _split_added(numbering, frees):
    """(hydro positions without the ylidene pair, [added-hydrogen position]) for a component whose ylidene
    atom and one sp3 atom were counted as 'dihydro' (P-29.3.4.1)."""
    hydro = tuple(_hydro(numbering))
    for atom, _ in frees:
        position = numbering.position_of.get(atom)
        if position is not None and position in hydro and len(hydro) == 2:
            return (), [next(p for p in hydro if p != position)]
    return hydro, []


def system_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, free=None, ring_suffix=None, within=None):
    """(count, ((-count,), name, parts)) for the parent form, (name, True) for the substituent form `free`
    = (root, coming_from); None when `mol` is not two identical fused systems joined by one bond."""
    systems = [s for s in _systems(mol) if within is None or set(s[1]) <= within]
    if len(systems) != 2 or all(len(rings) == 1 for rings, _ in systems) and not all(
        _hetero_mancude_ring(mol, atoms) for _, atoms in systems
    ):
        return None
    frees = [] if free is None else [free] if isinstance(free[0], int) else list(free)
    (rings_a, atoms_a), (rings_b, atoms_b) = systems
    monocycles = all(len(rings) == 1 for rings, _ in systems)
    if _bare_key(mol, atoms_a) != _bare_key(mol, atoms_b) and not (monocycles and _skeleton_key(mol, atoms_a) == _skeleton_key(mol, atoms_b)) and not (
        frees and all(mol.GetBondBetweenAtoms(*f).GetBondTypeAsDouble() == 2.0 for f in frees)
        and _skeleton_key(mol, atoms_a) == _skeleton_key(mol, atoms_b)
    ):
        return None
    if frees:
        if any(f[0] not in atoms_a | atoms_b for f in frees):
            return None
        orders = {mol.GetBondBetweenAtoms(*f).GetBondTypeAsDouble() for f in frees}
        if orders - {1.0, 2.0} or len(orders) != 1:
            return None
        free_ylidene = orders == {2.0}
    elif principal is not None and not occurrences:
        return None
    joins = [(a, b) for a in atoms_a for b in atoms_b if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1 or mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() not in (1.0, 2.0):
        return None
    ylidene = mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() == 2.0
    if ylidene and frees:
        return None
    if not all(
        any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms)
        or is_hydro_fusion_system(mol, atoms)
        or _hetero_mancude_ring(mol, atoms)
        for _, atoms in systems
    ):
        return None
    from ._polyfunctional import _require_mancude_system

    for _, atoms in systems:
        _require_mancude_system(mol, atoms)
    join = joins[0]
    token = SUFFIX_ATOMS.set(frozenset(join) | frozenset(o[1] for o in occurrences))
    try:
        numbered = [
            _monocycle_numberings(mol, graph, atoms) if len(rings) == 1 else system_numberings(mol, graph, rings, atoms)
            for rings, atoms in systems
        ]
    finally:
        SUFFIX_ATOMS.reset(token)
    owned = set().union(*(o[2] for o in occurrences)) if occurrences else set()
    marked = [o[1] for o in occurrences] if not frees else [f[0] for f in frees]
    atoms_all = atoms_a | atoms_b
    roots = [
        (r, n.GetIdx())
        for r in atoms_all
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in atoms_all
        and n.GetIdx() not in owned
        and (r, n.GetIdx()) not in frees
    ]
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    best = None
    for unprimed in (0, 1):
        for first in numbered[unprimed]:
            stem_first = _stem(first)
            for second in numbered[1 - unprimed]:
                stem_second = _stem(second)
                if monocycles:
                    hydro_first = hydro_second = ()
                    added_first = added_second = []
                else:
                    hydro_first, added_first = _split_added(first, frees) if frees and free_ylidene else (_hydro(first), [])
                    hydro_second, added_second = _split_added(second, frees) if frees and free_ylidene else (_hydro(second), [])
                if stem_first is None or stem_first != stem_second or hydro_first != hydro_second:
                    continue
                locants = {a: (0, p) for a, p in first.position_of.items()}
                locants.update({a: (1, p) for a, p in second.position_of.items()})
                if any(a not in locants for a in marked + [r for r, _, _ in entries] + list(join)):
                    continue
                ih = sorted(
                    [(0, p) for p in _lit(mol, first, added_first, monocycles)]
                    + [(1, p) for p in _lit(mol, second, added_second, monocycles)],
                    key=_order,
                )
                members = (atoms_a, atoms_b)
                junction = (locants[join[0 if join[0] in members[unprimed] else 1]][1], locants[join[1 if join[0] in members[unprimed] else 0]][1])
                ih_key = tuple(_order(x) for x in ih)
                marked_key = tuple(sorted(_order(locants[a]) for a in marked))
                key = (
                    junction,
                    (marked_key, ih_key) if frees and free_ylidene else (ih_key, marked_key),
                    tuple(sorted(_order(locants[r]) for r, _, _ in entries)),
                    tuple(loc for loc, _ in sorted(((_order(locants[r]), name) for r, name, _ in entries), key=lambda e: (e[1], e[0]))),
                )
                if best is None or key < best[0]:
                    best = (key, locants, stem_first, ih, unprimed, hydro_first, _fully_hydro(first))
    if best is None:
        raise UnsupportedStructure("this ring assembly of fused systems needs hydro or added hydrogen, not supported yet")
    _, locants, stem, ih, unprimed, hydro, fully_hydro = best
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(_cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(re.match(r"\d+", text).group()), text.count(chr(39))))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    ih_text = (",".join(f"{_primed(p, n)}H" for n, p in ih) + "-") if ih else ""
    hydro_prefix = ""
    if hydro:
        cited = sorted([(0, p) for p in hydro] + [(1, p) for p in hydro], key=_order)
        hydro_word = multiplied_word(len(cited), "hydro")
        hydro_text = hydro_word if fully_hydro else f"{','.join(_primed(p, n) for n, p in cited)}-{hydro_word}"
        ih_text = hydro_text + "-" + ih_text
        hydro_prefix = hydro_text + "-"
    unprimed_join = join[0] if join[0] in (atoms_a, atoms_b)[unprimed] else join[1]
    primed_join = join[1] if unprimed_join == join[0] else join[0]
    spots = f"{locants[unprimed_join][1]},{_primed(locants[primed_join][1], 1)}"
    def base(elide):
        text = stem[:-1] if (elide or ylidene) and stem.endswith("e") else stem
        if ylidene:
            text += "ylidene"
        inner = f"bi({text})" if re.search(r"[\d\[-]|cyclo", text) else f"bi{text}"
        return f"{spots}-{inner}"

    if frees:
        free_spots = ",".join(_cite(c) for c in sorted((locants[f[0]] for f in frees), key=_order))
        if free_ylidene:
            added = ",".join(f"{_primed(p, n)}H" for n, p in ih)
            tail = f"{free_spots}({added})" if added else free_spots
            core = f"{hydro_prefix}[{base(multiplied_word(len(frees), 'ylidene')[0] in 'aeiouy')}]-{tail}-{multiplied_word(len(frees), 'ylidene')}"
        else:
            core = f"{ih_text}[{base(multiplied_word(len(frees), 'yl')[0] in 'aeiouy')}]-{free_spots}-{multiplied_word(len(frees), 'yl')}"
        return (assembly_join(prefix, core)), True
    count = len(occurrences)
    if principal is None:
        core = f"{ih_text}{base(False)}"
    else:
        from ._polyfunctional import _RING_SUFFIX

        word = multiplied_word(count, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        spots_text = ",".join(_cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: _order(locants[o[1]])))
        core = f"{ih_text}[{base(word[0] in 'aeiouy')}]-{spots_text}-{word}"
    name = assembly_join(prefix, core)
    return count, ((-count,), name, (None, None, None, 0, {a: PrimedLocant(*loc) for a, loc in locants.items()}, True, len(name) - len(core)))
