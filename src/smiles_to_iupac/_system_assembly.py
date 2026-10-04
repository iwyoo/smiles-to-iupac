"""Ring assemblies of two identical fused ring systems joined by one single bond (P-28.2.1, P-29.3.5):
[1,1'-binaphthalene]-2,2'-diol, 2'H-1,2'-biindole-type names, [1,2'-binaphthalen]-4'-yl."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._multiplicative import _bare_key
from ._multiplicative_ring import _SUFFIX_WORDS
from ._common import multiplied_word
from ._multiplicative_text import PrimedLocant
from ._ring_diyl_numbering import SUFFIX_ATOMS, is_hydro_fusion_system, system_numberings
from ._substituents import format_substituent_prefixes, name_branch


def _order(locant):
    return float(locant[1]), locant[0]


def _cite(locant):
    return f"{locant[1]}{chr(39) * locant[0]}"


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
        return numbering.stem
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


def system_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, free=None, ring_suffix=None):
    """(count, ((-count,), name, parts)) for the parent form, (name, True) for the substituent form `free`
    = (root, coming_from); None when `mol` is not two identical fused systems joined by one bond."""
    systems = _systems(mol)
    if len(systems) != 2 or all(len(rings) == 1 for rings, _ in systems):
        return None
    (rings_a, atoms_a), (rings_b, atoms_b) = systems
    if _bare_key(mol, atoms_a) != _bare_key(mol, atoms_b):
        return None
    if free is not None:
        if free[0] not in atoms_a | atoms_b:
            return None
    elif principal is not None and not occurrences:
        return None
    joins = [(a, b) for a in atoms_a for b in atoms_b if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1 or mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() not in (1.0, 2.0):
        return None
    ylidene = mol.GetBondBetweenAtoms(*joins[0]).GetBondTypeAsDouble() == 2.0
    if not all(
        any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms) or is_hydro_fusion_system(mol, atoms) for _, atoms in systems
    ):
        return None
    from ._polyfunctional import _require_mancude_system

    for _, atoms in systems:
        _require_mancude_system(mol, atoms)
    join = joins[0]
    token = SUFFIX_ATOMS.set(frozenset(join))
    try:
        numbered = [system_numberings(mol, graph, rings, atoms) for rings, atoms in systems]
    finally:
        SUFFIX_ATOMS.reset(token)
    owned = set().union(*(o[2] for o in occurrences)) if occurrences else set()
    marked = [o[1] for o in occurrences] if free is None else [free[0]]
    atoms_all = atoms_a | atoms_b
    roots = [
        (r, n.GetIdx())
        for r in atoms_all
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in atoms_all
        and n.GetIdx() not in owned
        and not (free is not None and r == free[0] and n.GetIdx() == free[1])
    ]
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    best = None
    for unprimed in (0, 1):
        for first in numbered[unprimed]:
            stem_first = _stem(first)
            for second in numbered[1 - unprimed]:
                stem_second = _stem(second)
                if stem_first is None or stem_first != stem_second or _hydro(first) != _hydro(second):
                    continue
                locants = {a: (0, p) for a, p in first.position_of.items()}
                locants.update({a: (1, p) for a, p in second.position_of.items()})
                if any(a not in locants for a in marked + [r for r, _, _ in entries] + list(join)):
                    continue
                ih = sorted(
                    [(0, p) for p in _lit_hydrogens(mol, first)] + [(1, p) for p in _lit_hydrogens(mol, second)], key=_order
                )
                members = (atoms_a, atoms_b)
                junction = (locants[join[0 if join[0] in members[unprimed] else 1]][1], locants[join[1 if join[0] in members[unprimed] else 0]][1])
                key = (
                    junction,
                    tuple(_order(x) for x in ih),
                    tuple(sorted(_order(locants[a]) for a in marked)),
                    tuple(sorted(_order(locants[r]) for r, _, _ in entries)),
                    tuple(loc for loc, _ in sorted(((_order(locants[r]), name) for r, name, _ in entries), key=lambda e: (e[1], e[0]))),
                )
                if best is None or key < best[0]:
                    best = (key, locants, stem_first, ih, unprimed, _hydro(first), _fully_hydro(first))
    if best is None:
        raise UnsupportedStructure("this ring assembly of fused systems needs hydro or added hydrogen, not supported yet")
    _, locants, stem, ih, unprimed, hydro, fully_hydro = best
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(_cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(re.match(r"\d+", text).group()), text.count(chr(39))))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    ih_text = (",".join(f"{p}{chr(39) * n}H" for n, p in ih) + "-") if ih else ""
    if hydro:
        cited = sorted([(0, p) for p in hydro] + [(1, p) for p in hydro], key=_order)
        hydro_word = multiplied_word(len(cited), "hydro")
        hydro_text = hydro_word if fully_hydro else f"{','.join(f'{p}{chr(39) * n}' for n, p in cited)}-{hydro_word}"
        ih_text = hydro_text + "-" + ih_text
    unprimed_join = join[0] if join[0] in (atoms_a, atoms_b)[unprimed] else join[1]
    primed_join = join[1] if unprimed_join == join[0] else join[0]
    spots = f"{locants[unprimed_join][1]},{locants[primed_join][1]}'"
    def base(elide):
        text = stem[:-1] if (elide or ylidene) and stem.endswith("e") else stem
        if ylidene:
            text += "ylidene"
        inner = f"bi({text})" if re.search(r"[\d\[-]|cyclo", text) else f"bi{text}"
        return f"{spots}-{inner}"

    if free is not None:
        core = f"{ih_text}[{base(True)}]-{_cite(locants[free[0]])}-yl"
        return (f"{prefix}-{core}" if prefix else core), True
    count = len(occurrences)
    if principal is None:
        core = f"{ih_text}{base(False)}"
    else:
        from ._polyfunctional import _RING_SUFFIX

        word = multiplied_word(count, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        spots_text = ",".join(_cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: _order(locants[o[1]])))
        core = f"{ih_text}[{base(word[0] in 'aeiouy')}]-{spots_text}-{word}"
    name = f"{prefix}-{core}" if prefix else core
    return count, ((-count,), name, (None, None, None, 0, {a: PrimedLocant(*loc) for a, loc in locants.items()}, True))
