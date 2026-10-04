"""Ring parent structures of multiplicative names (P-15.3.1.3, P-15.3.2.4):
monocycle numbering with the junction to the linker ranked right after the
principal characteristic groups (P-15.3.2.4.2), ring-diyl linker components
(P-29.3), and bare fused/bridged/spiro units via a probe name.
"""

import re
from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    alpha_sort_key,
    heteroaromatic_monocycle_name,
    multiplied_word,
    ring_cycle,
    unsaturation_suffix,
)
from ._multiplicative_groups import SUFFIX_RANKS
from ._multiplicative_prefix import SIMPLE_PREFIXES, prefix_name, probe_name, subtree
from ._numerals import alkane_name, alkyl_name
from ._substituents import _ring_base_name, format_substituent_prefixes

class _SuffixWords(dict):
    def __missing__(self, key):
        if isinstance(key, str) and key.startswith("acid:"):
            from ._acid_lexicon import carbo_suffix, spec_from_key

            return carbo_suffix(spec_from_key(key), 1)
        raise KeyError(key)

    def __contains__(self, key):
        return dict.__contains__(self, key) or (isinstance(key, str) and key.startswith("acid:"))


_SUFFIX_WORDS = _SuffixWords({
    "carboxylic_acid": "carboxylic acid",
    "sulfonic_acid": "sulfonic acid",
    "amide": "carboxamide",
    "sulfonamide": "sulfonamide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
    "ketone": "one",
    "alcohol": "ol",
    "thiol": "thiol",
    "amine": "amine",
})
_RETAINED_BENZENE = {
    "acid:C:O:O!": "benzoate",
    "carboxylic_acid": "benzoic acid",
    "amide": "benzamide",
    "nitrile": "benzonitrile",
    "aldehyde": "benzaldehyde",
    "alcohol": "phenol",
    "amine": "aniline",
}
_PRIMARY_NITROGEN = {"amide", "sulfonamide", "amine"}
_HETERO_PARENTS = {"pyridine": "pyridine", "furan": "furan", "thiophene": "thiophene", "pyrrole": "1H-pyrrole"}
_VALENCE_WORDS = {2: "diyl", 3: "triyl", 4: "tetrayl"}


@dataclass
class RingSpec:
    cycle: list
    kind: str
    hetero: object
    multiple: tuple = ()

    @property
    def parent(self):
        if self.kind == "benzene":
            return "benzene"
        if self.kind in ("cycloalkane", "cycloalkene"):
            return "cyclo" + alkane_name(len(self.cycle))
        return _HETERO_PARENTS[self.kind]


@dataclass
class UnitText:
    text: str
    junction_locant: object
    substituted: bool
    has_locants: bool


def ene_spec(mol, ring_atoms):
    """A non-aromatic carbocycle with ring double or triple bonds, for the
    multiplicative unit, substituent and linker namers (cyclohex-2-en-1-yl)."""
    ring_atoms = list(ring_atoms)
    graph = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in ring_atoms] for a in ring_atoms}
    cycle = ring_cycle(graph, ring_atoms)
    atoms = [mol.GetAtomWithIdx(a) for a in cycle]
    if any(a.GetFormalCharge() or a.GetIsotope() or a.GetIsAromatic() or a.GetAtomicNum() != 6 for a in atoms):
        return None
    bonds = []
    for i in range(len(cycle)):
        a, b = cycle[i], cycle[(i + 1) % len(cycle)]
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order not in (1.0, 2.0, 3.0):
            return None
        if order > 1.0:
            bonds.append((a, b, order))
    if not bonds or any(order == 3.0 for _, _, order in bonds) and len(cycle) < 8:
        return None
    return RingSpec(cycle, "cycloalkene", None, tuple(bonds))


def spec_of(mol, ring_atoms):
    return monocycle_spec(mol, ring_atoms) or ene_spec(mol, ring_atoms)


def multiple_locants(spec, locants):
    ene, yne = [], []
    n = len(spec.cycle)
    for a, b, order in spec.multiple:
        low, high = sorted((locants[a], locants[b]))
        locant = n if (low, high) == (1, n) else low
        (ene if order == 2.0 else yne).append(locant)
    return sorted(ene), sorted(yne)


def parent_text(spec, locants):
    if spec.kind != "cycloalkene":
        return spec.parent
    ene, yne = multiple_locants(spec, locants)
    count = len(ene) + len(yne)
    stem = "cyclo" + alkane_name(len(spec.cycle))[:-3]
    body, needs_a = unsaturation_suffix(ene, yne)
    return f"{stem}{'a' if needs_a else ''}-{body}" if count else spec.parent


def monocycle_spec(mol, ring_atoms):
    ring_atoms = list(ring_atoms)
    graph = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in ring_atoms] for a in ring_atoms}
    cycle = ring_cycle(graph, ring_atoms)
    atoms = [mol.GetAtomWithIdx(a) for a in cycle]
    if any(a.GetFormalCharge() or a.GetIsotope() for a in atoms):
        return None
    if all(a.GetAtomicNum() == 6 for a in atoms):
        if len(cycle) == 6 and all(a.GetIsAromatic() for a in atoms):
            return RingSpec(cycle, "benzene", None)
        saturated = not any(a.GetIsAromatic() for a in atoms) and all(
            mol.GetBondBetweenAtoms(cycle[i], cycle[(i + 1) % len(cycle)]).GetBondTypeAsDouble() == 1
            for i in range(len(cycle))
        )
        return RingSpec(cycle, "cycloalkane", None) if saturated else None
    name = heteroaromatic_monocycle_name(mol, cycle)
    if name is None:
        return None
    (hetero,) = [a for a in cycle if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    return RingSpec(cycle, name, hetero)


def numberings(spec):
    n = len(spec.cycle)
    result = []
    for start in range(n):
        for step in (1, -1):
            locants = {spec.cycle[(start + step * k) % n]: k + 1 for k in range(n)}
            if spec.hetero is None or locants[spec.hetero] == 1:
                result.append(locants)
    return result


def principal_rank_of(groups, atoms):
    ranks = [g.rank for g in groups if g.anchor in atoms and g.rank <= SUFFIX_RANKS["amine"]]
    return min(ranks) if ranks else None


def _is_primary(mol, group):
    if group.name not in _PRIMARY_NITROGEN:
        return True
    nitrogen = [a for a in group.atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    return len(nitrogen) == 1 and mol.GetAtomWithIdx(nitrogen[0]).GetTotalNumHs() == 2


def _join(prefix_text, core):
    return f"{prefix_text}-{core}" if prefix_text and core[0].isdigit() else prefix_text + core


def _citation_key(entries):
    return tuple(loc for loc, _ in sorted(entries, key=lambda e: (alpha_sort_key(e[1]), e[0])))


def _prefix_text(entries, locants):
    grouped = {}
    for r, name, compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": compound})
        info["locants"].append(locants[r])
    return format_substituent_prefixes(grouped) if grouped else ""


def _ring_roots(mol, spec, skip=frozenset()):
    ring = set(spec.cycle)
    return [
        (r, n.GetIdx())
        for r in spec.cycle
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring and (r, n.GetIdx()) not in skip
    ]


def _group_at(mol, groups, ring_atom, root):
    atoms = subtree(mol, root, ring_atom)
    for g in groups:
        if g.name == "ketone" and g.ring_atom == g.anchor:
            anchor_root = next(a for a in g.atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 8)
            attached = g.anchor
        else:
            anchor_root, attached = g.anchor, g.ring_atom
        if (attached, anchor_root) == (ring_atom, root) and g.atoms - {ring_atom} == atoms and _is_primary(mol, g):
            return g
    return None


def ring_substituent_name(mol, ring_atoms, attach_atom, from_atom, groups, suffix_group, name_function=None):
    """(name, is_compound) of a monocyclic substituent group ('phenyl',
    '4-chlorophenyl', 'cyclohexyl', 'pyridin-2-yl') attached through
    `attach_atom`, or None when the ring isn't a supported monocycle."""
    if sum(1 for r in mol.GetRingInfo().AtomRings() if set(r) & set(ring_atoms)) != 1:
        return None
    spec = spec_of(mol, ring_atoms)
    if spec is None:
        return None
    entries = _prefix_entries(mol, _ring_roots(mol, spec, {(attach_atom, from_atom)}), groups, suffix_group, name_function)
    best = None
    for locants in numberings(spec):
        key = (
            locants[attach_atom],
            multiple_locants(spec, locants),
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    locants = best[1]
    if spec.kind == "benzene":
        core = "phenyl"
    elif spec.kind == "cycloalkane":
        core = "cyclo" + alkyl_name(len(spec.cycle))
    elif spec.kind == "cycloalkene":
        ene, yne = multiple_locants(spec, locants)
        core = _ring_base_name(len(spec.cycle), ene, yne, 1.0)
    else:
        parent = spec.parent
        core = f"{parent[:-1] if parent.endswith('e') else parent}-{locants[attach_atom]}-yl"
    prefix_text = _prefix_text(entries, locants)
    return _join(prefix_text, core), bool(prefix_text) or spec.hetero is not None or spec.kind == "cycloalkene"


def _prefix_entries(mol, roots, groups, suffix_group, name_function):
    entries = []
    for r, n in roots:
        group = _group_at(mol, groups, r, n)
        if group is not None and group.name in SIMPLE_PREFIXES:
            entries.append((r, SIMPLE_PREFIXES[group.name], False))
        else:
            name, compound = prefix_name(mol, n, r, suffix_group, name_function, groups)
            entries.append((r, name, compound))
    return entries


def _suffix_text(parent, suffix_name, locants, spec):
    count = len(locants)
    if count == 1 and spec.kind == "benzene" and suffix_name in _RETAINED_BENZENE:
        return _RETAINED_BENZENE[suffix_name], False
    suffix = multiplied_word(count, _SUFFIX_WORDS[suffix_name])
    stem = parent[:-1] if (parent.endswith("e") and suffix[0] in "aeiouy") else parent
    return f"{stem}-{','.join(str(x) for x in sorted(locants))}-{suffix}", True


def name_monocyclic_unit(mol, ring_atoms, junction, linker_atom, groups, unit_atoms, name_function=None):
    """UnitText for a monocyclic multiplied parent attached to its linker
    through ring atom `junction`; None when the unit's principal class can't
    be expressed as a suffix on this ring."""
    spec = spec_of(mol, ring_atoms)
    if spec is None:
        raise UnsupportedStructure("this monocyclic ring is not supported as a multiplied parent structure yet")
    principal = principal_rank_of(groups, unit_atoms)
    suffix_name = None
    suffix_atoms = []
    prefix_roots = []
    for r, n in _ring_roots(mol, spec, {(junction, linker_atom)}):
        group = _group_at(mol, groups, r, n)
        if group is not None and group.rank == principal and group.name in _SUFFIX_WORDS:
            if suffix_name not in (None, group.name):
                return None
            suffix_name = group.name
            suffix_atoms.append(r)
        else:
            prefix_roots.append((r, n))
    if principal is not None:
        principal_count = sum(1 for g in groups if g.anchor in unit_atoms and g.rank == principal)
        if suffix_name is None or principal_count != len(suffix_atoms):
            return None
    entries = _prefix_entries(mol, prefix_roots, groups, suffix_name, name_function)

    best = None
    for locants in numberings(spec):
        key = (
            tuple(sorted(locants[a] for a in suffix_atoms)),
            multiple_locants(spec, locants),
            locants[junction],
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    locants = best[1]

    if suffix_name is not None:
        core, has_locants = _suffix_text(parent_text(spec, locants), suffix_name, [locants[a] for a in suffix_atoms], spec)
    else:
        core, has_locants = parent_text(spec, locants), False
    prefix_text = _prefix_text(entries, locants)
    return UnitText(
        _join(prefix_text, core), locants[junction], bool(prefix_text), has_locants or any(ch.isdigit() for ch in core)
    )


def name_ring_component(mol, ring_atoms, attachments, groups, suffix_group, name_function=None, directed=None):
    """(name, has_prefix) of a monocyclic ring used as a linker component, or
    None. `attachments`: [(ring_atom, external_atom)] for each free valence;
    `directed`: (unit_side_atom, center_side_atom) for a concatenated arm,
    where the unit-side atom takes the lowest locant (P-15.3.1.2.2.4)."""
    spec = spec_of(mol, ring_atoms)
    if spec is None:
        return None
    free_atoms = [a for a, _ in attachments]
    roots = _ring_roots(mol, spec, set(attachments))
    entries = _prefix_entries(mol, roots, groups, suffix_group, name_function)
    best = None
    for locants in numberings(spec):
        if directed is not None:
            free_key = (locants[directed[0]], locants[directed[1]])
        else:
            free_key = tuple(sorted(locants[a] for a in free_atoms))
        key = (
            free_key,
            multiple_locants(spec, locants),
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    locants = best[1]
    if directed is not None:
        cited = [locants[directed[1]], locants[directed[0]]]
    else:
        cited = sorted(locants[a] for a in free_atoms)
    loc = ",".join(str(x) for x in cited)
    if spec.kind == "benzene" and len(cited) == 2:
        body = f"{loc}-phenylene"
    else:
        word = _VALENCE_WORDS.get(len(cited))
        if word is None:
            return None
        body = f"{parent_text(spec, locants)}-{loc}-{word}"
    prefix_text = _prefix_text(entries, locants)
    return _join(prefix_text, body), bool(prefix_text)


def bare_polycyclic_unit(mol, atoms, junction, name_function=None):
    """UnitText for an unsubstituted fused/bridged/spiro ring system, with
    the junction's locant read off the name of its iodo derivative (the
    lowest locant of the junction's orbit, P-15.3.1.3)."""
    rw = Chem.RWMol(mol)
    iodine = rw.AddAtom(Chem.Atom(53))
    rw.AddBond(junction, iodine, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        rw.RemoveAtom(idx)
    probe = rw.GetMol()
    Chem.SanitizeMol(probe)
    name = probe_name(Chem.MolToSmiles(probe), name_function)
    match = re.fullmatch(r"(\d+[a-z]?)-iodo-?(.+)", name)
    if match is None:
        raise UnsupportedStructure(f"could not read the junction locant from {name!r}")
    locant, parent = match.groups()
    return UnitText(parent, locant, False, any(ch.isdigit() for ch in parent))
