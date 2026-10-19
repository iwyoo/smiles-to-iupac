"""Several ocene units in one molecule (P-69.2.7, P-15.3, P-44.1.1): the units, the groups between them and the groups
on them form a tree. A characteristic group outside the rings outranks the ocenes, which are then ocenyl prefixes of the
senior chain; otherwise an ocene is the parent: identical units joined by one linker are multiplied, any other tree takes
the most senior ocene as parent and cites the rest as ocenyl groups."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._ocene import (
    _ORDER,
    _PREFIX_WORD,
    _PRIME,
    _SUFFIX_ORDER,
    OCENES,
    _component_smiles,
    _fmt,
    _name_contracted,
    _stand_in,
    _group_kinds,
    _grouped,
    _numberings,
    _suffix_kind,
    _suffix_text,
    swap_stand_ins,
)
from ._pin import mark
from ._prefix_groups import PrefixNamer
from ._substituents import format_substituent_prefixes

_METAL_ORDER = ["Ni", "Co", "Fe", "Ru", "Os", "Cr", "V"]


class _Assembly:
    def __init__(self, mol, graph, units):
        self.mol, self.graph = mol, graph
        self.metals = [m for m, _ in units]
        self.rings = [rings for _, rings in units]
        self.atoms = [{m} | set(rings[0]) | set(rings[1]) for m, rings in units]
        self.every = set().union(*self.atoms)
        self.owner = {a: k for k, atoms in enumerate(self.atoms) for a in atoms}
        self.parent = [OCENES[mol.GetAtomWithIdx(m).GetSymbol()] for m in self.metals]
        self.roots = [[] for _ in units]
        self.direct = [[] for _ in units]
        for k, rings in enumerate(self.rings):
            for a in (a for ring in rings for a in ring):
                for x in graph[a]:
                    if x in self.atoms[k]:
                        continue
                    if x in self.every:
                        self.direct[k].append((a, self.owner[x], x))
                    else:
                        self.roots[k].append((a, x))
        self.comp_of, self.comp_atoms = {}, []
        for start in range(mol.GetNumAtoms()):
            if start in self.every or start in self.comp_of:
                continue
            members, stack = {start}, [start]
            while stack:
                for n in graph[stack.pop()]:
                    if n not in self.every and n not in members:
                        members.add(n)
                        stack.append(n)
            for a in members:
                self.comp_of[a] = len(self.comp_atoms)
            self.comp_atoms.append(members)
        self.touch = {}
        for k, roots in enumerate(self.roots):
            for a, x in roots:
                self.touch.setdefault(self.comp_of[x], []).append((k, a, x))
        self.suffix = {(k, a, x): _suffix_kind(mol, graph, a, x) for k, roots in enumerate(self.roots) for a, x in roots}
        self.kinds = {
            cid: _group_kinds(mol, graph, self.comp_atoms[cid])
            for cid, touches in self.touch.items()
            if not any(self.suffix[t] for t in touches)
        }

    def is_tree(self):
        edges = sum(len(r) for r in self.roots) + sum(len(d) for d in self.direct) // 2
        return edges == len(self.rings) + len(self.touch) - 1

    def group_name(self, k, via):
        """Ocenyl group of unit k attached through `via` = ('root', component) or ('direct', unit)."""
        named, attach = self._named(k, via)
        positions, grouped, _ = self._numbering(k, named, None, attach)
        stem = self.parent[k][:-1] if self.parent[k].endswith("e") else self.parent[k]
        return f"{format_substituent_prefixes(grouped) if grouped else ''}{stem}-1-yl"

    def _named(self, k, via, skip=()):
        named, attach = {}, None
        for a, x in self.roots[k]:
            if via[0] == "root" and self.comp_of[x] == via[1]:
                if attach is not None:
                    raise UnsupportedStructure("a bridge between the rings of a metallocene beside further metallocenes is not supported here")
                attach = a
            elif (a, x) in skip:
                continue
            elif self.suffix[(k, a, x)]:
                named[(a, x)] = (_PREFIX_WORD[self.suffix[(k, a, x)]], False)
            else:
                named[(a, x)] = self._branch(k, a, x)
        for a, j, b in self.direct[k]:
            if via[0] == "direct" and via[1] == j:
                attach = a
            else:
                named[(a, ("unit", j))] = (self.group_name(j, ("direct", k)), True)
        return named, attach

    def _branch(self, k, a, x):
        editable = Chem.RWMol(self.mol)
        cid = self.comp_of[x]
        if sum(1 for j, _, _ in self.touch[cid] if j == k) > 1:
            raise UnsupportedStructure("a bridge between the rings of a metallocene beside further metallocenes is not supported here")
        stands = []
        for j, b, y in self.touch[cid]:
            if j == k:
                continue
            editable.RemoveBond(y, b)
            group = self.group_name(j, ("root", cid))
            if _stand_in(editable, y, group):
                stands.append(group)
        mol = editable.GetMol()
        mol.UpdatePropertyCache(strict=False)
        name, compound = PrefixNamer(mol, adjacency(mol)).name(x, a)
        return swap_stand_ins(name, stands), compound

    def _numbering(self, k, named, principal, attach=None):
        best = None
        for positions in _numberings(self.graph, self.rings[k]):
            if attach is not None and positions[attach] != (1, 0):
                continue
            suffix = sorted(
                positions[a] for a, x in self.roots[k] if principal and self.suffix[(k, a, x)] == principal
            )
            grouped = _grouped(positions, named)
            key = (suffix, sorted(positions[a] for a, _ in named), [grouped[n]["locants"] for n in sorted(grouped)])
            if best is None or key < best[0]:
                best = (key, positions, grouped, suffix)
        return best[1], best[2], best[3]

    def principal(self, k):
        return min(
            (kind for (j, _, _), kind in self.suffix.items() if j == k and kind), key=_SUFFIX_ORDER.index, default=None
        )

    def ocene_parent(self, k, principal):
        named, _ = self._named(k, ("root", None), skip=())
        for a, x in self.roots[k]:
            if principal and self.suffix[(k, a, x)] == principal:
                named.pop((a, x), None)
        positions, grouped, suffix = self._numbering(k, named, principal)
        prefixes = format_substituent_prefixes(grouped) if grouped else ""
        stem = _suffix_text(self.parent[k], principal, suffix) if principal else self.parent[k]
        return f"{prefixes}{stem}"

    def unit_text(self, k, link):
        named, _ = self._named(k, ("root", None), skip=(link,))
        principal = self.principal(k)
        for a, x in self.roots[k]:
            if principal and self.suffix[(k, a, x)] == principal:
                named.pop((a, x), None)
        named[link] = ("\0link", False)
        positions, grouped, suffix = self._numbering(k, named, principal)
        locant = _fmt(positions[link[0]])
        grouped.pop("\0link")
        prefixes = format_substituent_prefixes(grouped) if grouped else ""
        stem = _suffix_text(self.parent[k], principal, suffix) if principal else self.parent[k]
        return f"{prefixes}{stem}", locant, bool(prefixes or principal)

    def rank(self, k):
        symbol = self.mol.GetAtomWithIdx(self.metals[k]).GetSymbol()
        return _METAL_ORDER.index(symbol)

    def senior_unit(self, principal):
        def key(k):
            count = sum(1 for (j, _, _), kind in self.suffix.items() if j == k and principal and kind == principal)
            return (-count, self.rank(k), -(len(self.roots[k]) + len(self.direct[k])), k)

        return min(range(len(self.rings)), key=key)

    def multiplied(self):
        linkers = [cid for cid, touches in self.touch.items() if len({k for k, _, _ in touches}) > 1]
        count = len(self.rings)
        if len(linkers) != 1 or any(self.direct):
            return None
        touches = self.touch[linkers[0]]
        if sorted(k for k, _, _ in touches) != list(range(count)):
            return None
        arms, smiles = [], set()
        for k, a, x in sorted(touches):
            arm = self.atoms[k] | {y for a2, x2 in self.roots[k] if (a2, x2) != (a, x) for y in self.comp_atoms[self.comp_of[x2]]}
            arms.append((a, x, arm))
            smiles.add(_component_smiles(self.mol, arm, a))
        if len(smiles) != 1:
            return None
        text, locant, use_bis = self.unit_text(0, (arms[0][0], arms[0][1]))
        lead = ",".join(f"{locant}{_PRIME * (2 * i)}" for i in range(count)) + "-"
        from ._chain_multiplicative import _attempt

        return _attempt(self.mol, self.graph, [], arms, "ocene", ocene_unit=(text, lead, use_bis))

    def directly_joined(self):
        """Two identical units joined by one bond: the ring assembly 'biferrocene' (P-28.2)."""
        if len(self.rings) != 2 or sum(len(d) for d in self.direct) != 2 or self.parent[0] != self.parent[1]:
            return None
        cited = []
        for k in (0, 1):
            named, attach = self._named(k, ("direct", 1 - k))
            principal = self.principal(k)
            for a, x in self.roots[k]:
                if principal and self.suffix[(k, a, x)] == principal:
                    named.pop((a, x), None)
            positions, grouped, suffix = self._numbering(k, named, principal, attach)
            cited.append((_fmt(positions[attach]), {n: v["locants"] for n, v in grouped.items()}, suffix, principal, grouped))
        if any(cited[0][i] != cited[1][i] for i in range(4)):
            return None
        first, _, suffix, principal, grouped = cited[0]
        merged = {
            name: {"locants": sorted(entry["locants"] + [loc + _PRIME * 2 for loc in entry["locants"]]), "compound": entry["compound"]}
            for name, entry in grouped.items()
        }
        prefixes = format_substituent_prefixes(merged) + "-" if merged else ""
        stem = f"{first},{first}{_PRIME * 2}-bi{self.parent[0]}"
        if principal:
            return prefixes + _suffix_text(f"[{stem}]", principal, suffix + [(n, p + 2) for n, p in suffix])
        return prefixes + stem

    def chain_parent(self, kind):
        holders = [cid for cid, kinds in self.kinds.items() if kind in kinds]
        used = {}
        keep = set().union(*(self.comp_atoms[cid] for cid in holders))
        for cid in holders:
            for k, a, x in self.touch[cid]:
                used[k] = used.get(k, 0) + 1
                if used[k] > 1:
                    raise UnsupportedStructure("a metallocene joined to the senior group at several places is not supported here")
        marks = [(x, self.group_name(k, ("root", cid))) for cid in holders for k, a, x in self.touch[cid]]
        return _name_contracted(Chem.RWMol(self.mol), set(range(self.mol.GetNumAtoms())) - keep, marks)


def name_assembly(mol, units, graph):
    asm = _Assembly(mol, graph, units)
    if not asm.is_tree():
        raise UnsupportedStructure("metallocenes joined in a ring other than by the carbon nodes of a phane are not supported here")
    ring_kinds = [kind for kind in asm.suffix.values() if kind]
    chain_kinds = [kind for kinds in asm.kinds.values() for kind in kinds]
    ring_best = min(ring_kinds, key=_SUFFIX_ORDER.index, default=None)
    chain_best = min(chain_kinds, key=_ORDER.index, default=None)
    ring_rank = _ORDER.index(ring_best) if ring_best else len(_ORDER)
    chain_rank = _ORDER.index(chain_best) if chain_best else len(_ORDER)
    if chain_best and (
        chain_rank < ring_rank or (chain_rank == ring_rank and chain_kinds.count(chain_best) > ring_kinds.count(ring_best))
    ):
        try:
            return asm.chain_parent(chain_best)
        except UnsupportedStructure:
            return mark(
                asm.ocene_parent(asm.senior_unit(ring_best), ring_best),
                "the senior group of a substituent is cited as a prefix because the substituent could not be named as the parent (P-41)",
            )
    joined = asm.directly_joined()
    if joined:
        return joined
    multiplied = asm.multiplied()
    if multiplied:
        return multiplied
    parent = asm.senior_unit(ring_best)
    return asm.ocene_parent(parent, ring_best)
