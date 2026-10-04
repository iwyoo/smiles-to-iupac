"""Substituent-prefix names for groups hung on a parent skeleton (P-65, P-66,
P-68): halogens, hydroxy/alkoxy, sulfanyl, amino, phosphanyl, silyl-type
groups, nitro, cyano and carbonyl-derived groups, with alkyl/aryl parts named
through `name_branch` and every heteroatom inside a branch resolved first.
"""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._hetero_prefixes import CHALCOGEN_PREFIXES, require_plain_chalcogen_kids, require_senior_group
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_SIMPLE = {"amino", "hydroxy", "sulfanyl", "methoxy", "ethoxy", "propoxy", "butoxy", "phenoxy", "nitro", "cyano"}
_METALLOID = {14: "silyl", 32: "germyl", 50: "stannyl", 82: "plumbyl"}
_CONTRACTED = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy", "phenyl": "phenoxy"}
_ACYL = {"methyl": "acetyl", "phenyl": "benzoyl", "ethyl": "propanoyl", "propyl": "butanoyl"}


def _oxy(name, compound):
    if name in _CONTRACTED:
        return _CONTRACTED[name], False
    if name.endswith("phenyl") and not compound:
        return name[: -len("phenyl")] + "phenoxy", True
    if name.endswith(("methyl", "ethyl", "propyl", "butyl")) and "cyclo" not in name:
        return name[:-2] + "oxy", compound or any(ch.isdigit() for ch in name)
    return name + "oxy", compound or any(ch.isdigit() for ch in name)


def enclose(name: str) -> str:
    if name in _SIMPLE:
        return name
    if "{" in name:
        return f"[{name}]"
    if "[" in name:
        return "{" + name + "}"
    if "(" in name:
        return f"[{name}]"
    return f"({name})"


class PrefixNamer:
    def __init__(self, mol, graph, sigma=frozenset()):
        self.mol, self.graph, self.sigma = mol, graph, frozenset(sigma)
        self._cache: dict[tuple[int, int], tuple[str, bool]] = {}
        self._marked: set[int] = set()

    def _z(self, i):
        return self.mol.GetAtomWithIdx(i).GetAtomicNum()

    def placeholder(self, n, parent):
        name, compound = self.name(n, parent)
        return enclose(name) if compound else name

    def _rings_of(self, atom):
        return [frozenset(r) for r in self.mol.GetRingInfo().AtomRings() if atom in r]

    def _region_halogens(self, root, parent):
        seen, stack, hetero = {root}, [root], {}
        while stack:
            c = stack.pop()
            for y in self.graph[c]:
                if y == parent or y in seen:
                    continue
                if self._z(y) == 6:
                    ring_y = self._rings_of(y)
                    if ring_y and not any(c in r for r in ring_y) and y not in hetero:
                        hetero[y] = self.placeholder(y, c)
                        continue
                    seen.add(y)
                    stack.append(y)
                elif y not in hetero:
                    if self.mol.GetAtomWithIdx(y).IsInRing() and any(c in r for r in self._rings_of(y)):
                        raise UnsupportedStructure("a ring heteroatom needs the ring-group namer")
                    hetero[y] = self.placeholder(y, c)
        return hetero

    def name(self, n, parent):
        key = (n, parent)
        if key not in self._cache:
            name, compound = self._name(n, parent)
            if n in self.sigma:
                name = f"{name}-\u03ba{self.mol.GetAtomWithIdx(n).GetSymbol()}"
                compound = True
                self._marked.add(n)
            else:
                deep = (self.sigma & set(self._group_atoms(n, parent))) - self._marked
                if deep:
                    name = f"{name}{self._deep_kappa(n, parent, deep)}"
                    compound = True
                    self._marked |= deep
            self._cache[key] = (name, compound)
        return self._cache[key]

    def _deep_kappa(self, n, parent, deep):
        group = self._group_atoms(n, parent)
        symbols = sorted(self.mol.GetAtomWithIdx(i).GetSymbol() for i in deep)
        for symbol in set(symbols):
            total = sum(1 for i in group if self.mol.GetAtomWithIdx(i).GetSymbol() == symbol)
            if total != symbols.count(symbol):
                raise UnsupportedStructure("the kappa donor needs a locant inside this substituent group")
        seen: dict[str, int] = {}
        cited = []
        for symbol in symbols:
            cited.append(symbol + "'" * seen.get(symbol, 0))
            seen[symbol] = seen.get(symbol, 0) + 1
        return f"-\u03ba{len(symbols) if len(symbols) > 1 else ''}{','.join(cited)}"

    def _group_atoms(self, n, parent):
        seen, stack = {n}, [n]
        while stack:
            for v in self.graph[stack.pop()]:
                if v != parent and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    def _entries(self, n, parent):
        return [self.name(q, n) for q in self.graph[n] if q != parent]

    def _name(self, n, parent):
        mol = self.mol
        atom = mol.GetAtomWithIdx(n)
        z = atom.GetAtomicNum()
        if z in HALOGEN_PREFIXES:
            return HALOGEN_PREFIXES[z], False
        if atom.GetFormalCharge() and z != 7:
            raise UnsupportedStructure("a charged substituent group is not supported here")
        if z == 6:
            return self._carbon(n, parent)
        if z in (8, 16, 34, 52) and any(self._z(q) != 6 for q in self.graph[n] if q != parent):
            raise UnsupportedStructure("adjacent heteroatoms are not supported here")
        entries = self._entries(n, parent)
        if z == 8:
            if not entries:
                if mol.GetBondBetweenAtoms(n, parent).GetBondTypeAsDouble() == 2.0:
                    return "oxo", False
                return "hydroxy", False
            return _oxy(*entries[0])
        if z in CHALCOGEN_PREFIXES:
            word = CHALCOGEN_PREFIXES[z]
            require_plain_chalcogen_kids(mol, z, [q for q in self.graph[n] if q != parent])
            if not entries:
                require_senior_group(mol, z)
                return word, False
            return format_mononuclear_prefixes(entries) + word, True
        if z == 7:
            return self._nitrogen(n, parent, entries)
        if z == 15:
            return format_mononuclear_prefixes(entries) + "phosphanyl", True
        if z in _METALLOID:
            return (format_mononuclear_prefixes(entries) if entries else "") + _METALLOID[z], bool(entries)
        raise UnsupportedStructure("this substituent group is not supported here")

    def _nitrogen(self, n, parent, entries):
        atom = self.mol.GetAtomWithIdx(n)
        oxygens = [q for q in self.graph[n] if q != parent and self._z(q) == 8 and self.mol.GetAtomWithIdx(q).GetDegree() == 1]
        if len(oxygens) == 2 and atom.GetDegree() == 3:
            return "nitro", False
        if atom.GetFormalCharge() or any(self._z(q) not in (6,) for q in self.graph[n] if q != parent):
            raise UnsupportedStructure("this nitrogen substituent is not supported here")
        if not entries:
            return "amino", False
        return format_mononuclear_prefixes(entries) + "amino", True

    def _carbon(self, n, parent):
        mol, graph = self.mol, self.graph
        atom = mol.GetAtomWithIdx(n)
        bond = mol.GetBondBetweenAtoms(n, parent).GetBondTypeAsDouble()
        others = [q for q in graph[n] if q != parent]
        triple_n = [q for q in others if self._z(q) == 7 and mol.GetBondBetweenAtoms(n, q).GetBondTypeAsDouble() == 3.0]
        if triple_n and len(others) == 1 and bond == 1.0:
            return "cyano", False
        oxo = [q for q in others if self._z(q) == 8 and mol.GetBondBetweenAtoms(n, q).GetBondTypeAsDouble() == 2.0]
        if oxo and bond == 1.0 and not atom.GetIsAromatic():
            rest = [q for q in others if q != oxo[0]]
            if not rest:
                return "formyl", False
            if len(rest) == 1:
                q = rest[0]
                z = self._z(q)
                if z == 8:
                    inner = [x for x in graph[q] if x != n]
                    if not inner:
                        return "carboxy", False
                    return _oxy(*self.name(inner[0], q))[0] + "carbonyl", True
                if z == 7 and mol.GetAtomWithIdx(q).GetDegree() == 1:
                    return "carbamoyl", False
                if z == 6:
                    base = self.name(q, n)
                    if base[0] in _ACYL and not base[1]:
                        return _ACYL[base[0]], False
        ring = self._ring_group(n, parent) if atom.IsInRing() else None
        if ring is not None:
            return ring
        return name_branch(graph, n, parent, self._region_halogens(n, parent), mol=mol)

    def _ring_group(self, n, parent):
        """Monocyclic hetero-aromatic substituent group ('pyridin-2-yl',
        '5-methylfuran-2-yl'): the heteroatom is 1, the attachment atom takes
        the lowest locant, then the substituents."""
        rings = self._rings_of(n)
        if len(rings) != 1:
            return None
        ring = set(rings[0])
        if any(len(self._rings_of(a)) != 1 for a in ring):
            return None
        mol = self.mol
        hetero = [a for a in ring if self._z(a) != 6]
        if len(ring) not in (5, 6) or not all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring) or len(hetero) != 1:
            return None
        from ._hapto_ext import _HETERO5, _HETERO6, _ring_orders, _yl

        table = _HETERO5 if len(ring) == 5 else _HETERO6
        z = self._z(hetero[0])
        if z not in table:
            return None
        best = None
        for order in _ring_orders(self.graph, ring, hetero[0]):
            pos = {a: i + 1 for i, a in enumerate(order)}
            grouped: dict = {}
            for a in order:
                for q in self.graph[a]:
                    if q in ring or (a == n and q == parent):
                        continue
                    name, compound = self.name(q, a)
                    entry = grouped.setdefault(name, {"locants": [], "compound": compound})
                    entry["locants"].append(pos[a])
            key = (pos[n], sorted(l for e in grouped.values() for l in e["locants"]),
                   [grouped[k]["locants"] for k in sorted(grouped)])
            if best is None or key < best[0]:
                best = (key, grouped, pos[n])
        _, grouped, attach = best
        prefixes = format_substituent_prefixes(grouped)
        stem = _yl(table[z], "yl", attach)
        if table[z].endswith("e") is False:
            stem = table[z] + f"-{attach}-yl"
        return f"{prefixes}{stem}", True
