"""Compounds and substituent groups named on the retained parent structures of Appendix 3 / P-101.2.7
(`_appendix3_table.py`).

A molecule is named on a parent when it contains one skeleton whole, matched on connectivity alone: the rings of the
skeleton are rings of the molecule and no other ring touches them. The difference in multiple bonds between skeleton and molecule becomes 'ene'/'yne' endings in a saturated portion
(P-101.6.1) and 'hydro'/'dehydro' prefixes elsewhere (P-101.6.4 - P-101.6.7); the characteristic groups follow P-41 and
P-65/P-66 (P-101.7.1); a free valence gives a substituent group (P-101.7.3). P-10 identifies no preferred IUPAC names,
so every name is returned as valid but not preferred."""

import re
from collections import Counter
from itertools import product

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._appendix3_groups import SENIORITY, branch_counts, classify, reject_exotic
from ._appendix3_naming import N_CLASSES, Choice, assemble, n_roots
from ._appendix3_stereo import ParentStereo, describe, deviations
from ._appendix3_stereo_data import STEREO
from ._appendix3_table import SKELETONS
from ._common import UnsupportedStructure
from ._pin import mark
from ._substituents import BRANCH_STEREO

_NO_PIN = "P-101 identifies no preferred IUPAC names for natural-product parent structures"
_SUPERSCRIPT = {0: "", 1: "¹", 2: "²", 3: "³"}
_PRIMES = {0: "", 1: "′", 2: "″"}
_UNLABELED = 90000
_FORM_CAP = 256
_LABEL_PARTS = re.compile(r"(\d+)([a-c]?)([¹²³]?)([′″]?)")


def _label(map_number):
    if map_number >= _UNLABELED:
        return f"_{map_number - _UNLABELED}"
    number = map_number % 100
    prime = (map_number % 1000) // 100
    letter = (map_number % 10000) // 1000
    superscript = map_number // 10000
    return f"{number}{chr(96 + letter) if letter else ''}{_SUPERSCRIPT[superscript]}{_PRIMES[prime]}"


def sort_key(label):
    if label.startswith("_"):
        return (3, int(label[1:]), "", "")
    number, letter, superscript, prime = _LABEL_PARTS.fullmatch(label).groups()
    return (len(prime), int(number), letter, superscript)


def _flatten(mol):
    flat = Chem.RWMol(mol)
    for bond in flat.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    for atom in flat.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
    flat = flat.GetMol()
    flat.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(flat)
    return flat


def _perfect_matchings(atoms, edges):
    results = []

    def extend(free, pairs):
        if len(results) > _FORM_CAP:
            return
        if not free:
            results.append(tuple(pairs))
            return
        first = min(free)
        for other in edges[first]:
            if other in free:
                extend(free - {first, other}, pairs + [tuple(sorted((first, other)))])

    extend(frozenset(atoms), [])
    if not results or len(results) > _FORM_CAP:
        raise UnsupportedStructure("too many Kekule structures to compare with the parent")
    return results


def _kekule_forms(mol, restrict=None):
    """Every Kekule arrangement of the multiple bonds as a frozenset of (atom pair, order), limited to atoms in
    `restrict`: the non-aromatic multiple bonds as drawn, and each way of pairing the atoms of an aromatic system."""
    base = Chem.Mol(mol)
    Chem.Kekulize(base, clearAromaticFlags=False)
    inside = lambda a: restrict is None or a in restrict
    fixed = set()
    aromatic_edges = {}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = round(base.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble())
        if bond.GetIsAromatic():
            if order == 2:
                aromatic_edges.setdefault(a, set()).add(b)
                aromatic_edges.setdefault(b, set()).add(a)
            else:
                aromatic_edges.setdefault(a, set())
                aromatic_edges.setdefault(b, set())
            continue
        if order >= 2 and inside(a) and inside(b):
            fixed.add((tuple(sorted((a, b))), order))
    needing = {a for a in aromatic_edges if any(base.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() == 2.0 for n in aromatic_edges[a])}
    graph = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in needing and mol.GetBondBetweenAtoms(a, n.GetIdx()).GetIsAromatic()] for a in needing}
    components, seen = [], set()
    for atom in sorted(needing):
        if atom in seen:
            continue
        stack, comp = [atom], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(graph[x])
        seen |= comp
        if any(inside(a) for a in comp):
            components.append(comp)
    options = []
    for comp in components:
        options.append(
            {
                frozenset((pair, 2) for pair in matching if inside(pair[0]) and inside(pair[1]))
                for matching in _perfect_matchings(comp, {a: graph[a] for a in comp})
            }
        )
    forms = set()
    for combination in product(*options) if options else [()]:
        pairs = set(fixed)
        for part in combination:
            pairs |= part
        forms.add(frozenset(pairs))
    if len(forms) > _FORM_CAP:
        raise UnsupportedStructure("too many Kekule structures to compare with the parent")
    return forms


def _matching_forms(mol, restrict):
    """Every arrangement of double bonds that pairs up the atoms of `restrict` which carry a double bond in `mol`:
    a mancude system has no fixed Kekule structure (P-25.7.1.3)."""
    base = Chem.Mol(mol)
    Chem.Kekulize(base, clearAromaticFlags=True)
    atoms = sorted(
        a.GetIdx()
        for a in base.GetAtoms()
        if a.GetIdx() in restrict
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtomIdx(a.GetIdx()) in restrict for b in a.GetBonds())
    )
    edges = {a: [b.GetOtherAtomIdx(a) for b in base.GetAtomWithIdx(a).GetBonds() if b.GetOtherAtomIdx(a) in atoms] for a in atoms}
    return {frozenset((pair, 2) for pair in matching) for matching in _perfect_matchings(atoms, edges)}


def _systems(rings):
    """Sorted (atoms, rings) of every fused ring system among `rings` (tuples of atom indices)."""
    groups = []
    for ring in rings:
        atoms = set(ring)
        merged = [g for g in groups if g[0] & atoms]
        for g in merged:
            groups.remove(g)
            atoms |= g[0]
        groups.append((atoms, 1 + sum(g[1] for g in merged)))
    return sorted((len(atoms), count) for atoms, count in groups)


def _contains(big, small):
    big = list(big)
    for item in small:
        if item not in big:
            return False
        big.remove(item)
    return True


def _polyene_chain(name):
    """Main chain of a carotenoid or retinoid, in order, for the cis/trans citation of its double bonds."""
    if name == "retinal":
        return [str(n) for n in range(6, 16)]
    if name.endswith("carotene"):
        unprimed, primed = name[: -len("-carotene")].split(",")
        head = ["4", "5"] if unprimed == "ψ" else []
        tail = ["5′", "4′"] if primed == "ψ" else []
        return head + [str(n) for n in range(6, 16)] + [f"{n}′" for n in range(15, 5, -1)] + tail
    return []


class Skeleton:
    """One Appendix 3 parent: its graph, locants, mancude region, double-bond arrangements and implied configuration."""

    def __init__(self, name, smiles):
        self.name = name
        query = Chem.MolFromSmiles(smiles)
        self.labels = {a.GetIdx(): _label(a.GetAtomMapNum()) for a in query.GetAtoms()}
        for atom in query.GetAtoms():
            atom.SetAtomMapNum(0)
        self.query = query
        self.size = query.GetNumAtoms()
        self.flat = _flatten(query)
        self.ring_atoms = {a.GetIdx() for a in query.GetAtoms() if a.IsInRing()}
        self.ring_count = query.GetRingInfo().NumRings()
        self.systems = _systems(query.GetRingInfo().AtomRings())
        self.elements = Counter(a.GetAtomicNum() for a in query.GetAtoms())
        self.terminals = {self.labels[a.GetIdx()] for a in query.GetAtoms() if not a.IsInRing() and a.GetDegree() == 1}
        self.mancude = bool(re.search(r"\dH[,-]", name))
        self.forms = [
            {(self._pair(a, b), order) for (a, b), order in form}
            for form in (_matching_forms(query, set(self.labels)) if self.mancude else _kekule_forms(query))
        ]
        self.region_bonds = self._region_bonds()
        self.region = {l for form in self.forms for (pair, _) in form for l in pair if pair in self.region_bonds}
        self.outside_forms = [{(pair, o) for pair, o in form if pair not in self.region_bonds} for form in self.forms]
        self.unsaturated = frozenset(l for pair, _ in self.forms[0] for l in pair)
        self.polyene = _polyene_chain(name)
        self.polyene_bonds = self._polyene_bonds()
        self.stereo = None

    def _pair(self, a, b):
        return tuple(sorted((self.labels[a], self.labels[b]), key=sort_key))

    def _region_bonds(self):
        """The bonds of the mancude part, whose arrangement of double bonds is not fixed: the aromatic bonds, or every
        bond between doubly bonded atoms of a skeleton named with indicated hydrogen."""
        kekule = Chem.Mol(self.query)
        Chem.Kekulize(kekule, clearAromaticFlags=True)
        sp2 = {a.GetIdx() for a in kekule.GetAtoms() if any(b.GetBondTypeAsDouble() == 2.0 for b in a.GetBonds())}
        return {
            self._pair(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            for b in self.query.GetBonds()
            if (self.mancude and b.GetBeginAtomIdx() in sp2 and b.GetEndAtomIdx() in sp2) or b.GetIsAromatic()
        }

    def _polyene_bonds(self):
        """How many chain double bonds of a carotenoid or retinoid can be cited as cis or trans."""
        chain = self.polyene
        count = 0
        for pair, order in self.forms[0]:
            if order == 2 and all(l in chain for l in pair):
                low, high = sorted(chain.index(l) for l in pair)
                count += 0 < low and high + 1 < len(chain)
        return count


_SKELETONS = {name: Skeleton(name, smiles) for name, smiles in SKELETONS.items()}
for _name, (_smiles, _ref, _anticlockwise) in STEREO.items():
    _SKELETONS[_name].stereo = ParentStereo(_SKELETONS[_name], _smiles, _ref, _anticlockwise, _label, sort_key)


def _adjacent(pair):
    """True when the two locants of a double bond are consecutive plain numbers of one series, so the bond is cited by
    the lower alone (P-31.1.4.2.4)."""
    if any(label.startswith("_") for label in pair):
        return False
    parts = [_LABEL_PARTS.fullmatch(label).groups() for label in pair]
    (n1, l1, s1, p1), (n2, l2, s2, p2) = parts
    return not (l1 or l2 or s1 or s2) and p1 == p2 and int(n2) - int(n1) == 1


def _unsaturation(skeleton, mol, mapping):
    """(lost, gained, unsaturated, hydro): outside the mancude part of the skeleton, the double-bond units the molecule
    lacks (`lost`: label pairs) and adds (`gained`: label pair, old order, new order) under the Kekule pair of fewest
    changes, compound locants avoided, then lowest locants; inside it, the atoms that are no longer doubly bonded
    (`hydro`). `unsaturated`: the labels of the skeleton's own multiple bonds."""
    atom_label = {atom: label for label, atom in mapping.items()}
    pairs_of = lambda form: {
        (tuple(sorted((atom_label[a], atom_label[b]), key=sort_key)), order) for (a, b), order in form
    }
    mol_forms = [pairs_of(form) for form in _kekule_forms(mol, set(atom_label))]
    region_pairs = skeleton.region_bonds
    sp2_in_mol = {l for form in mol_forms for (pair, _) in form for l in pair}
    hydro = sorted((l for l in skeleton.region if l not in sp2_in_mol), key=sort_key)
    outside_mol = [{(pair, o) for pair, o in form if pair not in region_pairs} for form in mol_forms]
    best = None
    for ks, km in product(skeleton.outside_forms, outside_mol):
        s_orders, m_orders = dict(ks), dict(km)
        lost, gained = [], []
        for pair in set(s_orders) | set(m_orders):
            q, p = s_orders.get(pair, 1), m_orders.get(pair, 1)
            if q > p:
                lost.extend([pair] * (q - p))
            elif p > q:
                gained.append((pair, q, p))
        units = len(lost) + sum(p - q for _, q, p in gained)
        locants = (
            tuple(sorted(sort_key(l) for pair in lost for l in pair)),
            tuple(sorted(sort_key(l) for pair, _, _ in gained for l in pair)),
        )
        compound = sum(1 for pair, _, _ in gained if not _adjacent(pair)) + sum(1 for pair in lost if not _adjacent(pair))
        rank = (units, compound, locants)
        if best is None or rank < best[0]:
            best = (rank, lost, gained)
    _, lost, gained = best
    key = lambda p: tuple(sort_key(l) for l in p)
    return sorted(lost, key=key), sorted(gained, key=lambda g: key(g[0])), skeleton.unsaturated, hydro


def _candidates(skeleton, mol, flat):
    found = []
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    for match in flat.GetSubstructMatches(skeleton.flat, uniquify=False, maxMatches=5000):
        matched_ring = {match[q] for q in skeleton.ring_atoms}
        touching = [r for r in rings if r & matched_ring]
        if any(not r <= matched_ring for r in touching) or len(touching) != skeleton.ring_count:
            continue
        if not all(mol.GetAtomWithIdx(a).IsInRing() for a in matched_ring):
            continue
        if any(
            mol.GetAtomWithIdx(match[q.GetIdx()]).IsInRing()
            for q in skeleton.query.GetAtoms()
            if q.GetIdx() not in skeleton.ring_atoms
        ):
            continue
        if any(
            mol.GetAtomWithIdx(match[q.GetIdx()]).GetFormalCharge() != q.GetFormalCharge()
            for q in skeleton.query.GetAtoms()
        ):
            continue
        found.append({skeleton.labels[q]: match[q] for q in range(len(match))})
    return found


def _skeleton_choices(skeleton, mol, flat):
    choices = []
    for mapping in _candidates(skeleton, mol, flat):
        try:
            groups, branches, attach = classify(mol, mapping, skeleton.terminals)
            attach = tuple(sorted(attach, key=lambda item: sort_key(item[0]))) if attach else None
            lost, gained, unsaturated, hydro = _unsaturation(skeleton, mol, mapping)
        except UnsupportedStructure:
            continue
        choices.append(Choice(mapping, groups, branches, attach, lost, gained, unsaturated, hydro))
    return choices


def _choice_key(choice):
    principal = next((c for c in SENIORITY if any(g.cls == c for g in choice.groups)), None)
    return (
        len(choice.lost) + len(choice.hydro) + sum(p - q for _, q, p in choice.gained),
        tuple(sort_key(l) for l, _ in choice.attach) if choice.attach else (),
        tuple(sorted(sort_key(g.label) for g in choice.groups if g.cls == principal)),
        tuple(sorted(sort_key(l) for pair in choice.lost for l in pair)),
        tuple(sorted(sort_key(l) for pair, _, _ in choice.gained for l in pair)),
        tuple(sorted(sort_key(g.label) for g in choice.groups)),
        tuple(sorted(sort_key(l) for l, _ in choice.branches)),
    )


def _senior_counts(mol, choice):
    """(own, branch, arms): the members of each principal class on the skeleton, the most in one side branch, and
    the classes of each plain side branch as (skeleton atom, root, Counter)."""
    mapped = set(choice.mapping.values())
    own = {}
    for group in choice.groups:
        cls = "ester" if group.cls == "ester_o" else group.cls
        own[cls] = own.get(cls, 0) + 1
    arms = [
        (choice.mapping[label], root, branch_counts(mol, root, choice.mapping[label], mapped))
        for label, root in choice.branches
    ]
    reached = [(root, anchor) for anchor, root, _ in arms]
    for group in choice.groups:
        if group.cls in N_CLASSES:
            reached += [
                (root, nitrogen)
                for _, nitrogen, roots in n_roots(mol, group)
                for root in roots
                if mol.GetAtomWithIdx(root).GetAtomicNum() == 6
            ]
    branch = {}
    for root, anchor in reached:
        for cls, count in branch_counts(mol, root, anchor, mapped).items():
            branch[cls] = max(branch.get(cls, 0), count)
    return own, branch, arms


def _parent_adequate(mol, choice):
    """P-44.1.1: the skeleton is the parent only when it carries at least as many members of the most senior class
    present as any side branch does."""
    if choice.attach is not None:
        return True
    own, branch, _ = _senior_counts(mol, choice)
    top = next((c for c in SENIORITY if c in own or c in branch), None)
    return top is None or own.get(top, 0) >= max(1, branch.get(top, 0))


def _best_skeleton(mol, attach_allowed):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    systems = _systems(mol.GetRingInfo().AtomRings())
    counts = Counter(atom.GetAtomicNum() for atom in mol.GetAtoms())
    flat = None
    best = None
    for skeleton in _SKELETONS.values():
        if mol.GetNumAtoms() < skeleton.size or not _contains(systems, skeleton.systems):
            continue
        if any(counts[z] < n for z, n in skeleton.elements.items()):
            continue
        flat = flat or _flatten(mol)
        choices = [c for c in _skeleton_choices(skeleton, mol, flat) if attach_allowed or c.attach is None]
        if not choices:
            continue
        chosen = min(choices, key=_choice_key)
        changes = len(chosen.lost) + len(chosen.hydro) + sum(p - q for _, q, p in chosen.gained)
        if 2 * len(chosen.lost) + len(chosen.hydro) + 2 * sum(p - q for _, q, p in chosen.gained) > skeleton.size // 2:
            continue
        rank = (skeleton.size, -changes, -deviations(mol, chosen.mapping, skeleton.stereo))
        if best is None or rank > best[0]:
            best = (rank, skeleton, chosen)
    return best and best[1:]


def name_on_skeleton(mol, attach_allowed=False):
    """(name, valence, substituted) of `mol` on an Appendix 3 parent, or None."""
    best = _best_skeleton(mol, attach_allowed)
    if best is None:
        return None
    skeleton, choice = best
    reject_exotic(mol, set(choice.mapping.values()))
    if not _parent_adequate(mol, choice):
        return None
    stereo = describe(mol, skeleton, choice.mapping, sort_key, skeleton.stereo)
    valence = len(choice.attach) if choice.attach else 0
    return assemble(mol, skeleton, choice, stereo, sort_key), valence, len(mol.GetAtoms()) - valence > len(choice.mapping)


def name_appendix3_skeleton(mol):
    """Name of `mol` on an Appendix 3 retained parent, or None."""
    try:
        result = name_on_skeleton(mol)
    except Exception:
        return None
    return mark(result[0], _NO_PIN) if result else None


_MIN_SIZE = min(skeleton.size for skeleton in _SKELETONS.values())


def _branch_of(graph, roots, blocked):
    branch, stack = set(roots), list(roots)
    while stack:
        for neighbor in graph[stack.pop()]:
            if neighbor not in blocked and neighbor not in branch:
                branch.add(neighbor)
                stack.append(neighbor)
    return branch


def _holds_a_skeleton(mol, branch):
    if len(branch) < _MIN_SIZE:
        return False
    rings = [ring for ring in mol.GetRingInfo().AtomRings() if branch.issuperset(ring)]
    return bool(rings) and any(_contains(_systems(rings), skeleton.systems) for skeleton in _SKELETONS.values())


def _cite_stereo(mol, branch):
    context = BRANCH_STEREO.get()
    if context:
        context["used"].update(("atom", a) for a in branch)
        context["used"].update(
            ("bond", (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in branch and b.GetEndAtomIdx() in branch
        )


def appendix3_group(mol, graph, root, coming_from):
    """(name, True) of the substituent group rooted at `root` when it is an Appendix 3 parent with a free valence at
    `root` (P-101.7.3), else None. The atom it hangs from is kept as a dummy so that its centres keep their tags and
    the CIP labels of the whole molecule."""
    branch = _branch_of(graph, [root], {coming_from})
    if not _holds_a_skeleton(mol, branch):
        return None
    try:
        result = name_on_skeleton(_fragment(mol, branch, [coming_from]), attach_allowed=True)
    except Exception:
        return None
    if not result or result[1] != 1:
        return None
    _cite_stereo(mol, branch)
    return mark(result[0], _NO_PIN), True


def appendix3_multivalent_group(mol, graph, ring_atoms, attachments):
    """(name, substituted, atoms) of the multivalent group on an Appendix 3 parent that contains the ring system `ring_atoms`
    and is joined to the rest of the molecule by `attachments` ([(group atom, outside atom)]), or None (P-101.7.3,
    P-15.3)."""
    outside = {external for _, external in attachments}
    if len(outside) != len(attachments):
        return None
    branch = _branch_of(graph, ring_atoms, outside)
    if not _holds_a_skeleton(mol, branch) or any(external in branch for external in outside):
        return None
    try:
        result = name_on_skeleton(_fragment(mol, branch, sorted(outside)), attach_allowed=True)
    except Exception:
        return None
    if not result or result[1] != len(attachments):
        return None
    _cite_stereo(mol, branch)
    return mark(result[0], _NO_PIN), result[2], branch


def appendix3_central_group(mol, graph):
    """(name, attachments, atoms) of the multivalent group of an Appendix 3 parent whose side branches carry the most
    senior class present, while the parent carries none of it (P-15.3, P-44.1.1): `attachments` are the
    [(skeleton atom, branch root)] pairs, `atoms` the group without the branches. None when there is no such parent."""
    try:
        best = _best_skeleton(mol, False)
        if best is None:
            return None
        _, choice = best
        own, branch, arms = _senior_counts(mol, choice)
        top = next((c for c in SENIORITY if c in own or c in branch), None)
        if top is None or own.get(top):
            return None
        attachments = [(atom, root) for atom, root, counts in arms if counts.get(top)]
        if len(attachments) < 2:
            return None
        found = appendix3_multivalent_group(mol, graph, set(choice.mapping.values()), attachments)
    except UnsupportedStructure:
        return None
    return None if found is None else (found[0], attachments, found[2])


def _fragment(mol, branch, attached_to):
    labelled = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(labelled)
    keep = sorted(branch | set(attached_to))
    editable = Chem.RWMol(labelled)
    for index in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in keep), reverse=True):
        editable.RemoveAtom(index)
    fragment = editable.GetMol()
    for atom in attached_to:
        dummy = fragment.GetAtomWithIdx(keep.index(atom))
        dummy.SetAtomicNum(0)
        dummy.SetFormalCharge(0)
        dummy.SetIsotope(0)
        dummy.SetNumExplicitHs(0)
        dummy.SetNoImplicit(True)
        dummy.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        dummy.SetIsAromatic(False)
    fragment.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(fragment)
    fragment.SetBoolProp("cip_assigned", True)
    return fragment
