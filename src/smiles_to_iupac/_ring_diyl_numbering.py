"""Numbering candidates and parent names for ring systems cited as multivalent groups (P-29.2, P-31.1.4):
monocycles (carbocycle, benzene, mancude/hydro/saturated heterocycles by retained and Hantzsch-Widman names,
P-22.2.2), ortho-fused arene chains, ortho-fused mancude heterocycles, von Baeyer and monospiro carbocycles.
"""

from itertools import combinations

from rdkit import Chem

from . import _aromatic
from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import (
    UnsupportedStructure,
    ring_bond_locants,
    von_baeyer_unsaturation_citations,
)
from ._fusion_numbering_general import general_peripheral_numberings
from ._hetero_monocyclic import (
    _ROLE_SEQUENCES,
    saturated_five_membered_1_2_two_heteroatom_ring_name,
    saturated_five_membered_1_3_two_heteroatom_ring_name,
    saturated_ring_name,
    saturated_seven_membered_1_2_two_heteroatom_ring_name,
    saturated_seven_membered_1_3_two_heteroatom_ring_name,
    saturated_seven_membered_1_4_two_heteroatom_ring_name,
    saturated_two_heteroatom_1_4_ring_name,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix
from ._polycyclic import find_polycyclic_core, iter_polycyclic_candidates
from ._spiro import find_monospiro_atom, iter_monospiro_numberings
from ._unsaturated import _unsaturation_suffix_from_citations
from ._common import multiplied_word

_SENIORITY = ["F", "Cl", "Br", "I", "O", "S", "Se", "Te", "N", "P", "As", "Sb", "Bi", "Si", "Ge", "Sn", "Pb", "B", "Al", "Ga", "In", "Tl"]
_RANK = {e: i for i, e in enumerate(_SENIORITY)}
_PREFIX = {
    "O": "oxa", "S": "thia", "Se": "selena", "Te": "tellura", "N": "aza", "P": "phospha", "As": "arsa", "Sb": "stiba",
    "Bi": "bisma", "Si": "sila", "Ge": "germa", "Sn": "stanna", "Pb": "plumba", "B": "bora", "Al": "aluma",
    "Ga": "galla", "In": "indiga", "Tl": "thalla",
}
_CHALCOGENS = {"O", "S", "Se", "Te"}
_NO_DOUBLE_BOND = _CHALCOGENS | {"F", "Cl", "Br", "I"}
_SIX_A = {"O", "S", "Se", "Te", "Bi"}
_SIX_B = {"N", "Si", "Ge", "Sn", "Pb"}
_STEMS_UNSAT = {3: "irene", 4: "ete", 5: "ole", 7: "epine", 8: "ocine", 9: "onine", 10: "ecine"}
_STEMS_SAT = {3: "irane", 4: "etane", 5: "olane", 7: "epane", 8: "ocane", 9: "onane", 10: "ecane"}
_STEMS_SAT_N = {3: "iridine", 4: "etidine", 5: "olidine"}
_TWO_HETERO_SATURATED = {
    (6, 1, 4): saturated_two_heteroatom_1_4_ring_name,
    (5, 1, 3): saturated_five_membered_1_3_two_heteroatom_ring_name,
    (5, 1, 2): saturated_five_membered_1_2_two_heteroatom_ring_name,
    (7, 1, 4): saturated_seven_membered_1_4_two_heteroatom_ring_name,
    (7, 1, 3): saturated_seven_membered_1_3_two_heteroatom_ring_name,
    (7, 1, 2): saturated_seven_membered_1_2_two_heteroatom_ring_name,
}
_MANCUDE_RETAINED = {
    tuple(role[0] for role in roles): name[3:] if name.startswith("1H-") else name
    for name, roles in _ROLE_SEQUENCES.items()
    if name != "phosphinine"
}
_MANCUDE_RETAINED.update({(e, "C", "C", "C", "C", "C"): n for e, n in (("O", "pyran"), ("S", "thiopyran"), ("Se", "selenopyran"), ("Te", "telluropyran"))})
_LACKS_REGULAR_NUMBERING = ("acridine", "carbazole", "xanthene", "purine", "anthracene", "phenanthrene")


class Numbering:
    def __init__(self, position_of, text, pre_key=(), unsat_key=(), ih=()):
        self.position_of = position_of
        self.text = text
        self.pre_key = pre_key
        self.unsat_key = unsat_key
        self.ih = ih


def _walks(ring_order):
    n = len(ring_order)
    for base in (ring_order, list(reversed(ring_order))):
        for start in range(n):
            yield base[start:] + base[:start]


def _yl(valence):
    return f"{multiplying_prefix(valence) if valence > 1 else ''}yl"


def _tail(stem, locants, valence, suffix="yl"):
    if suffix == "carboxylate":
        return f"{stem}-{_locs(locants)}-carboxylate"
    if valence == 1:
        base = stem[:-1] if stem.endswith("e") else stem
        return f"{base}-{locants[0]}-yl"
    return f"{stem}-{_locs(locants)}-{_yl(valence)}"


def _ring_size(stem):
    from ._numerals import alkane_name as _an

    return next(n for n in range(3, 60) if _an(n) == stem[len("cyclo"):])


def _locs(locants):
    return ",".join(str(loc) for loc in sorted(locants))


def _hantzsch_widman(elements, saturated):
    size = len(elements)
    hetero = [e for e in elements if e != "C"]
    if size < 3 or size > 10 or not hetero or any(e not in _PREFIX for e in hetero):
        return None
    counts = {}
    for e in hetero:
        counts[e] = counts.get(e, 0) + 1
    ordered = sorted(counts, key=lambda e: _RANK[e])
    last = ordered[-1]
    if size == 6:
        group = "A" if last in _SIX_A else "B" if last in _SIX_B else "C"
        stem = ("ane" if group == "A" else "inane") if saturated else ("ine" if group != "C" else "inine")
    elif saturated:
        stem = _STEMS_SAT_N[size] if size in _STEMS_SAT_N and "N" in counts else _STEMS_SAT[size]
    else:
        stem = "irine" if size == 3 and set(counts) == {"N"} else _STEMS_UNSAT[size]
    pieces = []
    for e in ordered:
        mult = "" if counts[e] == 1 else multiplying_prefix(counts[e])
        pieces.append(mult + _PREFIX[e])
    text = ""
    for piece in pieces + [stem]:
        if text and piece[0] in "aeiou" and text.endswith("a"):
            text = text[:-1]
        text += piece
    return text


def _max_matching(atoms_can_hold, ring_order):
    n = len(ring_order)
    can = [a in atoms_can_hold for a in ring_order]
    edges = [(i, (i + 1) % n) for i in range(n) if can[i] and can[(i + 1) % n]]
    best = 0
    for r in range(len(edges), 0, -1):
        for combo in combinations(edges, r):
            used = [x for edge in combo for x in edge]
            if len(set(used)) == len(used):
                return r
    return best


def _hetero_monocycle(mol, ring_order, attached):
    atoms = [mol.GetAtomWithIdx(i) for i in ring_order]
    size = len(atoms)
    ring_set = set(ring_order)
    sym = {a.GetIdx(): a.GetSymbol() for a in atoms}
    if any(a.GetFormalCharge() != 0 for a in atoms):
        raise UnsupportedStructure("a charged ring atom is not supported as a diyl yet")
    aromatic = all(a.GetIsAromatic() for a in atoms)
    ring_double = {
        i
        for b in mol.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
        for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    }
    for a in atoms:
        for n in a.GetNeighbors():
            exocyclic_double = (
                n.GetIdx() not in ring_set and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
            )
            if exocyclic_double and (aromatic or ring_double):
                raise UnsupportedStructure("a ring atom with an exocyclic double bond is not supported as a diyl yet")
    if not aromatic and any(
        b.GetBondTypeAsDouble() not in (1.0, 2.0)
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
    ):
        raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
    can_hold = {i for i in ring_order if sym[i] not in _NO_DOUBLE_BOND}
    if aromatic:
        in_double = {
            i
            for i in can_hold
            if not (sym[i] != "C" and (mol.GetAtomWithIdx(i).GetTotalNumHs() > 0 or mol.GetAtomWithIdx(i).GetDegree() == 3))
        }
    else:
        in_double = ring_double
    saturated_atoms = can_hold - in_double
    matching = _max_matching(can_hold, ring_order)
    mancude_sp3 = len(can_hold) - 2 * matching
    fully_saturated = not ring_double and not aromatic
    if not fully_saturated:
        extra = len(saturated_atoms) - mancude_sp3
        if extra < 0 or extra % 2:
            raise UnsupportedStructure("this partly hydrogenated ring is not expressible with hydro prefixes yet")
    results = []
    best_pre = None
    senior = min(_RANK.get(sym[i], 99) for i in ring_order if sym[i] != "C")

    def starts_at_senior(walk):
        return _RANK.get(sym[walk[0]], 99) == senior

    for walk in _walks(ring_order):
        if not starts_at_senior(walk):
            continue
        hetero = [(i + 1, _RANK.get(sym[a], 99)) for i, a in enumerate(walk) if sym[a] != "C"]
        pre = (tuple(p for p, _ in hetero), tuple(r for _, r in hetero))
        if best_pre is None or pre < best_pre:
            best_pre = pre
    for walk in _walks(ring_order):
        if not starts_at_senior(walk):
            continue
        hetero = [(i + 1, _RANK.get(sym[a], 99)) for i, a in enumerate(walk) if sym[a] != "C"]
        pre = (tuple(p for p, _ in hetero), tuple(r for _, r in hetero))
        if pre != best_pre:
            continue
        position_of = {a: i + 1 for i, a in enumerate(walk)}
        elements = tuple(sym[a] for a in walk)
        if fully_saturated:
            stem = _saturated_name(elements, hetero)
            results.append((position_of, stem, (), ()))
            continue
        stem = _MANCUDE_RETAINED.get(elements)
        if stem is None:
            stem = _hantzsch_widman(elements, saturated=False)
            if stem is None:
                raise UnsupportedStructure("this heteromonocycle has no supported mancude parent name yet")
            stem = _with_hetero_locants(stem, elements, hetero)
        sat_pos = sorted(position_of[a] for a in saturated_atoms)
        ih = tuple(sat_pos[:mancude_sp3])
        hydro = tuple(sat_pos[mancude_sp3:])
        results.append((position_of, stem, ih, hydro))
    return best_pre, results


def _saturated_name(elements, hetero):
    positions = [p for p, _ in hetero]
    size = len(elements)
    names = [e for e in elements if e != "C"]
    if len(names) == 1:
        stem = saturated_ring_name(names[0], size)
    elif len(names) == 2 and (size, positions[0], positions[1]) in _TWO_HETERO_SATURATED:
        stem = _TWO_HETERO_SATURATED[(size, positions[0], positions[1])](names)
    else:
        stem = None
    if stem is None:
        stem = _hantzsch_widman(elements, saturated=True)
    if stem is None:
        raise UnsupportedStructure("this saturated heteromonocycle has no supported name as a diyl yet")
    return _with_hetero_locants(stem, elements, hetero) if stem[0].isalpha() else stem


def _with_hetero_locants(stem, elements, hetero):
    names = [e for e in elements if e != "C"]
    if len(names) < 2 or elements.count("C") <= 1:
        return stem
    pairs = sorted(zip(names, [p for p, _ in hetero]), key=lambda ep: (_RANK[ep[0]], ep[1]))
    return ",".join(str(p) for _, p in pairs) + "-" + stem


def monocycle_numberings(mol, ring_order, attached, valence, ene_bonds_getter=None):
    ring_set = set(ring_order)
    atoms = [mol.GetAtomWithIdx(i) for i in ring_order]
    size = len(ring_order)
    if all(a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in atoms) and size == 6:
        out = []
        for walk in _walks(ring_order):
            position_of = {a: i + 1 for i, a in enumerate(walk)}
            out.append(Numbering(position_of, _benzene_text))
        return out
    if all(a.GetAtomicNum() == 6 for a in atoms):
        if any(a.GetIsAromatic() for a in atoms):
            raise UnsupportedStructure("an aromatic carbocycle other than benzene is not supported as a diyl yet")
        ring_bonds = [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set and b.GetBondTypeAsDouble() != 1.0
        ]
        if any(order != 2.0 for _, _, order in ring_bonds):
            raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
        stem = "cyclo" + alkane_name(size)
        out = []
        for walk in _walks(ring_order):
            position_of = {a: i + 1 for i, a in enumerate(walk)}
            ene = tuple(ring_bond_locants(position_of, ring_bonds, size)[0])
            out.append(Numbering(position_of, _carbocycle_text(stem, ene), unsat_key=ene))
        return out
    best_pre, results = _hetero_monocycle(mol, ring_order, attached)
    out = []
    for position_of, stem, ih, hydro in results:
        out.append(
            Numbering(
                position_of,
                _hetero_text(stem, ih, hydro),
                pre_key=best_pre + (ih,),
                unsat_key=hydro,
                ih=ih,
            )
        )
    return out


def _benzene_text(locants, valence, substituted=frozenset(), suffix="yl"):
    loc = _locs(locants)
    if suffix == "carboxylate":
        return "benzoate"
    if valence == 1:
        return "phenyl"
    return f"{loc}-phenylene" if valence == 2 else f"benzene-{loc}-{_yl(valence)}"


def _carbocycle_text(stem, ene):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        if suffix == "carboxylate" and not ene and not substituted:
            return f"{stem}carboxylate"
        if not ene:
            return "cyclo" + alkyl_name(_ring_size(stem)) if valence == 1 and suffix == "yl" else _tail(stem, locants, valence, suffix)
        base = stem[:-3] + ("a" if len(ene) > 1 else "")
        return _tail(f"{base}-{','.join(map(str, ene))}-{multiplied_word(len(ene), 'ene')}", locants, valence, suffix)

    return text


def _hetero_text(stem, ih, hydro):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        shown_ih = [p for p in ih if hydro or not (p == 1 and p in substituted)]
        ih_text = ",".join(f"{p}H" for p in shown_ih) + "-" if shown_ih else ""
        hydro_text = f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" + ("-" if ih_text else "") if hydro else ""
        return hydro_text + ih_text + _tail(stem, locants, valence, suffix)

    return text


import re


def _bare_skeleton(mol, skeleton_atoms, mancude=False):
    """Skeleton-only copy of the ring system; `mancude` makes every atom/bond aromatic (the mancude parent)."""
    order = sorted(skeleton_atoms)
    rw = Chem.RWMol()
    new_of = {}
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        copy = Chem.Atom(atom.GetAtomicNum())
        copy.SetIsAromatic(True if mancude else atom.GetIsAromatic())
        copy.SetFormalCharge(atom.GetFormalCharge())
        new_of[old] = rw.AddAtom(copy)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in new_of and b in new_of:
            rw.AddBond(new_of[a], new_of[b], Chem.BondType.AROMATIC if mancude else bond.GetBondType())
    bare = rw.GetMol()
    for old in order:
        atom = mol.GetAtomWithIdx(old)
        if not mancude and atom.GetIsAromatic() and atom.GetAtomicNum() != 6:
            heavy_inside = sum(1 for n in atom.GetNeighbors() if n.GetIdx() in new_of)
            if atom.GetTotalNumHs() > 0 or atom.GetDegree() > heavy_inside or (atom.GetAtomicNum() == 7 and atom.GetDegree() == 3):
                bare.GetAtomWithIdx(new_of[old]).SetNumExplicitHs(1)
    if not mancude:
        Chem.SanitizeMol(bare)
    return bare, new_of, {v: k for k, v in new_of.items()}


def _sanitized_mancude(mol, skeleton_atoms, sp3):
    """Skeleton-only mancude parent: every atom aromatic, with one explicit [nH] or one sp3 CH2 where needed."""
    bare, new_of, old_of = _bare_skeleton(mol, skeleton_atoms, mancude=True)
    nitrogens = [a for a in sorted(skeleton_atoms) if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    attempts = [("none", None)] + [("nh", a) for a in nitrogens] + [("ch2", a) for a in sorted(sp3) if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    for mode, atom_idx in attempts:
        trial = Chem.RWMol(bare)
        if mode == "nh":
            trial.GetAtomWithIdx(new_of[atom_idx]).SetNumExplicitHs(1)
        elif mode == "ch2":
            target = new_of[atom_idx]
            trial.GetAtomWithIdx(target).SetIsAromatic(False)
            for bond in list(trial.GetAtomWithIdx(target).GetBonds()):
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
        try:
            Chem.SanitizeMol(trial)
        except Exception:
            continue
        return trial.GetMol(), new_of, old_of
    raise UnsupportedStructure("the mancude parent of this ring system could not be built")


def _hydro_text(hydro):
    return f"{_locs(hydro)}-{multiplied_word(len(hydro), 'hydro')}" if hydro else ""


def _vb_text(parent, ene_citations, hetero_prefix=""):
    def text(locants, valence, substituted=frozenset(), suffix="yl"):
        stem = parent
        if ene_citations:
            body, needs_a = _unsaturation_suffix_from_citations(ene_citations, [])
            stem = parent[:-3] + ("a" if needs_a else "") + "-" + body
        return _tail(hetero_prefix + stem, locants, valence, suffix)

    return text


def _replacement_prefix(bare, order_new):
    """('2-oxa-5-aza' style text, hetero locants, element ranks) for a von Baeyer numbering."""
    by_element = {}
    for i, atom_idx in enumerate(order_new):
        symbol = bare.GetAtomWithIdx(atom_idx).GetSymbol()
        if symbol != "C":
            by_element.setdefault(symbol, []).append(i + 1)
    if not by_element:
        return "", (), ()
    if any(e not in _PREFIX for e in by_element):
        raise UnsupportedStructure("an unsupported skeletal heteroatom in a von Baeyer/spiro skeleton")
    pieces = []
    locants = []
    ranks = []
    for element in sorted(by_element, key=lambda e: _RANK[e]):
        locs = sorted(by_element[element])
        mult = "" if len(locs) == 1 else multiplying_prefix(len(locs))
        pieces.append(f"{','.join(map(str, locs))}-{mult}{_PREFIX[element]}")
        locants.extend(locs)
        ranks.extend([_RANK[element]] * len(locs))
    pairs = sorted(zip(locants, ranks))
    return "-".join(pieces), tuple(p for p, _ in pairs), tuple(r for _, r in pairs)


def _von_baeyer(mol, skeleton_atoms):
    bare, new_of, old_of = _bare_skeleton(mol, skeleton_atoms)
    if any(a.GetIsAromatic() for a in bare.GetAtoms()):
        raise UnsupportedStructure("an aromatic von Baeyer/spiro skeleton is not supported as a diyl yet")
    bonds_new = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
        for b in bare.GetBonds()
        if b.GetBondTypeAsDouble() != 1.0
    ]
    if any(order not in (2.0, 3.0) for _, _, order in bonds_new):
        raise UnsupportedStructure("an unsupported ring bond order in a polycyclic diyl")
    ring_count = bare.GetNumBonds() - bare.GetNumAtoms() + 1
    candidates = []
    spiro = find_monospiro_atom(bare) if ring_count == 2 else None
    if spiro is not None:
        for parent, order in iter_monospiro_numberings(bare, spiro):
            candidates.append((parent, order, ()))
    elif ring_count == 2:
        core = find_bicyclic_core(bare)
        if core is None:
            raise UnsupportedStructure("this bicyclic skeleton is not supported as a diyl yet")
        parent = bicyclic_parent_name(core)
        for order in iter_bicyclic_numberings(core):
            candidates.append((parent, order, ()))
    else:
        core = find_polycyclic_core(bare, ring_count)
        if core is None:
            raise UnsupportedStructure("this polycyclic skeleton is not supported as a diyl yet")
        for order, parent, outer in iter_polycyclic_candidates(core, ring_count):
            candidates.append((parent, order, outer))
    if not candidates:
        raise UnsupportedStructure("this ring skeleton has no supported von Baeyer/spiro numbering as a diyl yet")
    out = []
    for parent, order, outer in candidates:
        position_new = {a: i + 1 for i, a in enumerate(order)}
        position_of = {old_of[a]: p for a, p in position_new.items()}
        prefix, hetero_locs, hetero_ranks = _replacement_prefix(bare, order)
        pre = tuple(outer) + (hetero_locs, hetero_ranks)
        if bonds_new:
            ene, yne, compound_count, primary, full = von_baeyer_unsaturation_citations(position_new, bonds_new)
            if yne:
                raise UnsupportedStructure("a ring triple bond is not supported as a diyl yet")
            unsat = (compound_count, tuple(sorted(primary)), tuple(sorted(full)))
            out.append(Numbering(position_of, _vb_text(parent, ene, prefix), pre_key=pre, unsat_key=unsat))
        else:
            out.append(Numbering(position_of, _vb_text(parent, [], prefix), pre_key=pre))
    return out


def _sp3_ring_atoms(mol, skeleton_atoms):
    """Skeleton atoms that are neither aromatic nor in a skeleton double bond."""
    sp3 = set()
    for a in skeleton_atoms:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetIsAromatic():
            continue
        if any(
            mol.GetBondBetweenAtoms(a, n.GetIdx()).GetBondTypeAsDouble() != 1.0 and n.GetIdx() in skeleton_atoms
            for n in atom.GetNeighbors()
        ):
            continue
        sp3.add(a)
    return sp3


def _arene_chain(mol, graph, rings, skeleton_atoms):
    atom_rings = [tuple(r) for r in rings]
    ring_atom_sets = [set(r) for r in atom_rings]
    fusion_bond_idxs = set()
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if sum(1 for s in ring_atom_sets if a in s and b in s) >= 2:
            fusion_bond_idxs.add(bond.GetIdx())
    membership = _aromatic._atom_ring_membership(atom_rings)
    if any(len(r) >= 3 for r in membership.values()):
        raise UnsupportedStructure("a peri-fused arene is not supported as a diyl yet")
    adj, pairs = _aromatic._ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    ring_order = _aromatic._ring_path_order(adj, len(atom_rings))
    shape = _aromatic._classify_shape(graph, atom_rings, ring_order, pairs)
    parent = _aromatic._retained_chain_name(len(atom_rings), shape)
    if parent is None:
        return None
    if parent == "anthracene":
        candidates = _aromatic._anthracene_candidates(graph, ring_atom_sets, pairs, ring_order)
    elif parent == "phenanthrene":
        candidates = _aromatic._phenanthrene_candidates(graph, ring_atom_sets, pairs, ring_order)
    else:
        candidates = _aromatic._straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)
    sp3 = _sp3_ring_atoms(mol, skeleton_atoms)
    if len(sp3) % 2:
        raise UnsupportedStructure("an odd number of saturated atoms needs indicated hydrogen, not supported yet")
    out = []
    for c in candidates:
        position_of = dict(c)
        if any(a not in position_of for a in sp3):
            raise UnsupportedStructure("a saturated ring-fusion atom needs lettered hydro locants, not supported yet")
        hydro = tuple(sorted(position_of[a] for a in sp3))

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro):
            return _hydro_text(hydro) + _tail(parent, locants, valence, suffix)

        out.append(Numbering(position_of, text, unsat_key=hydro))
    return out


def _fused_mancude(mol, skeleton_atoms):
    from .core import smiles_to_iupac

    sp3 = _sp3_ring_atoms(mol, skeleton_atoms)
    bare, new_of, old_of = _sanitized_mancude(mol, skeleton_atoms, sp3)
    parent = smiles_to_iupac(Chem.MolToSmiles(bare))
    if any(x in parent for x in _LACKS_REGULAR_NUMBERING):
        raise UnsupportedStructure("this fused parent has a retained non-peripheral numbering not supported yet")
    match = re.match(r"^(\d+H(?:,\d+H)*)-(.*)$", parent)
    ih_count = len(match.group(1).split(",")) if match else 0
    stem = match.group(2) if match else parent
    numberings = general_peripheral_numberings(bare, ignore_indicated=True)
    if not numberings:
        raise UnsupportedStructure("this fused skeleton has no supported peripheral numbering as a diyl yet")
    out = []
    for numbering in numberings:
        position_of = {old_of[new]: int(loc) for new, loc in numbering.items() if loc.isdigit()}
        if any(a not in position_of for a in sp3):
            continue
        saturated = sorted(
            position_of[a]
            for a in skeleton_atoms
            if a in position_of
            and (
                a in sp3
                or (
                    mol.GetAtomWithIdx(a).GetIsAromatic()
                    and mol.GetAtomWithIdx(a).GetAtomicNum() != 6
                    and (mol.GetAtomWithIdx(a).GetTotalNumHs() > 0 or mol.GetAtomWithIdx(a).GetDegree() == 3)
                )
            )
        )
        if len(saturated) < ih_count or (len(saturated) - ih_count) % 2:
            continue
        ih = tuple(saturated[:ih_count])
        hydro = tuple(saturated[ih_count:])
        ih_text = (",".join(f"{p}H" for p in ih) + "-") if ih else ""

        def text(locants, valence, substituted=frozenset(), suffix="yl", hydro=hydro, ih_text=ih_text):
            return _hydro_text(hydro) + ("-" if hydro and ih_text else "") + ih_text + _tail(stem, locants, valence, suffix)

        out.append(Numbering(position_of, text, pre_key=(ih,), unsat_key=hydro))
    if not out:
        raise UnsupportedStructure("no numbering of this partly hydrogenated fused system fits its hydro/indicated hydrogen")
    return out


def system_numberings(mol, graph, rings, skeleton_atoms):
    if all(
        len(r) == 6 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in r) for r in rings
    ) and any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        arene = _arene_chain(mol, graph, rings, skeleton_atoms)
        if arene is not None:
            return arene
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in skeleton_atoms):
        return _fused_mancude(mol, skeleton_atoms)
    return _von_baeyer(mol, skeleton_atoms)


def chain_numberings(mol, paths, valence):
    from ._common import name_from_substituents

    out = []
    for path in paths:
        path_set = set(path)
        for atom in path:
            for n in mol.GetAtomWithIdx(atom).GetNeighbors():
                bond = mol.GetBondBetweenAtoms(atom, n.GetIdx())
                if n.GetIdx() not in path_set and bond.GetBondTypeAsDouble() != 1.0 and n.GetAtomicNum() == 6:
                    raise UnsupportedStructure("an ylidene-type multiple bond off the polyol chain is not supported yet")
        for direction in (path, list(reversed(path))):
            position_of = {a: i + 1 for i, a in enumerate(direction)}
            ene, yne = [], []
            for i in range(len(direction) - 1):
                bond = mol.GetBondBetweenAtoms(direction[i], direction[i + 1])
                order = bond.GetBondTypeAsDouble()
                if order == 2.0:
                    ene.append(i + 1)
                elif order == 3.0:
                    yne.append(i + 1)
                elif order != 1.0:
                    raise UnsupportedStructure("an aromatic or unusual bond on the polyol chain")
            length = len(direction)

            def text(locants, valence, substituted=frozenset(), suffix="yl", ene=tuple(ene), yne=tuple(yne), length=length):
                if length == 1 and valence == 2:
                    return "methylene"
                if valence == 1 and not ene and not yne:
                    return alkyl_name(length) if locants[0] == 1 else name_from_substituents(length, [], [], "yl", own_locants=list(locants))
                return name_from_substituents(length, list(ene), list(yne), _yl(valence), own_locants=sorted(locants))

            out.append(Numbering(position_of, text, unsat_key=(tuple(ene), tuple(yne))))
    return out
