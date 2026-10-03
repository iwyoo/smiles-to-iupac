"""Substituents and principal characteristic groups on a metallacycle ring
(P-69.4 parents are skeletal-replacement hydrides, so -one/-ol/-amine are
cited as suffixes in that seniority order and the rest as prefixes).
"""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, plain_phenyl_substituent_atoms
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

KINDS = ("acid", "one", "ol", "amine")
_PREFIX_FORM = {"acid": "carboxy", "one": "oxo", "ol": "hydroxy", "amine": "amino"}
_CONTRACTED = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy", "phenyl": "phenoxy"}


def _is_carboxyl(mol, n, ring_atom):
    atom = mol.GetAtomWithIdx(n)
    if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetDegree() != 3 or atom.GetTotalNumHs():
        return False
    others = [x for x in atom.GetNeighbors() if x.GetIdx() != ring_atom]
    oxo = [x for x in others if x.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n, x.GetIdx()).GetBondTypeAsDouble() == 2.0]
    hydroxy = [x for x in others if x.GetAtomicNum() == 8 and x.GetTotalNumHs() == 1 and x.GetDegree() == 1]
    return len(oxo) == 1 and len(hydroxy) == 1


def is_carboxyl_bond(mol, a, b):
    """True for the C=O bond of a ring-attached carboxylic acid group."""
    carbon, oxygen = (a, b) if mol.GetAtomWithIdx(a).GetAtomicNum() == 6 else (b, a)
    return mol.GetAtomWithIdx(oxygen).GetAtomicNum() == 8 and any(
        _is_carboxyl(mol, carbon, x.GetIdx()) for x in mol.GetAtomWithIdx(carbon).GetNeighbors() if x.GetIdx() != oxygen
    )


def group_kind(mol, ring_atom, n):
    atom = mol.GetAtomWithIdx(n)
    z = atom.GetAtomicNum()
    if _is_carboxyl(mol, n, ring_atom):
        return "acid"
    if z == 8 and atom.GetDegree() == 1 and atom.GetFormalCharge() == 0:
        if mol.GetBondBetweenAtoms(ring_atom, n).GetBondTypeAsDouble() == 2.0:
            return "one"
        if atom.GetTotalNumHs() == 1:
            return "ol"
    if z == 7 and atom.GetFormalCharge() == 0 and not atom.GetIsAromatic() and atom.GetDegree() <= 3:
        if all(b.GetBondTypeAsDouble() == 1.0 for b in atom.GetBonds()) and all(
            x.GetAtomicNum() == 6 for x in atom.GetNeighbors()
        ):
            return "amine"
    return None


def principal_kind(mol, graph, ring_set, skip):
    found = {
        group_kind(mol, a, n)
        for a in ring_set - {skip}
        for n in graph[a]
        if n not in ring_set
    }
    return next((k for k in KINDS if k in found), None)


def _oxy(name, compound):
    if name in _CONTRACTED:
        return _CONTRACTED[name], False
    if name.endswith("phenyl") and not compound:
        return name[: -len("phenyl")] + "phenoxy", True
    if name.endswith(("methyl", "ethyl", "propyl", "butyl")) and "cyclo" not in name:
        return name[:-2] + "oxy", compound or any(ch.isdigit() for ch in name)
    return name + "oxy", compound or any(ch.isdigit() for ch in name)


def _amino(mol, graph, ring_atom, n):
    subs = [name_branch(graph, x, n, {}, mol=mol) for x in graph[n] if x != ring_atom]
    return (format_mononuclear_prefixes(subs) + "amino", True) if subs else ("amino", False)


def n_substituents(mol, graph, ring_atom, n):
    return [name_branch(graph, x, n, {}, mol=mol) for x in graph[n] if x != ring_atom]


def ring_substituents(mol, graph, ring_atom, ring_set, principal):
    """(prefix entries, n-entries, principal count) for one ring atom;
    n-entries are N-substituents of a principal amine (locant 'N')."""
    entries, n_entries, count = [], [], 0
    for n in graph[ring_atom]:
        if n in ring_set:
            continue
        atom = mol.GetAtomWithIdx(n)
        z = atom.GetAtomicNum()
        kind = group_kind(mol, ring_atom, n)
        if kind is not None and kind == principal:
            count += 1
            if kind == "amine":
                n_entries.extend(n_substituents(mol, graph, ring_atom, n))
        elif kind == "amine":
            entries.append(_amino(mol, graph, ring_atom, n))
        elif kind is not None:
            entries.append((_PREFIX_FORM[kind], False))
        elif z in HALOGEN_PREFIXES:
            entries.append((HALOGEN_PREFIXES[z], False))
        elif z == 6:
            if plain_phenyl_substituent_atoms(mol, graph, {n}):
                entries.append(("phenyl", False))
            else:
                entries.append(name_branch(graph, n, ring_atom, {}, mol=mol))
        elif z == 8 and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 0:
            carbon = next(x for x in graph[n] if x != ring_atom)
            if mol.GetAtomWithIdx(carbon).GetAtomicNum() != 6:
                raise UnsupportedStructure("this ring substituent is not supported here")
            entries.append(_oxy(*name_branch(graph, carbon, n, {}, mol=mol)))
        else:
            raise UnsupportedStructure("this ring substituent is not supported here")
    return entries, n_entries, count


def add_n_entries(grouped, n_entries):
    out = {k: {"locants": list(v["locants"]), "compound": v["compound"]} for k, v in grouped.items()}
    for name, compound in n_entries:
        out.setdefault(name, {"locants": [], "compound": compound})["locants"].append("N")
    return out


def with_suffix(stem, locants, kind):
    """`stem` is a saturated/unsaturated ring name ending in 'e'."""
    if not locants:
        return stem
    word = kind
    loc = ",".join(str(x) for x in sorted(locants))
    if kind == "acid":
        return f"{stem}-{loc}-{multiplying_prefix(len(locants)) if len(locants) > 1 else ''}carboxylic acid"
    if len(locants) == 1:
        return f"{stem[:-1]}-{loc}-{word}"
    return f"{stem}-{loc}-{multiplying_prefix(len(locants))}{word}"


def is_exocyclic_oxo(mol, a, b, ring_set):
    """True for a C=O bond whose carbon is in `ring_set` and oxygen is not."""
    ring_end, other = (a, b) if a in ring_set else (b, a)
    atom = mol.GetAtomWithIdx(other)
    return (
        other not in ring_set
        and atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and mol.GetAtomWithIdx(ring_end).GetAtomicNum() == 6
    )
