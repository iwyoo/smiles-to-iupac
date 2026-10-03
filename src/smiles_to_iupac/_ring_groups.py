"""Substituents and principal characteristic groups on a metallacycle ring
(P-69.4 parents are skeletal-replacement hydrides, so -one/-ol/-amine are
cited as suffixes in that seniority order and the rest as prefixes).
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, halogen_substituents, plain_phenyl_substituent_atoms
from ._numerals import multiplying_prefix
from ._prefix_groups import PrefixNamer, _oxy
from ._substituents import branch_atom_locant, format_mononuclear_prefixes, name_branch

KINDS = ("acid", "one", "ol", "amine")
_PREFIX_FORM = {"acid": "carboxy", "one": "oxo", "ol": "hydroxy", "amine": "amino"}


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


def _amino(mol, graph, ring_atom, n):
    subs = [name_branch(graph, x, n, {}, mol=mol) for x in graph[n] if x != ring_atom]
    return (format_mononuclear_prefixes(subs) + "amino", True) if subs else ("amino", False)


def n_substituents(mol, graph, ring_atom, n):
    return [name_branch(graph, x, n, {}, mol=mol) for x in graph[n] if x != ring_atom]


def ring_substituents(mol, graph, ring_atom, ring_set, principal):
    """(prefix entries, n-entries, principal count) for one ring atom;
    n-entries are N-substituents of a principal amine (locant 'N')."""
    entries, n_entries, count = [], [], 0
    namer = PrefixNamer(mol, graph)
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
                try:
                    entries.append(_branch_with_stereo(mol, graph, n, ring_atom))
                except UnsupportedStructure:
                    if _has_stereo(mol, graph, n, ring_atom):
                        raise
                    entries.append(namer.name(n, ring_atom))
        elif z == 8 and atom.GetDegree() == 2 and atom.GetTotalNumHs() == 0:
            carbon = next(x for x in graph[n] if x != ring_atom)
            if mol.GetAtomWithIdx(carbon).GetAtomicNum() != 6:
                raise UnsupportedStructure("this ring substituent is not supported here")
            entries.append(_oxy(*name_branch(graph, carbon, n, {}, mol=mol)))
        else:
            entries.append(namer.name(n, ring_atom))
    return entries, n_entries, count


def _has_stereo(mol, graph, n, ring_atom):
    group, stack = {n}, [n]
    while stack:
        for v in graph[stack.pop()]:
            if v != ring_atom and v not in group:
                group.add(v)
                stack.append(v)
    hetero_multiple = any(
        {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= group
        and b.GetBondTypeAsDouble() >= 2.0
        and any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in mol.GetBonds()
    )
    return hetero_multiple or any(mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for i in group)


def _branch_with_stereo(mol, graph, n, ring_atom):
    halogens = {**halogen_substituents(mol), **PrefixNamer(mol, graph)._region_halogens(n, ring_atom)}
    name, compound = name_branch(graph, n, ring_atom, halogens, mol=mol)
    group, stack = {n}, [n]
    while stack:
        for v in graph[stack.pop()]:
            if v != ring_atom and v not in group:
                group.add(v)
                stack.append(v)
    for bond in mol.GetBonds():
        if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} <= group and bond.GetBondTypeAsDouble() >= 2.0:
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())):
                raise UnsupportedStructure("a characteristic group inside a ring substituent is out of scope here")
    specified = [i for i in group if mol.GetAtomWithIdx(i).GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
    if not specified:
        return name, compound
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    cited = []
    for i in specified:
        atom = probe.GetAtomWithIdx(i)
        if not atom.HasProp("_CIPCode"):
            raise UnsupportedStructure("the stereocentre of this substituent has no R/S descriptor")
        cited.append((branch_atom_locant(graph, n, ring_atom, i, halogens, mol=mol), atom.GetProp("_CIPCode")))
    cited.sort()
    return f"({','.join(f'{l}{c}' for l, c in cited)})-{name}", True


def add_n_entries(grouped, n_entries):
    out = {k: {**v, "locants": list(v["locants"])} for k, v in grouped.items()}
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
    """True for an exocyclic C=O or C=C(H..) bond on a ring carbon (oxo /
    ylidene substituent); the other end must not be in `ring_set`."""
    ring_end, other = (a, b) if a in ring_set else (b, a)
    atom = mol.GetAtomWithIdx(other)
    if other in ring_set or mol.GetAtomWithIdx(ring_end).GetAtomicNum() != 6:
        return False
    if mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 2.0:
        return False
    return (atom.GetAtomicNum() == 8 and atom.GetDegree() == 1) or (atom.GetAtomicNum() == 6 and not atom.GetIsAromatic())
