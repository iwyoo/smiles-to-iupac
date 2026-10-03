"""Two phosphane groups joined by a carbon linker (P-15.3 multiplicative
nomenclature on the P-68 phosphane parent, `tmp/bluebook/P6a.txt` line
~8920): identical halves give 'methylenebis(dimethylphosphane)' and
'ethane-1,2-diylbis(diphenylphosphane)'; different halves are named
substitutively on the better-substituted phosphorus. Only alkyl/phenyl P
substituents and an unbranched alkane or benzene linker are supported.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds, plain_phenyl_substituent_atoms
from ._metal_pair import _brackets, _mark
from ._numerals import alkane_name, alkyl_name
from ._substituents import format_mononuclear_prefixes, name_branch


def has_polyphosphane_shape(mol) -> bool:
    phosphorus = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 15]
    if len(phosphorus) != 2:
        return False
    if any(n.GetAtomicNum() == 15 for p in phosphorus for n in p.GetNeighbors()):
        return False
    return all(a.GetAtomicNum() in (6, 15) for a in mol.GetAtoms())


def _linker(mol, graph, p1, p2):
    trails = []
    stack = [(p1, [p1])]
    while stack:
        cur, trail = stack.pop()
        if cur == p2:
            trails.append(trail)
            continue
        for nxt in graph[cur]:
            if nxt not in trail and (nxt == p2 or mol.GetAtomWithIdx(nxt).GetAtomicNum() == 6):
                stack.append((nxt, trail + [nxt]))
    if not trails:
        raise UnsupportedStructure("the two phosphorus atoms are not connected")
    trail = min(trails, key=len)
    chain = trail[1:-1]
    if all(not mol.GetAtomWithIdx(c).GetIsAromatic() for c in chain):
        for c in chain:
            extra = [n for n in graph[c] if n not in trail]
            if extra:
                raise UnsupportedStructure("a branched linker is not supported here")
        if any(mol.GetAtomWithIdx(c).IsInRing() for c in chain):
            raise UnsupportedStructure("a cyclic non-aromatic linker is not supported here")
        return chain, None
    ring = next((set(r) for r in mol.GetRingInfo().AtomRings() if set(chain) <= set(r)), None)
    if ring is None or len(ring) != 6 or len(chain) < 2 or len(chain) > 4:
        raise UnsupportedStructure("this aromatic linker is not supported here")
    if any(mol.GetRingInfo().NumAtomRings(i) != 1 for i in ring):
        raise UnsupportedStructure("a fused-ring linker is not supported here")
    if any(set(graph[i]) - ring - {p1, p2} for i in ring):
        raise UnsupportedStructure("a substituted aromatic linker is not supported here")
    return chain, ring


def _half(mol, graph, p, linker_atom, phenyl_atoms):
    entries = []
    for n in graph[p]:
        if n == linker_atom:
            continue
        if n in phenyl_atoms:
            entries.append(("phenyl", False))
        else:
            if mol.GetAtomWithIdx(n).GetIsAromatic():
                raise UnsupportedStructure("only plain phenyl and alkyl substituents are supported")
            entries.append(name_branch(graph, n, p, {}, mol=mol))
    return entries


def name_polyphosphane(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetFormalCharge() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    p1, p2 = (a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 15)
    graph = adjacency(mol)
    chain, ring = _linker(mol, graph, p1, p2)
    link1, link2 = chain[0], chain[-1]
    roots = {n for p in (p1, p2) for n in graph[p]}
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots - set(chain))
    covered = set(chain) | {p1, p2} | phenyl_atoms | (ring or set())
    halves = [_half(mol, graph, p, link, phenyl_atoms) for p, link in ((p1, link1), (p2, link2))]
    for entries, p in zip(halves, (p1, p2)):
        for n in graph[p]:
            if n not in chain and n not in phenyl_atoms:
                covered.update(_branch(graph, n, p))
    if covered != set(range(mol.GetNumAtoms())):
        raise UnsupportedStructure("this structure contains atoms outside the supported shapes")
    if any(
        not (mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(b).GetIsAromatic())
        for a, b, *_ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure("an unsaturated substituent is out of scope here")

    if ring is None:
        length = len(chain)
        linker_name = "methylene" if length == 1 else f"{alkane_name(length)}-1,{length}-diyl"
    else:
        span = len(chain) - 1
        linker_name = {1: "benzene-1,2-diyl", 2: "benzene-1,3-diyl", 3: "benzene-1,4-diyl"}[span]
    names = [format_mononuclear_prefixes(e) + "phosphane" if e else "phosphane" for e in halves]
    if names[0] == names[1]:
        return f"{linker_name}bis({names[0]})"

    order = sorted(range(2), key=lambda i: (-len(halves[i]), names[i]))
    parent, other = order
    other_group = format_mononuclear_prefixes(halves[other]) + "phosphanyl" if halves[other] else "phosphanyl"
    if ring is None:
        length = len(chain)
        text = f"({other_group})methyl" if length == 1 else f"{length}-({other_group}){alkyl_name(length)}"
    else:
        text = f"{len(chain)}-({other_group})phenyl"
    substituent = (_mark(text), True)
    return _brackets(format_mononuclear_prefixes(halves[parent] + [substituent]) + "phosphane")


def _branch(graph, root, boundary):
    atoms, stack = set(), [root]
    while stack:
        a = stack.pop()
        if a in atoms:
            continue
        atoms.add(a)
        stack.extend(x for x in graph[a] if x != boundary)
    return atoms
