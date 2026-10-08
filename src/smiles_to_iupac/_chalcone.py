"""The retained name chalcone (P-64.2.1.1).

Chalcone is (2E)-1,3-diphenylprop-2-en-1-one with ring substitution only, and only by groups cited as prefixes
(characteristic groups lower than ketone). The ring on the carbonyl carbon carries the primed locants, the ring on the
alkene carbon the unprimed ones; each ring is numbered from its point of attachment in the direction that gives its
substituents the lowest locants.
"""

from rdkit import Chem

from ._cited_group import cited_group, subtree
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, numerical_term

_CHALCONE = Chem.MolFromSmarts("O=[CX3]([c;R1])[CX3;H1]=[CX3;H1][c;R1]")
_COMPLEX = {2: "bis", 3: "tris", 4: "tetrakis", 5: "pentakis"}
_PRIME = "′"


def _plain_prefix_group(mol, graph, root, ring_atom):
    """True when the group is named as a prefix beside a ketone parent: no carbonyl-type, nitrile, sulfonyl or other
    group that outranks the ketone."""
    for a in subtree(graph, root, ring_atom):
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() not in (6, 7, 8, 9, 16, 17, 35, 53) or atom.GetIsotope() or atom.GetNumRadicalElectrons():
            return False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if atom.GetAtomicNum() == 6 and bond.GetBondTypeAsDouble() > 1.0 and not bond.GetIsAromatic() and other.GetAtomicNum() in (7, 8, 16):
                return False
        if atom.GetAtomicNum() == 16 and any(b.GetBondTypeAsDouble() > 1.0 for b in atom.GetBonds()):
            return False
        if atom.GetFormalCharge() and not any(n.GetFormalCharge() == -atom.GetFormalCharge() for n in atom.GetNeighbors()):
            return False
    return True


def _ring_locants(mol, ring, attach):
    """[(substituted ring atom, locant)] for the best of the two directions, with the substituent roots."""
    start = ring.index(attach)
    options = []
    for step in (1, -1):
        numbers = {ring[(start + step * k) % 6]: k + 1 for k in range(6)}
        options.append(numbers)
    subs = {}
    for atom in ring:
        for n in mol.GetAtomWithIdx(atom).GetNeighbors():
            if n.GetIdx() not in ring and atom != attach:
                subs[atom] = n.GetIdx()
    best = min(options, key=lambda numbers: sorted(numbers[a] for a in subs))
    return best, subs


def _prefix_segment(locants, name, compound):
    count = len(locants)
    ordered = sorted(locants, key=lambda item: (item[0], item[1]))
    text = ",".join(f"{n}{_PRIME if primed else ''}" for n, primed in ordered)
    multiplier = "" if count == 1 else (_COMPLEX[count] if compound else numerical_term(count))
    group = f"({name})" if compound else name
    return f"{text}-{multiplier}{group}"


def chalcone_name(mol):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() and a.GetAtomicNum() not in (7, 8) for a in mol.GetAtoms()):
        return None
    matches = mol.GetSubstructMatches(_CHALCONE)
    if len(matches) != 1:
        return None
    oxygen, carbonyl, ring_a_atom, alpha, beta, ring_b_atom = matches[0]
    double = mol.GetBondBetweenAtoms(alpha, beta)
    if double.GetStereo() not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS):
        return None
    rings = mol.GetRingInfo().AtomRings()
    ring_a = next((r for r in rings if ring_a_atom in r and len(r) == 6), None)
    ring_b = next((r for r in rings if ring_b_atom in r and len(r) == 6), None)
    if ring_a is None or ring_b is None or set(ring_a) & set(ring_b):
        return None
    if any(not all(mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring) for ring in (ring_a, ring_b)):
        return None
    graph = adjacency(mol)
    in_rings = set(ring_a) | set(ring_b)
    chalcone_core = in_rings | {oxygen, carbonyl, alpha, beta}
    if any(n.GetIdx() not in chalcone_core for a in (oxygen, carbonyl, alpha, beta) for n in mol.GetAtomWithIdx(a).GetNeighbors()):
        return None
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE and b.GetIdx() != double.GetIdx() for b in mol.GetBonds()):
        return None
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        return None
    groups = {}
    try:
        for ring, attach, primed in ((ring_a, ring_a_atom, True), (ring_b, ring_b_atom, False)):
            numbers, subs = _ring_locants(mol, list(ring), attach)
            for atom, root in subs.items():
                if not _plain_prefix_group(mol, graph, root, atom):
                    return None
                name, compound = cited_group(mol, graph, root, atom)
                groups.setdefault((name, compound), []).append((numbers[atom], primed))
    except UnsupportedStructure:
        return None
    segments = [
        _prefix_segment(locants, name, compound)
        for (name, compound), locants in sorted(groups.items(), key=lambda item: alpha_sort_key(item[0][0]))
    ]
    return "-".join(segments) + "chalcone"


def has_chalcone_shape(mol) -> bool:
    return chalcone_name(mol) is not None


def name_chalcone(mol) -> str:
    return chalcone_name(mol)
