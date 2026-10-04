"""Retained name 'sphinganine' and its unsaturated, N-substituted and O-substituted derivatives (P-107.4.3.1).

An unbranched C18 chain H2N-C(2), HO-C(1) and HO-C(3), (2S,3R), is 'sphinganine'; double bonds from C-4 onward
give 'sphing-4-enine', 'sphinga-4,14-dienine' with the double-bond descriptors cited first. Alkyl, glycosyl and
similar groups on N, O-1 or O-3 are cited as 'N-' and 'n-O-' prefixes; an acyl group makes an amide or ester that is
named systematically. Other chain lengths and other diastereoisomers are named systematically.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._numerals import multiplying_prefix
from ._substituents import BRANCH_STEREO, is_plain_stem_prefix, name_branch

_CHAIN_LENGTH = 18


def _carbon_neighbors(mol, graph, atom):
    return [n for n in graph[atom] if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]


def _hetero_neighbor(mol, graph, carbon, element):
    found = [n for n in graph[carbon] if mol.GetAtomWithIdx(n).GetAtomicNum() == element]
    return found[0] if len(found) == 1 else None


def _chain_from(mol, graph, c2):
    """C-1..C-18 of an unbranched chain through the amino-bearing carbon C-2, or None."""
    carbons = _carbon_neighbors(mol, graph, c2)
    if len(carbons) != 2:
        return None
    for c1, c3 in (carbons, carbons[::-1]):
        if mol.GetAtomWithIdx(c1).GetTotalNumHs() != 2 or len(_carbon_neighbors(mol, graph, c1)) != 1:
            continue
        chain = [c1, c2, c3]
        while len(chain) < _CHAIN_LENGTH:
            onward = [n for n in _carbon_neighbors(mol, graph, chain[-1]) if n != chain[-2]]
            if len(onward) != 1:
                break
            chain.append(onward[0])
        if len(chain) != _CHAIN_LENGTH or any(a.IsInRing() for a in map(mol.GetAtomWithIdx, chain)):
            continue
        if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for a in chain[3:] for n in graph[a]):
            continue
        if any(n not in chain for n in _carbon_neighbors(mol, graph, chain[-1])):
            continue
        return chain
    return None


def _skeleton(mol, graph):
    """(chain, nitrogen, O-1, O-3) of the sphingoid skeleton, or None."""
    for nitrogen in mol.GetAtoms():
        if nitrogen.GetAtomicNum() != 7 or nitrogen.IsInRing():
            continue
        for c2 in _carbon_neighbors(mol, graph, nitrogen.GetIdx()):
            chain = _chain_from(mol, graph, c2)
            if chain is None:
                continue
            o1 = _hetero_neighbor(mol, graph, chain[0], 8)
            o3 = _hetero_neighbor(mol, graph, chain[2], 8)
            if o1 is None or o3 is None or _hetero_neighbor(mol, graph, c2, 7) != nitrogen.GetIdx():
                continue
            if any(mol.GetAtomWithIdx(n).GetAtomicNum() not in (6, 7, 8) for a in chain[:3] for n in graph[a]):
                continue
            return chain, nitrogen.GetIdx(), o1, o3
    return None


def _acyl_like(mol, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    return atom.GetAtomicNum() in (6, 15, 16) and any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (8, 16) for b in atom.GetBonds()
    )


def _prefix_text(substituents):
    """'N,N-dimethyl', '1-O-methyl', '1,3-di-O-methyl' pieces in alphanumerical order (P-14.5.2, P-16.3)."""
    grouped = {}
    for locant, entries in substituents.items():
        for name, compound in entries:
            info = grouped.setdefault(name, {"compound": compound, "locants": []})
            info["locants"].append(locant)
    pieces = []
    for name in sorted(grouped, key=alpha_sort_key):
        info = grouped[name]
        locants = sorted(info["locants"], key=lambda loc: (loc != "N", loc))
        kinds = {loc == "N" for loc in locants}
        if len(kinds) > 1:
            raise UnsupportedStructure("the same group on both nitrogen and oxygen is not supported yet")
        compound = info["compound"] and not is_plain_stem_prefix(name)
        multiplier = multiplying_prefix(len(locants), compound=compound) if len(locants) > 1 else ""
        cited = f"({name})" if compound else name
        if locants[0] == "N":
            pieces.append(f"{','.join(locants)}-{multiplier}{cited}")
        else:
            pieces.append(f"{','.join(locants)}-{multiplier}{'-' if multiplier else ''}O-{cited}")
    return "-".join(pieces)


def _name(mol):
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    found = _skeleton(mol, graph)
    if found is None:
        return None
    chain, nitrogen, o1, o3 = found
    c1, c2, c3 = chain[:3]
    doubles = []
    for position in range(_CHAIN_LENGTH - 1):
        bond = mol.GetBondBetweenAtoms(chain[position], chain[position + 1])
        order = bond.GetBondTypeAsDouble()
        if order == 2.0 and position >= 3:
            doubles.append((position + 1, bond))
        elif order != 1.0:
            return None
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    if (
        probe.GetAtomWithIdx(c2).GetPropsAsDict().get("_CIPCode") != "S"
        or probe.GetAtomWithIdx(c3).GetPropsAsDict().get("_CIPCode") != "R"
    ):
        return None
    geometry = [
        f"{locant}{probe.GetBondWithIdx(bond.GetIdx()).GetPropsAsDict()['_CIPCode']}"
        for locant, bond in doubles
        if probe.GetBondWithIdx(bond.GetIdx()).HasProp("_CIPCode")
    ]

    skeleton = set(chain) | {nitrogen, o1, o3}
    context = {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.HasProp("_CIPCode") and a.GetIdx() not in skeleton},
        "bonds": {},
        "used": set(),
    }
    substituents = {}
    covered = set(skeleton)
    halogens = halogen_substituents(mol)
    token = BRANCH_STEREO.set(context)
    try:
        for locant, heteroatom, carbon in (("N", nitrogen, c2), ("1", o1, c1), ("3", o3, c3)):
            for neighbor in graph[heteroatom]:
                if neighbor == carbon:
                    continue
                if _acyl_like(mol, neighbor):
                    return None
                try:
                    entry = name_branch(graph, neighbor, heteroatom, halogens, mol=mol)
                except UnsupportedStructure:
                    return None
                substituents.setdefault(locant, []).append(entry)
                stack = [neighbor]
                while stack:
                    atom = stack.pop()
                    if atom in covered:
                        continue
                    covered.add(atom)
                    stack.extend(n for n in graph[atom] if n != heteroatom)
    finally:
        BRANCH_STEREO.reset(token)
    if len(covered) != mol.GetNumAtoms() or any(("atom", a) not in context["used"] for a in context["atoms"]):
        return None

    count = len(doubles)
    if not count:
        stem = "sphinganine"
    else:
        multiple = "en" if count == 1 else multiplying_prefix(count) + "en"
        stem = f"{'sphing' if count == 1 else 'sphinga'}-{','.join(str(n) for n, _ in doubles)}-{multiple}ine"
    prefix = _prefix_text(substituents)
    return (f"({','.join(geometry)})-" if geometry else "") + prefix + stem


def has_sphingoid_shape(mol) -> bool:
    return _name(mol) is not None


def name_sphingoid(mol) -> str:
    return _name(mol)
