"""Thiourea, selenourea and tellurourea (H2N-C(=E)-NH2, E = S/Se/Te) with any number of N-substituents
(P-66.1.4.1); unlike `_urea.py` no semicarbazide (amino-substituted nitrogen) handling. The N/N' citation helpers
are shared with `_urea.py`."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_nitro_nitrogen,
)
from ._substituents import alpha_sort_key, format_mononuclear_prefixes, format_substituent_prefixes, name_branch


def is_oxo_nitrogen(mol, atom):
    """A nitro or nitroso nitrogen, which is cited as a prefix on the nitrogen it is bonded to (P-61.5, P-61.11)."""
    if atom.GetAtomicNum() != 7:
        return False
    if is_nitro_nitrogen(mol, atom.GetIdx()):
        return True
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
    return (
        atom.GetDegree() == 2
        and len(oxygens) == 1
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), oxygens[0].GetIdx()).GetBondTypeAsDouble() == 2.0
    )


def is_core_substituent_root(mol, atom):
    """A carbon, a halogen or a nitro/nitroso nitrogen: the atoms a urea or guanidine nitrogen may carry."""
    return atom.GetAtomicNum() in (6, *HALOGEN_PREFIXES) or is_oxo_nitrogen(mol, atom)


_UREA_SUBSTITUENT_ELEMENTS = {6, 7, 8, 16, 34, 52, *HALOGEN_PREFIXES}


def is_urea_substituent_root(mol, atom):
    """An atom a urea nitrogen may carry: any element whose group is cited as a prefix, so that hydroxy, alkoxy,
    amino and sulfanyl groups on the nitrogen are ordinary substituents (P-66.1.6.1.1.2)."""
    return atom.GetAtomicNum() in _UREA_SUBSTITUENT_ELEMENTS


_SENIOR_TO_UREA = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,Se,Te])[OX2H1,OX1-]",
        "[CX3](=[O,S,Se,Te])[OX2][#6]",
        "[CX3](=[O,S,Se,Te])[SX2,SeX2,TeX2]",
        "[CX3](=[O,S,Se,Te])[F,Cl,Br,I]",
        "[CX3](=[O,S,Se,Te])[NX3]",
        "[CX3](=[NX2])[NX3]",
        "[SX4,SX3,SeX4,TeX4](=O)[OX2H1,OX1-,OX2,NX3]",
    )
]


def has_group_senior_to_urea(mol, core_atoms):
    """A carboxylic or sulfonic acid, ester, acid halide, amide or amidine outside the urea core, or a charge other
    than a nitro group's: the class that outranks urea as the parent (P-41, P-66.1.6.1.1.5)."""
    for query in _SENIOR_TO_UREA:
        if any(match[0] not in core_atoms for match in mol.GetSubstructMatches(query)):
            return True
    return any(
        a.GetFormalCharge() and a.GetIdx() not in core_atoms and not is_nitro_nitrogen(mol, a.GetIdx()) and not any(
            is_nitro_nitrogen(mol, n.GetIdx()) for n in a.GetNeighbors()
        )
        for a in mol.GetAtoms()
    )


def _n_substituent_roots(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx)


def _substituent_atoms(graph, roots, core_atoms):
    """Every atom reachable from `roots` without passing through the core, or
    None when a substituent loops back into it (a ring-fused urea)."""
    reached = set()
    stack = list(roots)
    while stack:
        idx = stack.pop()
        if idx in reached:
            continue
        reached.add(idx)
        stack.extend(n for n in graph[idx] if n not in core_atoms)
    for idx in reached:
        if sum(1 for n in graph[idx] if n in core_atoms) > (1 if idx in roots else 0):
            return None
    return reached


def n_substituent_names(mol, graph, core_atoms, nitrogens, carbon_idx, junior_groups=False):
    """One list of substituent names per nitrogen of `nitrogens`, for every atom outside `core_atoms`
    (P-66.1.4.1); raises when a ring fuses into the core or a non-halogen
    acyclic heteroatom is present, unless `junior_groups` lets groups junior to the amide be prefixes."""
    roots = [
        tuple(r for r in _n_substituent_roots(mol, n, carbon_idx) if r not in core_atoms) for n in nitrogens
    ]
    for group in roots:
        arms = [_substituent_atoms(graph, (r,), core_atoms) for r in group]
        if any(arm is None for arm in arms) or sum(len(arm) for arm in arms) != len(set().union(*arms)):
            raise UnsupportedStructure(
                "a ring through a nitrogen of the group is named as a ring parent, not as N-substituents"
            )
    outside = _substituent_atoms(graph, tuple(r for group in roots for r in group), core_atoms)
    if outside is None or len(outside) + len(core_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a ring-fused urea (e.g. hydantoin) or a characteristic group outside the urea core "
            "and its N-substituents is not supported yet"
        )
    halogens = halogen_substituents(mol)
    if junior_groups:
        if has_group_senior_to_urea(mol, core_atoms):
            raise UnsupportedStructure("a group senior to the urea is the parent, not an N-substituent")
        if any(
            mol.GetAtomWithIdx(idx).GetAtomicNum() not in _UREA_SUBSTITUENT_ELEMENTS and not mol.GetAtomWithIdx(idx).IsInRing()
            for idx in outside
        ):
            raise UnsupportedStructure("this element in an N-substituent of a urea is not supported yet")
        return [[name_branch(graph, c, n, halogens, mol=mol) for c in group] for n, group in zip(nitrogens, roots)]
    oxo_atoms = {
        n.GetIdx()
        for idx in outside
        if is_oxo_nitrogen(mol, mol.GetAtomWithIdx(idx))
        for n in [mol.GetAtomWithIdx(idx), *mol.GetAtomWithIdx(idx).GetNeighbors()]
        if n.GetAtomicNum() in (7, 8)
    }
    for idx in outside:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 and idx not in halogens and idx not in oxo_atoms and not atom.IsInRing():
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the urea core and its N-substituents "
                "is not supported yet"
            )
    return [[name_branch(graph, c, n, halogens, mol=mol) for c in group] for n, group in zip(nitrogens, roots)]


def n_prefix(n1_names, n2_names):
    """The 'N'/'N'' substituent prefix block: the nitrogen with more
    substituents takes the unprimed locant, then the one whose substituent
    comes first alphanumerically (P-14.3.5, P-14.5.2). Four identical substituents leave no substitutable hydrogen,
    so the locants are omitted (P-14.3.4.5)."""
    names = n1_names + n2_names
    if len(names) == 4 and len(set(names)) == 1:
        return format_mononuclear_prefixes(names)
    if len(n1_names) != len(n2_names):
        unprimed, primed = (n1_names, n2_names) if len(n1_names) > len(n2_names) else (n2_names, n1_names)
    else:
        key = lambda names: min((alpha_sort_key(name) for name, _ in names), default="")
        unprimed, primed = sorted((n1_names, n2_names), key=key)
    positions = {"N": unprimed, "N'": primed}
    return format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))


@dataclass(frozen=True)
class Chalcogenourea:
    atomic_num: int
    symbol: str
    word: str

    def _core(self, mol):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
                continue
            if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
                continue
            neighbors = atom.GetNeighbors()
            chalcogens = [n for n in neighbors if n.GetAtomicNum() == self.atomic_num]
            nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
            if len(chalcogens) != 1 or len(nitrogens) != 2:
                continue
            (chalcogen,) = chalcogens
            if (
                chalcogen.GetDegree() != 1
                or mol.GetBondBetweenAtoms(atom.GetIdx(), chalcogen.GetIdx()).GetBondTypeAsDouble() != 2.0
            ):
                continue
            if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
                continue
            if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
                continue
            if any(
                not is_urea_substituent_root(mol, nn)
                for n in nitrogens
                for nn in n.GetNeighbors()
                if nn.GetIdx() != atom.GetIdx()
            ):
                continue
            return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
        return None

    def has_shape(self, mol) -> bool:
        return self._core(mol) is not None

    def name(self, mol) -> str:
        word = self.word
        core = self._core(mol)
        if core is None:
            raise UnsupportedStructure(
                f"no {word} (H2N-C(={self.symbol})-NH2 or an N-substituted derivative) "
                f"shape found; this module only handles {word} and simple "
                f"N-substituted {word}s"
            )
        carbon_idx, (n1_idx, n2_idx) = core

        if len(Chem.GetMolFrags(mol)) > 1:
            raise UnsupportedStructure("multi-fragment structures are not supported yet")

        (chalcogen_idx,) = (
            n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == self.atomic_num
        )
        graph = adjacency(mol)
        n1_names, n2_names = n_substituent_names(
            mol, graph, {carbon_idx, chalcogen_idx, n1_idx, n2_idx}, (n1_idx, n2_idx), carbon_idx, junior_groups=True
        )
        return f"{n_prefix(n1_names, n2_names)}{word}"
