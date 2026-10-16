"""Thiourea, selenourea and tellurourea (H2N-C(=E)-NH2, E = S/Se/Te) with any number of N-substituents
(P-66.1.4.1); unlike `_urea.py` no semicarbazide (amino-substituted nitrogen) handling. The N/N' citation helpers
are shared with `_urea.py`."""

from dataclasses import dataclass

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    group_substituents,
    halogen_substituents,
    is_nitro_nitrogen,
)
from ._substituents import (
    alpha_sort_key,
    cited_branch_stereo,
    format_mononuclear_prefixes,
    format_substituent_prefixes,
    name_branch,
)


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
    """A carbon, a halogen, a nitro/nitroso nitrogen or a ring nitrogen: the atoms a urea or guanidine nitrogen may carry."""
    return atom.GetAtomicNum() in (6, *HALOGEN_PREFIXES) or is_oxo_nitrogen(mol, atom) or (
        atom.GetAtomicNum() == 7 and atom.IsInRing()
    )


_HYDRIDE_GROUP_CENTRES = {5, 13, 14, 15, 31, 32, 33, 49, 50, 51, 81, 82, 83}
_UREA_SUBSTITUENT_ELEMENTS = {6, 7, 8, 16, 34, 52, *_HYDRIDE_GROUP_CENTRES, *HALOGEN_PREFIXES}


def is_element_hydride_group(atom):
    """The centre of a boranyl, silyl, germyl, stannyl, phosphanyl, ... group on a nitrogen: a Group 13-15 atom with single
    bonds to hydrogen, carbon, halogen and that nitrogen only, so that no acid or oxo group makes it a senior class (P-68.2.1)."""
    return (
        atom.GetAtomicNum() in _HYDRIDE_GROUP_CENTRES
        and not atom.GetFormalCharge()
        and not atom.IsInRing()
        and all(
            b.GetBondTypeAsDouble() == 1.0 and b.GetOtherAtom(atom).GetAtomicNum() in (6, 7, *HALOGEN_PREFIXES)
            for b in atom.GetBonds()
        )
    )


def carries_only_element_hydride_groups(mol):
    """Whether the molecule holds a Group 13-15 atom and every such atom is the centre of a plain hydride group, none of
    them an oxoacid or amide centre that would outrank a urea or guanidine (P-41)."""
    centres = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _HYDRIDE_GROUP_CENTRES]
    return bool(centres) and all(is_element_hydride_group(a) for a in centres)


def is_urea_substituent_root(mol, atom):
    """An atom a urea nitrogen may carry: any element whose group is cited as a prefix, so that hydroxy, alkoxy,
    amino and sulfanyl groups on the nitrogen are ordinary substituents (P-66.1.6.1.1.2)."""
    return atom.GetAtomicNum() in _UREA_SUBSTITUENT_ELEMENTS


_SENIOR_TO_UREA = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,Se,Te])[OX2H1,OX1-]",
        "[CX3;!$([CX3;R]([OX2,NX3,SX2;R]))](=[O,S,Se,Te])[OX2][#6]",
        "[CX3;!$([CX3;R]([OX2,NX3,SX2;R]))](=[O,S,Se,Te])[SX2,SeX2,TeX2]",
        "[CX3](=[O,S,Se,Te])[F,Cl,Br,I;X1]",
        "[CX3;!$([CX3;R]([OX2,NX3,SX2;R]))](=[O,S,Se,Te])[NX3;!R]",
        "[SX4,SX3,SeX4,TeX4](=O)[OX2H1,OX1-,OX2,NX3;!R]",
    )
]


def has_group_senior_to_urea(mol, core_atoms, hydrazide_junior=False, junior_cores=frozenset()):
    """A carboxylic or sulfonic acid, ester, acid halide, amide or amidine outside the urea core, or a charge other
    than a nitro group's: the class that outranks urea as the parent (P-41, P-66.1.6.1.1.5). With `hydrazide_junior`
    an acyl or sulfonyl group on a nitrogen of the core makes a hydrazide, which ranks below an amidine (P-41).
    `junior_cores` are carbons of urea groups ranking equal to or below the parent urea (P-41, P-45.2)."""
    for query in _SENIOR_TO_UREA:
        if any(
            match[0] not in core_atoms
            and match[0] not in junior_cores
            and not (hydrazide_junior and match[-1] in core_atoms)
            for match in mol.GetSubstructMatches(query)
        ):
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


def _named_groups(mol, graph, halogens, nitrogens, roots, core_atoms):
    with cited_branch_stereo(mol, graph, core_atoms, [r for group in roots for r in group]):
        return [[name_branch(graph, c, n, halogens, mol=mol) for c in group] for n, group in zip(nitrogens, roots)]


def n_substituent_names(
    mol, graph, core_atoms, nitrogens, carbon_idx, junior_groups=False, hydrazide_junior=False, junior_cores=frozenset()
):
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
        if has_group_senior_to_urea(mol, core_atoms, hydrazide_junior, junior_cores):
            raise UnsupportedStructure("a group senior to the urea is the parent, not an N-substituent")
        if any(
            mol.GetAtomWithIdx(idx).GetAtomicNum() not in _UREA_SUBSTITUENT_ELEMENTS and not mol.GetAtomWithIdx(idx).IsInRing()
            for idx in outside
        ):
            raise UnsupportedStructure("this element in an N-substituent of a urea is not supported yet")
        return _named_groups(mol, graph, halogens, nitrogens, roots, core_atoms)
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
    return _named_groups(mol, graph, halogens, nitrogens, roots, core_atoms)


def n_prefix(n1_names, n2_names):
    """The 'N'/'N'' substituent prefix block: the nitrogen with more
    substituents takes the unprimed locant, then the one whose substituent
    comes first alphanumerically (P-14.3.5, P-14.5.2). A sole substituent sits on one of two equivalent nitrogens
    (P-14.3.4.3) and four identical ones leave no substitutable hydrogen (P-14.3.4.5), so no locant is cited."""
    names = n1_names + n2_names
    if names == [("carbamimidoyl", False)]:
        return "N-carbamimidoyl"
    if len(names) == 1 or (len(names) == 4 and len(set(names)) == 1):
        return format_mononuclear_prefixes(names)
    unprimed, primed = (n1_names, n2_names) if unprimed_first(n1_names, n2_names) else (n2_names, n1_names)
    positions = {"N": unprimed, "N'": primed}
    return format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))


def unprimed_first(n1_names, n2_names):
    """Whether the first of two amino nitrogens takes the unprimed locant: the one with more substituents, then the one
    holding the alphanumerically first substituent."""
    if len(n1_names) != len(n2_names):
        return len(n1_names) > len(n2_names)
    key = lambda names: min((alpha_sort_key(name) for name, _ in names), default="")
    return key(n1_names) <= key(n2_names)


def unprimed_first_labelled(n1_idx, n2_idx, n1_names, n2_names, labels):
    """`unprimed_first` for nitrogens that may carry nuclides: of two otherwise equivalent nitrogens the labelled one takes
    the unprimed locant (P-14.3.5)."""
    if n1_names != n2_names:
        return unprimed_first(n1_names, n2_names)
    if n1_idx in labels and n2_idx in labels and labels[n1_idx] != labels[n2_idx]:
        raise UnsupportedStructure("differently modified equivalent amino nitrogens are not supported yet")
    return n1_idx in labels or n2_idx not in labels


@dataclass(frozen=True)
class Chalcogenourea:
    atomic_num: int

    def has_shape(self, mol) -> bool:
        from ._urea import senior_chalcogen

        return senior_chalcogen(mol) == self.atomic_num

    def name(self, mol) -> str:
        from ._urea import name_cores

        return name_cores(mol, self.atomic_num)
