"""Heterones and nitriles of an acyclic polychalcogane chain, e.g. 1-methyl-2-phenyl-1λ6,2λ6-disulfane-1,1,2,2-tetrone
(P-64.1.2.2, P-64.4.1, P-64.7.2) and 1-hydroxy-1,1-diiodo-1λ4-disulfane-2-carbonitrile (P-67.3.1): the chain of S, Se
or Te atoms is the parent hydride, each doubly bonded oxygen is cited by the suffix 'one', a nitrile carbon on the chain
by 'carbonitrile', and each chalcogen of non-standard valence by its λ label (P-14.1). The numbering gives lowest
locants first to the λ atoms, then to the suffix, then to the prefixes, then to the alphabetically first prefix
(P-31.1.4).
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    alpha_sort_key,
    group_substituents,
    halogen_substituents,
    multiplied_word,
)
from ._numerals import numerical_term
from ._phosphanone import _SENIOR
from ._substituents import format_substituent_prefixes, name_branch

_STEMS = {16: "sulfane", 34: "selane", 52: "tellane"}
_OXO_VALENCES = (4, 6)


def _terminal_oxygens(atom):
    return [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 8
        and n.GetDegree() == 1
        and not n.GetFormalCharge()
        and atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]


def _chain(mol):
    """The atoms of the single chalcogen-chalcogen chain in order, or None."""
    chalcogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS]
    if len(chalcogens) < 2 or len({a.GetAtomicNum() for a in chalcogens}) != 1:
        return None
    ids = {a.GetIdx() for a in chalcogens}
    links = {
        a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors() if n.GetIdx() in ids] for a in chalcogens
    }
    ends = [i for i, ns in links.items() if len(ns) == 1]
    if len(ends) != 2 or any(len(ns) > 2 or not ns for ns in links.values()):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(ids):
        following = [n for n in links[order[-1]] if n != previous]
        if not following:
            return None
        previous = order[-1]
        order.append(following[0])
    return order


def _nitrile_carbons(mol, order):
    """{chain atom: nitrile carbon} for the C#N groups bonded to the chain."""
    found = {}
    for i in order:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if (
                n.GetAtomicNum() == 6
                and n.GetDegree() == 2
                and any(
                    m.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 3.0
                    for m in n.GetNeighbors()
                )
            ):
                found[i] = n.GetIdx()
    return found


def has_chalcogen_chain_heterone_shape(mol) -> bool:
    order = _chain(mol)
    if order is None:
        return False
    atoms = [mol.GetAtomWithIdx(i) for i in order]
    return any(_terminal_oxygens(a) for a in atoms) or any(a.GetTotalValence() != 2 for a in atoms)


def name_chalcogen_chain_heterone(mol) -> str:
    order = _chain(mol)
    if order is None or not has_chalcogen_chain_heterone_shape(mol):
        raise UnsupportedStructure("no chalcogen chain carrying doubly bonded oxygen found")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    nitriles = _nitrile_carbons(mol, order)
    if any(
        mol.HasSubstructMatch(query) for query in _SENIOR if Chem.MolToSmarts(query) != "[C&X2]#[N&X1]"
    ) or len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX2]#[NX1]"))) != len(nitriles):
        raise UnsupportedStructure("an acid, ester, amide, aldehyde or nitrile outside the chain outranks the chain")
    if any(mol.GetAtomWithIdx(i).IsInRing() for i in order):
        raise UnsupportedStructure("a ring chalcogen is not an acyclic heterone parent")
    graph = adjacency(mol)
    oxo = {i: [o.GetIdx() for o in _terminal_oxygens(mol.GetAtomWithIdx(i))] for i in order}
    lambdas = {}
    for i in order:
        valence = mol.GetAtomWithIdx(i).GetTotalValence()
        if valence not in (2, *_OXO_VALENCES) or (oxo[i] and valence == 2):
            raise UnsupportedStructure("this chalcogen valence is not supported yet")
        if valence != 2:
            lambdas[i] = valence
    chain_atoms = set(order) | {o for ids in oxo.values() for o in ids} | set(nitriles.values())
    chain_atoms |= {m.GetIdx() for c in nitriles.values() for m in mol.GetAtomWithIdx(c).GetNeighbors()}
    ketones = sum(
        1
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and any(n.GetIdx() not in chain_atoms and n.GetAtomicNum() == 8 and n.GetDegree() == 1 for n in a.GetNeighbors())
    )
    if ketones > len(chain_atoms) - len(order):
        raise UnsupportedStructure("more ketone groups than heterone oxygens outrank the chalcogen chain")
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    attached = []
    for i in order:
        for n in graph[i]:
            if n in chain_atoms:
                continue
            if mol.GetBondBetweenAtoms(i, n).GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("an unsaturated link to the chalcogen chain is not supported yet")
            attached.append((i, name_branch(graph, n, i, halogens, aromatic, mol=mol, unsaturated=True)))

    def numbering(sequence):
        position = {atom: k + 1 for k, atom in enumerate(sequence)}
        prefixes = sorted(
            ((alpha_sort_key(name), position[i]) for i, (name, _) in attached), key=lambda item: (item[0], item[1])
        )
        return (
            sorted(position[i] for i in lambdas),
            sorted(position[i] for i in sequence for _ in oxo[i] + ([nitriles[i]] if i in nitriles else [])),
            sorted(position[i] for i, _ in attached),
            [loc for _, loc in prefixes],
        ), position

    best_key, position = min(
        (numbering(sequence) for sequence in (order, order[::-1])), key=lambda candidate: candidate[0]
    )
    grouped = group_substituents(
        {loc: [entry for i, entry in attached if position[i] == loc] for loc in sorted(position.values())}
    )
    prefix_text = format_substituent_prefixes({n: v for n, v in grouped.items() if v["locants"]})
    lambda_text = ",".join(f"{position[i]}λ{v}" for i, v in sorted(lambdas.items(), key=lambda item: position[item[0]]))
    suffix_locants = ",".join(str(loc) for loc in best_key[1])
    stem = _STEMS[mol.GetAtomWithIdx(order[0]).GetAtomicNum()]
    parent = f"{numerical_term(len(order))}{stem}"
    head = f"{prefix_text}-" if prefix_text else ""
    lambda_head = f"{lambda_text}-" if lambda_text else ""
    if nitriles:
        if oxo and any(oxo.values()):
            raise UnsupportedStructure("a nitrile and a heterone suffix on one chalcogen chain are not named together yet")
        return f"{head}{lambda_head}{parent}-{suffix_locants}-{multiplied_word(len(best_key[1]), 'carbonitrile')}"
    if not best_key[1]:
        return f"{head}{lambda_head}{parent}"
    return f"{head}{lambda_head}{parent}-{suffix_locants}-{multiplied_word(len(best_key[1]), 'one')}"
