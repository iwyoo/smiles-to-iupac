"""Fullerene names (P-27.2): a closed cage of three-coordinate carbon atoms is named
`(C<n>-<point group>[(<ring surround>)])[<ring sizes>]fullerene` from its planar
embedding (`_fullerene_symmetry.py`); the 24 isolated-pentagon C84 isomers also cite
their atlas number (`_fullerene_spiral.py`). Silicon-for-carbon C60 is an exact match."""

from functools import lru_cache

import networkx as nx
from rdkit import Chem

from ._common import RuleUndefinedStructure, UnsupportedStructure
from ._fullerene_spiral import match_c84_isomer
from ._fullerene_symmetry import read_cage
from ._pin import mark

_FULLERENE_C60_SMILES = (
    "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1"
    "C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4="
    "C3C3=C5C(=C81)C%10=C23"
)
_FULLERENE_C70_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%10=C%10C8=C5C1=C%10C1=C%13C5=C8C1=C2"
    "C1=C3C2=C3C%10=C%13C%14=C3C1=C8C1=C3C5=C%12C5=C8C%11=C%11C9=C7C7=C9C6=C4C2=C2"
    "C%10=C4C(=C29)C2=C6C(=C8C8=C9C6=C4C%13=C9C(=C%141)C3=C85)C%11=C27"
)
_FULLERENE_C76_SMILES = (
    "C12=C3C4=C5C6=C1C7=C8C2=C9C1=C2C%10=C%11C(=C13)C1=C4C3=C4C%12=C%13C3=C5C3=C6"
    "C5=C7C6=C7C5=C5C3=C%13C3=C%13C%12=C%12C%14=C%15C%16=C%17C%18=C%19C%20=C%16C"
    "(=C%13%14)C(=C35)C7=C%20C3=C%19C5=C7C%18=C%13C%17=C%14C%15=C%15C%12=C4C1=C%11"
    "C%15=C%14C%10=C%13C2=C7C9=C5C8=C63"
)

_SILA_C60_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C(=C61)C1=C6C%10=C%11C(=C21)C1=C3C2=C4C3=C4C5=C7C5=C7C8="
    "C8C9=C6C6=C%10C9=C%10C%11=C1C1=C2C2=C3C3=C%11C%12=C2C1=C%10C1=C9C2=C9C(=C7C(="
    "C%11[Si]9=C1%12)C5=C43)C8=C62"
)

_C84_D6H_ISOMER_24_SMILES = (
    "C12=C3C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15C%16=C%17C%18=C%19C%20=C%21C"
    "%22=C(C1=C1C(=C36)C(=C7%10)C(=C%11%14)C(=C%18%15)C1=C%22%19)C1=C2C2=C4C3=C4"
    "C5=C8C5=C6C9=C%12C7=C8C%13=C%16C9=C%10C%17=C%20C%11=C%12C%21=C1C1=C2C2=C3C3="
    "C%13C%14=C2C1=C%12C1=C%14C2=C(C%10=C%111)C9=C8C1=C2C%13=C(C5=C43)C6=C71"
)

_FULLERENE_NAMES = {Chem.CanonSmiles(_SILA_C60_SMILES): "sila(C60-Ih)[5,6]fullerene"}

_NO_RING_SURROUND = (
    "P-27.2.1 gives no single ring-surround subscript for the rings on the principal axis of this cage"
)


def _is_carbon_cage(mol) -> bool:
    if mol.GetNumAtoms() < 20 or len(Chem.GetMolFrags(mol)) != 1:
        return False
    return all(
        a.GetAtomicNum() == 6
        and a.GetDegree() == 3
        and a.GetTotalNumHs() == 0
        and a.GetFormalCharge() == 0
        and a.GetIsotope() == 0
        and a.GetNumRadicalElectrons() == 0
        for a in mol.GetAtoms()
    )


def _systematic_name(mol):
    if not _is_carbon_cage(mol):
        return None
    cage = read_cage(_cage_graph(mol, range(mol.GetNumAtoms())))
    if cage is None:
        return None
    subscript = f"({cage.surround})" if cage.surround else ""
    sizes = ",".join(map(str, cage.sizes))
    name = f"(C{mol.GetNumAtoms()}-{cage.symbol}{subscript})[{sizes}]fullerene"
    return name, _NO_RING_SURROUND if cage.ambiguous else None


def _name_and_reason(mol):
    exact = _FULLERENE_NAMES.get(Chem.MolToSmiles(mol)) or match_c84_isomer(mol)
    return (exact, None) if exact else _systematic_name(mol)


def has_fullerene_name(mol) -> bool:
    return _name_and_reason(mol) is not None


def name_fullerene(mol) -> str:
    name, reason = _name_and_reason(mol)
    return mark(name, reason) if reason else name


def is_fullerene_cage(mol, atoms) -> bool:
    """`atoms` form a closed carbon cage of twelve five-membered rings and only five- or six-membered rings."""
    atoms = set(atoms)
    if len(atoms) < 20 or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in atoms):
        return False
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= atoms]
    return sum(len(r) == 5 for r in rings) == 12 and all(len(r) in (5, 6) for r in rings)


def has_substituted_fullerene_cage(mol) -> bool:
    ring_info = mol.GetRingInfo()
    if sum(len(r) == 5 for r in ring_info.AtomRings()) < 12:
        return False
    cage = {a for r in ring_info.AtomRings() for a in r}
    return mol.GetNumAtoms() > len(cage) and is_fullerene_cage(mol, cage)


def _cage_graph(mol, atoms):
    atoms = set(atoms)
    return nx.Graph(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
    )


@lru_cache(maxsize=None)
def _reference_graph(smiles):
    return _cage_graph(Chem.MolFromSmiles(smiles), range(Chem.MolFromSmiles(smiles).GetNumAtoms()))


_NUMBERED_CAGE_STEMS = (
    (_FULLERENE_C60_SMILES, "(C60-Ih)[5,6]fullerene"),
    (_FULLERENE_C70_SMILES, "(C70-D5h(6))[5,6]fullerene"),
)


def numbered_cage_stem(mol, atoms):
    """The substitutive stem of the C60-Ih or C70-D5h(6) cage formed by `atoms`, or None for any other cage."""
    graph = _cage_graph(mol, atoms)
    for smiles, stem in _NUMBERED_CAGE_STEMS:
        reference = _reference_graph(smiles)
        if graph.number_of_nodes() == reference.number_of_nodes() and nx.is_isomorphic(graph, reference):
            return stem
    return None


UNDEFINED_NUMBERING = (
    "the locants of a fullerene other than (C60-Ih)[5,6]fullerene and (C70-D5h(6))[5,6]fullerene are not given by the "
    "rules (P-27.0, P-27.3, P-27.7)"
)


def require_defined_fullerene_numbering(mol, atoms):
    if is_fullerene_cage(mol, atoms) and numbered_cage_stem(mol, atoms) is None:
        raise RuleUndefinedStructure(UNDEFINED_NUMBERING)
