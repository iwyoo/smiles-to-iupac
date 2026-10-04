"""Benzo-fused heteromonocycle parents named benzo + component (P-25.2.2.4): 1,3-benzoxazole, 1H-benzimidazole, 2H-1,3-benzodioxole, 2H-1-benzopyran, 1H-indazole."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_numbering_general import _HETERO_RANK, general_peripheral_numberings

_HETEROATOMS = {"O", "S", "Se", "Te", "N"}
_RETAINED_COMPONENTS = {"pyridine", "pyridazine", "pyrimidine", "pyrazine"}
_UNAMBIGUOUS_COMPONENTS = {"imidazole", "triazole"}


def _split_rings(mol):
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings) != 2 or len(rings[0] & rings[1]) != 2:
        return None
    benzene = [r for r in rings if len(r) == 6 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and mol.GetAtomWithIdx(a).GetIsAromatic() for a in r)]
    if len(benzene) != 1:
        return None
    other = rings[0] if rings[1] is benzene[0] else rings[1]
    if len(other) not in (5, 6, 7):
        return None
    return benzene[0], other


def _saturated(mol, ring_atoms):
    sat = []
    for a in ring_atoms:
        atom = mol.GetAtomWithIdx(a)
        symbol = atom.GetSymbol()
        if symbol in ("O", "S", "Se", "Te"):
            continue
        has_double = atom.GetIsAromatic() or any(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())
        if symbol == "N" and atom.GetTotalNumHs() > 0 and (atom.GetIsAromatic() or not has_double):
            sat.append(a)
        elif not has_double:
            sat.append(a)
    return sat


def _mancude_requirement(mol, ring_atoms):
    can = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetSymbol() not in ("O", "S", "Se", "Te")]
    adj = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in can] for a in can}

    def best(remaining):
        if not remaining:
            return 0
        first = min(remaining)
        options = [best(remaining - {first})]
        options += [1 + best(remaining - {first, o}) for o in adj[first] if o in remaining]
        return max(options)

    return len(can) - 2 * best(frozenset(can))


def _component(mol, other):
    from .core import smiles_to_iupac
    from ._ring_diyl_numbering import _sanitized_mancude

    sp3 = {a for a in other if mol.GetAtomWithIdx(a).GetAtomicNum() == 6}
    bare, _, _ = _sanitized_mancude(mol, other, sp3)
    name = smiles_to_iupac(Chem.MolToSmiles(bare))
    return re.sub(r"^(?:\d+H-)?(?:[\d,]+-)?", "", name)


def _match(mol):
    if mol.GetRingInfo().NumRings() != 2 or any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms()):
        return None
    split = _split_rings(mol)
    if split is None:
        return None
    benzene, other = split
    if mol.GetNumAtoms() != len(benzene | other):
        return None
    hetero = [a for a in other if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if not hetero or any(mol.GetAtomWithIdx(a).GetSymbol() not in _HETEROATOMS for a in hetero):
        return None
    if any(a in benzene for a in hetero):
        return None
    sat = _saturated(mol, other | benzene)
    if len(sat) != _mancude_requirement(mol, other | benzene):
        return None
    return benzene, other, hetero, sat


def has_benzo_heterocycle_name(mol) -> bool:
    match = _match(mol)
    return match is not None and _component(mol, match[1]) not in _RETAINED_COMPONENTS


def name_benzo_heterocycle(mol) -> str:
    benzene, other, hetero, sat = _match(mol)
    numberings = general_peripheral_numberings(mol, ignore_indicated=True)
    if not numberings:
        raise UnsupportedStructure("this benzo-fused ring has no supported numbering")
    locant = lambda numbering, a: int(numbering[a].rstrip("abcdefgh"))
    numbering = min(numberings, key=lambda n: sorted(locant(n, a) for a in sat))
    ih = sorted(locant(numbering, a) for a in sat)
    ordered = sorted(hetero, key=lambda a: (_HETERO_RANK[mol.GetAtomWithIdx(a).GetSymbol()], locant(numbering, a)))
    hetero_locants = [locant(numbering, a) for a in ordered]
    component = _component(mol, other)
    ih_text = ",".join(f"{p}H" for p in ih) + "-" if ih else ""
    if component == "pyrrole":
        return f"{ih_text}{'indole' if hetero_locants == [1] else 'isoindole'}"
    if component == "pyrazole":
        return f"{ih_text}indazole"
    stem = ("benz" if component[0] in "aeiou" else "benzo") + component
    if component in _UNAMBIGUOUS_COMPONENTS:
        return f"{ih_text}{stem}"
    return f"{ih_text}{','.join(map(str, hetero_locants))}-{stem}"
