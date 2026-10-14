"""Nitro and nitroso groups as named monovalent prefixes (P-61.5.1, P-61.5.2).
`contract_nitro_groups` turns each such group into an iodine placeholder whose atom map number names the prefix, so
every halogen-aware module cites 'nitro' and 'nitroso' without learning them and the placeholder survives a SMILES
round trip.
"""

from rdkit import Chem

from ._common import NITRO_MAP_NUMBER, NITROSO_MAP_NUMBER, is_nitro_nitrogen

_MAP_OF = {"nitro": NITRO_MAP_NUMBER, "nitroso": NITROSO_MAP_NUMBER}


def _group(mol, atom):
    """(prefix, oxygen atoms) of the nitro or nitroso group whose nitrogen is `atom`, or None."""
    if atom.GetAtomicNum() != 7 or atom.IsInRing() or atom.GetIsotope():
        return None
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    hosts = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 8 or n.GetDegree() != 1]
    if len(hosts) != 1 or hosts[0].GetAtomicNum() not in (6, 7):
        return None
    if is_nitro_nitrogen(mol, atom.GetIdx()):
        return "nitro", oxygens
    if (
        len(oxygens) == 1
        and atom.GetDegree() == 2
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), oxygens[0].GetIdx()).GetBondTypeAsDouble() == 2.0
    ):
        return "nitroso", oxygens
    return None


def contract_nitro_groups(mol):
    """`mol` with each nitro and nitroso group collapsed to a placeholder, or None. The nitrogen itself becomes the
    placeholder, so the neighbours of the host keep their order and its stereodescriptor stays valid."""
    if any(a.GetAtomicNum() == 53 or a.GetAtomMapNum() for a in mol.GetAtoms()):
        return None
    groups = [(a.GetIdx(), *found) for a in mol.GetAtoms() if (found := _group(mol, a))]
    if not groups:
        return None
    rw = Chem.RWMol(mol)
    drop = set()
    for nitrogen, prefix, oxygens in groups:
        atom = rw.GetAtomWithIdx(nitrogen)
        atom.SetAtomicNum(53)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
        atom.SetAtomMapNum(_MAP_OF[prefix])
        drop |= {o.GetIdx() for o in oxygens}
    for idx in sorted(drop, reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out
