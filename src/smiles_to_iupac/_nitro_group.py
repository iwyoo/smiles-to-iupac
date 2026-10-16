"""Nitro, nitroso and azide groups as named monovalent prefixes (P-61.5.1, P-61.5.2, P-61.7).
`contract_nitro_groups` turns each such group into an iodine placeholder whose atom map number names the prefix, so
every halogen-aware module cites 'nitro', 'nitroso' and 'azido' without learning them and the placeholder survives a
SMILES round trip.
"""

from rdkit import Chem

from ._common import AZIDO_MAP_NUMBER, NITRO_MAP_NUMBER, NITROSO_MAP_NUMBER, is_nitro_nitrogen

_MAP_OF = {"nitro": NITRO_MAP_NUMBER, "nitroso": NITROSO_MAP_NUMBER, "azido": AZIDO_MAP_NUMBER}


def _azide_group(mol, atom):
    """('azido', the two further nitrogens) of the azide group -N=[N+]=[N-] whose first nitrogen is `atom`, or None."""
    if atom.GetFormalCharge() or atom.GetDegree() != 2 or atom.GetTotalNumHs():
        return None
    for middle in atom.GetNeighbors():
        ends = [n for n in middle.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        host = [n for n in atom.GetNeighbors() if n.GetIdx() != middle.GetIdx()]
        if (
            middle.GetAtomicNum() == 7
            and middle.GetFormalCharge() == 1
            and middle.GetDegree() == 2
            and len(ends) == 1
            and ends[0].GetAtomicNum() == 7
            and ends[0].GetFormalCharge() == -1
            and ends[0].GetDegree() == 1
            and len(host) == 1
            and host[0].GetAtomicNum() in (6, 7)
            and mol.GetBondBetweenAtoms(atom.GetIdx(), middle.GetIdx()).GetBondTypeAsDouble() == 2.0
            and mol.GetBondBetweenAtoms(middle.GetIdx(), ends[0].GetIdx()).GetBondTypeAsDouble() == 2.0
            and mol.GetBondBetweenAtoms(atom.GetIdx(), host[0].GetIdx()).GetBondTypeAsDouble() == 1.0
        ):
            return "azido", [middle, ends[0]]
    return None


def nitro_or_nitroso_group(mol, atom):
    """(prefix, atoms to drop) of the nitro, nitroso or azide group whose nitrogen is `atom`, or None."""
    if atom.GetAtomicNum() != 7 or atom.IsInRing() or atom.GetIsotope():
        return None
    azide = _azide_group(mol, atom)
    if azide is not None:
        return azide
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
    groups = [(a.GetIdx(), *found) for a in mol.GetAtoms() if (found := nitro_or_nitroso_group(mol, a))]
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
