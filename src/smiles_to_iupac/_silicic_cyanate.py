"""Silicic acid with hydroxy groups replaced by cyanato-type groups (P-67.1.2.4.1.3): the group -OCN would make an anhydride,
which ranks below an acid, so the acid is named with the prefix 'cyanato' ('cyanatosilicic acid')."""

from rdkit import Chem

from ._numerals import numerical_term

_GROUPS = (
    ("[O;D2]-C#N", "cyanato"),
    ("[N;D2]=C=O", "isocyanato"),
    ("[S;D2]-C#N", "thiocyanato"),
    ("[N;D2]=C=S", "isothiocyanato"),
    ("[Se;D2]-C#N", "selenocyanato"),
)
_PATTERNS = [(Chem.MolFromSmarts(smarts), prefix) for smarts, prefix in _GROUPS]


def silicic_cyanate_name(mol):
    silicons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 14]
    if len(silicons) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    silicon = silicons[0]
    if silicon.GetDegree() != 4 or silicon.GetFormalCharge() or silicon.GetIsotope():
        return None
    groups, hydroxy, covered = [], 0, {silicon.GetIdx()}
    for neighbor in silicon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 1:
            hydroxy += 1
            covered.add(neighbor.GetIdx())
            continue
        for pattern, prefix in _PATTERNS:
            match = next((m for m in mol.GetSubstructMatches(pattern) if m[0] == neighbor.GetIdx()), None)
            if match is not None:
                groups.append(prefix)
                covered |= set(match)
                break
        else:
            return None
    if not groups or len(set(groups)) != 1 or covered != set(range(mol.GetNumAtoms())):
        return None
    count = len(groups)
    return f"{numerical_term(count) if count > 1 else ''}{groups[0]}silicic acid"
