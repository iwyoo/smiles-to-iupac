"""Fusion-locant-letter naming (P-25.3.1.3) for a plain benzo ring ortho-
fused onto quinoline or isoquinoline, mirroring
`_polycyclic_component_fusion.py`'s mechanism for indole/1-benzofuran:
walk the base's own fixed retained-name numbering to compute the letter.
Both bases support plain ortho-fusion at the same three periphery bonds
on their carbocyclic ring ('f'/'g'/'h'); the pyridine-ring bonds and
fusion-atom-adjacent bonds are out of scope. Verified against PubChem
CIDs 6796/520238/9191/123043/601692/160447 (see tests).
"""

from rdkit import Chem

from ._common import UnsupportedStructure

_QUINOLINE_REF = Chem.MolFromSmiles("c1ccc2ncccc2c1")
_ISOQUINOLINE_REF = Chem.MolFromSmiles("c1ccc2cnccc2c1")

# Role -> ref atom index, from each base's own fixed retained-name
# numbering (see module docstring for the periphery walk).
_QUINOLINE_ROLE_TO_IDX = {"C5": 9, "C6": 0, "C7": 1, "C8": 2}
_ISOQUINOLINE_ROLE_TO_IDX = {"C5": 9, "C6": 0, "C7": 1, "C8": 2}

_SUPPORTED_FUSION_BONDS = {
    frozenset(("C5", "C6")): "f",
    frozenset(("C6", "C7")): "g",
    frozenset(("C7", "C8")): "h",
}

_BASES = {
    "quinoline": (_QUINOLINE_REF, _QUINOLINE_ROLE_TO_IDX, "benzo[{letter}]quinoline"),
    "isoquinoline": (_ISOQUINOLINE_REF, _ISOQUINOLINE_ROLE_TO_IDX, "benzo[{letter}]isoquinoline"),
}


def _find_fusion_letter(mol, reference, role_to_idx):
    matches = mol.GetSubstructMatches(reference, useChirality=False)
    if len(matches) != 1:
        return None
    (match,) = matches
    idx_to_role = {v: k for k, v in role_to_idx.items()}
    target_to_role = {match[i]: role for i, role in idx_to_role.items()}
    core_atoms = set(match)
    extra_atoms = set(range(mol.GetNumAtoms())) - core_atoms
    if len(extra_atoms) != 4:
        return None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(a).GetIsAromatic() for a in extra_atoms):
        return None

    fusion_roles = set()
    for a in extra_atoms:
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            if n.GetIdx() in target_to_role:
                fusion_roles.add(target_to_role[n.GetIdx()])
            elif n.GetIdx() not in extra_atoms:
                return None
    if len(fusion_roles) != 2:
        return None
    return _SUPPORTED_FUSION_BONDS.get(frozenset(fusion_roles))


def _find_core(mol):
    if mol.GetNumAtoms() != 14:
        return None
    if mol.GetRingInfo().NumRings() != 3:
        return None
    if any(not atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return None
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if atom.GetAtomicNum() not in (6, 7):
            return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7) != 1:
        return None

    for reference, role_to_idx, template in _BASES.values():
        letter = _find_fusion_letter(mol, reference, role_to_idx)
        if letter is not None:
            return template.format(letter=letter)
    return None


def has_pyridine_bicyclic_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_pyridine_bicyclic_fusion(mol) -> str:
    name = _find_core(mol)
    if name is None:
        raise UnsupportedStructure(
            "this tricyclic system is not a supported benzo-fused "
            "quinoline/isoquinoline shape (see P-25.3.1.3)"
        )
    return name
