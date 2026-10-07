"""Azines (R2C=N-N=CR2, P-68.3.1.2.3): hydrazine with an ylidene prefix on each nitrogen, named by `_hydrazine.py`."""



def _find_cn_bonds(mol):
    return [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
    ]


def _find_azine_nitrogens(mol, cn_bonds):
    """Given exactly two C=N double bonds, return (c1, n1, c2, n2) if their
    two nitrogens are joined by a single bond and each nitrogen has no
    other neighbor (the R2C=N-N=CR2 azine shape); else None."""
    if len(cn_bonds) != 2:
        return None
    carbons, nitrogens = [], []
    for bond in cn_bonds:
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        carbons.append(a if a.GetAtomicNum() == 6 else b)
        nitrogens.append(a if a.GetAtomicNum() == 7 else b)
    n1, n2 = nitrogens
    if n1.GetIdx() == n2.GetIdx():
        return None
    nn_bond = mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx())
    if nn_bond is None or nn_bond.GetBondTypeAsDouble() != 1.0:
        return None
    if n1.GetDegree() != 2 or n2.GetDegree() != 2:
        return None
    return carbons[0].GetIdx(), n1.GetIdx(), carbons[1].GetIdx(), n2.GetIdx()


def has_azine_shape(mol) -> bool:
    return _find_azine_nitrogens(mol, _find_cn_bonds(mol)) is not None


def name_azine(mol) -> str:
    """P-68.3.1.2.3: an azine is hydrazine carrying an ylidene prefix on each nitrogen."""
    from ._hydrazine import name_hydrazine

    return name_hydrazine(mol)
