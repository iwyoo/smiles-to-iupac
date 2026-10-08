"""Carbon monoxide, the one binary oxide whose multiplying prefix 'mono' is kept (P-16.4.2, P-13.3.5)."""


def has_carbon_monoxide_shape(mol) -> bool:
    if mol.GetNumAtoms() != 2 or mol.GetNumBonds() != 1:
        return False
    atoms = {a.GetSymbol(): a for a in mol.GetAtoms()}
    if set(atoms) != {"C", "O"} or any(a.GetIsotope() for a in mol.GetAtoms()):
        return False
    carbon, oxygen = atoms["C"], atoms["O"]
    bond = mol.GetBondWithIdx(0).GetBondTypeAsDouble()
    dipolar = bond == 3.0 and carbon.GetFormalCharge() == -1 and oxygen.GetFormalCharge() == 1
    carbene = bond == 2.0 and not carbon.GetFormalCharge() and not oxygen.GetFormalCharge() and carbon.GetNumRadicalElectrons() == 2
    return dipolar or carbene


def name_carbon_monoxide(mol) -> str:
    return "carbon monoxide"
