from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._alcohol import name_alcohol
from ._aldehyde import name_aldehyde
from ._amine import name_amine
from ._aromatic import find_aromatic_fused_core, name_aromatic_fused
from ._bicyclic import find_bicyclic_core, name_bicycloalkane
from ._carboxylic_acid import has_carboxylic_acid_shape, name_carboxylic_acid
from ._common import UnsupportedStructure, non_single_bonds
from ._cyclic import name_cycloalkane
from ._ester import has_ester_shape, name_ester
from ._ether import has_ether_shape, name_ether
from ._ketone import name_ketone
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._polyspiro import find_linear_polyspiro_chain, name_linear_polyspiro
from ._spiro import find_monospiro_atom, name_monospiro
from ._tricyclic import find_propellane_core, name_propellane
from ._unsaturated import name_acyclic_unsaturated


def _is_aldehyde_shaped(carbonyl_oxygen):
    (carbon,) = carbonyl_oxygen.GetNeighbors()
    return carbon.GetAtomicNum() == 6 and sum(1 for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6) == 1


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        # A plain -O- ether (P-63.2.1) has no suffix, so it must be routed
        # here before the carbonyl/alcohol checks below, none of which
        # accept a degree-2 oxygen at all.
        if has_ether_shape(mol):
            return name_ether(mol)
        # A carbon bearing both a carbonyl oxygen and a second, carbon-bonded
        # oxygen is an ester (-COO-), which must be routed before the
        # carboxylic-acid/aldehyde/ketone checks below: its carbonyl half
        # would otherwise look aldehyde/ketone-shaped, and (for a rejected,
        # out-of-scope case) its non-carbonyl oxygen would never satisfy the
        # carboxylic acid module's hydroxyl (O-H) requirement anyway.
        if has_ester_shape(mol):
            return name_ester(mol)
        # A carbon bearing both a carbonyl and a hydroxyl oxygen is a -COOH
        # group (Table 3.3's most senior suffix here) and must be routed
        # before the aldehyde/ketone/alcohol checks below, which would
        # otherwise misread its carbonyl or hydroxyl half in isolation.
        if has_carboxylic_acid_shape(mol):
            return name_carboxylic_acid(mol)
        # A doubly-bonded, monovalent oxygen is carbonyl-shaped (aldehyde or
        # ketone, depending on how many carbon neighbors its carbon has);
        # anything else falls to the alcohol module, which itself rejects a
        # coexisting carbonyl oxygen it finds among otherwise hydroxyl-only
        # atoms (Table 3.3 seniority between 'ol'/'one'/'al' isn't handled).
        carbonyl_oxygens = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetBonds()[0].GetBondTypeAsDouble() == 2.0
        ]
        if carbonyl_oxygens:
            if any(_is_aldehyde_shaped(o) for o in carbonyl_oxygens):
                return name_aldehyde(mol)
            return name_ketone(mol)
        return name_alcohol(mol)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return name_amine(mol)

    num_rings = mol.GetRingInfo().NumRings()
    # Aromatic rings carry non-single (order 1.5) bonds, which every other
    # ring module's non_single_bonds check rejects; an aromatic ring
    # system's carbon skeleton can also be graph-isomorphic to a *saturated*
    # bicyclic through pentacyclic core (e.g. naphthalene <-> decahydro-
    # naphthalene), so this check must run, and must succeed for any
    # in-scope aromatic shape, before num_rings==1 or any saturated
    # find_*_core below gets a chance to misdetect it and raise the wrong
    # ("unsaturated ... not supported yet") error (see _aromatic.py).
    if num_rings >= 1:
        aromatic_core = find_aromatic_fused_core(mol)
        if aromatic_core is not None:
            return name_aromatic_fused(mol, aromatic_core)
    if num_rings == 0:
        bonds = non_single_bonds(mol)
        if not bonds:
            return name_acyclic_alkane(mol)
        if all(order in (2.0, 3.0) for _, _, order in bonds):
            return name_acyclic_unsaturated(mol)
        raise UnsupportedStructure(
            "a bond order other than double or triple is not supported yet "
            "(see P-31.1.1.1)"
        )
    if num_rings == 1:
        return name_cycloalkane(mol)

    # num_rings >= 2 from here on. RDKit's SSSR can overcount rings for
    # symmetric bridged bicyclics (see _bicyclic.py's find_bicyclic_core
    # docstring), so bicyclic detection isn't gated on num_rings == 2 either.
    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro(mol, spiro_atom)
    polyspiro_chain = find_linear_polyspiro_chain(mol)
    if polyspiro_chain is not None:
        return name_linear_polyspiro(mol, polyspiro_chain)
    bicyclic_core = find_bicyclic_core(mol)
    if bicyclic_core is not None:
        return name_bicycloalkane(mol, bicyclic_core)
    for ring_count in (3, 4, 5):
        core = find_polycyclic_core(mol, ring_count)
        if core is not None:
            return name_polycycloalkane(mol, core, ring_count)
    propellane_core = find_propellane_core(mol)
    if propellane_core is not None:
        return name_propellane(mol, propellane_core)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
