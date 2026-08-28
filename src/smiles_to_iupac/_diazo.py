"""Naming of diazo compounds (the 'diazo' substituent prefix, =N2) on
acyclic saturated hydrocarbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-61.4 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  '=N2' (=N+=N-) group attached to a single carbon atom (replacing two of
  its hydrogens) is named with the simple substituent prefix 'diazo' --
  confirmed by the worked example 'diazomethane (PIN)'.
- Unlike `_nitro.py`/`_azide.py`, the group's first nitrogen is attached to
  its carbon by a *double* bond rather than a single one -- but
  `_acyclic.py`'s `name_from_carbon_graph` 'terminals' dict mechanism is
  bond-order-agnostic (it works purely off graph connectivity, from
  `_common.adjacency`), so this doesn't require any special handling: the
  substituted carbon stays a completely ordinary skeletal chain atom
  (confirmed for both a chain-terminal position, 'diazomethane'/
  '1-diazopropane', and an internal one, '2-diazopropane', where the
  carbon keeps two other skeletal-carbon neighbors). Unlike
  `_isocyanate.py`/`_isocyanide.py`, no atom belonging to the diazo group
  itself is a carbon, so there's no "exclude the group's own carbon from
  the chain search" wrinkle here either.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any other heteroatom (including halogen substituents, unverified for
  this group) besides a diazo group's own two nitrogens.
- More than one diazo group, any unsaturation elsewhere in the molecule,
  any ring, or aromatic rings.
- A diazo group coexisting with a characteristic-group suffix elsewhere in
  the molecule (e.g. the Blue Book's own '1-diazo-...-2-one' example) --
  no Table 3.3 seniority handling here.
"""

from ._common import UnsupportedStructure, adjacency, carbon_adjacency, non_single_bonds
from ._acyclic import name_from_carbon_graph

_DIAZO_ALLOWED_ATOMIC_NUMS = {6, 7}


def _diazo_root_nitrogens(mol):
    """The carbon-attached nitrogen of each diazo group (=N+=N-): N1
    (degree 2, +1 charge) double-bonded to exactly one carbon and to a
    terminal nitrogen N2 (degree 1, -1 charge, also double-bonded)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() != 1:
            continue
        neighbors = atom.GetNeighbors()
        double_bond_carbons = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        double_bond_nitrogens = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(double_bond_carbons) != 1 or len(double_bond_nitrogens) != 1:
            continue
        (n2,) = double_bond_nitrogens
        if n2.GetDegree() != 1 or n2.GetFormalCharge() != -1:
            continue
        matches.append(atom)
    return matches


def has_diazo_shape(mol) -> bool:
    return bool(_diazo_root_nitrogens(mol))


def name_diazo(mol) -> str:
    root_nitrogens = _diazo_root_nitrogens(mol)
    if not root_nitrogens:
        raise UnsupportedStructure("no diazo (=N2) group found; this module only handles diazo compounds")
    if len(root_nitrogens) > 1:
        raise UnsupportedStructure("more than one diazo group is out of scope for this module")

    diazo_atom_idxs = set()
    for n1 in root_nitrogens:
        (n2,) = (
            n
            for n in n1.GetNeighbors()
            if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(n1.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        )
        diazo_atom_idxs.update((n1.GetIdx(), n2.GetIdx()))

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _DIAZO_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a diazo group's own two nitrogens "
                "(P-61.4) are not supported yet"
            )
        if atomic_num == 7 and atom.GetIdx() not in diazo_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain diazo group is "
                "out of scope for this module"
            )
        if atom.GetFormalCharge() not in (0, 1, -1) or atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")
    if any(a not in diazo_atom_idxs and b not in diazo_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.4's scope "
            "here is limited to a saturated chain)"
        )

    (n1,) = root_nitrogens
    terminals = {n1.GetIdx(): "diazo"}

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals)
