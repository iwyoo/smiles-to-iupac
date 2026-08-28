"""Naming of isothiocyanates (the 'isothiocyanato' substituent prefix,
-NCS) on acyclic saturated hydrocarbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-61.8 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): the
  sulfur chalcogen analogue of the isocyanate group (`_isocyanate.py`,
  P-61.8's own 'thio' functional-replacement infix), named the same way
  with '-N=C=S' in place of '-N=C=O' and the substituent prefix
  'isothiocyanato' in place of 'isocyanato'. `_isocyanate.py` itself
  explicitly scoped this out as a chalcogen-analogue follow-up; this
  module mirrors its structure directly. Confirmed via PubChem PUG REST:
  CID 10966 (`CC[N]=C=S`) -> "isothiocyanatoethane", CID 69403
  (`CCC[N]=C=S`) -> "1-isothiocyanatopropane" -- both agree with the
  mechanical prefix pattern. The mononuclear case (`C[N]=C=S`) has the
  same PubChem-autoname-vs-PIN mismatch `_isocyanate.py` already
  documents for its own oxygen analogue (CID 11167's auto-generated name
  is "methylimino(sulfanylidene)methane", not a prefix form) -- by the
  same direct analogy that already settled `_isocyanate.py`'s
  'isocyanatomethane' from P-61.8's own 'isocyanatoborane (PIN)' worked
  example, 'isothiocyanatomethane' is taken as the PIN here too.

Structural wrinkle (shared with `_isocyanate.py`/`_isocyanide.py`, see
their module docstrings for the general explanation): the atom attached to
the parent chain is the group's nitrogen, but the group's own central atom
('-C=S' of '-N=C=S') is itself a carbon that `_common.carbon_adjacency`
would otherwise mistake for a skeletal chain atom. This module builds its
own carbon-only adjacency graph excluding that atom, mirroring
`_isocyanate.py` exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- The other chalcogen analogues (isoselenocyanato/isotellurocyanato).
- Halogen substituents (unverified for this group, same as
  `_isocyanate.py`).
- More than one isothiocyanate group, any other heteroatom, any
  unsaturation elsewhere in the molecule, any ring, or aromatic rings.
"""

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._acyclic import name_from_carbon_graph

_ISOTHIOCYANATE_ALLOWED_ATOMIC_NUMS = {6, 7, 16}


def _carbon_adjacency_excluding(mol, excluded_idxs):
    """Like `_common.carbon_adjacency`, but omitting `excluded_idxs` from
    the node set entirely -- needed here since an isothiocyanate's own
    central carbon must not be mistaken for a skeletal chain atom (module
    docstring). Mirrors `_isocyanate.py`'s identical helper."""
    graph = {
        atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in excluded_idxs
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def _isothiocyanate_nitrogens(mol):
    """The carbon-attached nitrogen of each isothiocyanate group
    (-N=C=S): N (degree 2, neutral) bonded to exactly one skeletal carbon
    (single bond) and one central carbon (double bond, degree 2,
    neutral), that central carbon itself double-bonded to a terminal,
    degree-1 sulfur."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        single_bond_carbons = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        double_bond_carbons = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(single_bond_carbons) != 1 or len(double_bond_carbons) != 1:
            continue
        (cc,) = double_bond_carbons
        if cc.GetDegree() != 2 or cc.GetFormalCharge() != 0:
            continue
        sulfurs = [
            n
            for n in cc.GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
            and n.GetAtomicNum() == 16
            and mol.GetBondBetweenAtoms(cc.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(sulfurs) != 1:
            continue
        (s,) = sulfurs
        if s.GetDegree() != 1 or s.GetFormalCharge() != 0:
            continue
        matches.append(atom)
    return matches


def has_isothiocyanate_shape(mol) -> bool:
    return bool(_isothiocyanate_nitrogens(mol))


def _isothiocyanate_group_atoms(mol, nitrogens):
    idxs = set()
    for n in nitrogens:
        idxs.add(n.GetIdx())
        (cc,) = (
            nbr
            for nbr in n.GetNeighbors()
            if mol.GetBondBetweenAtoms(n.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 2.0
        )
        idxs.add(cc.GetIdx())
        idxs.update(s.GetIdx() for s in cc.GetNeighbors() if s.GetAtomicNum() == 16)
    return idxs


def name_isothiocyanate(mol) -> str:
    nitrogens = _isothiocyanate_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure(
            "no isothiocyanate (-NCS) group found; this module only handles isothiocyanates"
        )
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isothiocyanate group is out of scope for this module")

    group_atom_idxs = _isothiocyanate_group_atoms(mol, nitrogens)

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOTHIOCYANATE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isothiocyanate group's own N/C/S "
                "(P-61.8) are not supported yet"
            )
        if atomic_num in (7, 16) and atom.GetIdx() not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen or sulfur atom not shaped like a plain "
                "isothiocyanate group is out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.8's scope "
            "here is limited to a saturated chain)"
        )

    (n1,) = nitrogens
    terminals = {n1.GetIdx(): "isothiocyanato"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals)
