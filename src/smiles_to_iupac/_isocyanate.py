"""Naming of isocyanates (the 'isocyanato' substituent prefix, -NCO) on
acyclic saturated hydrocarbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-61.8 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  '-N=C=O' group attached to a parent hydride is named with the simple
  substituent prefix 'isocyanato' -- the same structural role as
  'nitro'/'azido' (`_nitro.py`/`_azide.py`), confirmed by the worked
  example 'isocyanatocyclohexane (PIN)'.
- The prefix applies even to a very small parent hydride -- confirmed by
  P-61.8's own 'isocyanatoborane (PIN)' example (BH2-NCO) -- so a
  mononuclear alkane parent (e.g. methyl isocyanate) is expected to follow
  the same pattern ('isocyanatomethane'), even though PubChem's own
  auto-generated name for that exact structure (CID 12228) is
  'methylimino(oxo)methane' instead: a different, non-prefix parent
  selection this project's own worked-example citation above overrides.
  This is treated as another instance of the PubChem-autoname-vs-PIN
  mismatch documented elsewhere in this project (see e.g.
  `_hydroxylamine.py`, `_sulfoxide.py`, `_nitro.py`'s own two-carbon-chain
  case) rather than a real exception to the prefix rule.
- Two- and three-carbon chains are confirmed directly against PubChem
  ('isocyanatoethane' CID 8022, '2-isocyanatopropane' CID 61277), which
  agree with the mononuclear case's own mechanical pattern.

Structural wrinkle (shared with `_isocyanide.py`, see its module docstring
for the general explanation): the atom attached to the parent chain is the
group's nitrogen, but the group's own central atom ('-C=O' of '-N=C=O') is
itself a carbon that `_common.carbon_adjacency` would otherwise mistake
for a skeletal chain atom. This module builds its own carbon-only
adjacency graph excluding that atom, mirroring `_isocyanide.py` exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- Chalcogen analogues (isothiocyanato/isoselenocyanato/isotellurocyanato,
  P-61.8's own 'thio'/'seleno'/'telluro' functional-replacement infixes).
- Halogen substituents (unverified for this group, unlike `_nitro.py`/
  `_azide.py`).
- More than one isocyanate group, any other heteroatom, any unsaturation
  elsewhere in the molecule, any ring, or aromatic rings.
"""

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._acyclic import name_from_carbon_graph

_ISOCYANATE_ALLOWED_ATOMIC_NUMS = {6, 7, 8}


def _carbon_adjacency_excluding(mol, excluded_idxs):
    """Like `_common.carbon_adjacency`, but omitting `excluded_idxs` from
    the node set entirely -- needed here since an isocyanate's own central
    carbon must not be mistaken for a skeletal chain atom (module
    docstring). Mirrors `_isocyanide.py`'s identical helper."""
    graph = {
        atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in excluded_idxs
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def _isocyanate_nitrogens(mol):
    """The carbon-attached nitrogen of each isocyanate group (-N=C=O): N
    (degree 2, neutral) bonded to exactly one skeletal carbon (single bond)
    and one central carbon (double bond, degree 2, neutral), that central
    carbon itself double-bonded to a terminal, degree-1 oxygen."""
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
        oxygens = [
            n
            for n in cc.GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
            and n.GetAtomicNum() == 8
            and mol.GetBondBetweenAtoms(cc.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(oxygens) != 1:
            continue
        (o,) = oxygens
        if o.GetDegree() != 1 or o.GetFormalCharge() != 0:
            continue
        matches.append(atom)
    return matches


def has_isocyanate_shape(mol) -> bool:
    return bool(_isocyanate_nitrogens(mol))


def _isocyanate_group_atoms(mol, nitrogens):
    idxs = set()
    for n in nitrogens:
        idxs.add(n.GetIdx())
        (cc,) = (
            nbr
            for nbr in n.GetNeighbors()
            if mol.GetBondBetweenAtoms(n.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 2.0
        )
        idxs.add(cc.GetIdx())
        idxs.update(o.GetIdx() for o in cc.GetNeighbors() if o.GetAtomicNum() == 8)
    return idxs


def name_isocyanate(mol) -> str:
    nitrogens = _isocyanate_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure("no isocyanate (-NCO) group found; this module only handles isocyanates")
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isocyanate group is out of scope for this module")

    group_atom_idxs = _isocyanate_group_atoms(mol, nitrogens)

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOCYANATE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isocyanate group's own N/C/O "
                "(P-61.8) are not supported yet"
            )
        if atomic_num in (7, 8) and atom.GetIdx() not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen or oxygen atom not shaped like a plain "
                "isocyanate group is out of scope for this module"
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
    terminals = {n1.GetIdx(): "isocyanato"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals)
