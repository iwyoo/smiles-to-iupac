"""Naming of isocyanides (the 'isocyano' substituent prefix, -NC) on
acyclic saturated hydrocarbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-61.9 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  '-NC' (-N+#C-) group attached to a parent hydride is named with the
  simple substituent prefix 'isocyano' -- the same structural role as
  'nitro'/'azido' (`_nitro.py`/`_azide.py`), which this module mirrors.
- Confirmed via PubChem for both the mononuclear case ('isocyanomethane',
  CID 11646) and a two-carbon chain ('isocyanoethane', CID 12226) -- unlike
  `_isocyanate.py`'s sibling module, there's no PubChem-autoname
  discrepancy to work around here.

Structural wrinkle this module has to handle that `_nitro.py`/`_azide.py`
don't: the atom attached directly to the parent chain is the group's
nitrogen (correctly excluded from the carbon-only chain search, same as
any heteroatom), but the isocyanide's *own* terminal atom is itself a
carbon (the '-C-' of '-N+#C-'). `_common.carbon_adjacency` includes every
carbon atom in the molecule regardless of context, so passing it unchanged
would let that isocyanide carbon appear as its own disconnected
single-atom "chain" candidate -- for the mononuclear case (e.g. methyl
isocyanide, one parent carbon vs. one isocyanide carbon, both length-1),
`_acyclic.py`'s chain search has no length-based way to prefer the real
parent carbon over the isocyanide's own carbon. This module therefore
builds its own carbon-only adjacency graph that excludes each isocyanide
group's carbon up front, verified against the mononuclear case above.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any other heteroatom (including halogen substituents, unverified for
  this group) besides an isocyanide group's own N/C.
- More than one isocyanide group, any unsaturation elsewhere in the
  molecule, any ring, or aromatic rings.
"""

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._acyclic import name_from_carbon_graph

_ISOCYANIDE_ALLOWED_ATOMIC_NUMS = {6, 7}


def _isocyanide_nitrogens(mol):
    """The carbon-attached nitrogen of each isocyanide group (-N+#C-): N
    (degree 2, +1 charge) bonded to exactly one skeletal carbon (single
    bond) and one terminal carbon (triple bond, degree 1, -1 charge)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() != 1:
            continue
        neighbors = atom.GetNeighbors()
        single_bond_carbons = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        triple_bond_carbons = [
            n
            for n in neighbors
            if n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 3.0
        ]
        if len(single_bond_carbons) != 1 or len(triple_bond_carbons) != 1:
            continue
        (cx,) = triple_bond_carbons
        if cx.GetDegree() != 1 or cx.GetFormalCharge() != -1:
            continue
        matches.append(atom)
    return matches


def has_isocyanide_shape(mol) -> bool:
    return bool(_isocyanide_nitrogens(mol))


def _isocyanide_group_atoms(mol, nitrogens):
    idxs = set()
    for n in nitrogens:
        idxs.add(n.GetIdx())
        (cx,) = (
            nbr
            for nbr in n.GetNeighbors()
            if mol.GetBondBetweenAtoms(n.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 3.0
        )
        idxs.add(cx.GetIdx())
    return idxs


def _carbon_adjacency_excluding(mol, excluded_idxs):
    """Like `_common.carbon_adjacency`, but omitting `excluded_idxs` from
    the node set entirely -- needed here since an isocyanide's own terminal
    carbon must not be mistaken for a skeletal chain atom (module
    docstring)."""
    graph = {
        atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in excluded_idxs
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def name_isocyanide(mol) -> str:
    nitrogens = _isocyanide_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure("no isocyanide (-NC) group found; this module only handles isocyanides")
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isocyanide group is out of scope for this module")

    group_atom_idxs = _isocyanide_group_atoms(mol, nitrogens)

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOCYANIDE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isocyanide group's own N/C "
                "(P-61.9) are not supported yet"
            )
        if atomic_num == 7 and atom.GetIdx() not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain isocyanide group "
                "is out of scope for this module"
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
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.9's scope "
            "here is limited to a saturated chain)"
        )

    (n1,) = nitrogens
    terminals = {n1.GetIdx(): "isocyano"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals)
