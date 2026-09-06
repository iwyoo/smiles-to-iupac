"""Naming of isoselenocyanates (the 'isoselenocyanato' substituent
prefix, -N=C=Se) on acyclic saturated hydrocarbon chains, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-61.8 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): the
  selenium chalcogen analogue of the isocyanate group (`_isocyanate.py`,
  `_isothiocyanate.py`, P-61.8's own 'seleno' functional-replacement
  infix), named the same way with '-N=C=Se' in place of '-N=C=O'/
  '-N=C=S' and the substituent prefix 'isoselenocyanato'. Confirmed via
  PubChem PUG REST: CID 12618729 (`CC[N]=C=[Se]`) ->
  "isoselenocyanatoethane", CID 134989633 (`CCC[N]=C=[Se]`) ->
  "1-isoselenocyanatopropane" -- both agree with the mechanical prefix
  pattern. The mononuclear case (`C[N]=C=[Se]`) has the same
  PubChem-autoname-vs-PIN mismatch `_isocyanate.py`/`_isothiocyanate.py`
  already document for their own analogues (CID 138232's auto-generated
  name is "methylimino(selanylidene)methane", not a prefix form) -- by
  the same direct analogy, 'isoselenocyanatomethane' is taken as the PIN
  here too.

Structural wrinkle (shared with `_isocyanate.py`/`_isothiocyanate.py`/
`_isocyanide.py`, see their module docstrings for the general
explanation): the atom attached to the parent chain is the group's
nitrogen, but the group's own central atom ('-C=Se' of '-N=C=Se') is
itself a carbon that `_common.carbon_adjacency` would otherwise mistake
for a skeletal chain atom. This module builds its own carbon-only
adjacency graph excluding that atom, mirroring `_isothiocyanate.py`
exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- The other chalcogen analogue (isotellurocyanato).
- Halogen substituents (unverified for this group, same as
  `_isocyanate.py`/`_isothiocyanate.py`).
- More than one isoselenocyanate group, any other heteroatom, any
  unsaturation elsewhere in the molecule, or any non-benzene ring.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_isoselenocyanate_chain`), mirroring
`_isocyanate.py`'s/`_isothiocyanate.py`'s equivalent: since
'isoselenocyanato' has no suffix form, P-44.1.2.2 rule (1) makes the ring
the parent regardless of the chain's length. Confirmed via PubChem PUG
REST (CID 555335 'c1ccccc1N=C=[Se]' -> 'isoselenocyanatobenzene', CID
12506035 'c1ccccc1CN=C=[Se]' -> 'isoselenocyanatomethylbenzene').
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._acyclic import name_from_carbon_graph
from ._substituents import name_branch

_SELENIUM = 34
_ISOSELENOCYANATE_ALLOWED_ATOMIC_NUMS = {6, 7, _SELENIUM}


def _carbon_adjacency_excluding(mol, excluded_idxs):
    """Like `_common.carbon_adjacency`, but omitting `excluded_idxs` from
    the node set entirely -- needed here since an isoselenocyanate's own
    central carbon must not be mistaken for a skeletal chain atom (module
    docstring). Mirrors `_isothiocyanate.py`'s identical helper."""
    graph = {
        atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in excluded_idxs
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def _isoselenocyanate_nitrogens(mol):
    """The carbon-attached nitrogen of each isoselenocyanate group
    (-N=C=Se): N (degree 2, neutral) bonded to exactly one skeletal
    carbon (single bond) and one central carbon (double bond, degree 2,
    neutral), that central carbon itself double-bonded to a terminal,
    degree-1 selenium."""
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
        seleniums = [
            n
            for n in cc.GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
            and n.GetAtomicNum() == _SELENIUM
            and mol.GetBondBetweenAtoms(cc.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(seleniums) != 1:
            continue
        (se,) = seleniums
        if se.GetDegree() != 1 or se.GetFormalCharge() != 0:
            continue
        matches.append(atom)
    return matches


def has_isoselenocyanate_shape(mol) -> bool:
    return bool(_isoselenocyanate_nitrogens(mol))


def _isoselenocyanate_group_atoms(mol, nitrogens):
    idxs = set()
    for n in nitrogens:
        idxs.add(n.GetIdx())
        (cc,) = (
            nbr
            for nbr in n.GetNeighbors()
            if mol.GetBondBetweenAtoms(n.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 2.0
        )
        idxs.add(cc.GetIdx())
        idxs.update(se.GetIdx() for se in cc.GetNeighbors() if se.GetAtomicNum() == _SELENIUM)
    return idxs


def _validate_isoselenocyanate_atoms(mol, group_atom_idxs, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_isothiocyanate.py`'s
    `_validate_isothiocyanate_atoms`): every atom outside `ring_atoms`
    (empty for the plain-chain path) must be a non-aromatic chain carbon
    or one of the isoselenocyanate group's own N/C/Se atoms."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOSELENOCYANATE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isoselenocyanate group's own "
                "N/C/Se (P-61.8) are not supported yet"
            )
        if atomic_num in (7, _SELENIUM) and idx not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen or selenium atom not shaped like a plain "
                "isoselenocyanate group is out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )


def _name_benzene_ring_isoselenocyanate_chain(mol, n1, group_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'isoselenocyanato' has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the chain's length --
    mirrors `_isothiocyanate.py`'s
    `_name_benzene_ring_isothiocyanate_chain`."""
    _validate_isoselenocyanate_atoms(mol, group_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an isoselenocyanate chain is not supported yet"
        )
    ring_atom, chain_root = attachment

    non_ring_unsaturation = [
        (a, b)
        for a, b, _ in non_single_bonds(mol)
        if a not in group_atom_idxs and b not in group_atom_idxs and a not in ring_atoms and b not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "isoselenocyanate chain is not supported yet"
        )

    terminals = {n1.GetIdx(): "isoselenocyanato"}
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals, mol=mol)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_isoselenocyanate(mol) -> str:
    nitrogens = _isoselenocyanate_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure(
            "no isoselenocyanate (-NCSe) group found; this module only handles isoselenocyanates"
        )
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isoselenocyanate group is out of scope for this module")

    (n1,) = nitrogens
    group_atom_idxs = _isoselenocyanate_group_atoms(mol, nitrogens)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_isoselenocyanate_chain(mol, n1, group_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_isoselenocyanate_atoms(mol, group_atom_idxs)
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.8's scope "
            "here is limited to a saturated chain)"
        )

    terminals = {n1.GetIdx(): "isoselenocyanato"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals, mol=mol)
