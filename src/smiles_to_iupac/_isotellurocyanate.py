"""Naming of isotellurocyanates (the 'isotellurocyanato' substituent
prefix, -N=C=Te) on acyclic saturated hydrocarbon chains, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-61.8 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): the
  tellurium chalcogen analogue of the isocyanate group (`_isocyanate.py`,
  `_isothiocyanate.py`, `_isoselenocyanate.py`, P-61.8's own 'telluro'
  functional-replacement infix), named the same way with '-N=C=Te' in
  place of '-N=C=O'/'-N=C=S'/'-N=C=Se' and the substituent prefix
  'isotellurocyanato'. Unlike the other three chalcogen analogues, no
  PubChem-listed compound was found for any of `CC[N]=C=[Te]`,
  `CCC[N]=C=[Te]`, or `C[N]=C=[Te]` (all CID 0) -- this is included as a
  reviewed (eyeballed), not independently verified, mechanical extension
  of the identical pattern already confirmed twice over for each of the
  oxygen/sulfur/selenium analogues.

Structural wrinkle (shared with `_isocyanate.py`/`_isothiocyanate.py`/
`_isoselenocyanate.py`/`_isocyanide.py`, see their module docstrings for
the general explanation): the atom attached to the parent chain is the
group's nitrogen, but the group's own central atom ('-C=Te' of
'-N=C=Te') is itself a carbon that `_common.carbon_adjacency` would
otherwise mistake for a skeletal chain atom. This module builds its own
carbon-only adjacency graph excluding that atom, mirroring
`_isoselenocyanate.py` exactly.

Explicitly out of scope (raise `UnsupportedStructure`):
- Halogen substituents (unverified for this group, same as the other
  isocyanate-family modules).
- More than one isotellurocyanate group, any other heteroatom, any
  unsaturation elsewhere in the molecule, or any non-benzene ring.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_isotellurocyanate_chain`), mirroring the other three
chalcogen analogues' identical path: since 'isotellurocyanato' has no
suffix form, P-44.1.2.2 rule (1) makes the ring the parent regardless of
the chain's length. Not independently PubChem-verified here (same
sparse-tellurium-data gap noted above for `c1ccccc1N=C=[Te]`/
`c1ccccc1CN=C=[Te]`, both CID 0) -- inherited unchanged from the
identical mechanism already confirmed for oxygen/sulfur/selenium.
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

_TELLURIUM = 52
_ISOTELLUROCYANATE_ALLOWED_ATOMIC_NUMS = {6, 7, _TELLURIUM}


def _carbon_adjacency_excluding(mol, excluded_idxs):
    """Like `_common.carbon_adjacency`, but omitting `excluded_idxs` from
    the node set entirely -- needed here since an isotellurocyanate's own
    central carbon must not be mistaken for a skeletal chain atom (module
    docstring). Mirrors `_isoselenocyanate.py`'s identical helper."""
    graph = {
        atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in excluded_idxs
    }
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def _isotellurocyanate_nitrogens(mol):
    """The carbon-attached nitrogen of each isotellurocyanate group
    (-N=C=Te): N (degree 2, neutral) bonded to exactly one skeletal
    carbon (single bond) and one central carbon (double bond, degree 2,
    neutral), that central carbon itself double-bonded to a terminal,
    degree-1 tellurium."""
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
        telluriums = [
            n
            for n in cc.GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
            and n.GetAtomicNum() == _TELLURIUM
            and mol.GetBondBetweenAtoms(cc.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(telluriums) != 1:
            continue
        (te,) = telluriums
        if te.GetDegree() != 1 or te.GetFormalCharge() != 0:
            continue
        matches.append(atom)
    return matches


def has_isotellurocyanate_shape(mol) -> bool:
    return bool(_isotellurocyanate_nitrogens(mol))


def _isotellurocyanate_group_atoms(mol, nitrogens):
    idxs = set()
    for n in nitrogens:
        idxs.add(n.GetIdx())
        (cc,) = (
            nbr
            for nbr in n.GetNeighbors()
            if mol.GetBondBetweenAtoms(n.GetIdx(), nbr.GetIdx()).GetBondTypeAsDouble() == 2.0
        )
        idxs.add(cc.GetIdx())
        idxs.update(te.GetIdx() for te in cc.GetNeighbors() if te.GetAtomicNum() == _TELLURIUM)
    return idxs


def _validate_isotellurocyanate_atoms(mol, group_atom_idxs, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_isoselenocyanate.py`'s
    `_validate_isoselenocyanate_atoms`): every atom outside `ring_atoms`
    (empty for the plain-chain path) must be a non-aromatic chain carbon
    or one of the isotellurocyanate group's own N/C/Te atoms."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOTELLUROCYANATE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isotellurocyanate group's own "
                "N/C/Te (P-61.8) are not supported yet"
            )
        if atomic_num in (7, _TELLURIUM) and idx not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen or tellurium atom not shaped like a plain "
                "isotellurocyanate group is out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )


def _name_benzene_ring_isotellurocyanate_chain(mol, n1, group_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'isotellurocyanato' has no suffix form
    (module docstring), a single, otherwise-unsubstituted benzene ring is
    always the parent hydride, regardless of the chain's length --
    mirrors `_isoselenocyanate.py`'s
    `_name_benzene_ring_isoselenocyanate_chain`."""
    _validate_isotellurocyanate_atoms(mol, group_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an isotellurocyanate chain is not supported yet"
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
            "isotellurocyanate chain is not supported yet"
        )

    terminals = {n1.GetIdx(): "isotellurocyanato"}
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_isotellurocyanate(mol) -> str:
    nitrogens = _isotellurocyanate_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure(
            "no isotellurocyanate (-NCTe) group found; this module only handles isotellurocyanates"
        )
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isotellurocyanate group is out of scope for this module")

    (n1,) = nitrogens
    group_atom_idxs = _isotellurocyanate_group_atoms(mol, nitrogens)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_isotellurocyanate_chain(mol, n1, group_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_isotellurocyanate_atoms(mol, group_atom_idxs)
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.8's scope "
            "here is limited to a saturated chain)"
        )

    terminals = {n1.GetIdx(): "isotellurocyanato"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals)
