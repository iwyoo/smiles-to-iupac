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
  elsewhere in the molecule, or any non-benzene ring.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_isocyanate_chain`), mirroring `_nitro.py`'s
equivalent: since 'isocyanato' has no suffix form, P-44.1.2.2 rule (1)
makes the ring the parent regardless of the chain's length, and the
chain (with its single isocyanato substituent) is named as one
substituent prefix on it via `name_branch` -- unlike the plain-chain path
above, no separate carbon-graph-excluding trick is needed here, since
`name_branch` already stops descending into a `terminals`-listed atom.
Confirmed via PubChem PUG REST (CID 7672 'c1ccccc1N=C=O' ->
'isocyanatobenzene', 'c1ccccc1CN=C=O' -> 'isocyanatomethylbenzene',
'c1ccccc1CCN=C=O' -> '2-isocyanatoethylbenzene').
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._common import (
    UnsupportedStructure,
    adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._acyclic import name_from_carbon_graph
from ._substituents import name_branch

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


def _validate_isocyanate_atoms(mol, group_atom_idxs, ring_atoms=frozenset()):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path (mirrors `_nitro.py`'s
    `_validate_nitro_atoms`): every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a non-aromatic chain carbon or one of
    the isocyanate group's own N/C/O atoms."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOCYANATE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isocyanate group's own N/C/O "
                "(P-61.8) are not supported yet"
            )
        if atomic_num in (7, 8) and idx not in group_atom_idxs:
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


def _name_benzene_ring_isocyanate_chain(mol, n1, group_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): since 'isocyanato' has no suffix form (module
    docstring), a single, otherwise-unsubstituted benzene ring is always
    the parent hydride, regardless of the chain's length -- mirrors
    `_nitro.py`'s `_name_benzene_ring_nitro_chain`. The chain (with its
    single isocyanato substituent) is named as one substituent prefix on
    the ring via `name_branch`."""
    _validate_isocyanate_atoms(mol, group_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an isocyanate chain is not supported yet"
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
            "isocyanate chain is not supported yet"
        )

    terminals = {n1.GetIdx(): "isocyanato"}
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals, mol=mol)
    display = enclose(branch_name) if is_compound else branch_name
    return f"{display}benzene"


def name_isocyanate(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    nitrogens = _isocyanate_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure("no isocyanate (-NCO) group found; this module only handles isocyanates")
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isocyanate group is out of scope for this module")

    (n1,) = nitrogens
    group_atom_idxs = _isocyanate_group_atoms(mol, nitrogens)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_isocyanate_chain(mol, n1, group_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_isocyanate_atoms(mol, group_atom_idxs)
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.8's scope "
            "here is limited to a saturated chain)"
        )

    terminals = {n1.GetIdx(): "isocyanato"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals, mol=mol)
