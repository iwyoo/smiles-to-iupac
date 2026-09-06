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
  molecule, or any non-benzene ring.

A single, otherwise-unsubstituted benzene ring gets a dedicated path
(`_name_benzene_ring_isocyanide_chain`), mirroring `_nitro.py`'s/
`_azide.py`'s identical benzene-ring path (P-44.1.2.2 rule (1): since
'isocyano' has no suffix form, a plain benzene ring is always senior to a
chain of the same hydrocarbon class regardless of chain length). The
isocyanide's own carbon-adjacency wrinkle (see above) doesn't arise here,
since `name_branch` treats the nitrogen as a leaf substituent (via the
same terminals-dict injection `_nitro.py` uses) and never recurses into
its own carbon. Confirmed via PubChem PUG REST: CID 13606
'c1ccccc1[N+]#[C-]' -> 'isocyanobenzene' (direct ring-nitrogen bond), CID
82558 'c1ccccc1C[N+]#[C-]' -> 'isocyanomethylbenzene', CID 427710
'c1ccccc1CC[N+]#[C-]' -> '2-isocyanoethylbenzene'.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._acyclic import name_from_carbon_graph
from ._substituents import name_branch

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


def _validate_isocyanide_atoms(mol, group_atom_idxs, ring_atoms):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path: every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a chain carbon or the isocyanide
    group's own N/C. Mirrors `_nitro.py`'s `_validate_nitro_atoms`."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ISOCYANIDE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an isocyanide group's own N/C "
                "(P-61.9) are not supported yet"
            )
        if atomic_num == 7 and idx not in group_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain isocyanide group "
                "is out of scope for this module"
            )
        if atom.GetFormalCharge() not in (0, 1, -1) or atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "more than one aromatic ring is not supported yet"
                if ring_atoms
                else "aromatic rings are out of scope for this module (see "
                "the separate aromatic-ring module)"
            )


def _name_benzene_ring_isocyanide_chain(mol, nitrogen, group_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): a single, otherwise-unsubstituted benzene ring
    is senior to a chain of the same (plain-hydrocarbon) class regardless
    of the chain's length -- since 'isocyano' has no suffix form, mirrors
    `_nitro.py`'s `_name_benzene_ring_nitro_chain` exactly. So the ring is
    always the parent hydride, and the chain (with its isocyanide
    substituent) is named as a single substituent prefix on it via
    `name_branch`."""
    _validate_isocyanide_atoms(mol, group_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an isocyanide chain is not supported yet"
        )
    ring_atom, chain_root = attachment

    non_ring_unsaturation = [
        (a, b)
        for a, b, _ in non_single_bonds(mol)
        if a not in group_atom_idxs
        and b not in group_atom_idxs
        and a not in ring_atoms
        and b not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "isocyanide chain is not supported yet"
        )

    terminals = {nitrogen.GetIdx(): "isocyano"}
    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals, mol=mol)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_isocyanide(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    nitrogens = _isocyanide_nitrogens(mol)
    if not nitrogens:
        raise UnsupportedStructure("no isocyanide (-NC) group found; this module only handles isocyanides")
    if len(nitrogens) > 1:
        raise UnsupportedStructure("more than one isocyanide group is out of scope for this module")
    (n1,) = nitrogens

    group_atom_idxs = _isocyanide_group_atoms(mol, nitrogens)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_isocyanide_chain(mol, n1, group_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_isocyanide_atoms(mol, group_atom_idxs, frozenset())
    if any(a not in group_atom_idxs and b not in group_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.9's scope "
            "here is limited to a saturated chain)"
        )

    terminals = {n1.GetIdx(): "isocyano"}
    excluded_carbons = group_atom_idxs - {n1.GetIdx()}
    carbon_graph = _carbon_adjacency_excluding(mol, excluded_carbons)

    return name_from_carbon_graph(adjacency(mol), carbon_graph, terminals, mol=mol)
