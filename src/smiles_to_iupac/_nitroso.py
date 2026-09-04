"""Naming of nitroso compounds (the 'nitroso' substituent prefix, -N=O)
on acyclic saturated hydrocarbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-61.5.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  '-N=O' is named with the simple substituent prefix 'nitroso', the same
  "no parent-hydride-selection role, no seniority of its own" pattern as
  its sibling group 'nitro' (`_nitro.py`, -NO2) -- confirmed via PubChem
  PUG REST: CID 70075 (`CN=O`) -> "nitrosomethane", CID 79124 (`CCN=O`)
  -> "nitrosoethane" (the 2-carbon case omits its own locant too, same
  P-14.3.4.2(b) rule `_nitro.py`'s own docstring already confirms for
  'nitroethane'), CID 21528903 (`ClCCN=O`) ->
  "1-chloro-2-nitrosoethane" (coexists freely with a halogen prefix), CID
  13116308 (`O=NCCN=O`) -> "1,2-dinitrosoethane" (multiple nitroso groups
  combine with the ordinary multiplying prefix, same mechanism
  `_nitro.py` already uses for 'dinitro').
- This module reuses `_acyclic.py`'s `name_from_carbon_graph` exactly the
  way `_nitro.py`/`_ether.py`/`_peroxide.py` do: each nitroso group's
  nitrogen is passed in as a precomputed 'terminals' entry (a leaf
  substituent excluded from the carbon-only chain search), merged with
  the ordinary halogen terminals dict.

Explicitly out of scope (raise `UnsupportedStructure`), matching
`_nitro.py`'s own scope:
- Any other heteroatom besides a nitroso group's own N/O and (optionally)
  halogen substituents -- in particular, this module doesn't attempt any
  Table 3.3 seniority coexistence with a characteristic-group suffix.
- Any unsaturation, any ring, or aromatic rings (a separate module's
  territory).
- aci-Nitroso tautomers or oxime-shaped nitrogen (a structurally
  different group this module doesn't recognize at all -- see
  `_imine.py`'s oxime handling for the =N-OH shape instead).
"""

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    halogen_substituents,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._acyclic import name_from_carbon_graph
from ._substituents import name_branch

_NITROSO_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _nitroso_nitrogen_atoms(mol):
    """Nitrogen atoms shaped like a nitroso group: bonded to exactly one
    carbon (single bond) and one oxygen (double bond, terminal,
    monovalent), formal charge 0."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 1:
            continue
        (carbon,) = carbons
        (oxygen,) = oxygens
        if oxygen.GetDegree() != 1 or oxygen.GetFormalCharge() != 0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        matches.append(atom)
    return matches


def has_nitroso_shape(mol) -> bool:
    return bool(_nitroso_nitrogen_atoms(mol))


def _validate_nitroso_atoms(mol, nitroso_atom_idxs, ring_atoms):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path: every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a chain carbon, a halogen, or one of a
    nitroso group's own N/O atoms."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _NITROSO_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitroso group's own N/O "
                "(P-61.5.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atomic_num == 7 and idx not in nitroso_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain nitroso group is "
                "out of scope for this module (see e.g. nitro, P-61.5.1)"
            )
        if atomic_num == 8 and idx not in nitroso_atom_idxs:
            raise UnsupportedStructure(
                "an oxygen atom not belonging to a nitroso group is out of "
                "scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "more than one aromatic ring is not supported yet"
                if ring_atoms
                else "aromatic rings are out of scope for this module (see "
                "the separate aromatic-ring module)"
            )


def _name_benzene_ring_nitroso_chain(mol, nitroso_nitrogens, nitroso_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): a single, otherwise-unsubstituted benzene ring
    is senior to a chain of the same (plain-hydrocarbon) class regardless
    of the chain's length -- since 'nitroso' has no suffix form, there is
    no characteristic group here to force the chain to be parent instead
    (mirrors `_nitro.py`'s `_name_benzene_ring_nitro_chain`, PR #321).
    Confirmed via PubChem: CID 11473 (nitrosobenzene, direct bond), CID
    12267972/21470433 (chain spacer cases -- see module's own
    'nitroso'/'nitro' sibling-group docstring cross-reference)."""
    _validate_nitroso_atoms(mol, nitroso_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a nitroso chain is not supported yet"
        )
    ring_atom, chain_root = attachment

    non_ring_unsaturation = [
        (a, b)
        for a, b, _ in non_single_bonds(mol)
        if a not in nitroso_atom_idxs
        and b not in nitroso_atom_idxs
        and a not in ring_atoms
        and b not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "nitroso chain is not supported yet"
        )

    terminals = dict(halogen_substituents(mol))
    for n in nitroso_nitrogens:
        terminals[n.GetIdx()] = "nitroso"

    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_nitroso(mol) -> str:
    nitroso_nitrogens = _nitroso_nitrogen_atoms(mol)
    if not nitroso_nitrogens:
        raise UnsupportedStructure("no nitroso (-N=O) group found; this module only handles nitroso compounds")

    nitroso_atom_idxs = set()
    for n in nitroso_nitrogens:
        nitroso_atom_idxs.add(n.GetIdx())
        nitroso_atom_idxs.update(o.GetIdx() for o in n.GetNeighbors() if o.GetAtomicNum() == 8)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_nitroso_chain(mol, nitroso_nitrogens, nitroso_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_nitroso_atoms(mol, nitroso_atom_idxs, frozenset())
    if any(a not in nitroso_atom_idxs and b not in nitroso_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.5.1's "
            "scope here is limited to a saturated chain)"
        )

    terminals = dict(halogen_substituents(mol))
    for n in nitroso_nitrogens:
        terminals[n.GetIdx()] = "nitroso"

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals)
