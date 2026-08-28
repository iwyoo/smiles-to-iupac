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
    non_single_bonds,
)
from ._acyclic import name_from_carbon_graph

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


def name_nitroso(mol) -> str:
    nitroso_nitrogens = _nitroso_nitrogen_atoms(mol)
    if not nitroso_nitrogens:
        raise UnsupportedStructure("no nitroso (-N=O) group found; this module only handles nitroso compounds")

    nitroso_atom_idxs = set()
    for n in nitroso_nitrogens:
        nitroso_atom_idxs.add(n.GetIdx())
        nitroso_atom_idxs.update(o.GetIdx() for o in n.GetNeighbors() if o.GetAtomicNum() == 8)

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _NITROSO_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitroso group's own N/O "
                "(P-61.5.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atomic_num == 7 and atom.GetIdx() not in nitroso_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain nitroso group is "
                "out of scope for this module (see e.g. nitro, P-61.5.1)"
            )
        if atomic_num == 8 and atom.GetIdx() not in nitroso_atom_idxs:
            raise UnsupportedStructure(
                "an oxygen atom not belonging to a nitroso group is out of "
                "scope for this module"
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
    if any(a not in nitroso_atom_idxs and b not in nitroso_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.5.1's "
            "scope here is limited to a saturated chain)"
        )

    terminals = dict(halogen_substituents(mol))
    for n in nitroso_nitrogens:
        terminals[n.GetIdx()] = "nitroso"

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals)
