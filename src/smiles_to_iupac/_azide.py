"""Naming of azides (the 'azido' substituent prefix, -N3) on acyclic
saturated hydrocarbon chains, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-61.7 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  '-N3' (-N=N+=N-) group attached to a parent hydride is named with the
  simple substituent prefix 'azido' -- the same structural role as 'nitro'
  (`_nitro.py`, P-61.5.1), which this module mirrors closely. Confirmed by
  the worked example '(2-azidoethyl)benzene (PIN)'.
- Like nitro, 'azido' coexists freely with halogen substituents (P-35.2.1)
  as an ordinary prefix -- confirmed here for acyclic chains via PubChem,
  e.g. '1-azido-2-chloroethane'.
- Reuses `_acyclic.py`'s `name_from_carbon_graph` exactly the way
  `_nitro.py` does: each azide group's first (carbon-attached) nitrogen is
  passed in as a precomputed 'terminals' entry, merged with the ordinary
  halogen terminals dict.
- Unlike `_nitro.py`'s own two-carbon-chain finding, the two-carbon case
  here ('azidoethane', PubChem CID 79118) *does* match the general,
  unconditional P-14.3.4.2(b) locant-omission rule -- no discrepancy to
  document for this group.

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_nitro.py`'s own scope limits:
- Any other heteroatom besides an azide group's own three nitrogens and
  (optionally) halogen substituents -- no Table 3.3 seniority coexistence
  with a characteristic-group suffix.
- Any unsaturation elsewhere in the molecule, any ring, or aromatic rings
  (a separate module's territory, e.g. the Blue Book's own
  '(2-azidoethyl)benzene' example).
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

_AZIDE_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _azide_root_nitrogens(mol):
    """The carbon-attached nitrogen of each azide group (-N=N+=N-): a
    three-nitrogen chain N1-N2-N3 where N1 is bonded to exactly one carbon
    (single bond) and N2 (non-single bond), N2 is bonded to only N1 and N3
    (both non-single bonds), N3 is terminal (bonded only to N2, non-single
    bond), and the canonical +1/-1 charge separation is on N2/N3 -- the
    representation RDKit gives '-N=[N+]=[N-]'. The bond-order and charge
    checks (not just connectivity) matter here: a plain triazane chain
    (-NH-NH-NH2, all single bonds, all neutral) has the identical
    connectivity shape and must NOT be mistaken for an azide."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        n2_candidates = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(n2_candidates) != 1:
            continue
        (n2,) = n2_candidates
        if n2.GetDegree() != 2 or n2.GetFormalCharge() != 1:
            continue
        n3_candidates = [n for n in n2.GetNeighbors() if n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() == 7]
        if len(n3_candidates) != 1:
            continue
        (n3,) = n3_candidates
        if n3.GetDegree() != 1 or n3.GetFormalCharge() != -1:
            continue
        bond_n1_n2 = mol.GetBondBetweenAtoms(atom.GetIdx(), n2.GetIdx())
        bond_n2_n3 = mol.GetBondBetweenAtoms(n2.GetIdx(), n3.GetIdx())
        if bond_n1_n2.GetBondTypeAsDouble() != 2.0 or bond_n2_n3.GetBondTypeAsDouble() != 2.0:
            continue
        matches.append(atom)
    return matches


def has_azide_shape(mol) -> bool:
    return bool(_azide_root_nitrogens(mol))


def name_azide(mol) -> str:
    root_nitrogens = _azide_root_nitrogens(mol)
    if not root_nitrogens:
        raise UnsupportedStructure("no azide (-N3) group found; this module only handles azides")

    azide_atom_idxs = set()
    for n1 in root_nitrogens:
        (n2,) = (n for n in n1.GetNeighbors() if n.GetAtomicNum() == 7)
        (n3,) = (n for n in n2.GetNeighbors() if n.GetIdx() != n1.GetIdx() and n.GetAtomicNum() == 7)
        azide_atom_idxs.update((n1.GetIdx(), n2.GetIdx(), n3.GetIdx()))

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _AZIDE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an azide group's own three "
                "nitrogens (P-61.7) and halogen substituents (P-35.2.1) "
                "are not supported yet"
            )
        if atomic_num == 7 and atom.GetIdx() not in azide_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain azide group is "
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
    if any(a not in azide_atom_idxs and b not in azide_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.7's scope "
            "here is limited to a saturated chain)"
        )

    terminals = dict(halogen_substituents(mol))
    for n1 in root_nitrogens:
        terminals[n1.GetIdx()] = "azido"

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals)
