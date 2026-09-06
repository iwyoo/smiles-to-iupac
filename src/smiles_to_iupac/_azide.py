"""Naming of azides (the 'azido' substituent prefix, -N3) on acyclic
saturated hydrocarbon chains, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-61.7 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): a
  '-N3' (-N=N+=N-) group attached to a parent hydride is named with the
  simple substituent prefix 'azido' -- the same structural role as 'nitro'
  (`_nitro.py`, P-61.5.1), which this module mirrors closely. Confirmed by
  the worked example '(2-azidoethyl)benzene (PIN)'.
- A single, otherwise-unsubstituted benzene ring with one exocyclic chain
  carrying the azide group(s) is also supported (that same worked
  example): 'azido' has no suffix form, so no characteristic-group
  seniority ever forces the chain to be the parent (P-41) -- the ordinary
  P-44.1.2.2 ring-vs-chain rule applies, and rule (1) there makes the ring
  senior to a chain of the same (plain-hydrocarbon) class *regardless of
  the chain's length* ('heptylbenzene (PIN)' is the worked example for the
  unsubstituted case -- not '1-phenylheptane', even though the chain has
  more skeletal atoms). So, unlike the benzene-ring slice already done for
  -COOH/ketone/etc. (where the chain's suffix forces it to stay parent),
  there is no "chain wins when longer" sub-case to support here -- the
  ring always wins.
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
- Any unsaturation elsewhere in the molecule, or any ring other than the
  single-benzene-ring case above (more than one ring, a non-benzene ring,
  a substituted/fused benzene ring, more than one exocyclic ring
  attachment, or a branched chain hanging off the ring).
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


def _azide_atom_idxs(root_nitrogens):
    azide_atom_idxs = set()
    for n1 in root_nitrogens:
        (n2,) = (n for n in n1.GetNeighbors() if n.GetAtomicNum() == 7)
        (n3,) = (n for n in n2.GetNeighbors() if n.GetIdx() != n1.GetIdx() and n.GetAtomicNum() == 7)
        azide_atom_idxs.update((n1.GetIdx(), n2.GetIdx(), n3.GetIdx()))
    return azide_atom_idxs


def _validate_azide_atoms(mol, azide_atom_idxs, ring_atoms):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path: every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a chain carbon, a halogen, or one of the
    azide group's own three nitrogens."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _AZIDE_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an azide group's own three "
                "nitrogens (P-61.7) and halogen substituents (P-35.2.1) "
                "are not supported yet"
            )
        if atomic_num == 7 and idx not in azide_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain azide group is "
                "out of scope for this module"
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


def _name_benzene_ring_azide_chain(mol, root_nitrogens, azide_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): a single, otherwise-unsubstituted benzene ring
    is senior to a chain of the same (plain-hydrocarbon) class regardless
    of the chain's length -- since 'azido' has no suffix form, there is no
    characteristic group here to force the chain to be parent instead
    (contrast the suffix modules' 'benzene ring is always a phenyl
    substituent' pattern). So the ring is always the parent hydride, and
    the chain (with its azide/halogen substituents) is named as a single
    substituent prefix on it via `name_branch` -- mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`, but with
    the ring/chain roles reversed."""
    _validate_azide_atoms(mol, azide_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside an azide chain is not supported yet"
        )
    ring_atom, chain_root = attachment

    non_ring_unsaturation = [
        (a, b)
        for a, b, _ in non_single_bonds(mol)
        if a not in azide_atom_idxs
        and b not in azide_atom_idxs
        and a not in ring_atoms
        and b not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "azide chain is not supported yet"
        )

    terminals = dict(halogen_substituents(mol))
    for n1 in root_nitrogens:
        terminals[n1.GetIdx()] = "azido"

    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals, mol=mol)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_azide(mol) -> str:
    root_nitrogens = _azide_root_nitrogens(mol)
    if not root_nitrogens:
        raise UnsupportedStructure("no azide (-N3) group found; this module only handles azides")
    azide_atom_idxs = _azide_atom_idxs(root_nitrogens)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_azide_chain(mol, root_nitrogens, azide_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_azide_atoms(mol, azide_atom_idxs, frozenset())
    if any(a not in azide_atom_idxs and b not in azide_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.7's scope "
            "here is limited to a saturated chain)"
        )

    terminals = dict(halogen_substituents(mol))
    for n1 in root_nitrogens:
        terminals[n1.GetIdx()] = "azido"

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals, mol=mol)
