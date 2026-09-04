"""Naming of nitro compounds (the 'nitro' substituent prefix, -NO2) on
acyclic saturated hydrocarbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-61.5.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  '-NO2' is named with the simple substituent prefix 'nitro', the same way
  a halogen is (P-35.2.1) -- confirmed by the worked example 'nitromethane
  (PIN)' for CH3-NO2. Unlike a characteristic-group suffix, 'nitro' has no
  parent-hydride-selection role and no seniority of its own: it coexists
  freely with halogen substituents as an ordinary prefix (P-61.5.1's own
  worked examples show this for aromatic rings; confirmed here for acyclic
  chains via PubChem, e.g. '1-chloro-2-nitroethane').
- This module reuses `_acyclic.py`'s `name_from_carbon_graph` exactly the
  way `_ether.py`/`_peroxide.py` do: each nitro group's nitrogen is passed
  in as a precomputed 'terminals' entry (a leaf substituent excluded from
  the carbon-only chain search), merged with the ordinary halogen
  terminals dict so the two combine through the same multiplying-prefix/
  alphabetization machinery (P-14.5.2) that already handles 'chloro' +
  'nitro' or multiple 'nitro' groups ('dinitro', ...).

Explicitly out of scope (raise `UnsupportedStructure`), matching how
`_peroxide.py`/`_ether.py` scope their own single new heteroatom-group
extension:
- Nitroso (-NO, P-61.5.1's sibling prefix) -- a different, unimplemented
  group shape.
- Any other heteroatom besides a nitro group's own N/O and (optionally)
  halogen substituents -- in particular, this module doesn't attempt any
  Table 3.3 seniority coexistence with a characteristic-group suffix
  (P-61.5.2's own scope, which this module doesn't implement).
- Any unsaturation, any ring, or aromatic rings (a separate module's
  territory, e.g. '2-nitronaphthalene').
- aci-Nitro tautomers (P-61.5.3) -- a structurally different (=N(O)OH)
  group this module doesn't recognize at all.
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

_NITRO_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _nitro_nitrogen_atoms(mol):
    """Nitrogen atoms shaped like a nitro group: bonded to exactly one
    carbon and two oxygens, one double-bonded (terminal, neutral) and one
    single-bonded (terminal, -1 formally charged) -- the representation
    RDKit gives '[N+](=O)[O-]' -- or, equivalently, two double-bonded
    oxygens (a neutral N(=O)=O tautomer some SMILES writers use)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        if any(o.GetDegree() != 1 for o in oxygens):
            continue
        bond_orders = sorted(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() for o in oxygens
        )
        if bond_orders == [2.0, 2.0]:
            if atom.GetFormalCharge() != 0:
                continue
        elif bond_orders == [1.0, 2.0]:
            charges = sorted(o.GetFormalCharge() for o in oxygens)
            if charges != [-1, 0] or atom.GetFormalCharge() != 1:
                continue
        else:
            continue
        matches.append(atom)
    return matches


def has_nitro_shape(mol) -> bool:
    return bool(_nitro_nitrogen_atoms(mol))


def _validate_nitro_atoms(mol, nitro_atom_idxs, ring_atoms):
    """Shared per-atom validation for both the plain-chain path and the
    single-benzene-ring path: every atom outside `ring_atoms` (empty for
    the plain-chain path) must be a chain carbon, a halogen, or one of a
    nitro group's own N/O atoms."""
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _NITRO_ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a nitro group's own N/O "
                "(P-61.5.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atomic_num == 7 and idx not in nitro_atom_idxs:
            raise UnsupportedStructure(
                "a nitrogen atom not shaped like a plain nitro group is "
                "out of scope for this module (see e.g. nitroso, P-61.5.1)"
            )
        if atomic_num == 8 and idx not in nitro_atom_idxs:
            raise UnsupportedStructure(
                "an oxygen atom not belonging to a nitro group is out of "
                "scope for this module"
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


def _name_benzene_ring_nitro_chain(mol, nitro_nitrogens, nitro_atom_idxs, ring_atoms) -> str:
    """P-44.1.2.2 rule (1): a single, otherwise-unsubstituted benzene ring
    is senior to a chain of the same (plain-hydrocarbon) class regardless
    of the chain's length -- since 'nitro' has no suffix form, there is no
    characteristic group here to force the chain to be parent instead
    (mirrors `_azide.py`'s `_name_benzene_ring_azide_chain`, and confirmed
    directly here too: PubChem CID 54143589 '7-nitroheptylbenzene', not
    '1-nitroheptylbenzene' or any chain-parent form). So the ring is
    always the parent hydride, and the chain (with its nitro/halogen
    substituents) is named as a single substituent prefix on it via
    `name_branch`."""
    _validate_nitro_atoms(mol, nitro_atom_idxs, ring_atoms)

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a nitro chain is not supported yet"
        )
    ring_atom, chain_root = attachment

    non_ring_unsaturation = [
        (a, b)
        for a, b, _ in non_single_bonds(mol)
        if a not in nitro_atom_idxs
        and b not in nitro_atom_idxs
        and a not in ring_atoms
        and b not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "nitro chain is not supported yet"
        )

    terminals = dict(halogen_substituents(mol))
    for n in nitro_nitrogens:
        terminals[n.GetIdx()] = "nitro"

    branch_name, is_compound = name_branch(graph, chain_root, ring_atom, terminals)
    display = f"({branch_name})" if is_compound else branch_name
    return f"{display}benzene"


def name_nitro(mol) -> str:
    nitro_nitrogens = _nitro_nitrogen_atoms(mol)
    if not nitro_nitrogens:
        raise UnsupportedStructure("no nitro (-NO2) group found; this module only handles nitro compounds")

    nitro_atom_idxs = set()
    for n in nitro_nitrogens:
        nitro_atom_idxs.add(n.GetIdx())
        nitro_atom_idxs.update(o.GetIdx() for o in n.GetNeighbors() if o.GetAtomicNum() == 8)

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_benzene_ring_nitro_chain(mol, nitro_nitrogens, nitro_atom_idxs, ring_atoms)
        raise UnsupportedStructure("rings are not supported by this module yet")
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    _validate_nitro_atoms(mol, nitro_atom_idxs, frozenset())
    if any(a not in nitro_atom_idxs and b not in nitro_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-61.5.1's "
            "scope here is limited to a saturated chain)"
        )

    terminals = dict(halogen_substituents(mol))
    for n in nitro_nitrogens:
        terminals[n.GetIdx()] = "nitro"

    return name_from_carbon_graph(adjacency(mol), carbon_adjacency(mol), terminals)
