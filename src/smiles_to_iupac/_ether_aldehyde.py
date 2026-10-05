"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate aldehyde (-CHO) on an acyclic saturated carbon skeleton, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py`/`_ether_thiol.py`/`_ether_ketone.py`
  document in full: an ether oxygen is *always* the 'R-oxy' substituent
  prefix (P-63.2.2.1.1), never the parent, so there is no seniority
  competition to resolve -- still wired through
  `_coexisting_groups.name_via_senior_acyclic` for its formal
  `_seniority.senior_class` assertion, mirroring the other migrated
  pairwise modules. This module mirrors `_ether_ketone.py`'s structure,
  swapping in `_aldehyde.py`'s own chain-naming machinery.
- `COCC=O` -> PubChem's own '2-methoxyacetaldehyde' confirms the aldehyde
  is always the suffix parent, the ether always the 'alkoxy' prefix.
- `_aldehyde.py` gained the identical `carbon_graph` parameter
  `_thiol.py`/`_ketone.py` gained (PR #428/#429), for the identical
  reason: an ether oxygen isn't itself a carbon, so its alkoxy branch
  would otherwise form a separate component that could wrongly outrank
  the real aldehyde-bearing chain in the global-diameter longest-chain
  search.

Scope, deliberately narrow (mirrors `_ether_ketone.py`): exactly one
plain ether oxygen (both sides acyclic saturated carbon) plus exactly one
aldehyde (-CHO, its carbonyl carbon has exactly one other carbon
neighbor, single-bonded), on one acyclic *saturated* skeleton, halogens
allowed. Explicitly out of scope (raise `UnsupportedStructure`): more
than one ether oxygen or aldehyde, any other heteroatom (including a
coexisting hydroxyl -- `_aldehyde.py` already supports that pairing on
its own, but the three-way combination is deferred), any ring, any chain
unsaturation (ene/yne) besides the aldehyde's own C=O, and any specified
stereocenter.
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._aldehyde import _name_acyclic_aldehyde
from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    find_ether_oxygens,
    non_single_bonds,
    specified_stereocenters,
    validate_allowed_atoms,
)
from ._ether import _oxy_prefix
from ._substituents import name_branch


def _find_aldehydes(mol):
    aldehydes = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6 or carbon.GetIsAromatic():
            continue
        if carbon.GetDegree() != 2:
            # A real aldehyde carbon has exactly two heavy-atom neighbors
            # (the carbonyl oxygen and one carbon, leaving room for its
            # one implicit H) -- degree 3 means a third heavy substituent
            # (e.g. an ether oxygen directly on this same carbon), which
            # is really an ester's acyl carbon, not an aldehyde.
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        if mol.GetBondBetweenAtoms(carbon.GetIdx(), carbon_neighbors[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        aldehydes.add(atom.GetIdx())
    return aldehydes


def has_ether_aldehyde_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ethers = find_ether_oxygens(mol)
    if len(ethers) != 1:
        return False
    aldehydes = _find_aldehydes(mol)
    if len(aldehydes) != 1:
        return False
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    return total_oxygens == 2


def _validate(mol, ether_oxygen, aldehyde_oxygen):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a plain ether oxygen (P-63.2.1) and an "
        "aldehyde carbonyl oxygen (P-33.3) are not supported yet",
        [
            (
                8,
                {ether_oxygen, aldehyde_oxygen},
                "an oxygen that isn't the single plain ether oxygen or the "
                "single aldehyde carbonyl is out of scope for this module",
            ),
        ],
    )


def name_ether_aldehyde(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/aldehyde combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ethers = find_ether_oxygens(mol)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with an aldehyde "
            "is supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    aldehydes = _find_aldehydes(mol)
    if len(aldehydes) != 1:
        raise UnsupportedStructure(
            "exactly one aldehyde (-CHO) coexisting with the single ether "
            "is supported here (see _aldehyde.py for a plain aldehyde)"
        )
    (aldehyde_oxygen,) = aldehydes
    _validate(mol, ether_oxygen_idx, aldehyde_oxygen)

    other_unsaturation = [b for b in non_single_bonds(mol) if aldehyde_oxygen not in (b[0], b[1])]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/aldehyde "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/aldehyde "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    (aldehyde_carbon,) = full_graph[aldehyde_oxygen]

    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    aldehyde_component, _ = bfs(full_carbon_graph, aldehyde_carbon)
    main_side = [c for c in ether_carbons if c in aldehyde_component]
    other_side = [c for c in ether_carbons if c not in aldehyde_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the aldehyde's own carbon "
            "skeleton and a separate alkoxy branch for this module"
        )
    (r_prime_carbon,) = other_side

    r_prime_component, _ = bfs(full_carbon_graph, r_prime_carbon)
    main_carbon_graph = {
        k: [n for n in v if n not in r_prime_component]
        for k, v in full_carbon_graph.items()
        if k not in r_prime_component
    }

    sub_name, sub_compound = name_branch(full_graph, r_prime_carbon, ether_oxygen_idx, {}, mol=mol)
    oxy_term = _oxy_prefix(sub_name)
    if sub_compound:
        oxy_term = enclose(oxy_term)
    extra_names = {ether_oxygen_idx: oxy_term}

    return name_via_senior_acyclic(
        _name_acyclic_aldehyde,
        "aldehyde",
        "ether",
        (mol, aldehydes, set(), ()),
        extra_names,
        carbon_graph=main_carbon_graph,
    )
