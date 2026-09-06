"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate ketone (>C=O) on an acyclic saturated carbon skeleton, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py`/`_ether_thiol.py` document in full: an
  ether oxygen is *always* the 'R-oxy' substituent prefix (P-63.2.2.1.1),
  never the parent, so there is no seniority competition to resolve. This
  module mirrors `_ether_thiol.py`'s structure, swapping in `_ketone.py`'s
  own chain-naming machinery.
- `COCC(C)=O` -> PubChem's own '1-methoxypropan-2-one' confirms the
  ketone is always the suffix parent, the ether always the 'alkoxy'
  prefix.
- `_alcohol.py` already supports a coexisting ether natively (pre-dating
  this pattern); `_thiol.py` was the first pairwise-module pilot (PR
  #428, adding one new optional `carbon_graph` parameter). `_ketone.py`
  gained the identical parameter here, for the identical reason: an
  ether oxygen isn't itself a carbon, so its alkoxy branch would
  otherwise form a separate component that could wrongly outrank the
  real ketone-bearing chain in the global-diameter longest-chain search.

Scope, deliberately narrow (mirrors `_ether_thiol.py`): exactly one plain
ether oxygen (both sides acyclic saturated carbon) plus exactly one
ketone (>C=O, both carbonyl-carbon neighbors themselves carbon), on one
acyclic *saturated* skeleton, halogens allowed. Explicitly out of scope
(raise `UnsupportedStructure`): more than one ether oxygen or ketone, any
other heteroatom (including a coexisting hydroxyl -- `_ketone.py` already
supports that pairing on its own, but the three-way ether+hydroxyl+ketone
combination is deferred), any ring, any chain unsaturation (ene/yne), and
any specified stereocenter.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    non_single_bonds,
    specified_stereocenters,
)
from ._ether import _oxy_prefix
from ._ketone import _name_acyclic_ketone
from ._substituents import name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _find_ether_oxygens(mol):
    """The molecule's sole oxygen, if it's a plain ether -- see
    `_ether_amine.py`'s identical helper for why requiring the *only*
    oxygen avoids misfiring on an ester/carbamate/lactone's own bridging
    oxygen. Since a ketone molecule may have a *second* oxygen (its own
    carbonyl), this module instead requires exactly one *ether-shaped*
    oxygen among possibly several, with the other(s) validated as the
    ketone's own carbonyl in `_validate` below."""
    return [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 8
        and atom.GetDegree() == 2
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    ]


def _find_ketones(mol):
    ketones = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6 or carbon.GetIsAromatic():
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            ketones.add(atom.GetIdx())
    return ketones


def has_ether_ketone_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ethers = _find_ether_oxygens(mol)
    if len(ethers) != 1:
        return False
    ketones = _find_ketones(mol)
    if len(ketones) != 1:
        return False
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    return total_oxygens == 2


def _validate(mol, ether_oxygen, ketone_oxygen):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a plain ether oxygen (P-63.2.1) "
                "and a ketone carbonyl oxygen (P-33.4) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in (ether_oxygen, ketone_oxygen):
                raise UnsupportedStructure(
                    "an oxygen that isn't the single plain ether oxygen or "
                    "the single ketone carbonyl is out of scope for this "
                    "module"
                )
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")


def name_ether_ketone(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/ketone combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ethers = _find_ether_oxygens(mol)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with a ketone is "
            "supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    ketones = _find_ketones(mol)
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "exactly one ketone (>C=O) coexisting with the single ether is "
            "supported here (see _ketone.py for a plain ketone)"
        )
    (ketone_oxygen,) = ketones
    _validate(mol, ether_oxygen_idx, ketone_oxygen)

    other_unsaturation = [b for b in non_single_bonds(mol) if ketone_oxygen not in (b[0], b[1])]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/ketone "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/ketone "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    (ketone_carbon,) = full_graph[ketone_oxygen]

    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    ketone_component, _ = bfs(full_carbon_graph, ketone_carbon)
    main_side = [c for c in ether_carbons if c in ketone_component]
    other_side = [c for c in ether_carbons if c not in ketone_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the ketone's own carbon "
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
    if sub_compound:
        sub_name = f"({sub_name})"
    oxy_term = _oxy_prefix(sub_name)
    extra_names = {ether_oxygen_idx: oxy_term}

    return _name_acyclic_ketone(mol, ketones, set(), (), carbon_graph=main_carbon_graph, extra_names=extra_names)
