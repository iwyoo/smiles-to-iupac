"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate thiol (-SH) on an acyclic saturated carbon skeleton, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py` documents in full: an ether oxygen is
  *always* the 'R-oxy' substituent prefix (P-63.2.2.1.1), never the
  parent, so there is no seniority competition to resolve. This module
  mirrors `_ether_amine.py`'s structure, swapping in `_thiol.py`'s own
  chain-naming machinery.
- `COCCS` -> PubChem's own '2-methoxyethanethiol' confirms the thiol is
  always the suffix parent, the ether always the 'alkoxy' prefix.
- `_alcohol.py` already supports a coexisting ether natively (its own
  `_ether_oxygens`/`_alkoxy_name` helpers, pre-dating this pattern) --
  this module is the second pilot extending the same idea, this time as a
  separate pairwise module since `_thiol.py` itself is left unmodified
  except for one new optional parameter (see below).

`_thiol.py`'s `_name_acyclic_thiol` already had `extra_names`/
`required_atoms` extension points (from the earlier amine work), but it
computed its own carbon-only chain-search graph from the whole molecule
internally -- an ether oxygen (not itself a carbon) leaves its alkoxy
branch as a separate connected component in that graph, which could wrongly
outrank the real thiol-bearing chain in the global-diameter longest-chain
search. `_thiol.py` gained one new optional `carbon_graph` parameter
(default `None`, meaning "compute it the old way" -- every existing
caller is unaffected) so this module can pass a graph with the alkoxy
branch's component already removed, mirroring how `_amide.py` already
excludes its N-alkyl substituent's component before its own chain search.

Scope, deliberately narrow (mirrors `_ether_amine.py`): exactly one plain
ether oxygen (both sides acyclic saturated carbon) plus exactly one thiol
(-SH), on one acyclic *saturated* skeleton, halogens allowed. Explicitly
out of scope (raise `UnsupportedStructure`): more than one ether oxygen or
thiol, any other heteroatom, any ring, any chain unsaturation (ene/yne),
and any specified stereocenter.
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
from ._substituents import name_branch
from ._thiol import _name_acyclic_thiol

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}


def _find_ether_oxygens(mol):
    """The molecule's sole oxygen, if it's a plain ether -- see
    `_ether_amine.py`'s identical helper for why requiring the *only*
    oxygen avoids misfiring on an ester/carbamate/lactone's own bridging
    oxygen."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxygens) != 1:
        return []
    (oxygen,) = oxygens
    if oxygen.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in oxygen.GetNeighbors()):
        return [oxygen]
    return []


def _find_thiols(mol):
    thiols = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 1:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            thiols.add(atom.GetIdx())
    return thiols


def has_ether_thiol_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    if len(_find_ether_oxygens(mol)) != 1:
        return False
    return len(_find_thiols(mol)) == 1


def _validate(mol, ether_oxygen, thiol_sulfur):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a plain ether oxygen (P-63.2.1) "
                "and a thiol sulfur (P-63.1.1) are not supported yet"
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
            if atom.GetIdx() != ether_oxygen:
                raise UnsupportedStructure(
                    "an oxygen that isn't the single plain ether oxygen is "
                    "out of scope for this module"
                )
        elif atomic_num == 16:
            if atom.GetIdx() != thiol_sulfur:
                raise UnsupportedStructure(
                    "a sulfur that isn't a single plain thiol (-SH) is out "
                    "of scope for this module"
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


def name_ether_thiol(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/thiol combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ethers = _find_ether_oxygens(mol)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with a thiol is "
            "supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    thiols = _find_thiols(mol)
    if len(thiols) != 1:
        raise UnsupportedStructure(
            "exactly one thiol (-SH) coexisting with the single ether is "
            "supported here (see _thiol.py for a plain thiol)"
        )
    (thiol_sulfur,) = thiols
    _validate(mol, ether_oxygen_idx, thiol_sulfur)

    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/thiol "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/thiol "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    (thiol_carbon,) = full_graph[thiol_sulfur]

    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    thiol_component, _ = bfs(full_carbon_graph, thiol_carbon)
    main_side = [c for c in ether_carbons if c in thiol_component]
    other_side = [c for c in ether_carbons if c not in thiol_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the thiol's own carbon "
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

    return _name_acyclic_thiol(mol, thiols, (), carbon_graph=main_carbon_graph, extra_names=extra_names)
