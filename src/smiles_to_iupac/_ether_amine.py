"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate primary amine (-NH2) on an acyclic saturated carbon skeleton, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- unlike every
  other X+amine pairwise module in this project (`_alcohol_amine.py` etc.,
  which each demote the *junior* class via `_seniority.senior_class`),
  there is no seniority competition here: an ether oxygen is *always* the
  'R-oxy' substituent prefix (P-63.2.2.1.1), never the parent, regardless
  of anything else in the molecule. So this module reuses `_ether.py`'s
  existing `_oxy_prefix` naming helper directly and calls straight into
  `_amine.py`'s own `_best_chain_name`/`_amine_locants` (which already
  accept an arbitrary {atom_idx -> prefix name} dict for exactly this
  purpose -- the same mechanism other single-group modules expose as
  `extra_names`), rather than going through `_coexisting_groups.py`'s
  senior/junior assertion (which doesn't apply here).
- `NCCOO` sibling comparison: PubChem PUG REST confirms 'COCCN' ->
  '2-methoxyethanamine' and 'CCOCCN' -> '2-ethoxyethanamine' -- the amine
  is always the suffix parent, the ether always the 'alkoxy' prefix.
- The ether oxygen's two carbon neighbors sit in different connected
  components of the carbon-only graph (an ether oxygen is not itself a
  carbon, so `carbon_adjacency` naturally cuts the skeleton there); this
  module identifies which side holds the amine (the main chain) and
  excludes the other side's whole component from the chain search before
  calling `_amine.py`'s machinery, the same "remove the substituent's own
  component before picking the principal chain" pattern `_amide.py` uses
  for its N-alkyl substituents.

Scope, deliberately narrow (first pilot of the general "ether coexists
with any suffix module" fix):
exactly one plain ether oxygen (both sides acyclic saturated carbon) plus
exactly one primary amine, on one acyclic *saturated* skeleton, halogens
allowed. Explicitly out of scope (raise `UnsupportedStructure`): more than
one ether oxygen, any other heteroatom, any ring, any chain unsaturation
(ene/yne), any specified stereocenter, a secondary/tertiary amine, and an
ether oxygen bonded to the amine's own nitrogen-bearing carbon in a way
that leaves it off the main chain entirely (not expected for this shape,
defensive only).
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._amine import _best_chain_name
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    find_ether_oxygens,
    find_primary_amines,
    halogen_substituents,
    non_single_bonds,
    specified_stereocenters,
    validate_allowed_atoms,
)
from ._ether import _oxy_prefix
from ._substituents import name_branch


def has_ether_amine_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ethers = find_ether_oxygens(mol, require_sole=True)
    if len(ethers) != 1:
        return False
    amines = find_primary_amines(mol, set())
    return len(amines) == 1


def _validate(mol, ether_oxygen, amine_nitrogen):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a plain ether oxygen (P-63.2.1) and a "
        "primary amine nitrogen (P-41, Table 4.1) are not supported yet",
        [
            (
                8,
                {ether_oxygen},
                "an oxygen that isn't the single plain ether oxygen is out "
                "of scope for this module",
            ),
            (
                7,
                {amine_nitrogen},
                "a nitrogen that isn't a single plain primary amine (-NH2) "
                "is out of scope for this module (secondary/tertiary "
                "amines, imines, and nitriles are not supported)",
            ),
        ],
    )


def name_ether_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ethers = find_ether_oxygens(mol, require_sole=True)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with a primary "
            "amine is supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    amines = find_primary_amines(mol, set())
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "ether is supported here (see _amine.py for a plain amine)"
        )
    (amine_nitrogen,) = amines
    _validate(mol, ether_oxygen_idx, amine_nitrogen)

    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/amine "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/amine "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    (amine_carbon,) = full_graph[amine_nitrogen]

    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    amine_component, _ = bfs(full_carbon_graph, amine_carbon)
    main_side = [c for c in ether_carbons if c in amine_component]
    other_side = [c for c in ether_carbons if c not in amine_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the amine's own carbon "
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
    halogens = {**halogen_substituents(mol), ether_oxygen_idx: oxy_term}

    best_name, _ = _best_chain_name(main_carbon_graph, full_graph, halogens, {amine_nitrogen}, ())
    return best_name
