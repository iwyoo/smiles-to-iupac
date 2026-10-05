"""Naming of a molecule combining a plain ether (-O-R') with an ester
(R-CO-O-R''), the ether sitting on the ester's acyl (R) chain, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py`/`_ether_thiol.py`/`_ether_ketone.py`/
  `_ether_aldehyde.py`/`_ether_amide.py`/`_ether_hydroperoxide.py`
  document in full: an ether oxygen is *always* the 'R-oxy' substituent
  prefix (P-63.2.2.1.1), never the parent, so there is no seniority
  competition to resolve.
- `COCC(=O)OC` -> PubChem's own 'methyl 2-methoxyacetate' confirms the
  ester is always the suffix parent, the ether always the 'alkoxy'
  prefix on the acyl chain.
- Wired through `_coexisting_groups.name_via_senior_acyclic` like every
  other module in this milestone, via `_ester.py`'s own
  `_name_acyclic_ester` entry point (a thin join of its existing
  `_name_alcohol_part`/`_name_acyl_part` calls). Unlike every earlier
  pilot in this series, `_ester.py`'s `_name_acyl_part` needed **no**
  `carbon_graph` parameter addition: it already isolates the acyl
  carbon's own reachable component via
  `_component_subgraph(carbon_graph, acyl_carbon_idx)` before searching
  for the principal chain, and an ether oxygen (not itself a carbon)
  already breaks that reachability at the ether, so the alkoxy branch's
  component is automatically excluded with zero changes to `_ester.py`.
  Only the existing `extra_names` parameter (added in PR #425) is used
  here.

Scope, deliberately narrow: exactly one plain ether oxygen on the acyl
(R) chain, alongside an otherwise-ordinary ester whose alcohol part (R')
is a plain unsubstituted acyclic alkyl (the same restriction
`_ester.py`'s own acyclic path already applies). Explicitly out of scope
(raise `UnsupportedStructure`): an ether on the alcohol part instead,
more than one ether oxygen, any other heteroatom, any ring, any chain
unsaturation (ene/yne), and any specified stereocenter.
"""

from rdkit import Chem

from ._multiplicative_text import enclose
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
from ._ester import _name_acyclic_ester
from ._ether import _oxy_prefix
from ._substituents import name_branch


def _find_ester_group(mol):
    """Locate a single ester group's (acyl_carbon, carbonyl_oxygen,
    ester_oxygen, alcohol_carbon) atoms, or None -- like `_ester.py`'s own
    `_find_ester_group`, but without that function's "exactly two oxygens
    in the whole molecule" restriction (which would reject this module's
    own coexisting third, ether-shaped oxygen)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if len(carbonyls) == 1 and len(ester_oxygens) == 1:
            acyl_carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(acyl_carbon_neighbors) > 1:
                continue
            matches.append((atom, carbonyls[0], ester_oxygens[0]))
    if len(matches) != 1:
        return None
    acyl_carbon, carbonyl_oxygen, ester_oxygen = matches[0]
    alcohol_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != acyl_carbon.GetIdx())
    return acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon

def has_ether_ester_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ester = _find_ester_group(mol)
    if ester is None:
        return False
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = ester
    ethers = find_ether_oxygens(mol, exclude={ester_oxygen.GetIdx()})
    if len(ethers) != 1:
        return False
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    return total_oxygens == 3


def _validate(mol, ether_oxygen, carbonyl_oxygen, ester_oxygen):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a plain ether oxygen (P-63.2.1) and the "
        "ester's own carbonyl/ester-oxygen pair (P-65.6.3) are not "
        "supported yet",
        [
            (
                8,
                {ether_oxygen, carbonyl_oxygen, ester_oxygen},
                "an oxygen that isn't the single plain ether oxygen or the "
                "ester's own carbonyl/ester-oxygen pair is out of scope "
                "for this module",
            ),
        ],
    )


def name_ether_ester(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/ester combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ester = _find_ester_group(mol)
    if ester is None:
        raise UnsupportedStructure(
            "exactly one ester group coexisting with a separate ether is "
            "supported here (see _ester.py for a plain ester)"
        )
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = ester

    ethers = find_ether_oxygens(mol, exclude={ester_oxygen.GetIdx()})
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with an ester is "
            "supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    _validate(mol, ether_oxygen_idx, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx())

    excluded = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    other_unsaturation = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/ester "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/ester "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    acyl_component, _ = bfs(full_carbon_graph, acyl_carbon.GetIdx())
    main_side = [c for c in ether_carbons if c in acyl_component]
    other_side = [c for c in ether_carbons if c not in acyl_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the ester's own acyl "
            "chain and a separate alkoxy branch for this module"
        )
    (r_prime_carbon,) = other_side
    if alcohol_carbon.GetIdx() in bfs(full_carbon_graph, r_prime_carbon)[0]:
        raise UnsupportedStructure(
            "an ether on the alcohol part (R') rather than the acyl chain "
            "is not supported yet"
        )

    sub_name, sub_compound = name_branch(full_graph, r_prime_carbon, ether_oxygen_idx, {}, mol=mol)
    oxy_term = _oxy_prefix(sub_name)
    if sub_compound:
        oxy_term = enclose(oxy_term)
    extra_names = {ether_oxygen_idx: oxy_term}

    return name_via_senior_acyclic(
        _name_acyclic_ester,
        "ester",
        "ether",
        (mol, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon),
        extra_names,
    )
