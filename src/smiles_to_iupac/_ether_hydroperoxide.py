"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate hydroperoxide (-OOH) on an acyclic saturated carbon skeleton, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py`/`_ether_thiol.py`/`_ether_ketone.py`/
  `_ether_aldehyde.py`/`_ether_amide.py` document in full: an ether
  oxygen is *always* the 'R-oxy' substituent prefix (P-63.2.2.1.1), never
  the parent, so there is no seniority competition to resolve -- still
  wired through `_coexisting_groups.name_via_senior_acyclic` for its
  formal `_seniority.senior_class` assertion, mirroring the other
  migrated pairwise modules. This module mirrors `_ether_thiol.py`'s
  structure, swapping in `_hydroperoxide.py`'s own chain-naming
  machinery.
- **Note on PubChem divergence** (same situation `_hydroperoxide_amine.py`
  already documents, PR #426): PubChem's own auto-generated IUPACName for
  'COCCOO' is "1-hydroperoxy-2-methoxyethane" -- naming the molecule as a
  substituted ethane with *both* groups as prefixes, not using the
  'peroxol' suffix at all. This contradicts P-41/P-43's own text (a
  hydroperoxide is always cited via the 'peroxol' suffix when it's the
  senior/only suffix-eligible group present, per `_hydroperoxide.py`'s
  own module docstring and the plain 'ethaneperoxol' PIN). This module
  follows the Blue Book primary source directly, the same policy
  `_hydroperoxide_amine.py` already established for this exact class's
  PubChem-divergence pattern.
- `_hydroperoxide.py` gained the identical `carbon_graph` parameter
  `_thiol.py`/`_ketone.py`/`_aldehyde.py` gained (PR #428-430), for the
  identical reason.

Scope, deliberately narrow (mirrors `_ether_thiol.py`): exactly one plain
ether oxygen (both sides acyclic saturated carbon) plus exactly one
hydroperoxide (-OOH), on one acyclic *saturated* skeleton, halogens
allowed. Explicitly out of scope (raise `UnsupportedStructure`): more
than one ether oxygen or hydroperoxide, any other heteroatom, any ring,
any chain unsaturation (ene/yne), and any specified stereocenter.
"""

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
from ._hydroperoxide import _name_acyclic_hydroperoxide
from ._substituents import name_branch

def _find_hydroperoxides(mol):
    """(attach, terminal) oxygen atoms of a plain -O-O-H hydroperoxide,
    or None -- like `_hydroperoxide.py`'s own `_hydroperoxide_oxygens`,
    but without that function's "exactly two oxygens in the whole
    molecule" restriction (which would reject this module's own
    coexisting third, ether-shaped oxygen)."""
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    pairs = []
    for o1 in oxygens:
        for o2 in o1.GetNeighbors():
            if o2.GetAtomicNum() != 8 or o2.GetIdx() <= o1.GetIdx():
                continue
            bond = mol.GetBondBetweenAtoms(o1.GetIdx(), o2.GetIdx())
            if bond.GetBondTypeAsDouble() != 1.0:
                continue
            degrees = sorted((o1.GetDegree(), o2.GetDegree()))
            if degrees != [1, 2]:
                continue
            attach, terminal = (o1, o2) if o1.GetDegree() == 2 else (o2, o1)
            if terminal.GetTotalNumHs() != 1:
                continue
            (other,) = (n for n in attach.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
            if other.GetAtomicNum() != 6:
                continue
            pairs.append((attach, terminal))
    if len(pairs) != 1:
        return None
    return pairs[0]


def has_ether_hydroperoxide_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ethers = find_ether_oxygens(mol)
    if len(ethers) != 1:
        return False
    if _find_hydroperoxides(mol) is None:
        return False
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    return total_oxygens == 3


def _validate(mol, ether_oxygen, hydroperoxide_atoms):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a plain ether oxygen (P-63.2.1) and a "
        "hydroperoxide's own two oxygens (P-56.1) are not supported yet",
        [
            (
                8,
                {ether_oxygen} | set(hydroperoxide_atoms),
                "an oxygen that isn't the single plain ether oxygen or "
                "part of the single hydroperoxide's -O-O-H pair is out of "
                "scope for this module",
            ),
        ],
    )


def name_ether_hydroperoxide(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/hydroperoxide combination on/in a ring uses a "
            "different naming construction, out of scope for this "
            "acyclic-only module"
        )
    ethers = find_ether_oxygens(mol)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with a "
            "hydroperoxide is supported here (see _ether.py for a plain "
            "ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    oxygens = _find_hydroperoxides(mol)
    if oxygens is None:
        raise UnsupportedStructure(
            "exactly one hydroperoxide (-OOH) coexisting with the single "
            "ether is supported here (see _hydroperoxide.py for a plain "
            "hydroperoxide)"
        )
    attach, terminal = oxygens
    hydroperoxide_atoms = {attach.GetIdx(), terminal.GetIdx()}
    _validate(mol, ether_oxygen_idx, hydroperoxide_atoms)

    other_unsaturation = [b for b in non_single_bonds(mol) if b[0] not in hydroperoxide_atoms and b[1] not in hydroperoxide_atoms]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/"
            "hydroperoxide coexistence is out of scope for this module "
            "(1st-pass scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/hydroperoxide "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    (site,) = (n for n in attach.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
    site_idx = site.GetIdx()

    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    site_component, _ = bfs(full_carbon_graph, site_idx)
    main_side = [c for c in ether_carbons if c in site_component]
    other_side = [c for c in ether_carbons if c not in site_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the hydroperoxide's own "
            "carbon skeleton and a separate alkoxy branch for this module"
        )
    (r_prime_carbon,) = other_side

    r_prime_component, _ = bfs(full_carbon_graph, r_prime_carbon)
    main_carbon_graph = {
        k: [n for n in v if n not in r_prime_component]
        for k, v in full_carbon_graph.items()
        if k not in r_prime_component
    }

    sub_name, sub_compound = name_branch(full_graph, r_prime_carbon, ether_oxygen_idx, {}, mol=mol)
    oxy_term = _oxy_prefix(sub_name, sub_compound)
    extra_names = {ether_oxygen_idx: oxy_term}

    return name_via_senior_acyclic(
        _name_acyclic_hydroperoxide,
        "hydroperoxide",
        "ether",
        (mol, site_idx, hydroperoxide_atoms),
        extra_names,
        carbon_graph=main_carbon_graph,
    )
