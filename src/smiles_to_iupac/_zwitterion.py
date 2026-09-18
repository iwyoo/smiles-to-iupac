"""Naming of an amino-acid/betaine-type zwitterion -- a single-fragment,
net-neutral molecule carrying exactly one ammonium-shaped (+1) nitrogen
and exactly one carboxylate-shaped (-1) oxygen pair on the same acyclic
carbon chain -- per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-74.1.3: "anionic centers are preferred for lower locants and become
  the parent structure, into which the cationic part is substituted [as a
  prefix]." This project's mechanism: the carboxylate is always the
  parent (`_carboxylate.py`), the ammonium nitrogen is always cited as a
  substituent prefix.
- The ammonium-nitrogen prefix itself is built by isolating the nitrogen
  and its own substituents into a standalone fragment (replacing the bond
  toward the rest of the chain with an extra hydrogen -- the same "cut a
  bond, add an H" trick this project's other demoted-group modules use
  for their own coexisting-group carbons), naming that isolated fragment
  with `_ammonium.py`'s own `name_ammonium` (already the established PIN
  derivation for every substitution level, per that module's own
  docstring -- explicitly NOT the alternative direct 'azanium'-
  multiplicative substitutive style), then appending 'yl' -- e.g.
  glycine's (`C(C(=O)[O-])[NH3+]`) isolated nitrogen fragment is plain
  NH4+ -> 'azanium' -> 'azaniumyl'; a mono-substituted nitrogen isolates
  to R-NH3+ -> '<R>aminium' -> '<R>aminiumyl'; betaine's
  (`C[N+](C)(C)CC(=O)[O-]`) trimethyl-substituted nitrogen isolates to
  (CH3)3NH+ (trimethylammonium) -> `name_ammonium`'s own
  'N,N-dimethylmethanaminium' -> 'N,N-dimethylmethanaminiumyl'.
- Reuses `_carboxylate.py`'s `_name_acyclic_carboxylate` directly (its
  `extra_names` extension point, mirroring `_carboxylic_acid_amine.py`'s
  identical 'amino' injection) rather than `_coexisting_groups.py`'s
  `name_via_senior_acyclic`: an ammonium prefix isn't a suffix-vs-suffix
  seniority demotion at all (ammonium has no entry in
  `_seniority.SUFFIX_CLASS_RANK`) -- P-74's zwitterion citation order is
  its own, separate mechanism (anion always the parent, full stop, no
  seniority comparison to make).

Scope, deliberately narrow (M1 of this project's zwitterion coverage):
exactly one +1 ammonium-shaped nitrogen and exactly one -1 carboxylate
oxygen pair, both on one connected, saturated, acyclic all-carbon-chain-
plus-one-nitrogen fragment, net formal charge 0, no other heteroatom,
halogen, or charge.
Explicitly out of scope (raise `UnsupportedStructure`): the ionic center
in a ring, more than one of either ionic center, the ammonium nitrogen
bonded directly to the carboxylate carbon itself (P-74.1.1/1.2's
same-parent case -- a different mechanism, not handled here), any chain
unsaturation, and any specified stereocenter.
"""

from rdkit import Chem

from ._ammonium import has_ammonium_shape, name_ammonium
from ._carboxylate import _find_carboxylate_group, _name_acyclic_carboxylate
from ._common import UnsupportedStructure, adjacency, bfs, specified_stereocenters


def has_zwitterion_shape(mol) -> bool:
    """True if `mol` looks like an amino-acid/betaine-type zwitterion --
    used by `core.py` to route here before `has_salt_shape`, since
    `_salt.py`'s cation loop would otherwise crash trying (and failing) to
    name the whole single-fragment molecule as a bare ammonium cation."""
    if len(Chem.GetMolFrags(mol)) > 1:
        return False
    if not has_ammonium_shape(mol):
        return False
    try:
        _find_carboxylate_group(mol)
    except UnsupportedStructure:
        return False
    return True


def _chain_neighbor(mol, nitrogen_idx, carboxylate_carbon_idx):
    """The one carbon neighbor of the nitrogen that lies on the same side
    as the carboxylate carbon (as opposed to a plain N-alkyl substituent
    that goes nowhere else) -- the molecule is a tree once the nitrogen
    itself is removed, so exactly one neighbor's component contains the
    carboxylate carbon."""
    graph = adjacency(mol)
    graph_without_nitrogen = {
        atom: [n for n in neighbors if n != nitrogen_idx]
        for atom, neighbors in graph.items()
        if atom != nitrogen_idx
    }
    for neighbor in graph[nitrogen_idx]:
        dist, _ = bfs(graph_without_nitrogen, neighbor)
        if carboxylate_carbon_idx in dist:
            return neighbor
    return None


def _ammonium_prefix(mol, nitrogen_idx, chain_neighbor_idx):
    """Isolate the ammonium nitrogen and its own substituents (replacing
    the bond toward the rest of the chain with an extra hydrogen) and name
    that standalone fragment via `_ammonium.py`'s own `name_ammonium`,
    then append 'yl' (see module docstring)."""
    rw = Chem.RWMol(mol)
    nitrogen = rw.GetAtomWithIdx(nitrogen_idx)
    total_hs = nitrogen.GetTotalNumHs()
    rw.RemoveBond(nitrogen_idx, chain_neighbor_idx)
    nitrogen.SetNoImplicit(True)
    nitrogen.SetNumExplicitHs(total_hs + 1)
    isolated_mol = rw.GetMol()
    Chem.SanitizeMol(isolated_mol)

    mapping = []
    fragments = Chem.GetMolFrags(isolated_mol, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping)
    (nitrogen_fragment,) = (
        fragment for fragment, atom_indices in zip(fragments, mapping) if nitrogen_idx in atom_indices
    )
    prefix = name_ammonium(nitrogen_fragment) + "yl"
    # A substituted nitrogen's own name carries its own locants (e.g.
    # 'N,N-dimethylmethanaminiumyl') -- `_substituents.name_branch` always
    # reports an injected `extra_names` entry as non-compound (it has no
    # way to know otherwise), so this module wraps it in parentheses
    # itself before injection, the same way `_ether.py`'s callers wrap a
    # compound alkoxy substituent (P-29.4).
    if "-" in prefix:
        prefix = f"({prefix})"
    return prefix


def name_zwitterion(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ionic center in a ring uses a different naming "
            "construction, out of scope for this acyclic-only module"
        )

    (nitrogen,) = (atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1)
    nitrogen_idx = nitrogen.GetIdx()
    carboxylate_carbon, carbonyl_oxygen, anion_oxygen = _find_carboxylate_group(mol)
    carboxylate_carbon_idx = carboxylate_carbon.GetIdx()
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), anion_oxygen.GetIdx()}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx == nitrogen_idx or idx in excluded_oxygens:
            continue
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms other than the single ammonium nitrogen and "
                "the single carboxylate's own oxygens are not supported "
                "yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "ammonium nitrogen and carboxylate anion oxygen is not "
                "supported yet"
            )
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")

    if any(bond.GetBondTypeAsDouble() not in (1.0, 2.0) for bond in mol.GetBonds()):
        raise UnsupportedStructure("a bond order other than single or the carbonyl double bond is not supported")
    non_carbonyl_unsaturation = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} != {carboxylate_carbon_idx, carbonyl_oxygen.GetIdx()}
    ]
    if non_carbonyl_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a zwitterion's own "
            "carboxylate/ammonium pair is out of scope for this module "
            "(1st-pass scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure("a specified stereocenter alongside a zwitterion is not supported yet")

    chain_neighbor_idx = _chain_neighbor(mol, nitrogen_idx, carboxylate_carbon_idx)
    if chain_neighbor_idx is None or chain_neighbor_idx == carboxylate_carbon_idx:
        raise UnsupportedStructure(
            "an ammonium nitrogen bonded directly to the carboxylate "
            "carbon itself uses a different naming construction (P-74.1.1/"
            "1.2's same-parent case), out of scope for this module"
        )

    prefix = _ammonium_prefix(mol, nitrogen_idx, chain_neighbor_idx)
    return _name_acyclic_carboxylate(
        mol,
        carboxylate_carbon_idx,
        excluded_oxygens,
        (),
        extra_names={nitrogen_idx: prefix},
        required_atoms={chain_neighbor_idx},
    )
