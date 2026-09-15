"""Naming of a molecule combining exactly one ester (R-CO-O-R') with one or
more internal ketone (=O) carbonyls on the same acyclic saturated acyl
chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'oate' (`_ester.py`) outranks 'one' (`_ketone.py`), so a
  coexisting ketone on the acyl chain is demoted to the 'oxo' substituent
  prefix instead of its own '-one' suffix, e.g.
  'CC(=O)CC(=O)OC' (methyl acetoacetate) -> 'methyl 3-oxobutanoate'
  (confirmed against this well-known worked example). This mirrors
  `_aldehyde_ketone.py`'s aldehyde+ketone demotion, reusing the same
  {atom_idx -> 'oxo'} injection into the acyl chain's substituent-prefix
  machinery.
- Otherwise mirrors `_ester.py` exactly: the acyl carbon is always a chain
  terminus and always becomes C1 of the acyl chain (P-14.3.3); a ketone
  carbon's locant, by contrast, is always cited.

Scope, deliberately narrow (same acyl-chain restriction as `_ester.py`, plus
one or more ketones on that chain): a single ester group whose acyl part (R)
carries one or more ketone carbonyls, with halogen substituents allowed on
the acyl chain. The alcohol part (R') is restricted the same way
`_ester.py` restricts it (plain unbranched saturated alkyl). Explicitly out
of scope (raise `UnsupportedStructure`): any chain unsaturation (ene/yne) on
the acyl chain, any ring, more than one ester group, a coexisting
hydroxyl/ether/other heteroatom, an aldehyde-shaped (rather than
ketone-shaped) extra carbonyl, and any ketone not captured by the acyl
chain's single longest path.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    component_subgraph,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    linear_branch,
    longest_branched_chain,
    longest_chains,
    non_single_bonds,
    ring_chain_attachment,
    substituent_locant_set_and_citation,
    validate_allowed_atoms,
)
from ._numerals import alkane_name, alkyl_name
from ._substituents import format_substituent_prefixes, name_branch



def _find_ester_group(mol):
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
        if carbonyls and ester_oxygens:
            matches.append((atom, carbonyls, ester_oxygens))
    if len(matches) != 1:
        return None
    acyl_carbon, carbonyls, ester_oxygens = matches[0]
    if len(carbonyls) != 1 or len(ester_oxygens) != 1:
        return None
    acyl_carbon_neighbors = [n for n in acyl_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(acyl_carbon_neighbors) > 1:
        return None
    carbonyl_oxygen = carbonyls[0]
    ester_oxygen = ester_oxygens[0]
    alcohol_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != acyl_carbon.GetIdx())
    return acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon


def _extra_ketones(mol, excluded_oxygens):
    ketones = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in excluded_oxygens:
            continue
        if atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6:
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            ketones.add(atom.GetIdx())
    return ketones


def has_ketone_ester_shape(mol) -> bool:
    group = _find_ester_group(mol)
    if group is None:
        return False
    _, carbonyl_oxygen, ester_oxygen, _ = group
    excluded = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    return bool(_extra_ketones(mol, excluded))


def _validate(mol, acyl_carbon, carbonyl_oxygen, ester_oxygen, ketones, aromatic_ring_atoms=frozenset()):
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()} | ketones
    validate_allowed_atoms(
        mol,
        "heteroatoms other than the ester's own oxygens and ketone "
        "carbonyl oxygens (P-41, Table 3.3) and halogen substituents "
        "(P-35.2.1) are not supported yet",
        [
            (
                8,
                excluded_oxygens,
                "an oxygen that isn't part of the single ester group or a "
                "ketone-shaped carbonyl is out of scope for this module "
                "(e.g. a coexisting hydroxyl or ether)",
            ),
        ],
        aromatic_ring_atoms=aromatic_ring_atoms,
    )


def _alcohol_component(full_graph, alcohol_carbon_idx, ester_oxygen_idx):
    component = set()
    stack = [alcohol_carbon_idx]
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        for neighbor in full_graph[node]:
            if neighbor != ester_oxygen_idx and neighbor not in component:
                stack.append(neighbor)
    return component


def _name_alcohol_part(mol, alcohol_carbon, ester_oxygen_idx):
    full_graph = adjacency(mol)
    component = _alcohol_component(full_graph, alcohol_carbon.GetIdx(), ester_oxygen_idx)
    for idx in component:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "the alcohol part (R') must be a plain alkyl group; "
                "heteroatoms/halogens there are not supported yet (P-65.6.3)"
            )
    for a, b, _ in non_single_bonds(mol):
        if a in component and b in component:
            raise UnsupportedStructure(
                "unsaturation in the alcohol part (R') is not supported yet "
                "(P-65.6.3)"
            )
    carbon_graph = carbon_adjacency(mol)
    length = linear_branch(carbon_graph, alcohol_carbon.GetIdx(), None)
    if length is None:
        raise UnsupportedStructure(
            "a branched alcohol part (R') is not supported yet (P-65.6.3)"
        )
    return alkyl_name(length)


def _name_from_substituents(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "oate"


def _candidate_key(chain_length, grouped):
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def _substituents_for_chain(graph, chain, names, excluded_oxygens, ring_atoms=frozenset(), mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_oxygens]
        if branch_roots:
            substituents[position] = [
                name_branch(graph, root, atom, names, ring_atoms, mol=mol) for root in branch_roots
            ]
    return substituents


def _name_acyl_part(mol, acyl_carbon, carbonyl_oxygen_idx, ester_oxygen_idx, ketones):
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    acyl_carbon_idx = acyl_carbon.GetIdx()
    excluded_oxygens = {carbonyl_oxygen_idx, ester_oxygen_idx}

    acyl_graph = component_subgraph(carbon_graph, acyl_carbon_idx)
    chains = longest_chains(acyl_graph)
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if any(full_graph[o][0] not in chain_set for o in ketones):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every ketone-bearing carbon lies on the acyl chain's single "
            "longest carbon chain; a shorter principal chain is not "
            "supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acyl_carbon_idx:
                # The ester carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            substituents = _substituents_for_chain(full_graph, candidate, names, excluded_oxygens, mol=mol)
            grouped = group_substituents(substituents)
            key, name = _candidate_key(chain_length, grouped)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    if best_name is None:
        raise UnsupportedStructure(
            "the ester's acyl carbon does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )
    return best_name


def _name_phenyl_chain_acyl_part(mol, acyl_carbon, carbonyl_oxygen_idx, ester_oxygen_idx, ketones, ring_atoms):
    """Name the acyl part (R) of a ketone+ester combination whose acyl
    chain hangs off a single unbranched chain from one atom of an
    otherwise-plain, unsubstituted benzene ring -- e.g.
    'methyl 3-oxo-4-phenylbutanoate'. Mirrors `_name_acyl_part` above but
    routes the acyl chain through `ring_chain_attachment` instead of
    `longest_chains`, the same way `_aldehyde_ketone.py`'s
    `_name_phenyl_chain_aldehyde_ketone` extends its own chain-terminus
    acyl search."""
    full_graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    acyl_carbon_idx = acyl_carbon.GetIdx()
    excluded_oxygens = {carbonyl_oxygen_idx, ester_oxygen_idx}

    attachment = ring_chain_attachment(full_graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain ketone/ester combination is not supported "
            "yet"
        )
    chain, _ = longest_branched_chain(full_graph, acyl_carbon_idx, ring_atoms, excluded_oxygens | set(names), halogens=halogen_substituents(mol))
    if any(full_graph[o][0] not in chain for o in ketones):
        raise UnsupportedStructure(
            "not every ketone-bearing carbon lies on the chain hanging "
            "off the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an ester group directly attached to the benzene ring uses a "
            "separate naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = _substituents_for_chain(full_graph, chain, names, excluded_oxygens, ring_atoms, mol=mol)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def _name_phenyl_chain_ketone_ester(mol, ring_atoms):
    """Name a ketone+ester combination whose acyl chain (R) hangs off a
    plain, unsubstituted benzene ring -- e.g.
    'methyl 3-oxo-4-phenylbutanoate'. The alcohol part (R') is unaffected
    (unbranched-alkyl-only, same as the acyclic path)."""
    group = _find_ester_group(mol)
    if group is None:
        raise UnsupportedStructure(
            "exactly one ester group is required; zero or multiple ester "
            "groups (e.g. a diester) are not supported yet (P-65.6.3)"
        )
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = group
    excluded = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    ketones = _extra_ketones(mol, excluded)
    if not ketones:
        raise UnsupportedStructure(
            "no coexisting ketone found; this module only handles an ester "
            "combined with at least one ketone on the acyl chain (see "
            "_ester.py for a plain ester)"
        )
    _validate(mol, acyl_carbon, carbonyl_oxygen, ester_oxygen, ketones, aromatic_ring_atoms=ring_atoms)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [
        b for b in all_non_single if b[0] in (excluded | ketones) or b[1] in (excluded | ketones)
    ]
    ring_bonds = [b for b in all_non_single if b[0] in ring_atoms and b[1] in ring_atoms]
    if len(carbonyl_bonds) + len(ring_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a benzene-ring-"
            "substituent ketone/ester combination is out of scope for "
            "this module"
        )

    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_phenyl_chain_acyl_part(
        mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx(), ketones, ring_atoms
    )
    return f"{alcohol_name} {acyl_name}"


def name_ketone_ester(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_ketone_ester(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached ester or ketone is out of scope for this "
            "acyclic-only module"
        )
    group = _find_ester_group(mol)
    if group is None:
        raise UnsupportedStructure(
            "exactly one ester group is required; zero or multiple ester "
            "groups (e.g. a diester) are not supported yet (P-65.6.3)"
        )
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = group
    excluded = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    ketones = _extra_ketones(mol, excluded)
    if not ketones:
        raise UnsupportedStructure(
            "no coexisting ketone found; this module only handles an ester "
            "combined with at least one ketone on the acyl chain (see "
            "_ester.py for a plain ester)"
        )
    _validate(mol, acyl_carbon, carbonyl_oxygen, ester_oxygen, ketones)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [
        b for b in all_non_single if b[0] in (excluded | ketones) or b[1] in (excluded | ketones)
    ]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a ketone-on-ester "
            "demotion is out of scope for this module"
        )

    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_acyl_part(mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx(), ketones)
    return f"{alcohol_name} {acyl_name}"
