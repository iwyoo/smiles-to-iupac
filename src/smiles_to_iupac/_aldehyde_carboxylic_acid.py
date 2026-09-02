"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more internal aldehyde (-CHO) carbonyls on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'oic acid' (`_carboxylic_acid.py`) outranks 'al'
  (`_aldehyde.py`), so a coexisting aldehyde is demoted to the 'oxo'
  substituent prefix instead of its own '-al' suffix, e.g.
  'O=CCC(=O)O' (malonaldehydic acid) -> '3-oxopropanoic acid' (a
  well-known worked example). This mirrors `_aldehyde_ketone.py`'s
  aldehyde+ketone demotion, reusing the same {atom_idx -> 'oxo'}
  injection into the chain's
  substituent-prefix machinery -- an aldehyde carbon that is part of the
  parent chain (rather than a branch) is cited as 'oxo', not 'formyl',
  the same way a demoted ketone is.
- Otherwise mirrors `_carboxylic_acid.py` exactly: the -COOH carbon is
  always a chain terminus and always becomes C1 with its own locant never
  cited (P-14.3.3); an aldehyde carbon's locant, by contrast, is always
  cited.

Scope, deliberately narrow (first extension of the general suffix-vs-suffix
demotion problem beyond aldehyde+ketone): a single carboxylic acid plus one
or more aldehydes, all on one acyclic saturated chain, with halogen
substituents allowed. Explicitly out of scope (raise `UnsupportedStructure`):
any chain unsaturation (ene/yne), any ring, more than one carboxylic acid, a
coexisting hydroxyl/ether/other heteroatom, and any aldehyde not captured by
a single longest chain.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    path_between,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _find_carboxylic_acid_carbon(mol):
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
        hydroxyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
        ]
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbonyls) == 1 and len(hydroxyls) == 1 and len(carbon_neighbors) <= 1:
            return atom, carbonyls[0], hydroxyls[0]
    return None


def _extra_aldehydes(mol, excluded_oxygens):
    aldehydes = set()
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
        # A genuine aldehyde carbon has exactly one oxygen neighbor (its own
        # carbonyl); a second oxygen neighbor means this is really a -COOH
        # (or other acid-shaped) carbon, which `_find_carboxylic_acid_carbon`
        # already accounts for and must not be double-counted here.
        oxygen_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 1:
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            aldehydes.add(atom.GetIdx())
    return aldehydes


def has_aldehyde_carboxylic_acid_shape(mol) -> bool:
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return False
    _, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    return bool(_extra_aldehydes(mol, excluded))


def _validate(mol, excluded_oxygens):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the acid's own oxygens and aldehyde "
                "carbonyl oxygens (P-41, Table 3.3) and halogen substituents "
                "(P-35.2.1) are not supported yet"
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
            if atom.GetIdx() in excluded_oxygens:
                continue
            raise UnsupportedStructure(
                "an oxygen that isn't part of the single carboxylic acid "
                "group or an aldehyde-shaped carbonyl is out of scope for "
                "this module (e.g. a coexisting hydroxyl or ether)"
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


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "oic acid"


def _candidate_key(chain_length, grouped):
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def _longest_chains(graph):
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def _substituents_for_chain(graph, chain, names, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names) for root in branch_roots]
    return substituents


def name_aldehyde_carboxylic_acid(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH group on/in a ring uses a different naming "
            "construction, out of scope for this acyclic-only module"
        )
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    aldehydes = _extra_aldehydes(mol, excluded_acid_oxygens)
    if not aldehydes:
        raise UnsupportedStructure(
            "no coexisting aldehyde found; this module only handles a "
            "carboxylic acid combined with at least one aldehyde (see "
            "_carboxylic_acid.py for a plain carboxylic acid)"
        )
    excluded = excluded_acid_oxygens | aldehydes
    _validate(mol, excluded)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an aldehyde-on-acid "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in aldehydes}}
    acid_carbon_idx = acid_carbon.GetIdx()
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if acid_carbon_idx not in chain:
            continue
        chain_set = set(chain)
        if any(graph[o][0] not in chain_set for o in aldehydes):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every acid/aldehyde-bearing carbon lies on a single "
            "longest carbon chain; a shorter principal chain is not "
            "supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acid_carbon_idx:
                # The -COOH carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            substituents = _substituents_for_chain(graph, candidate, names, excluded_acid_oxygens)
            grouped = _group(substituents)
            key, name = _candidate_key(chain_length, grouped)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
