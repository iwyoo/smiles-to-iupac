"""Naming of a molecule combining exactly one primary amide (-C(=O)NH2) with
one or more internal ketone (=O) carbonyls on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'amide' (`_amide.py`) outranks 'one' (`_ketone.py`), so a
  coexisting ketone is demoted to the 'oxo' substituent prefix instead of
  its own '-one' suffix, e.g. 'CC(=O)CC(N)=O' (acetoacetamide) ->
  '3-oxobutanamide' (a well-known worked example). This mirrors
  `_aldehyde_ketone.py`'s aldehyde+ketone demotion, reusing the same
  {atom_idx -> 'oxo'} injection into the chain's substituent-prefix
  machinery.
- Otherwise mirrors `_amide.py` exactly: the amide carbon is always a chain
  terminus and always becomes C1 with its own locant never cited
  (P-14.3.3); a ketone carbon's locant, by contrast, is always cited.

Scope, deliberately narrow (first extension of the general suffix-vs-suffix
demotion problem beyond aldehyde+ketone): a single primary amide plus one
or more ketones, all on one acyclic saturated chain, with halogen substituents
allowed. Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring, an N-substituted amide, more than one
amide, a coexisting hydroxyl/ether/other heteroatom, and any ketone not
captured by a single longest chain.
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

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _find_amide_carbon(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        carbonyls = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        amide_nitrogens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 7
            and n.GetDegree() == 1
            and n.GetTotalNumHs() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbonyls) == 1 and len(amide_nitrogens) == 1 and len(carbon_neighbors) <= 1:
            return atom, carbonyls[0], amide_nitrogens[0]
    return None


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


def has_ketone_amide_shape(mol) -> bool:
    found = _find_amide_carbon(mol)
    if found is None:
        return False
    _, amide_oxygen, _ = found
    return bool(_extra_ketones(mol, {amide_oxygen.GetIdx()}))


def _validate(mol, excluded_oxygens, amide_nitrogen):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the amide's own oxygen/nitrogen and "
                "ketone carbonyl oxygens (P-41, Table 3.3) and halogen "
                "substituents (P-35.2.1) are not supported yet"
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
                "an oxygen that isn't the amide's own carbonyl oxygen or a "
                "ketone-shaped carbonyl is out of scope for this module "
                "(e.g. a coexisting hydroxyl or ether)"
            )
        elif atomic_num == 7:
            if atom.GetIdx() != amide_nitrogen.GetIdx():
                raise UnsupportedStructure(
                    "more than one nitrogen, or a nitrogen that isn't the "
                    "amide's own, is out of scope for this module"
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
    return prefix + stem + "amide"


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


def name_ketone_amide(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amide on/in a ring (a lactam) is out of scope for this "
            "acyclic-only module"
        )
    found = _find_amide_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no primary amide (-CONH2) group found; this module only "
            "handles primary amides"
        )
    amide_carbon, amide_oxygen, amide_nitrogen = found
    ketones = _extra_ketones(mol, {amide_oxygen.GetIdx()})
    if not ketones:
        raise UnsupportedStructure(
            "no coexisting ketone found; this module only handles an amide "
            "combined with at least one ketone (see _amide.py for a plain "
            "amide)"
        )
    own_excluded = {amide_oxygen.GetIdx(), amide_nitrogen.GetIdx()}
    excluded = own_excluded | ketones
    _validate(mol, excluded, amide_nitrogen)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a ketone-on-amide "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    amide_carbon_idx = amide_carbon.GetIdx()
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if amide_carbon_idx not in chain:
            continue
        chain_set = set(chain)
        if any(graph[o][0] not in chain_set for o in ketones):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every amide/ketone-bearing carbon lies on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != amide_carbon_idx:
                # The amide carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            substituents = _substituents_for_chain(graph, candidate, names, own_excluded)
            grouped = _group(substituents)
            key, name = _candidate_key(chain_length, grouped)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
