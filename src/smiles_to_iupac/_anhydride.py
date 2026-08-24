"""Naming of symmetric acyclic carboxylic acid anhydrides (R-CO-O-CO-R',
R = R'), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-65.2, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.2; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): an anhydride's name is
  built from the name(s) of its constituent carboxylic acid(s), with the
  word 'acid' replaced by 'anhydride'. For a *symmetric* anhydride
  (R = R', this module's whole scope), that's simply the shared acid's own
  name, e.g. 'propanoic acid' -> 'propanoic anhydride'. An *unsymmetric*
  anhydride (R != R') instead cites both acid names (each with 'acid'
  dropped) in alphabetical order followed by a single 'anhydride' -- out of
  scope here, see below. Table 3.3 ranks anhydride senior to ester (junior
  only to carboxylic acid itself); this module doesn't implement that
  seniority (it only ever names an anhydride in isolation).
- Each acyl carbon (-C(=O)-O-) is a chain terminus exactly like
  `_carboxylic_acid.py`'s -COOH carbon: after its carbonyl oxygen and the
  bridging anhydride oxygen, it has room for at most one more substituent
  (another chain carbon, or nothing, for the formic/methanoic branch), so
  it always becomes C1 of its own branch, fixed rather than a
  locant-minimization choice.
- P-14.3.3: the acid stem's own locant is never cited, exactly as in
  `_carboxylic_acid.py` (e.g. 'propanoic anhydride', not
  'propan-1-oic anhydride').
- P-35.2.1: halogen substituents on either R chain are prefix-only and
  coexist freely.

Scope, deliberately narrow (first pass at this functional group): only a
*symmetric* anhydride (both acyl branches name identically) built from two
simple, saturated acyclic chains -- no chain unsaturation (ene/yne), no
other heteroatom anywhere, no ring, and no coexisting carbonyl-bearing
group. An unsymmetric anhydride, a cyclic anhydride (e.g. succinic
anhydride's ring, a structurally different two-carbonyl-in-one-ring
shape), and Table 3.3 seniority coexistence are future work. Explicitly
out of scope (raise `UnsupportedStructure`): rings, any chain
unsaturation, and two branches that name differently (unsymmetric).
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


def _anhydride_cores(mol):
    """List of (bridging_o_idx, acyl_carbon1, carbonyl_o1, acyl_carbon2,
    carbonyl_o2) for every -C(=O)-O-C(=O)- pattern: a degree-2, uncharged
    oxygen bridging two carbons that each also carry exactly one
    doubly-bonded, monovalent (terminal) oxygen."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 2 or atom.GetTotalNumHs() != 0:
            continue
        neighbors = atom.GetNeighbors()
        if any(n.GetAtomicNum() != 6 for n in neighbors):
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbors[1].GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue

        def _carbonyl_oxygen(carbon):
            carbonyls = [
                o
                for o in carbon.GetNeighbors()
                if o.GetAtomicNum() == 8
                and o.GetIdx() != atom.GetIdx()
                and o.GetDegree() == 1
                and mol.GetBondBetweenAtoms(carbon.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            ]
            return carbonyls[0].GetIdx() if len(carbonyls) == 1 else None

        c1, c2 = neighbors
        o1, o2 = _carbonyl_oxygen(c1), _carbonyl_oxygen(c2)
        if o1 is None or o2 is None:
            continue
        cores.append((atom.GetIdx(), c1.GetIdx(), o1, c2.GetIdx(), o2))
    return cores


def has_anhydride_shape(mol) -> bool:
    return bool(_anhydride_cores(mol))


def _validate_and_collect_anhydride(mol):
    cores = _anhydride_cores(mol)
    if not cores:
        raise UnsupportedStructure(
            "no acid anhydride (-C(=O)-O-C(=O)-) group found; this module "
            "only handles acid anhydrides"
        )
    if len(cores) > 1:
        raise UnsupportedStructure(
            "more than one acid anhydride group is out of scope for this "
            "module"
        )
    bridging_o, acyl1, carbonyl_o1, acyl2, carbonyl_o2 = cores[0]
    anhydride_atoms = {bridging_o, carbonyl_o1, carbonyl_o2}

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an anhydride's own oxygens "
                "(P-65.2) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in anhydride_atoms:
                raise UnsupportedStructure(
                    "an oxygen other than the anhydride's own three oxygens "
                    "(e.g. a coexisting ketone/aldehyde/carboxylic "
                    "acid/ester/amide) needs Table 3.3 seniority handling "
                    "not yet implemented here"
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
    for acyl in (acyl1, acyl2):
        carbon = mol.GetAtomWithIdx(acyl)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "an acyl carbon with more than one carbon neighbor is not "
                "a valid anhydride acyl carbon"
            )
    if len(non_single_bonds(mol)) != 2:
        # Exactly the anhydride's own two C=O bonds are always present;
        # anything else is chain unsaturation, out of scope for this module.
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) is out of scope for this module"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a cyclic anhydride (e.g. succinic anhydride's ring) is a "
            "structurally different shape and out of scope for this module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return acyl1, acyl2


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


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _acid_stem_name(chain_length, grouped):
    """'propanoic', '2-methylpropanoic', etc. -- the 'oic' word always
    begins with a vowel (no ene/yne segment is possible in this module's
    scope), so the parent stem's trailing 'e' is always elided, exactly as
    in `_carboxylic_acid.py`."""
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "oic"


def _candidate_key(substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return locant_set, citation_locants


def _branch_acid_name(chains, graph, halogens, acyl_carbon, carbonyl_oxygen, excluded):
    eligible = [c for c in chains if acyl_carbon in c]
    if not eligible:
        raise UnsupportedStructure(
            "the anhydride's acyl carbon does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )
    chain_length = len(eligible[0])

    best_key = None
    best_name = None
    for chain in eligible:
        # `acyl_carbon` has at most one carbon neighbor (see
        # `_validate_and_collect_anhydride`), so it can only ever be a
        # chain endpoint, never interior.
        candidate = chain if chain[0] == acyl_carbon else list(reversed(chain))
        substituents = _substituents_for_chain(graph, candidate, halogens, excluded | {carbonyl_oxygen})
        key = _candidate_key(substituents)
        if best_key is None or key < best_key:
            grouped = _group(substituents)
            best_key, best_name = key, _acid_stem_name(chain_length, grouped)
    return best_name


def name_anhydride(mol) -> str:
    acyl1, acyl2 = _validate_and_collect_anhydride(mol)
    cores = _anhydride_cores(mol)
    bridging_o, _, carbonyl_o1, _, carbonyl_o2 = cores[0]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {bridging_o}
    chains = _longest_chains(carbon_adjacency(mol))
    name1 = _branch_acid_name(chains, graph, halogens, acyl1, carbonyl_o1, excluded)
    name2 = _branch_acid_name(chains, graph, halogens, acyl2, carbonyl_o2, excluded)
    if name1 != name2:
        raise UnsupportedStructure(
            "an unsymmetric anhydride (the two acyl groups name "
            "differently) is out of scope for this module"
        )
    return name1 + " anhydride"
