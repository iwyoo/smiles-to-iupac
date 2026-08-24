"""Naming of symmetric acyclic imides (R-CO-NH-CO-R', R = R'), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-66.6.3, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-66.6.3; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): an acyclic imide is named
  substitutively, as the N-acyl derivative of a primary amide -- one acyl
  branch supplies the '-amide' parent (`_amide.py`'s own suffix
  construction) and the other supplies an 'N-...oyl' acyl substituent
  prefix (the same '-oyl' stem `_acyl_halide.py` builds), e.g.
  'N-ethanoylethanamide' for (CH3CO)2NH (a cyclic imide such as
  succinimide instead uses a completely different heterocyclic-dione axis,
  out of scope here, see the module the roadmap tracks it under). Since
  this module only ever handles a *symmetric* imide, which acyl branch
  plays which role is arbitrary; both are computed and their names
  compared to enforce symmetry (see below).
- Each acyl carbon (-C(=O)-N(H)-) is a chain terminus exactly like
  `_amide.py`'s own amide carbon: after its carbonyl oxygen and the
  bridging imide nitrogen, it has room for at most one more substituent
  (another chain carbon, or nothing), so it always becomes C1 of its own
  branch.
- P-14.3.3: neither acyl branch's own locant is ever cited, exactly as in
  `_amide.py`/`_carboxylic_acid.py`.
- An enclosing mark (parentheses) is needed around the N-acyl substituent
  only when its own name carries an internal locant/substituent (e.g.
  'N-(2-methylpropanoyl)-2-methylpropanamide'); a bare, unsubstituted acyl
  group needs none (e.g. 'N-ethanoylethanamide'). A hyphen separates the
  N-acyl prefix from the amide parent only when the amide name itself
  starts with a locant digit (the same case that needs the acyl
  substituent's own parentheses, by symmetry) -- otherwise the two words
  are written solid, mirroring 'N-acetylacetamide'-style examples.
- P-35.2.1: halogen substituents on either chain are prefix-only and
  coexist freely.

Scope, deliberately narrow (first pass at this functional group): only a
*symmetric* imide (both acyl branches name identically) built from two
simple, saturated acyclic chains -- no chain unsaturation (ene/yne), no
other heteroatom anywhere, no ring, and no coexisting carbonyl-bearing
group. An unsymmetric imide, N-substituted imides (an alkyl group on the
imide nitrogen instead of H), a cyclic imide (succinimide's ring, a
structurally different heterocyclic-dione shape), and Table 3.3 seniority
coexistence are future work. Explicitly out of scope (raise
`UnsupportedStructure`): rings, any chain unsaturation, an N-substituted
imide nitrogen, and two branches that name differently (unsymmetric).
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


def _imide_cores(mol):
    """List of (nitrogen_idx, acyl_carbon1, carbonyl_o1, acyl_carbon2,
    carbonyl_o2) for every -C(=O)-N(H)-C(=O)- pattern: a one-H nitrogen
    bonded to exactly two carbons that each also carry exactly one
    doubly-bonded, monovalent (terminal) oxygen."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetTotalNumHs() != 1:
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


def has_imide_shape(mol) -> bool:
    return bool(_imide_cores(mol))


def _validate_and_collect_imide(mol):
    cores = _imide_cores(mol)
    if not cores:
        raise UnsupportedStructure(
            "no imide (-C(=O)-NH-C(=O)-) group found; this module only "
            "handles imides"
        )
    if len(cores) > 1:
        raise UnsupportedStructure(
            "more than one imide group is out of scope for this module"
        )
    imide_n, acyl1, carbonyl_o1, acyl2, carbonyl_o2 = cores[0]
    imide_atoms = {imide_n, carbonyl_o1, carbonyl_o2}

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an imide's own nitrogen/oxygens "
                "(P-66.6.3) and halogen substituents (P-35.2.1) are not "
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
        elif atomic_num == 7:
            if atom.GetIdx() != imide_n:
                raise UnsupportedStructure(
                    "a nitrogen other than the imide's own bridging "
                    "nitrogen needs Table 3.3 seniority handling not yet "
                    "implemented here"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in imide_atoms:
                raise UnsupportedStructure(
                    "an oxygen other than the imide's own two carbonyl "
                    "oxygens (e.g. a coexisting ketone/aldehyde/carboxylic "
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
                "a valid imide acyl carbon"
            )
    if len(non_single_bonds(mol)) != 2:
        # Exactly the imide's own two C=O bonds are always present;
        # anything else is chain unsaturation, out of scope for this module.
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) is out of scope for this module"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a cyclic imide (e.g. succinimide's ring) is a structurally "
            "different shape and out of scope for this module"
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


def _candidate_key(substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return locant_set, citation_locants


def _best_branch(chains, graph, halogens, acyl_carbon, carbonyl_oxygen, excluded):
    """Return (chain_length, grouped substituents) for the best (lowest
    locant set) numbering of `acyl_carbon`'s own branch."""
    eligible = [c for c in chains if acyl_carbon in c]
    if not eligible:
        raise UnsupportedStructure(
            "the imide's acyl carbon does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )
    chain_length = len(eligible[0])

    best_key = None
    best_grouped = None
    for chain in eligible:
        # `acyl_carbon` has at most one carbon neighbor (see
        # `_validate_and_collect_imide`), so it can only ever be a chain
        # endpoint, never interior.
        candidate = chain if chain[0] == acyl_carbon else list(reversed(chain))
        substituents = _substituents_for_chain(graph, candidate, halogens, excluded | {carbonyl_oxygen})
        key = _candidate_key(substituents)
        if best_key is None or key < best_key:
            best_key, best_grouped = key, _group(substituents)
    return chain_length, best_grouped


def _acid_stem_name(chain_length, grouped):
    """'propanoic'-style stem, used only to compare the two branches for
    structural symmetry (see module docstring)."""
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "oic"


def _amide_name(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "amide"


def _acyl_prefix_name(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    is_compound = bool(prefix)
    name = prefix + stem + "oyl"
    return name, is_compound


def name_imide(mol) -> str:
    acyl1, acyl2 = _validate_and_collect_imide(mol)
    cores = _imide_cores(mol)
    imide_n, _, carbonyl_o1, _, carbonyl_o2 = cores[0]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {imide_n}
    chains = _longest_chains(carbon_adjacency(mol))
    length1, grouped1 = _best_branch(chains, graph, halogens, acyl1, carbonyl_o1, excluded)
    length2, grouped2 = _best_branch(chains, graph, halogens, acyl2, carbonyl_o2, excluded)

    if _acid_stem_name(length1, grouped1) != _acid_stem_name(length2, grouped2):
        raise UnsupportedStructure(
            "an unsymmetric imide (the two acyl groups name differently) "
            "is out of scope for this module"
        )

    amide_name = _amide_name(length1, grouped1)
    acyl_name, acyl_is_compound = _acyl_prefix_name(length2, grouped2)
    acyl_part = f"({acyl_name})" if acyl_is_compound else acyl_name
    separator = "-" if amide_name[0].isdigit() else ""
    return f"N-{acyl_part}{separator}{amide_name}"
