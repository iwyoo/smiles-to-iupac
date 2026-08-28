"""Naming of thiols (the '-thiol' suffix, -SH) on acyclic saturated or
unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-63.1.1 / P-65.1.1, Table 3.3 (Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'thiol' is the preselected
  suffix for -SH, the sulfur analogue of 'ol' — structurally the same
  substitutive-suffix construction, just with S instead of O, and (per
  Table 3.3's overall seniority order: acids > esters > amides > nitriles >
  aldehydes > ketones > alcohols > thiols > amines) ranked directly below
  alcohols and above amines. Unlike 'ol', 'thiol' begins with a consonant,
  so the parent hydride's final 'e' is never elided before it (P-16.3.3):
  'methane' + 'thiol' -> 'methanethiol', not 'methanthiol'.
- P-44.1.1 / P-44.4.1 / P-45.2 (Chapter P-4): same principal-chain and
  numbering machinery as `_alcohol.py` (maximum number of -SH groups first,
  then -SH locants ahead of ene/yne locants ahead of substituent-prefix
  locants) — this module deliberately mirrors that module's structure so
  the two stay easy to compare.
- P-14.3.4.2(a)/(b) (Chapter P-1): the same locant-omission rules as
  `_alcohol.py` apply here too (mononuclear parent, or a homogeneous
  two-carbon chain with exactly one substituent in total), e.g.
  'ethanethiol'.
- P-35.2.1 (Chapter P-3): halogen substituents are prefix-only and coexist
  freely with the -SH suffix, same as in `_alcohol.py`.

Scope, deliberately narrow (first pass at this functional group, mirroring
how `_amide.py`/`_nitrile.py`/etc. each started in isolation before any
cross-suffix seniority work): one or more -SH groups on an acyclic chain
(P-63.1.1's dithiol/trithiol/... multiplication, mirroring `_alcohol.py`'s
polyol support), with no other heteroatom (in particular no -OH or amine
nitrogen) anywhere in the molecule -- Table 3.3's alcohol/thiol/amine
seniority coexistence is future work, tracked as a separate roadmap item,
same as the analogous `multi-carbonyl-seniority.md` split for aldehyde/
ketone. Explicitly out of scope (raise `UnsupportedStructure`):
monocyclic/polycyclic/spiro rings, a sulfide (-S- ether-analogue) or any
other sulfur-oxidation-state group (sulfonic acid, etc.), and any oxygen
or nitrogen atom at all.
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
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 16, *HALOGEN_PREFIXES}


def _validate_and_collect_thiols(mol):
    thiols = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a thiol sulfur (P-63.1.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting -OH or amine nitrogen needs "
                "Table 3.3 seniority-coexistence handling not yet "
                "implemented for thiols"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num == 16:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a sulfur bonded to more than one heavy atom (e.g. a "
                    "sulfide) is out of scope; only an isolated thiol "
                    "(-SH) is supported (Table 3.3, P-63.1.1)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "a sulfur double-bonded to carbon is not a thiol"
                )
            if atom.GetTotalNumHs() != 1:
                raise UnsupportedStructure(
                    "an -S- atom that isn't a simple thiol (-SH) is out of "
                    "scope for this module"
                )
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() != 6:
                raise UnsupportedStructure("a thiol must be attached to a carbon atom")
            thiols.add(atom.GetIdx())
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
    if not thiols:
        raise UnsupportedStructure("no thiol (-SH) group found; this module only handles thiols")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return thiols


def _reject_enethiol_carbons(graph, thiols, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for s_idx in thiols:
        (carbon,) = graph[s_idx]
        if carbon in unsaturated_atoms:
            raise UnsupportedStructure(
                "a thiol on a carbon that is also part of a C=C/C#C bond "
                "is out of scope for this module"
            )


def _multiplied_word(count, base):
    """P-16.3.3: a multiplying prefix's terminal 'a' is elided before a
    suffix beginning with 'a' or 'o' (see `_common.py`'s `multiplied_word`
    docstring for the confirmed examples this mirrors)."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def _suffix_body(ene_locants, yne_locants, sh_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append((sorted(sh_locants), _multiplied_word(len(sh_locants), "thiol")))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, sh_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited.
        thiol_word = _multiplied_word(len(sh_locants), "thiol")
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + thiol_word

    if chain_length == 2 and not has_unsaturation and total_subs == 0 and len(sh_locants) == 1:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (here, the sole -SH) in total omits the locant, e.g.
        # 'ethanethiol'.
        return alkane_name(2) + "thiol"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, sh_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, sh_locants, ene_locants, yne_locants, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    sh_locant_set = lowest_locant_set(sh_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, sh_locants, ene_locants, yne_locants, grouped)
    return (
        (
            sh_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


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


def _bond_locant(chain, bond_atoms):
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def _bond_locants(chain, bonds):
    ene, yne = [], []
    for a, b, order in bonds:
        locant = _bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == _ENE_ORDER else yne).append(locant)
    return ene, yne


def _sh_locants(position_of, thiols, graph):
    locants = []
    for s in thiols:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, thiols):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in thiols]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def name_thiol(mol) -> str:
    thiols = _validate_and_collect_thiols(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enethiol_carbons(graph, thiols, bonds)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "cyclic thiols are not supported yet (this module only "
            "handles acyclic chains)"
        )

    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _sh_locants(position_of, thiols, graph) is None:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every thiol-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -SH groups, or an -SH expressed as a "
            "'sulfanyl' substituent prefix, is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sh_locants = _sh_locants(position_of, thiols, graph)
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, thiols)
            key, name = _candidate_key(chain_length, sh_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def has_thiol_shape(mol) -> bool:
    return any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
