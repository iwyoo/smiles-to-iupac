"""Naming of acyl halides (the '-oyl halide' suffix, -C(=O)X for X = F, Cl,
Br, I) on acyclic saturated or unsaturated carbon chains, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-65.1.5, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.1.5; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): an acyl halide's -OH (of
  the parent -COOH) is replaced by a halogen; the suffix is the parent
  hydride's '-oic acid' stem with 'oic acid' replaced by 'oyl' plus a
  separate halide word ('fluoride'/'chloride'/'bromide'/'iodide'), e.g.
  'propanoyl chloride'. Table 3.3 ranks acyl halide junior to carboxylic
  acid/ester/amide but senior to nitrile/aldehyde/ketone/alcohol/amine;
  this module doesn't implement that seniority (see scope note below).
- A -C(=O)X carbon is a chain terminus (its remaining two bonds after the
  carbonyl oxygen and the halogen leave room for at most one more carbon
  neighbor), exactly like `_carboxylic_acid.py`'s -COOH carbon: it always
  becomes C1, never a locant-minimization choice.
- P-14.3.3: the acyl halide's own locant is never cited, e.g. 'propanoyl
  chloride' (not 'propan-1-oyl chloride'), mirroring '-oic acid'.
- P-35.2.1: any *other* halogen (not the one forming the acyl halide) is
  an ordinary prefix substituent and coexists freely, e.g.
  '2-chloropropanoyl chloride'.
- P-91.3/P-92: a molecule
  with one or more *specified* tetrahedral stereocenters on the principal
  chain gets a "(<locant><R/S>,...)-" prefix, ascending locant order, e.g.
  '(2R)-2-chloropropanoyl chloride' -- same mechanism as
  `_carboxylic_acid.py`, since the acyl halide carbon's own fixed C1
  position already decides numbering before stereo is considered. A
  stereocenter on a substituent branch, mixed with an unspecified one, or
  alongside E/Z double-bond stereo remains out of scope (raises
  `UnsupportedStructure` via `_common.specified_stereocenters`).

Scope, deliberately narrow (mirrors `_carboxylic_acid.py`'s own first
pass): only a single acyl halide on an acyclic chain, with no other
carbonyl-bearing characteristic group (ketone/aldehyde/carboxylic
acid/ester/amide) anywhere in the molecule -- that Table 3.3 seniority
coexistence is future work, tracked under `multi-carbonyl-seniority.md`.
Explicitly out of scope (raise `UnsupportedStructure`): rings and two or
more acyl halide groups, except for one narrow ring case: an acyl
halide's chain hanging off a single plain, unsubstituted benzene ring
with no other substituent on the ring (`_name_phenyl_chain_acyl_halide`,
e.g. '3-phenylpropanoyl chloride'), mirroring `_aldehyde.py`/`_amide.py`/
`_nitrile.py`'s identical benzene-ring-substituent path. Narrower than
the acyclic path: no chain unsaturation, no specified stereocenter, and
no substituted benzene/naphthalene -- each a separate follow-up.
"""

from rdkit import Chem

from ._common import (
    elides_before,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    lowest_locant_set,
    non_single_bonds,
    path_between,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name, numerical_term
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0
_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}

_HALIDE_WORDS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}


def _acyl_halide_carbons(mol):
    """{carbon_idx -> (carbonyl_oxygen_idx, halogen_idx)} for every carbon
    shaped like an acyl halide: exactly one monovalent, doubly-bonded
    oxygen and exactly one monovalent halogen."""
    matches = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        halogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES]
        if len(oxygens) != 1 or len(halogens) != 1:
            continue
        (oxygen,) = oxygens
        (halogen,) = halogens
        if oxygen.GetDegree() != 1 or halogen.GetDegree() != 1:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        matches[atom.GetIdx()] = (oxygen.GetIdx(), halogen.GetIdx())
    return matches


def has_acyl_halide_shape(mol) -> bool:
    return bool(_acyl_halide_carbons(mol))


def _validate_and_collect_acyl_halides(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_acyl_halide`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_acyl_halide`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged."""
    matches = _acyl_halide_carbons(mol)
    if not matches:
        raise UnsupportedStructure(
            "no acyl halide (-C(=O)X) group found; this module only "
            "handles acyl halides"
        )
    if len(matches) > 1:
        raise UnsupportedStructure(
            "more than one acyl halide group is out of scope for this "
            "module"
        )
    (acyl_carbon, (carbonyl_oxygen, acyl_halogen)) = next(iter(matches.items()))

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an acyl halide's own oxygen/halogen "
                "(P-65.1.5) and other halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num == 8:
            if atom.GetIdx() != carbonyl_oxygen:
                raise UnsupportedStructure(
                    "an oxygen other than the acyl halide's own carbonyl "
                    "oxygen (e.g. a coexisting ketone/aldehyde/carboxylic "
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
    carbon = mol.GetAtomWithIdx(acyl_carbon)
    carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "an acyl halide carbon with more than one carbon neighbor is "
            "not a valid acyl halide carbon"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return acyl_carbon, carbonyl_oxygen, acyl_halogen


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


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'oyl' endings
    (e.g. '2-enoyl'); the acyl halide's own locant is never cited (P-14.3.3,
    see module docstring)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    oyl_word = "oyl"

    words = [word for _, word in segments] + [oyl_word]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    if segments:
        locant_parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words[:-1])
        ]
        body = "-".join(locant_parts) + words[-1]
    else:
        body = words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, ene_locants, yne_locants, halide_word, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if (ene_locants or yne_locants) else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body + " " + halide_word


def _candidate_key(chain_length, ene_locants, yne_locants, halide_word, substituents):
    grouped = _group(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, halide_word, grouped)
    return (
        (
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


def _substituents_for_chain(graph, chain, halogens, excluded, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_phenyl_chain_acyl_halide(mol, ring_atoms):
    """Name an acyl halide whose -C(=O)X lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropanoyl chloride. The ring
    is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than
    the acyclic path above: no chain unsaturation and no specified
    stereocenter -- each is a separate follow-up (see
    tasks/phenyl-substituent-on-acyl-halide-chain.md's scope note)."""
    acyl_carbon, carbonyl_oxygen, acyl_halogen = _validate_and_collect_acyl_halides(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "acyl halide chain is not supported yet"
        )
    excluded = {carbonyl_oxygen, acyl_halogen}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent acyl "
            "halide chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain acyl halide is not supported yet"
        )
    chain, branches = longest_branched_chain(graph, acyl_carbon, ring_atoms, excluded, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an acyl halide directly attached to the benzene ring uses a "
            "separate naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    halide_word = _HALIDE_WORDS[mol.GetAtomWithIdx(acyl_halogen).GetAtomicNum()]
    halogens = halogen_substituents(mol)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = _group(substituents)
    return _name_from_substituents(chain_length, [], [], halide_word, grouped)


def name_acyl_halide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_acyl_halide(mol, ring_atoms)
    acyl_carbon, carbonyl_oxygen, acyl_halogen = _validate_and_collect_acyl_halides(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an acyl halide group on/in a ring is out of scope for this "
            "acyclic-only module"
        )

    excluded = {carbonyl_oxygen, acyl_halogen}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    graph = adjacency(mol)
    halide_word = _HALIDE_WORDS[mol.GetAtomWithIdx(acyl_halogen).GetAtomicNum()]
    halogens = halogen_substituents(mol)
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo = specified_stereocenters(mol)
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if acyl_carbon not in chain:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo_atoms and any(
            acyl_carbon in chain and (not bonds or _bond_locants(chain, bonds) is not None)
            for chain in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the acyl-halide-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter principal "
            "chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acyl_carbon:
                # The acyl halide carbon must sit at C1 (see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, halide_word, substituents)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of

    if stereo:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
