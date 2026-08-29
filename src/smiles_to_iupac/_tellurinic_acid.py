"""Naming of tellurinic acids (the '-tellurinic acid' suffix, -Te(=O)OH) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-65.3.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf), the
  tellurium analogue of the sulfinic acid entry in the same table (the same
  chalcogen-pair pattern as `_seleninic_acid.py`'s own relationship to
  `_sulfinic_acid.py`): 'tellurinic acid' is the preselected suffix for
  -Te(=O)-OH, one oxidation state below telluronic acid (`_telluronic_acid.py`)
  -- the tellurium carries one carbon, one double-bonded oxygen, and one
  hydroxyl oxygen (degree 3, not 4), structurally identical to
  `_sulfinic_acid.py`'s -S(=O)-OH with tellurium instead of sulfur.
- Like 'sulfinic acid'/'seleninic acid', 'tellurinic acid' is cited as the
  parent hydride name followed directly by the two-word suffix with no
  elision -- 'methane' + 'tellurinic acid' -> 'methanetellurinic acid'
  (PubChem structure match -- the only tellurinic acid PubChem has
  registered at all; ethane/propane candidates come back as CID 0, the
  same sparse-data gap `_telluronic_acid.py` already documented). The
  chain locant-omission rules (P-14.3.4.2(a)/(b)) themselves are NOT
  independently re-verified for tellurium here -- inherited, unchanged,
  from the identical mechanism already independently confirmed by
  `_thiol.py`/`_alcohol.py`/`_sulfonic_acid.py`/`_selenonic_acid.py`/
  `_sulfinic_acid.py`/`_seleninic_acid.py`/`_telluronic_acid.py`, the same
  reduced-evidence bar `_selenol.py` first used for its own
  halogen-coexistence claim.
- P-44.4.1.8 / P-45.2: the -Te(=O)OH locant is minimized before ene/yne
  locants, which are minimized before substituent-prefix locants -- same
  ordering as `_sulfinic_acid.py`/`_telluronic_acid.py`.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -Te(=O)OH suffix.

Scope, deliberately narrow, mirroring `_telluronic_acid.py`'s own
chain-only scope (no monocyclic tellurinic acid has been found registered
on PubChem to verify that shape here): a single -Te(=O)OH on an acyclic
chain, with no other heteroatom anywhere in the molecule except the
tellurinic acid group's own two oxygens. Explicitly out of scope (raise
`UnsupportedStructure`): any ring, two or more -Te(=O)OH groups, and a
tellurinic acid on a carbon that is also part of a C=C/C#C bond.
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
_TELLURIUM = 52


def _tellurinic_tellurium_atoms(mol):
    """Tellurium atoms shaped like a tellurinic acid group: bonded to exactly
    one carbon, one double-bonded (terminal) oxygen, and one single-bonded
    hydroxyl oxygen (terminal, one H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _TELLURIUM or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        hydroxyl_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def has_tellurinic_acid_shape(mol) -> bool:
    return bool(_tellurinic_tellurium_atoms(mol))


def _validate_and_collect_tellurinic_acids(mol):
    tellurium_atoms = _tellurinic_tellurium_atoms(mol)
    if not tellurium_atoms:
        raise UnsupportedStructure(
            "no tellurinic acid (-Te(=O)OH) group found; this module only "
            "handles tellurinic acids"
        )
    if len(tellurium_atoms) > 1:
        raise UnsupportedStructure(
            "more than one tellurinic acid group is out of scope for this "
            "module"
        )
    tellurinic_atom_idxs = set()
    for te in tellurium_atoms:
        tellurinic_atom_idxs.add(te.GetIdx())
        tellurinic_atom_idxs.update(n.GetIdx() for n in te.GetNeighbors() if n.GetAtomicNum() == 8)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in tellurinic_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a tellurinic acid group (P-65.3.1) "
                "and halogen substituents (P-35.2.1) are not supported "
                "yet -- in particular a coexisting carboxylic/sulfonic/"
                "sulfinic/selenonic acid or other characteristic group "
                "needs acid-vs-acid Table 3.3 seniority handling not yet "
                "implemented here"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported yet (this module only handles "
            "acyclic chains)"
        )

    (tellurium,) = tellurium_atoms
    (carbon,) = (n for n in tellurium.GetNeighbors() if n.GetAtomicNum() == 6)
    return tellurium.GetIdx(), carbon.GetIdx()


def _reject_enetellurinic_carbon(graph, teoh_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if teoh_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a tellurinic acid on a carbon that is also part of a C=C/C#C "
            "bond is out of scope for this module"
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


def _suffix_body(ene_locants, yne_locants, teoh_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    segments.append(([teoh_locant], "tellurinic acid"))

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


def _name_from_substituents(chain_length, teoh_locant, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "tellurinic acid"

    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -Te(=O)OH) in total omits the locant, e.g.
        # 'ethanetellurinic acid'.
        return alkane_name(2) + "tellurinic acid"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, teoh_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, teoh_locant, ene_locants, yne_locants, substituents):
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
    name = _name_from_substituents(chain_length, teoh_locant, ene_locants, yne_locants, grouped)
    return (
        (
            teoh_locant,
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def name_tellurinic_acid(mol) -> str:
    tellurium_idx, teoh_carbon = _validate_and_collect_tellurinic_acids(mol)
    graph = adjacency(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and tellurium_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 1:
        # The one Se=O double bond is always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enetellurinic_carbon(graph, teoh_carbon, bonds)

    halogens = halogen_substituents(mol)
    excluded = {tellurium_idx}
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if teoh_carbon not in chain:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the tellurinic-acid-bearing carbon (and/or a multiple bond) "
            "does not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            teoh_locant = position_of[teoh_carbon]
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, teoh_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
