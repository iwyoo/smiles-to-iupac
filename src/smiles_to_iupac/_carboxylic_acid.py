"""Naming of carboxylic acids (the '-oic acid' suffix, -COOH) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.1.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.1.1; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'oic acid' is the
  preselected suffix for -COOH, ranked senior to every other suffix this
  project handles (ester, amide, nitrile, aldehyde, ketone, alcohol, amine).
  Since it outranks all of them, a coexisting carbonyl (aldehyde/ketone-
  shaped) is still rejected as an unresolved competition; but a coexisting
  standalone hydroxyl (-OH), being junior to '-oic acid', is cited as the
  'hydroxy' substituent prefix instead (P-41), e.g. 'OC(=O)CCO' ->
  '3-hydroxypropanoic acid'. Any other heteroatom, or any oxygen not part of
  a full -COOH pattern or a standalone hydroxyl, is simply rejected.
- P-65.1.1.2: a -COOH substituent on a ring is named with the separate
  'carboxylic acid' suffix (ring name + 'carboxylic acid'), a different
  construction from the chain-terminal 'oic acid' suffix this module
  builds. That case, and any ring anywhere in the molecule, is out of
  scope here (acyclic only).
- A -COOH carbon is always a chain terminus: after its carbonyl (=O) and
  hydroxyl (-OH) oxygens, it has room for at most one more substituent,
  which must be another chain carbon (or nothing, for formic acid,
  H-COOH). Unlike `_alcohol.py`/`_ketone.py`, this means the parent chain's
  numbering is never a locant-minimization choice for the -COOH group
  itself: whichever chain end carries a -COOH carbon simply becomes C1
  (P-65.1.1, mirroring the aldehyde '-al' suffix's own terminal-carbon
  fixing).
- P-14.3.3: the suffix locant for -COOH (and 'dioic acid' for two -COOH
  groups) is never cited, since it is always fully determined by the -COOH
  carbon(s) being chain terminus/termini with no other possible position,
  e.g. 'propanoic acid' (not 'propan-1-oic acid'), 'hexanedioic acid' (not
  'hexane-1,6-dioic acid'). Confirmed against PubChem. Only 'ene'/'yne'
  locants (P-31.1.1.1-.2, same mechanics as `_alcohol.py`) are ever cited in
  this module's suffix, e.g. 'but-2-enoic acid', 'but-2-enedioic acid'
  (fumaric/maleic acid's PIN).
- P-44.1.1 (Chapter P-4): the senior parent structure captures the maximum
  number of principal characteristic groups (here, -COOH). As in
  `_alcohol.py`, this module only accepts a candidate chain that is
  simultaneously one of the longest carbon chains (P-44.3.2) and carries
  every -COOH-bearing carbon in the molecule; anything requiring a shorter
  principal chain to capture more -COOH groups is out of scope.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with
  the -oic acid suffix, reusing `halogen_substituents`/
  `format_substituent_prefixes` unchanged.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule (P-65.1.1.2's territory; acyclic only,
  per this module's scope).
- Any oxygen that isn't part of a full -COOH pattern on some carbon, or a
  standalone hydroxyl on a carbon with no carbonyl (a lone carbonyl with no
  matching hydroxyl indicates an unresolved aldehyde/ketone competition; an
  ether or any other oxygen shape).
- A -COOH carbon with more than one carbon neighbor (impossible for a
  genuine carboxyl carbon, since its remaining two bonds are already the
  carbonyl and hydroxyl oxygens; guarded defensively).
- A standalone hydroxyl on a carbon that is also part of a C=C/C#C bond (an
  enol, tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- Any other heteroatom (N, S, ...).
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
_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def has_carboxylic_acid_shape(mol) -> bool:
    """True if some carbon carries both a doubly-bonded, monovalent oxygen
    and a singly-bonded, monovalent, one-H oxygen (a -COOH pattern),
    regardless of whether the rest of the molecule is in scope. Used by
    `core.py` to route ahead of the ketone/alcohol dispatch, since a -COOH
    carbon would otherwise look aldehyde- or alcohol-shaped to those
    modules."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        has_carbonyl = any(
            o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for o in oxygens
        )
        has_hydroxyl = any(
            o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
            for o in oxygens
        )
        if has_carbonyl and has_hydroxyl:
            return True
    return False


def _validate_and_collect_carboxyls(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls): the set
    of -COOH carbon atom indices, the set of all their carbonyl/hydroxyl
    oxygen atom indices (excluded from substituent detection), and the set of
    any coexisting standalone alcohol hydroxyl-oxygen atom indices. A
    standalone hydroxyl is junior to 'oic acid' in Table 3.3's suffix
    seniority order, so it is cited as the 'hydroxy' substituent prefix
    instead of competing for the suffix (P-41)."""
    has_carbon = False
    carbonyls_by_carbon = {}
    hydroxyls_by_carbon = {}
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than carboxylic-acid oxygens (P-65.1.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
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
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether or ester) is out of scope; only oxygens forming "
                    "an isolated -COOH group are supported (P-65.1.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a -COOH oxygen must be attached to a carbon atom")
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 2.0:
                carbonyls_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
            elif bond_order == 1.0 and atom.GetTotalNumHs() == 1:
                hydroxyls_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
            else:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "half of a -COOH group is out of scope for this module"
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

    carboxyl_carbons = set()
    carboxyl_oxygens = set()
    extra_hydroxyls = set()
    for carbon_idx in set(carbonyls_by_carbon) | set(hydroxyls_by_carbon):
        carbonyls = carbonyls_by_carbon.get(carbon_idx, [])
        hydroxyls = hydroxyls_by_carbon.get(carbon_idx, [])
        if not carbonyls and len(hydroxyls) == 1:
            # A hydroxyl-only carbon is a standalone alcohol (-OH), junior to
            # -COOH in Table 3.3 (see docstring) - cited as a 'hydroxy'
            # prefix rather than rejected outright.
            extra_hydroxyls.add(hydroxyls[0])
            continue
        if len(carbonyls) != 1 or len(hydroxyls) != 1:
            raise UnsupportedStructure(
                "an oxygen pattern that isn't exactly one carbonyl and one "
                "hydroxyl oxygen on the same carbon, or a standalone "
                "hydroxyl, is a more/less senior characteristic group than a "
                "plain carboxylic acid (Table 3.3), which this module does "
                "not attempt to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "a -COOH carbon with more than one carbon neighbor is not a "
                "valid carboxyl carbon"
            )
        carboxyl_carbons.add(carbon_idx)
        carboxyl_oxygens.add(carbonyls[0])
        carboxyl_oxygens.add(hydroxyls[0])

    if not carboxyl_carbons:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls


def _multiplied_word(count, base):
    if count == 0:
        return ""
    if count == 1:
        return base
    return numerical_term(count) + base


def _suffix_body(ene_locants, yne_locants, acid_count):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'oic' endings
    (e.g. '2-enoic', '2-enedioic'); the acid group's own locant is never
    cited (P-14.3.3, see module docstring). Also returns whether the stem's
    trailing 'e' should be elided at the stem/first-segment boundary (only
    relevant when there is no 'ene'/'yne', see `_name_from_substituents`)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), _multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), _multiplied_word(len(yne_locants), "yne")))
    acid_word = _multiplied_word(acid_count, "oic")

    words = [word for _, word in segments] + [acid_word]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
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


def _name_from_substituents(chain_length, acid_count, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, acid_count)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if (ene_locants or yne_locants) else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body + " acid"


def _candidate_key(chain_length, acid_count, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The -COOH group's
    own locant isn't part of this key: candidates are pre-filtered so a
    -COOH carbon always sits at C1 (see `_name_acyclic_carboxylic_acid`)."""
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
    name = _name_from_substituents(chain_length, acid_count, ene_locants, yne_locants, grouped)
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


def _substituents_for_chain(graph, chain, halogens, carboxyl_oxygens):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in carboxyl_oxygens]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_carboxylic_acid(mol, carboxyl_carbons, carboxyl_oxygens, hydroxyls, bonds):
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = _longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    acid_count = len(carboxyl_carbons)

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if not carboxyl_carbons <= chain_set:
            continue
        if bonds and _bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every -COOH-bearing carbon (and/or multiple bond) lies on "
            "a single longest carbon chain; a shorter principal chain "
            "capturing more -COOH groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] not in carboxyl_carbons:
                # A -COOH carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = _bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, carboxyl_oxygens)
            key, name = _candidate_key(chain_length, acid_count, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_carboxylic_acid(mol) -> str:
    carboxyl_carbons, carboxyl_oxygens, hydroxyls = _validate_and_collect_carboxyls(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH group on/in a ring uses the separate 'carboxylic acid' "
            "suffix construction (P-65.1.1.2), out of scope for this "
            "acyclic-only module"
        )

    all_non_single = [
        b for b in non_single_bonds(mol) if b[0] not in carboxyl_oxygens and b[1] not in carboxyl_oxygens
    ]
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    graph = adjacency(mol)
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior carbonyl form and "
                "is out of scope for this module (P-31.1.4.2.4)"
            )

    return _name_acyclic_carboxylic_acid(mol, carboxyl_carbons, carboxyl_oxygens, hydroxyls, bonds)
