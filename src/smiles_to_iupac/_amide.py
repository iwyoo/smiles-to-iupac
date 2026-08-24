"""Naming of primary amides (the '-amide' suffix, terminal -C(=O)NH2) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-66.1, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-66.1; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'amide' is the preselected
  suffix for -CONH2, ranked junior to '-oic acid'/ester but senior to
  nitrile/aldehyde/ketone/alcohol/amine. A coexisting standalone hydroxyl
  (-OH), being junior to 'amide', is cited as the 'hydroxy' substituent
  prefix instead (P-41), mirroring `_aldehyde.py`/`_carboxylic_acid.py`. Any
  other coexisting carbonyl (aldehyde/ketone/-COOH-shaped) is rejected as an
  unresolved suffix-vs-suffix seniority competition, same as those modules.
- Like a -COOH carbon (`_carboxylic_acid.py`), an amide carbon is always a
  chain terminus: after its carbonyl (=O) and amide nitrogen, it has room
  for at most one more substituent, which must be another chain carbon (or
  nothing, for methanamide/formamide, HCONH2). So the parent chain's
  numbering is never a locant-minimization choice for the amide group
  itself: whichever chain end carries the amide carbon simply becomes C1.
- P-14.3.3: the suffix locant for -CONH2 is never cited, since it is always
  fully determined by the amide carbon being a chain terminus with no other
  possible position, e.g. 'ethanamide' (not 'ethan-1-amide').
- P-31.0/P-31.1.1.1-.2: construction of the 'ene'/'yne' portion of a combined
  unsaturated-amide name reuses the same mechanics as `_aldehyde.py`/
  `_carboxylic_acid.py` (ending replaces 'ane' entirely, multiplying
  prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before 'yne', euphonic
  stem 'a' before a multiplied ending); 'amide' is appended directly onto
  the last such segment (eliding a trailing vowel the same way 'ene'+'al'
  does), since it never carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'amide' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.

Explicitly out of scope (raise `UnsupportedStructure`), per the task's
first-pass scope:
- An amide nitrogen with any substituent other than its two hydrogens
  (N-substituted amide, e.g. N-methylacetamide) - deferred entirely; only a
  primary amide (-CONH2) is supported here.
- An amide on/in a ring (a lactam) - a separate module's territory.
- More than one amide group in the same molecule (a diamide) - deferred
  entirely, along with any other multiple-principal-characteristic-group
  combination.
- Any oxygen that isn't the amide's own carbonyl oxygen or a standalone
  hydroxyl (ethers, esters, -COOH, and any oxygen bonded to more than one
  heavy atom).
- A carbonyl carbon with other than exactly one nitrogen neighbor, or more
  than one carbon neighbor (a ketone-shaped carbon is not an amide carbon).
- An aromatic carbonyl carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory.
- Any other heteroatom (S, ...).
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def has_amide_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded, monovalent carbonyl
    oxygen and a singly-bonded, monovalent, two-H nitrogen (a primary
    -CONH2 pattern), regardless of whether the rest of the molecule is in
    scope. Used by `core.py` to route ahead of the aldehyde/ketone dispatch,
    since an amide carbon would otherwise look aldehyde-shaped to those
    modules (both have exactly one carbon neighbor besides the carbonyl)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        has_carbonyl = any(
            n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
        has_primary_amide_n = any(
            n.GetAtomicNum() == 7
            and n.GetDegree() == 1
            and n.GetTotalNumHs() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in atom.GetNeighbors()
        )
        if has_carbonyl and has_primary_amide_n:
            return True
    return False


def _validate_and_collect_amide(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (amide_carbon, amide_oxygen, amide_nitrogen, hydroxyls): the
    single -CONH2 carbon/oxygen/nitrogen atom indices, and the set of any
    coexisting standalone hydroxyl-oxygen atom indices."""
    has_carbon = False
    amide_carbons = set()
    amide_oxygen_by_carbon = {}
    amide_nitrogen_by_carbon = {}
    hydroxyls = set()
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a primary amide's oxygen/nitrogen "
                "(P-66.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
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
                    "ether or ester) is out of scope; only an isolated "
                    "amide carbonyl or a standalone hydroxyl is supported "
                    "(P-66.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an amide/hydroxyl oxygen must be attached to a carbon atom")
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 1.0 and atom.GetTotalNumHs() == 1:
                hydroxyls.add(atom.GetIdx())
                continue
            if bond_order != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for "
                    "this module"
                )
            amide_oxygen_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
        elif atomic_num == 7:
            if atom.GetDegree() != 1 or atom.GetTotalNumHs() != 2:
                raise UnsupportedStructure(
                    "an amide nitrogen with any substituent other than its "
                    "two hydrogens (N-substituted amide) is out of scope "
                    "for this module"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure("an amide nitrogen must be singly bonded to its carbonyl carbon")
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an amide nitrogen must be attached to a carbon atom")
            amide_nitrogen_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
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

    for carbon_idx in set(amide_oxygen_by_carbon) | set(amide_nitrogen_by_carbon):
        oxygens = amide_oxygen_by_carbon.get(carbon_idx, [])
        nitrogens = amide_nitrogen_by_carbon.get(carbon_idx, [])
        if len(oxygens) != 1 or len(nitrogens) != 1:
            raise UnsupportedStructure(
                "an oxygen/nitrogen pattern that isn't exactly one carbonyl "
                "oxygen and one primary-amide nitrogen on the same carbon "
                "is a more/less senior characteristic group than a plain "
                "primary amide (Table 3.3), which this module does not "
                "attempt to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "an amide carbon with more than one carbon neighbor is not "
                "a valid terminal amide carbon (a ketone-shaped carbon is "
                "out of scope for this module)"
            )
        amide_carbons.add(carbon_idx)

    if not amide_carbons:
        raise UnsupportedStructure(
            "no primary amide (-CONH2) group found; this module only "
            "handles primary amides"
        )
    if len(amide_carbons) > 1:
        raise UnsupportedStructure(
            "more than one amide group (a diamide) is out of scope for "
            "this module"
        )
    (amide_carbon,) = amide_carbons
    (amide_oxygen,) = amide_oxygen_by_carbon[amide_carbon]
    (amide_nitrogen,) = amide_nitrogen_by_carbon[amide_carbon]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return amide_carbon, amide_oxygen, amide_nitrogen, hydroxyls


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'amide' ending
    (e.g. '2-enamide'); the amide group's own locant is never cited
    (P-14.3.3, see module docstring)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["amide"]
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


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
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
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The amide group's
    own locant isn't part of this key: candidates are pre-filtered so the
    amide carbon always sits at C1 (see `_name_acyclic_amide`)."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_amide(mol, amide_carbon, excluded, hydroxyls, bonds):
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if amide_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the amide-bearing carbon (and/or a multiple bond) does not lie "
            "on a single longest carbon chain; a shorter principal chain is "
            "not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != amide_carbon:
                # The amide carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_amide(mol) -> str:
    amide_carbon, amide_oxygen, amide_nitrogen, hydroxyls = _validate_and_collect_amide(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amide on/in a ring (a lactam) is out of scope for this "
            "acyclic-only module"
        )

    excluded = {amide_oxygen, amide_nitrogen}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
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

    return _name_acyclic_amide(mol, amide_carbon, excluded, hydroxyls, bonds)
