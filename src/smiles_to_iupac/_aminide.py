"""Naming of aminide anions (the '-aminide' suffix, -NH(-)), restricted to
a single such group on an acyclic saturated or unsaturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-72.2.2.2.3 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html): an
  anion formed by removing a hydron from a primary amine's nitrogen is
  named by adding 'ide' to the '-amine' suffix -- e.g. CH3-NH(-) ->
  'methanaminide (PIN)' (confirmed worked example). This project's
  `_amine.py` already names CH3-NH2 as 'methanamine', so this module
  mirrors that module's acyclic chain/locant logic exactly, with the
  suffix word changed from 'amine' to 'aminide'.
- Structure-verified via PubChem: `C[NH-]` (CID 21952893) and `CC[NH-]`
  (CID 187872) -- PubChem's own generated names use a different naming
  system ('methylazanide'/'ethylazanide', an azanide-parent-hydride style)
  rather than the Blue Book's '-aminide' suffix convention, so only the
  structures (not the names) are cross-checked, same as this project's
  other anion/isotope modules whose PubChem-generated names diverge from
  the Blue Book's own PIN.
- P-14.3.4.2(a)/(b), P-35.2.1: identical locant-omission and halogen-
  coexistence rules as `_amine.py`'s own '-amine' suffix.

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring anywhere in the molecule (acyclic-only, mirrors the other anion
  modules; `_amine.py` itself supports a simple ring, but this first
  aminide pass doesn't yet).
- Aromatic -NH(-) (an anilinide-type anion) -- a separate follow-up.
- Iminide anions (a deprotonated imine, P-72.2.2.2.3's other suffix) --
  a different characteristic group entirely.
- More than one -NH(-) group, or a doubly-charged nitrogen ('aminediide').
- Secondary/tertiary amine anions -- this project's cross-cutting
  secondary/tertiary amine blocker applies here too (see `_amine.py`).
- Any other heteroatom (O, S, ...), charge, or isotopic modification.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ENE_BOND_ORDER,
    YNE_BOND_ORDER,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _find_aminide_nitrogens(mol):
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if atom.GetFormalCharge() != -1 or atom.GetDegree() != 1:
            continue
        if atom.GetTotalNumHs() != 1:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() != 6:
            continue
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
        if bond.GetBondTypeAsDouble() != 1.0:
            continue
        matches.append(atom)
    return matches


def has_aminide_shape(mol) -> bool:
    return bool(_find_aminide_nitrogens(mol))


def _find_aminide_group(mol):
    matches = _find_aminide_nitrogens(mol)
    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one aminide (-NH(-)) group is required; zero or "
            "multiple such groups are not supported yet (P-72.2.2.2.3)"
        )
    return matches[0]


def _suffix_body(ene_locants, yne_locants, n_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append(([n_locant], "aminide"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    elide_stem = words[0][0] in "aeiouy"
    return "-".join(parts), elide_stem


def _name_from_substituents(chain_length, n_locant, ene_locants, yne_locants, grouped):
    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locants are always '1' and
        # never cited, however many substituents there are.
        stem = alkane_name(1)[:-1]  # 'aminide' starts with a vowel
        return format_substituent_prefixes(grouped, omit_locants=True) + stem + "aminide"

    has_unsaturation = bool(ene_locants or yne_locants)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and not has_unsaturation and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly one
        # substituent (here, the sole -NH(-)) in total has only one
        # possible structure, so the locant is omittable.
        return alkane_name(2)[:-1] + "aminide"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, n_locant)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, n_locant, ene_locants, yne_locants, substituents):
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
    name = _name_from_substituents(chain_length, n_locant, ene_locants, yne_locants, grouped)
    return (
        (
            (n_locant,),
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded_atoms):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _reject_enamine_carbon(graph, nitrogen_idx, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    (carbon,) = graph[nitrogen_idx]
    if carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "an aminide on a carbon that is also part of a C=C/C#C bond "
            "(an enamine-type structure) is out of scope for this module"
        )


def _name_acyclic_aminide(mol, nitrogen_idx, excluded_atoms, bonds):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    (nitrogen_carbon,) = [n.GetIdx() for n in mol.GetAtomWithIdx(nitrogen_idx).GetNeighbors()]

    eligible = []
    for chain in chains:
        if nitrogen_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "the aminide carbon (and/or multiple bonds) does not lie on a "
            "single longest carbon chain"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            n_locant = candidate.index(nitrogen_carbon) + 1
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded_atoms)
            key, name = _candidate_key(chain_length, n_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_aminide(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring is out of scope for this acyclic-only module"
        )

    nitrogen = _find_aminide_group(mol)
    excluded_atoms = {nitrogen.GetIdx()}

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atom.GetIdx() in excluded_atoms:
            continue
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the aminide's own nitrogen "
                "(P-72.2.2.2.3) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "single aminide anion nitrogen is not supported yet"
            )
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 7:
            raise UnsupportedStructure(
                "more than one nitrogen is not supported yet (P-72.2.2.2.3)"
            )
        elif atom.GetDegree() != 1:
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

    all_non_single = [
        b for b in non_single_bonds(mol) if b[0] not in excluded_atoms and b[1] not in excluded_atoms
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enamine_carbon(adjacency(mol), nitrogen.GetIdx(), bonds)

    return _name_acyclic_aminide(mol, nitrogen.GetIdx(), excluded_atoms, bonds)
