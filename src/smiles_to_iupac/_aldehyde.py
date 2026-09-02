"""Naming of aldehydes (the '-al' suffix, terminal -CHO) on acyclic saturated
or unsaturated carbon chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-33.3, Table 3.3 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  'al' is the preselected suffix for -CHO, cited in a combined chain/suffix
  name the same way 'ol'/'amine'/'one' are for -OH/-NH2/=O (`_alcohol.py`/
  `_amine.py`/`_ketone.py`), e.g. 'propanal', 'pentanedial'. Table 3.3 ranks
  'al' senior to 'one'/'ol'/'amine'; a carbonyl carbon shaped like a ketone
  (two carbon neighbors) or a carboxylic acid/ester/amide (a second oxygen on
  the same carbon) is rejected outright rather than silently named as if it
  were an aldehyde, and this module never attempts suffix-vs-suffix
  seniority competition against a coexisting senior group either (carboxylic
  acid/ester/amide/nitrile) since only C, halogen, and aldehyde/hydroxyl
  oxygen atoms are accepted at all. A coexisting hydroxyl (-OH), being junior
  to 'al', is *not* rejected: it is cited as the 'hydroxy' substituent prefix
  instead (P-41), e.g. 'OCC=O' -> '2-hydroxyethanal'.
- Unlike -OH/=O, a -CHO carbon is always a chain terminus (it has exactly one
  carbon neighbor, being otherwise saturated by =O and one H), so it is
  never a genuine locant choice: whichever end of the principal chain bears
  it is numbered C1 (P-44.4.1.8, suffix locants minimized first, same
  ordering as `_alcohol.py`/`_amine.py`/`_ketone.py`), and that locant is
  always '1' (or '1' and the last position, for a chain with -CHO at both
  ends) - so unlike those modules, this one never actually cites an 'al'
  locant in the name (P-14.3.3's general omission-when-unambiguous applies
  uniformly here, not just to a mononuclear/two-carbon parent): 'propanal',
  not 'propan-1-al'; 'pentanedial', not 'pentane-1,5-dial'. An 'ene'/'yne'
  locant elsewhere on the chain is still cited as usual, e.g. 'hex-4-enal'.
- P-31.0 / P-31.1.1.1-.2 (Chapter P-3): construction of the 'ene'/'yne'
  portion of a combined unsaturated-aldehyde name reuses the same mechanics
  as `_unsaturated.py`/`_ketone.py` (ending replaces 'ane' entirely,
  multiplying prefixes 'di'/'tri' for >=2 bonds of a kind, 'ene' before
  'yne', euphonic stem 'a' before a multiplied ending); the 'al'/'dial'
  suffix word is appended directly onto the last such segment (eliding a
  trailing vowel the same way 'ene'+'ol' does in `_alcohol.py`), since it
  never carries its own locant to separate it with a hyphen.
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  'al' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, e.g. '(2R)-2-chloropropanal',
  '(2R,3S)-2,3-dichlorobutanal', same pattern as `_carboxylic_acid.py`
  (CIP computation delegated entirely to RDKit).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any oxygen that isn't a doubly-bonded, isolated aldehyde carbonyl oxygen or
  a singly-bonded hydroxyl (ethers, and any oxygen bonded to more than one
  heavy atom).
- A carbonyl carbon with other than exactly one carbon neighbor: zero (a
  carbon-less carbonyl, e.g. formaldehyde) or two-or-more (a ketone) -
  neither is this module's territory.
- An aromatic carbonyl carbon, or any aromatic ring elsewhere in the
  molecule - a separate module's territory.
- Any other heteroatom (N, S, ...).
- A hydroxyl on a carbon that is also part of a C=C/C#C bond (an enol,
  tautomeric with a more senior carbonyl form) — same restriction as
  `_alcohol.py`'s own enol check.
- -CHO on a ring (the 'carbaldehyde' suffix, P-33.3.1.2, a substituent-style
  name rather than this module's parent-chain suffix) - deferred entirely;
  only an acyclic terminal -CHO is supported here.
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
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _validate_and_collect_aldehydes(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (aldehydes, hydroxyls): the set of aldehyde carbonyl-oxygen
    atom indices, and the set of any coexisting hydroxyl-oxygen atom indices.
    A hydroxyl is junior to 'al' in Table 3.3's suffix seniority order, so it
    is cited as the 'hydroxy' substituent prefix instead of competing for the
    suffix (P-41)."""
    aldehydes = set()
    hydroxyls = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than an aldehyde carbonyl oxygen (P-33.3) "
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
                    "ether) is out of scope; only an isolated aldehyde "
                    "carbonyl or hydroxyl is supported (Table 3.3, P-33.3)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("an aldehyde/hydroxyl oxygen must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() == 1.0:
                if atom.GetTotalNumHs() != 1:
                    raise UnsupportedStructure(
                        "an oxygen that isn't a carbonyl (=O) or hydroxyl "
                        "(-OH) is out of scope for this module"
                    )
                hydroxyls.add(atom.GetIdx())
                continue
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for this "
                    "module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                raise UnsupportedStructure(
                    "a carbonyl carbon with other than exactly one carbon "
                    "neighbor (a ketone, or a carbon-less carbonyl such as "
                    "formaldehyde) is not an aldehyde and is out of scope "
                    "for this module (Table 3.3)"
                )
            if mol.GetBondBetweenAtoms(carbon.GetIdx(), carbon_neighbors[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
                raise UnsupportedStructure(
                    "a carbonyl carbon that is itself doubly bonded to its "
                    "carbon neighbor (a cumulated double bond, e.g. a "
                    "ketene's C=C=O) is not an aldehyde and is out of scope "
                    "for this module"
                )
            aldehydes.add(atom.GetIdx())
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
    if not aldehydes:
        raise UnsupportedStructure("no aldehyde (-CHO) group found; this module only handles aldehydes")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return aldehydes, hydroxyls


def _suffix_body(ene_locants, yne_locants, al_count):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'al' ending
    (e.g. '4-enal'). Unlike `_ketone.py`'s '-one', 'al'/'dial' never carries
    its own locant (see module docstring), so it is glued directly onto the
    preceding word instead of getting a hyphenated locant segment of its
    own."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    al_word = multiplied_word(al_count, "al")
    words = [word for _, word in segments] + [al_word]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
            words[i] = words[i][:-1]

    locanted_parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words[:-1])
    ]
    body = "-".join(locanted_parts) + words[-1] if locanted_parts else words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _name_from_substituents(chain_length, al_count, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants, al_count)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.8 (suffix locants) ahead of
    P-44.4.1.10 (ene/yne locants) ahead of P-45.2 (substituent-prefix
    locants), most-preferred first. The 'al' locant set still drives
    orientation choice even though it is never printed (see module
    docstring)."""
    grouped = group_substituents(substituents)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    al_locant_set = lowest_locant_set(al_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, len(al_locants), ene_locants, yne_locants, grouped)
    return (
        (
            al_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _al_locants(position_of, aldehydes, graph):
    locants = []
    for o in aldehydes:
        (carbon,) = graph[o]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, aldehydes):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in aldehydes]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_aldehyde(mol, aldehydes, hydroxyls, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`'s identical treatment), and the
    winning candidate's own locants are used to format a
    "(<locant><R/S>,...)-" prefix onto the name, ascending locant order
    (P-91.3)."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _al_locants(position_of, aldehydes, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _al_locants({a: i + 1 for i, a in enumerate(c)}, aldehydes, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every aldehyde-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more -CHO groups (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            al_locants = _al_locants(position_of, aldehydes, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, aldehydes)
            key, name = _candidate_key(chain_length, al_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def name_aldehyde(mol) -> str:
    aldehydes, hydroxyls = _validate_and_collect_aldehydes(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    # Exclude each C=O carbonyl bond itself: `non_single_bonds` reports it as
    # order 2.0 same as a C=C, but it isn't a chain 'ene' bond (one endpoint
    # is the aldehyde oxygen, never part of any carbon chain).
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in aldehydes and b[1] not in aldehydes]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior carbonyl form and "
                "is out of scope for this module (P-31.1.4.2.4)"
            )

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "an aldehyde on a ring (the 'carbaldehyde' suffix, P-33.3.1.2) "
            "is out of scope for this module; only an acyclic terminal "
            "-CHO is supported"
        )
    return _name_acyclic_aldehyde(mol, aldehydes, hydroxyls, bonds, stereo)
