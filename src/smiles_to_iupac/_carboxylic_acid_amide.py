"""Naming of an amic acid (a molecule combining exactly one free carboxylic
acid (-COOH) with one or more coexisting primary amide (-CO-NH2) groups on
the same acyclic saturated carbon chain), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-66.1.1.3.3 ("Amic acids"): 'oic acid' (`_carboxylic_acid.py`) outranks
  'amide' (`_amide.py`), so a coexisting primary amide is demoted to
  *two* substituent prefixes at once -- 'amino' for its nitrogen, 'oxo'
  for its own carbon -- instead of its own '-amide' suffix, e.g.
  'H2N-CO-[CH2]8-COOH' -> '10-amino-10-oxodecanoic acid' (the Blue Book's
  own worked example) and 'H2N-CO-CH2-COOH' (malonamic acid) ->
  '3-amino-3-oxopropanoic acid'. This mirrors
  `_aldehyde_carboxylic_acid.py`'s single-'oxo' carbonyl demotion, just
  with a second, simultaneous 'amino' prefix at the same chain position
  for the nitrogen -- the amide carbon is a chain atom (like the aldehyde
  case), not a plain substituent atom (like `_carboxylic_acid_amine.py`'s
  amine nitrogen).
- Otherwise mirrors `_carboxylic_acid.py`/`_aldehyde_carboxylic_acid.py`
  exactly: the -COOH carbon is always a chain terminus and always
  becomes C1 with its own locant never cited (P-14.3.3); an amide
  carbon's locant, by contrast, is always cited.

Scope, deliberately narrow: a single carboxylic acid plus one or more
*plain, unsubstituted* primary amides (-CO-NH2, no N-substitution), all on
one acyclic saturated chain, with halogen substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring, more than one carboxylic acid, an
N-substituted amide, a coexisting hydroxyl/ether/other heteroatom, and any
acid/amide not captured by a single longest chain.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    name_from_substituents,
    non_single_bonds,
    substituent_locant_set_and_citation,
    validate_allowed_atoms,
)
from ._substituents import format_substituent_prefixes, substituents_for_chain


def _find_carboxylic_acid_carbon(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
        ]
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbonyls) == 1 and len(hydroxyls) == 1 and len(carbon_neighbors) <= 1:
            return atom, carbonyls[0], hydroxyls[0]
    return None


def _extra_amides(mol, excluded_oxygens):
    """Every plain, unsubstituted primary amide (-CO-NH2) other than the
    acid's own carbon, as (carbonyl_oxygen_idxs, nitrogen_idxs) -- mirrors
    `_aldehyde_carboxylic_acid.py`'s `_extra_aldehydes`, but an amide
    carbon's nitrogen substituent is also collected for its own 'amino'
    demotion, alongside the carbon's 'oxo' one."""
    oxygens = set()
    nitrogens = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        carbonyl_oxygens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetIdx() not in excluded_oxygens
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        amide_nitrogens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 7
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(carbonyl_oxygens) != 1 or len(amide_nitrogens) != 1:
            continue
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        oxygens.add(carbonyl_oxygens[0].GetIdx())
        nitrogens.add(amide_nitrogens[0].GetIdx())
    return oxygens, nitrogens


def has_carboxylic_acid_amide_shape(mol) -> bool:
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return False
    _, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    oxygens, _ = _extra_amides(mol, excluded)
    return bool(oxygens)


def _validate(mol, excluded_oxygens, amide_nitrogens):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than the acid's own oxygens, plain primary "
        "amide oxygens/nitrogens (P-66.1.1.3.3), and halogen substituents "
        "(P-35.2.1) are not supported yet",
        [
            (
                8,
                excluded_oxygens,
                "an oxygen that isn't part of the single carboxylic acid "
                "group or an amide-shaped carbonyl is out of scope for "
                "this module (e.g. a coexisting hydroxyl or ether)",
            ),
            (
                7,
                amide_nitrogens,
                "a nitrogen that isn't a plain, unsubstituted primary "
                "amide (-CO-NH2) is out of scope for this module "
                "(N-substituted amides are not supported)",
            ),
        ],
    )


def _name_from_substituents(chain_length, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "oic acid")


def _candidate_key(chain_length, grouped):
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def name_carboxylic_acid_amide(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH group on/in a ring uses a different naming "
            "construction, out of scope for this acyclic-only module"
        )
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    amide_oxygens, amide_nitrogens = _extra_amides(mol, excluded_acid_oxygens)
    if not amide_oxygens:
        raise UnsupportedStructure(
            "no coexisting primary amide found; this module only handles "
            "a carboxylic acid combined with at least one amide (see "
            "_carboxylic_acid.py for a plain carboxylic acid)"
        )
    excluded = excluded_acid_oxygens | amide_oxygens
    _validate(mol, excluded, amide_nitrogens)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an amide-on-acid "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {
        **halogen_substituents(mol),
        **{o: "oxo" for o in amide_oxygens},
        **{n: "amino" for n in amide_nitrogens},
    }
    acid_carbon_idx = acid_carbon.GetIdx()
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    demoted_carbons = {graph[o][0] for o in amide_oxygens}
    eligible = []
    for chain in chains:
        if acid_carbon_idx not in chain:
            continue
        chain_set = set(chain)
        if not demoted_carbons <= chain_set:
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every acid/amide-bearing carbon lies on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acid_carbon_idx:
                # The -COOH carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            substituents = substituents_for_chain(graph, candidate, names, excluded_acid_oxygens, mol=mol)
            grouped = group_substituents(substituents)
            key, name = _candidate_key(chain_length, grouped)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
