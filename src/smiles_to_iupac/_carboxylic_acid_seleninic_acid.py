"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more seleninic acid (-Se(=O)OH) groups on the same acyclic
saturated carbon chain, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-41/P-43 (`_seniority.py`): carboxylic acid (Table 4.4 class 1)
  outranks seleninic acid (class 13), so a coexisting seleninic acid is
  demoted to the 'selenino' substituent prefix (P-65.3.1) instead of its
  own '-seleninic acid' suffix -- e.g. 'OC(=O)C[Se](=O)O' ->
  '2-seleninoacetic acid' (PubChem CID 428823 confirms the systematic
  '2-seleninoacetic acid' form; this module always uses the systematic
  '...anoic acid' stem rather than a retained name, mirroring
  `_carboxylic_acid.py`'s and `_carboxylic_acid_sulfinic_acid.py`'s own
  established convention once any substituent is present). Mirrors
  `_carboxylic_acid_sulfinic_acid.py` (PR #315), with the demoted group
  swapped from sulfinic acid/'sulfino' (sulfur) to seleninic
  acid/'selenino' (selenium) -- structurally identical shape, one row
  down the chalcogen table (P-65.3.1's own sulfur/selenium/tellurium
  triad, same pattern as `_seleninic_acid.py`'s relationship to
  `_sulfinic_acid.py`).
- Otherwise mirrors `_carboxylic_acid.py`'s acyclic-chain path: P-14.3.3's
  locant-omission rule for the -COOH suffix itself (always C1, never
  cited), P-14.3.4.2(a)'s locant-omission rule for a mononuclear parent's
  own substituent locant, and 'selenino' is injected into the same
  {atom_idx -> name} map `_carboxylic_acid.py` already uses for halogens/
  a demoted hydroxyl, so it is formatted, alphabetized, and multiplied by
  the same general substituent-prefix machinery.

Scope, deliberately narrow (mirrors `_carboxylic_acid_sulfinic_acid.py`):
a single carboxylic acid plus one or more seleninic acids, all on one
acyclic saturated chain, with halogen substituents allowed. Explicitly
out of scope (raise `UnsupportedStructure`): any chain unsaturation
(ene/yne), any ring, more than one carboxylic acid, a coexisting
standalone hydroxyl/ether/other heteroatom, a specified stereocenter, and
any carboxylic acid/seleninic acid not captured by a single longest
chain. The tellurium analogue (tellurinic acid) has no registered PubChem
example to verify this shape against, so it is left out of this module's
scope.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
)
from ._numerals import alkane_name
from ._seniority import senior_class
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_SELENIUM = 34
_ALLOWED_ATOMIC_NUMS = {6, 8, _SELENIUM, *HALOGEN_PREFIXES}

assert senior_class("carboxylic_acid", "seleninic_acid") == "carboxylic_acid"


def _seleninic_selenium_atoms(mol):
    """Selenium atoms shaped like a seleninic acid group: bonded to
    exactly one carbon, one double-bonded (terminal) oxygen, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_seleninic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SELENIUM or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 2:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def _carboxyl_carbons(mol):
    """Carbon atoms shaped like a -COOH group: one double-bonded (terminal)
    oxygen and one single-bonded, one-H hydroxyl oxygen, both otherwise
    unaccounted-for. Mirrors the core shape check `_carboxylic_acid.py`
    uses, narrowed to plain -COOH only (no standalone-hydroxyl carve-out,
    kept out of scope here)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
        ]
        if len(carbonyls) != 1 or len(hydroxyls) != 1:
            continue
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            continue
        matches.append(atom)
    return matches


def has_carboxylic_acid_seleninic_acid_shape(mol) -> bool:
    return bool(_carboxyl_carbons(mol)) and bool(_seleninic_selenium_atoms(mol))


def _validate_and_collect(mol):
    carboxyl_carbons = _carboxyl_carbons(mol)
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "exactly one carboxylic acid is required; this module only "
            "handles a single carboxylic acid combined with one or more "
            "seleninic acids"
        )
    (carboxyl_carbon,) = carboxyl_carbons
    carboxyl_oxygens = {n.GetIdx() for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 8}

    seleninics = _seleninic_selenium_atoms(mol)
    if not seleninics:
        raise UnsupportedStructure(
            "no seleninic acid found; this module only handles a "
            "carboxylic acid combined with at least one seleninic acid "
            "(see _carboxylic_acid.py for a plain carboxylic acid)"
        )
    seleninic_idxs = {s.GetIdx() for s in seleninics}
    seleninic_oxygens = {
        o.GetIdx() for s in seleninics for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }

    accounted_oxygen_idxs = carboxyl_oxygens | seleninic_oxygens

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the carboxylic/seleninic acid "
                "groups' own oxygens (P-65.1.1/P-65.3.1) and halogen "
                "substituents (P-35.2.1) are not supported yet"
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
            if atom.GetIdx() not in accounted_oxygen_idxs:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the carboxylic/seleninic "
                    "acid groups is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == _SELENIUM:
            if atom.GetIdx() not in seleninic_idxs:
                raise UnsupportedStructure(
                    "a selenium atom not shaped like a seleninic acid "
                    "group is out of scope for this module"
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
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return carboxyl_carbon.GetIdx(), carboxyl_oxygens, seleninic_idxs


def _name_from_substituents(chain_length, grouped):
    # P-14.3.4.2(a): a mononuclear parent's substituent locant is always
    # '1' and never cited.
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    stem = alkane_name(chain_length)
    # P-14.3.3: the -COOH suffix locant is never cited, since a -COOH
    # carbon is always the chain-terminal C1.
    return f"{prefix}{stem[:-1]}oic acid"


def _candidate_key(chain_length, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def _substituents_for_chain(graph, chain, names, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names) for root in branch_roots]
    return substituents


def name_carboxylic_acid_seleninic_acid(mol) -> str:
    carboxyl_carbon, carboxyl_oxygens, seleninic_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a carboxylic acid/seleninic acid combination on a ring is "
            "out of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    excluded_from_unsaturation_check = {carboxyl_carbon} | seleninic_idxs
    if any(
        a not in excluded_from_unsaturation_check and b not in excluded_from_unsaturation_check
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a carboxylic acid/"
            "seleninic acid combination is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{s: "selenino" for s in seleninic_idxs}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    seleninic_carbons = {
        n.GetIdx()
        for s in seleninic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    eligible = [
        chain
        for chain in chains
        if carboxyl_carbon in chain and seleninic_carbons.issubset(set(chain))
    ]
    if not eligible:
        raise UnsupportedStructure(
            "not every carboxylic acid/seleninic acid-bearing carbon lies "
            "on a single longest carbon chain; a shorter principal chain "
            "is not supported yet"
        )

    excluded = carboxyl_oxygens
    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != carboxyl_carbon:
                # A -COOH carbon must sit at C1 (P-65.1.1); a direction
                # that doesn't start there is never valid.
                continue
            substituents = _substituents_for_chain(graph, candidate, names, excluded)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
