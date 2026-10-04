"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more internal aldehyde (-CHO) carbonyls on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'oic acid' (`_carboxylic_acid.py`) outranks 'al'
  (`_aldehyde.py`), so a coexisting aldehyde is demoted to the 'oxo'
  substituent prefix instead of its own '-al' suffix, e.g.
  'O=CCC(=O)O' (malonaldehydic acid) -> '3-oxopropanoic acid' (a
  well-known worked example). This mirrors `_aldehyde_ketone.py`'s
  aldehyde+ketone demotion, reusing the same {atom_idx -> 'oxo'}
  injection into the chain's
  substituent-prefix machinery -- an aldehyde carbon that is part of the
  parent chain (rather than a branch) is cited as 'oxo', not 'formyl',
  the same way a demoted ketone is.
- Otherwise mirrors `_carboxylic_acid.py` exactly: the -COOH carbon is
  always a chain terminus and always becomes C1 with its own locant never
  cited (P-14.3.3); an aldehyde carbon's locant, by contrast, is always
  cited.

Scope, deliberately narrow (first extension of the general suffix-vs-suffix
demotion problem beyond aldehyde+ketone): a single carboxylic acid plus one
or more aldehydes, all on one acyclic saturated chain, with halogen
substituents allowed. Explicitly out of scope (raise `UnsupportedStructure`):
any chain unsaturation (ene/yne), any ring, more than one carboxylic acid, a
coexisting hydroxyl/ether/other heteroatom, and any aldehyde not captured by
a single longest chain.
"""

from rdkit import Chem

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
from ._retained_acids import retained_chain_acid
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


def _extra_aldehydes(mol, excluded_oxygens):
    aldehydes = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in excluded_oxygens:
            continue
        if atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6:
            continue
        # A genuine aldehyde carbon has exactly one oxygen neighbor (its own
        # carbonyl); a second oxygen neighbor means this is really a -COOH
        # (or other acid-shaped) carbon, which `_find_carboxylic_acid_carbon`
        # already accounts for and must not be double-counted here.
        oxygen_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 1:
            continue
        # A genuine aldehyde carbon has no nitrogen neighbor either -- one
        # means this is really an amide carbon (an amic acid's coexisting
        # -CO-NH2, #820), a separate shape this module doesn't handle,
        # not an aldehyde miscounted by only checking carbon neighbors.
        if any(n.GetAtomicNum() == 7 for n in carbon.GetNeighbors()):
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            aldehydes.add(atom.GetIdx())
    return aldehydes


def has_aldehyde_carboxylic_acid_shape(mol) -> bool:
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return False
    _, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    return bool(_extra_aldehydes(mol, excluded))


def _validate(mol, excluded_oxygens):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than the acid's own oxygens and aldehyde "
        "carbonyl oxygens (P-41, Table 3.3) and halogen substituents "
        "(P-35.2.1) are not supported yet",
        [
            (
                8,
                excluded_oxygens,
                "an oxygen that isn't part of the single carboxylic acid "
                "group or an aldehyde-shaped carbonyl is out of scope for "
                "this module (e.g. a coexisting hydroxyl or ether)",
            ),
        ],
    )


def _name_from_substituents(chain_length, grouped):
    retained = retained_chain_acid(grouped, chain_length, [], [], 1, "acid")
    if retained is not None:
        return retained
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "oic acid")


def _candidate_key(chain_length, grouped):
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name




def name_aldehyde_carboxylic_acid(mol) -> str:
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
    aldehydes = _extra_aldehydes(mol, excluded_acid_oxygens)
    if not aldehydes:
        raise UnsupportedStructure(
            "no coexisting aldehyde found; this module only handles a "
            "carboxylic acid combined with at least one aldehyde (see "
            "_carboxylic_acid.py for a plain carboxylic acid)"
        )
    excluded = excluded_acid_oxygens | aldehydes
    _validate(mol, excluded)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an aldehyde-on-acid "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in aldehydes}}
    acid_carbon_idx = acid_carbon.GetIdx()
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if acid_carbon_idx not in chain:
            continue
        chain_set = set(chain)
        if any(graph[o][0] not in chain_set for o in aldehydes):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every acid/aldehyde-bearing carbon lies on a single "
            "longest carbon chain; a shorter principal chain is not "
            "supported yet"
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
