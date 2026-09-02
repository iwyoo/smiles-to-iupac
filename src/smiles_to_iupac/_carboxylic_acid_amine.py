"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
exactly one primary amine (-NH2) on the same acyclic saturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'oic acid' (`_carboxylic_acid.py`) far outranks 'amine'
  (`_amine.py`), so a coexisting primary amine is demoted to the 'amino'
  substituent prefix instead of its own '-amine' suffix, e.g.
  'NCCCC(=O)O' (4-aminobutanoic acid, GABA) -> '4-aminobutanoic acid'
  (PubChem CID 119's own IUPACName). This mirrors
  `_aldehyde_carboxylic_acid.py`'s carbonyl demotion, except the demoted
  group here is a substituent atom (the amine nitrogen itself), not a
  chain carbon,
  reusing `_carboxylic_acid.py`'s own `{atom_idx -> 'hydroxy'}` injection
  pattern for a coexisting standalone hydroxyl, just with 'amino' instead.
- Otherwise mirrors `_carboxylic_acid.py` exactly: the -COOH carbon is
  always a chain terminus and always becomes C1 with its own locant never
  cited (P-14.3.3); the amine nitrogen's locant, by contrast, is always
  cited via the 'amino' prefix.
- Glycine (H2N-CH2-COOH) is a real-world worked example, but its actual PIN
  ('2-aminoacetic acid', PubChem CID 750) uses the retained acid name
  'acetic acid' rather than the systematic 'ethanoic acid' this project's
  own `_carboxylic_acid.py` always produces (see that module's own tests,
  e.g. 'CC(=O)O' -> 'ethanoic acid', not 'acetic acid') -- this module
  follows that same already-established project convention, so it produces
  '2-aminoethanoic acid' for glycine's structure instead, a pre-existing
  divergence from PubChem inherited from `_carboxylic_acid.py`, not a new
  one introduced here.

Scope, deliberately narrow: a single carboxylic acid plus a single primary
amine, both on one acyclic *saturated* chain, with halogen substituents
allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring, more than one carboxylic acid or amine, a
secondary/tertiary amine, a coexisting hydroxyl/ether/other heteroatom, and
any acid/amine not captured by a single longest chain.

Any *specified* tetrahedral stereocenter (e.g. the alpha carbon of an amino
acid such as alanine/valine) is labeled via `_common.specified_stereocenters`
using the same P-91.3/P-92 mechanism as `_carboxylic_acid.py` (P-92: a
stereocenter on a substituent branch rather than the principal chain is out
of scope).
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    non_single_bonds,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


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


def _find_primary_amines(mol):
    amines = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0:
            continue
        if atom.GetTotalNumHs() != 2:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            amines.add(atom.GetIdx())
    return amines


def _amine_on_a_different_carbon(mol, acid_carbon_idx, amines):
    """False if some amine nitrogen is bonded directly to the acid carbon
    itself (H2N-COOH, carbamic acid -- a distinct retained functional class,
    `_carbamate.py`'s territory, not this module's "amine elsewhere on the
    chain" shape)."""
    return all(next(iter(mol.GetAtomWithIdx(n).GetNeighbors())).GetIdx() != acid_carbon_idx for n in amines)


def has_carboxylic_acid_amine_shape(mol) -> bool:
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return False
    amines = _find_primary_amines(mol)
    return bool(amines) and _amine_on_a_different_carbon(mol, found[0].GetIdx(), amines)


def _validate(mol, excluded_oxygens, amines):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the acid's own oxygens, a primary "
                "amine nitrogen (P-41, Table 3.3), and halogen substituents "
                "(P-35.2.1) are not supported yet"
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
            if atom.GetIdx() not in excluded_oxygens:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the single carboxylic "
                    "acid group is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 7:
            if atom.GetIdx() not in amines:
                raise UnsupportedStructure(
                    "a nitrogen that isn't a plain primary amine (-NH2) is "
                    "out of scope for this module (secondary/tertiary "
                    "amines, imines, and nitriles are not supported)"
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


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, grouped):
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    return prefix + stem + "oic acid"


def _candidate_key(chain_length, grouped):
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


def name_carboxylic_acid_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH and/or -NH2 group on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    amines = _find_primary_amines(mol)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "carboxylic acid is supported here (see _carboxylic_acid.py "
            "for a plain carboxylic acid, _amine.py for a plain amine)"
        )
    if not _amine_on_a_different_carbon(mol, acid_carbon.GetIdx(), amines):
        raise UnsupportedStructure(
            "a primary amine bonded directly to the acid carbon itself "
            "(H2N-COOH, carbamic acid) is a distinct retained functional "
            "class, out of scope here (see _carbamate.py)"
        )
    _validate(mol, excluded_acid_oxygens, amines)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded_acid_oxygens or b[1] in excluded_acid_oxygens]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a carboxylic-acid/"
            "amine seniority demotion is out of scope for this module "
            "(1st-pass scope: saturated only)"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{n: "amino" for n in amines}}
    acid_carbon_idx = acid_carbon.GetIdx()
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo = specified_stereocenters(mol)
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    def _carries_acid_and_amines(chain_set):
        return acid_carbon_idx in chain_set and all(graph[n][0] in chain_set for n in amines)

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if not _carries_acid_and_amines(chain_set):
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(_carries_acid_and_amines(set(c)) for c in chains):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every acid/amine-bearing carbon lies on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acid_carbon_idx:
                # The -COOH carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            substituents = _substituents_for_chain(graph, candidate, names, excluded_acid_oxygens)
            grouped = _group(substituents)
            key, name = _candidate_key(chain_length, grouped)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        # P-91.3: the -COOH carbon's own fixed C1 position (see module
        # docstring) already decides numbering before stereo is considered,
        # mirroring `_carboxylic_acid.py`'s identical treatment.
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name
