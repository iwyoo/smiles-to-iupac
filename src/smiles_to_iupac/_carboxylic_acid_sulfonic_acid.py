"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
one or more sulfonic acid (-SO3H) groups on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41/P-43 (`_seniority.py`): carboxylic acid (Table 4.4 class 1)
  outranks sulfonic acid (class 4), so a coexisting sulfonic acid is
  demoted to the 'sulfo' substituent prefix (P-65.3.1) instead of its own
  '-sulfonic acid' suffix -- e.g. 'OC(=O)CS(=O)(=O)O' -> '2-sulfoacetic
  acid' (PubChem CID 31257 confirms the systematic '2-sulfoacetic acid'
  form; this module always uses the systematic '...anoic acid' stem
  rather than a retained name like 'acetic acid', mirroring
  `_carboxylic_acid.py`'s own established convention once any substituent
  is present). Unlike the three existing `_seniority.py` consumers
  (`_sulfonic_acid_thiol.py`/`_sulfonic_acid_sulfinic_acid.py`/
  `_sulfonic_acid_sulfonamide.py`), which all demote *to* a prefix
  *around* a sulfonic-acid-suffixed parent, this is the first pairwise
  module where the winning suffix is carboxylic acid instead -- so the
  parent-chain naming mirrors `_carboxylic_acid.py` (the -COOH carbon is
  always the chain-terminal C1, P-65.1.1) rather than
  `_sulfonic_acid.py`.
- The acyclic-chain path reuses `_carboxylic_acid.py`'s own naming
  function directly (`_name_acyclic_carboxylic_acid`, via
  `_coexisting_groups.name_via_senior_acyclic`) rather than duplicating
  its chain-search/numbering logic, mirroring the
  `_carboxylic_acid_sulfinic_acid.py`/`_carboxylic_acid_sulfonamide.py`
  migrations (PR #510/#511) -- 'sulfo' is injected as an `extra_names`
  prefix, and the sulfonic-bearing carbon(s) are passed as
  `required_atoms`. The benzene-ring-substituent path below is untouched
  by this (the dispatcher is acyclic-only, same as every other pairwise
  module's phenyl-chain path).

Scope, deliberately narrow (mirrors `_sulfonic_acid_sulfinic_acid.py`): a
single carboxylic acid plus one or more sulfonic acids, all on one
acyclic saturated chain, with halogen substituents allowed. One narrow
*aromatic*-ring exception: `_name_phenyl_chain_carboxylic_acid_sulfonic_acid`
names a carboxylic acid/sulfonic acid chain hanging off a single plain,
unsubstituted benzene ring (e.g. '3-phenyl-3-sulfopropanoic acid',
PubChem CID 20265801), mirroring `_carboxylic_acid.py`'s/
`_sulfonic_acid_thiol.py`'s identical benzene-ring-substituent path.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring other than the single benzene-
substituent exception above, more than one carboxylic acid, a coexisting
standalone hydroxyl/ether/other heteroatom, a specified stereocenter, and
any carboxylic acid/sulfonic acid not captured by a single longest chain.
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._carboxylic_acid import _name_acyclic_carboxylic_acid
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}


def _sulfonic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonic acid group: bonded to exactly
    one carbon, two double-bonded (terminal) oxygens, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfonic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
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
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
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


def has_carboxylic_acid_sulfonic_acid_shape(mol) -> bool:
    return bool(_carboxyl_carbons(mol)) and bool(_sulfonic_sulfur_atoms(mol))


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    carboxyl_carbons = _carboxyl_carbons(mol)
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "exactly one carboxylic acid is required; this module only "
            "handles a single carboxylic acid combined with one or more "
            "sulfonic acids"
        )
    (carboxyl_carbon,) = carboxyl_carbons
    carboxyl_oxygens = {n.GetIdx() for n in carboxyl_carbon.GetNeighbors() if n.GetAtomicNum() == 8}

    sulfonics = _sulfonic_sulfur_atoms(mol)
    if not sulfonics:
        raise UnsupportedStructure(
            "no sulfonic acid found; this module only handles a "
            "carboxylic acid combined with at least one sulfonic acid "
            "(see _carboxylic_acid.py for a plain carboxylic acid)"
        )
    sulfonic_idxs = {s.GetIdx() for s in sulfonics}
    sulfonic_oxygens = {
        o.GetIdx() for s in sulfonics for o in s.GetNeighbors() if o.GetAtomicNum() == 8
    }

    accounted_oxygen_idxs = carboxyl_oxygens | sulfonic_oxygens

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the carboxylic/sulfonic acid "
                "groups' own oxygens (P-65.1.1/P-65.3.1) and halogen "
                "substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in accounted_oxygen_idxs:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the carboxylic/sulfonic "
                    "acid groups is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 16:
            if atom.GetIdx() not in sulfonic_idxs:
                raise UnsupportedStructure(
                    "a sulfur atom not shaped like a sulfonic acid group "
                    "is out of scope for this module"
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

    return carboxyl_carbon.GetIdx(), carboxyl_oxygens, sulfonic_idxs


def _name_from_substituents(chain_length, grouped):
    # P-14.3.4.2(a): a mononuclear parent's substituent locant is always
    # '1' and never cited (mirrors the same rule other `_seniority.py`
    # consumers apply to their own single-carbon case).
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    stem = alkane_name(chain_length)
    # P-14.3.3: the -COOH suffix locant is never cited, since a -COOH
    # carbon is always the chain-terminal C1.
    return f"{prefix}{stem[:-1]}oic acid"


def _substituents_for_chain(graph, chain, names, excluded, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names, mol=mol) for root in branch_roots]
    return substituents


def _name_phenyl_chain_carboxylic_acid_sulfonic_acid(mol, ring_atoms):
    """Name a carboxylic acid plus one or more sulfonic acids, all lying
    on a single unbranched chain hanging off one atom of an otherwise-
    plain, unsubstituted benzene ring -- e.g. '3-phenyl-3-sulfopropanoic
    acid' (PubChem CID 20265801). The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid` (the
    -COOH carbon is always the chain's far terminus from the ring, P-
    65.1.1) combined with `_sulfonic_acid_thiol.py`'s
    `_name_phenyl_chain_sulfonic_acid_thiol` (a demoted heteroatom
    injected into the shared {atom_idx -> name} substituent map).
    Narrower than the acyclic path above: no chain unsaturation and no
    specified stereocenter."""
    carboxyl_carbon, carboxyl_oxygens, sulfonic_idxs = _validate_and_collect(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "carboxylic acid/sulfonic acid chain is not supported yet"
        )
    excluded_from_unsaturation_check = {carboxyl_carbon} | sulfonic_idxs
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_from_unsaturation_check
        and b[1] not in excluded_from_unsaturation_check
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "carboxylic acid/sulfonic acid chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain carboxylic acid/sulfonic acid is not "
            "supported yet"
        )
    ring_atom, chain_root = attachment
    chain, _ = longest_branched_chain(graph, carboxyl_carbon, ring_atoms, carboxyl_oxygens | sulfonic_idxs, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -COOH group directly attached to the benzene ring (no "
            "intervening chain carbon) uses the separate 'carboxylic "
            "acid' suffix construction (P-65.1.1.2), out of scope for "
            "this acyclic-chain-parent module"
        )
    chain_set = set(chain)
    sulfonic_carbons = {
        n.GetIdx()
        for s in sulfonic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    if not sulfonic_carbons.issubset(chain_set):
        raise UnsupportedStructure(
            "a sulfonic acid outside the single unbranched chain hanging "
            "off the benzene ring is not supported yet"
        )

    chain_length = len(chain)
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    names = {**halogen_substituents(mol), **{s: "sulfo" for s in sulfonic_idxs}}
    substituents = _substituents_for_chain(graph, chain, names, carboxyl_oxygens | {ring_atom}, mol=mol)
    ring_entry = name_branch(graph, ring_atom, chain_root, names, ring_atoms, mol=mol)
    substituents.setdefault(position_of[chain_root], []).append(ring_entry)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def name_carboxylic_acid_sulfonic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_carboxylic_acid_sulfonic_acid(mol, ring_atoms)

    carboxyl_carbon, carboxyl_oxygens, sulfonic_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a carboxylic acid/sulfonic acid combination on a ring is out "
            "of scope for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    excluded_from_unsaturation_check = {carboxyl_carbon} | sulfonic_idxs
    if any(
        a not in excluded_from_unsaturation_check and b not in excluded_from_unsaturation_check
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a carboxylic acid/"
            "sulfonic acid combination is out of scope for this module"
        )

    sulfonic_carbons = {
        n.GetIdx()
        for s in sulfonic_idxs
        for n in mol.GetAtomWithIdx(s).GetNeighbors()
        if n.GetAtomicNum() == 6
    }
    return name_via_senior_acyclic(
        _name_acyclic_carboxylic_acid,
        "carboxylic_acid",
        "sulfonic_acid",
        (mol, {carboxyl_carbon}, carboxyl_oxygens, set(), []),
        {s: "sulfo" for s in sulfonic_idxs},
        required_atoms=sulfonic_carbons,
    )
