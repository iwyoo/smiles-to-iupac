"""Naming of carboxylate anions (the '-oate' suffix, R-COO-) on acyclic
saturated or unsaturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-72.2.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html): the
  preferred IUPAC name of an anion formed by removing a hydron from the
  chalcogen atom of an acid is formed by replacing the acid name's 'ic acid'
  ending with 'ate' -- e.g. CH3-CO-O(-) -> 'acetate' (PIN); the retained
  stems come from `_retained_acids.py`.
- The carboxylate carbon's shape mirrors `_carboxylic_acid.py`'s -COOH
  carbon exactly, except the hydroxyl oxygen (-OH, one H, neutral) is
  replaced by an anionic oxygen (no H, formal charge -1): a doubly-bonded,
  monovalent, neutral carbonyl oxygen plus a singly-bonded, monovalent,
  formal-charge -1 oxygen, both on the same carbon. Reuses the same
  chain-numbering and 'ene'/'yne' mechanics as `_carboxylic_acid.py` (the
  carboxylate carbon is always C1, its own locant never cited, P-14.3.3).
- `_name_acyclic_carboxylate` takes an optional `extra_names` map (mirrors
  `_carboxylic_acid.py`'s identical extension point) so `_zwitterion.py`
  can inject a coexisting ammonium nitrogen's prefix name without
  reimplementing this module's own chain search.
- Unlike `_carboxylic_acid.py`, this first pass only supports exactly one
  -COO- group (no 'dioate'), mirroring `_ester.py`'s own single-group scope
  and reusing its `_suffix_body` construction (no ' acid' word, 'oate'
  never multiplied).
- P-35.2.1: halogen substituents are prefix-only and coexist freely with the
  '-oate' suffix, reusing `halogen_substituents`/`format_substituent_prefixes`
  unchanged.
- P-91.3/P-92: a molecule with
  one or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself, no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
  ascending locant order, same pattern as `_carboxylic_acid.py` (the
  carboxylate carbon is always C1, so it never affects numbering). The
  carboxylate carbon itself is never a potential stereocenter (its two
  oxygens are both terminal and the carbon is otherwise sp2-shaped
  relative to the anion resonance), confirmed via RDKit
  `FindPotentialStereo` on `CC[C@@H](C)C(=O)[O-]`.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule (mirrors `_carboxylic_acid.py`'s own
  acyclic-only scope; P-65.1.1.2's ring-attached construction differs).
- More than one -COO- group, or any oxygen that isn't part of the single
  carboxylate's carbonyl/anion pair (an ether, alcohol, second carbonyl,
  or a second -COO-).
- Any charged or radical atom other than the single carboxylate oxygen's
  formal charge -1 (P-73's cations, P-71's radicals, P-74's zwitterions
  are all out of scope).
- Any other heteroatom (N, S, ...).
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    all_chains,
    carbon_adjacency,
    chain_bond_locants,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    lowest_locant_set,
    most_multiple_bonds,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    ring_branch_attachments,
    separate_aromatic_monocycles,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._retained_acids import retained_chain_acid
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    substituents_for_chain,
)

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def has_carboxylate_shape(mol) -> bool:
    """True if some carbon carries both a doubly-bonded, monovalent, neutral
    carbonyl oxygen and a singly-bonded, monovalent, formal-charge -1 oxygen
    (a -COO- pattern), regardless of whether the rest of the molecule is in
    scope. Used by `core.py` to route ahead of the carboxylic-acid/ketone/
    alcohol dispatch, since a carboxylate carbon would otherwise look
    carbonyl-shaped (and its own charged oxygen would be rejected by every
    other module, none of which expect a charged atom)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        has_carbonyl = any(
            o.GetDegree() == 1
            and o.GetFormalCharge() == 0
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for o in oxygens
        )
        has_anion = any(
            o.GetDegree() == 1
            and o.GetFormalCharge() == -1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            for o in oxygens
        )
        if has_carbonyl and has_anion:
            return True
    return False


def _find_carboxylate_group(mol):
    """Locate the molecule's single -COO- group and return
    (carboxylate_carbon, carbonyl_oxygen, anion_oxygen) atoms, after checking
    the molecule has exactly one such group and no other oxygens (see module
    docstring)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and o.GetFormalCharge() == 0
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        anions = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and o.GetFormalCharge() == -1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if carbonyls and anions:
            matches.append((atom, carbonyls, anions))

    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one carboxylate (-COO-) group is required; zero or "
            "multiple such groups are not supported yet (P-72.2.2.2.1.1)"
        )
    carboxylate_carbon, carbonyls, anions = matches[0]
    if len(carbonyls) != 1 or len(anions) != 1:
        raise UnsupportedStructure(
            "a carbon with more than one carbonyl or anionic oxygen does "
            "not match a simple carboxylate group"
        )
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 2:
        raise UnsupportedStructure(
            "an oxygen outside the single carboxylate group's carbonyl/"
            "anion pair (e.g. an ether or a hydroxyl) is out of scope for "
            "this module"
        )

    carbon_neighbors = [n for n in carboxylate_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "a carboxylate carbon with more than one carbon neighbor is "
            "not a valid carboxylate group"
        )

    return carboxylate_carbon, carbonyls[0], anions[0]


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    retained = retained_chain_acid(grouped, chain_length, ene_locants, yne_locants, 1, "anion")
    if retained is not None:
        return retained
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, ene_locants, yne_locants, "oate")


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
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




def _name_acyclic_carboxylate(
    mol, carboxylate_carbon_idx, excluded_oxygens, bonds, stereo=None, extra_names=None, required_atoms=frozenset()
):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`), and the winning candidate's own
    locants are used to format a "(<locant><R/S>,...)-" prefix onto the
    final name.

    `extra_names`: optional {atom_idx -> prefix name} for a coexisting
    characteristic group cited as a substituent prefix (e.g. a
    zwitterion's ammonium-nitrogen 'azaniumyl'-family prefix, see
    `_zwitterion.py`), mirroring `_carboxylic_acid.py`'s identical
    extension point. `None` keeps the original halogen-only behavior
    unchanged. `required_atoms`: additional atoms a candidate chain must
    also carry -- empty by default so existing callers are unaffected."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    chains = all_chains(carbon_adjacency(mol))
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if carboxylate_carbon_idx not in chain:
            continue
        chain_set = set(chain)
        if not required_atoms <= chain_set:
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            carboxylate_carbon_idx in c
            and required_atoms <= set(c)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the carboxylate carbon (and/or multiple bonds) does not lie "
            "on a single longest carbon chain"
        )

    best_key = None
    best_name = None
    best_position_of = None
    chain_length = max(len(c) for c in eligible)
    eligible = most_multiple_bonds([c for c in eligible if len(c) == chain_length], bonds)
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != carboxylate_carbon_idx:
                # The carboxylate carbon must sit at C1 (see module
                # docstring); a direction that doesn't start there is never
                # valid.
                continue
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            ene_locants, yne_locants = chain_bond_locants(candidate, bonds)
            substituents = substituents_for_chain(graph, candidate, halogens, excluded_oxygens, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _validate_and_prepare_carboxylate(mol, aromatic_ring_atoms=frozenset()):
    """(carboxylate_carbon, excluded_oxygens, bonds, stereo) after checking
    the molecule fits this module's scope (see module docstring).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_carboxylate`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_carboxylate`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors
    `_carboxylic_acid.py`'s equivalent aromatic-exemption pattern."""
    carboxylate_carbon, carbonyl_oxygen, anion_oxygen = _find_carboxylate_group(mol)
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), anion_oxygen.GetIdx()}
    stereo = specified_stereocenters(mol)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the carboxylate's own oxygens "
                "(P-72.2.2.2.1.1) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetIdx() in excluded_oxygens:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "single carboxylate oxygen is not supported yet"
            )
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num != 8 and atom.GetDegree() != 1:
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
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_oxygens
        and b[1] not in excluded_oxygens
        and (b[0] not in aromatic_ring_atoms or b[1] not in aromatic_ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    return carboxylate_carbon, excluded_oxygens, bonds, stereo


def _name_phenyl_chain_carboxylate(mol, ring_atoms):
    """Name a carboxylate whose -COO- lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropanoate. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`. Narrower
    than the acyclic path above: no chain unsaturation and no specified
    stereocenter."""
    carboxylate_carbon, excluded_oxygens, bonds, stereo = _validate_and_prepare_carboxylate(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if stereo:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "carboxylate chain is not supported yet"
        )
    non_ring_unsaturation = [b for b in bonds if b[0] not in ring_atoms or b[1] not in ring_atoms]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "carboxylate chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain carboxylate is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain(graph, carboxylate_carbon.GetIdx(), ring_atoms, excluded_oxygens, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -COO- group directly attached to the benzene ring "
            "(benzoate-type naming) uses a separate construction, out of "
            "scope for this acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def name_carboxylate(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        return _name_phenyl_chain_carboxylate(mol, set().union(*aromatic_rings))
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_carboxylate(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COO- group on/in a ring is out of scope for this "
            "acyclic-only module"
        )

    carboxylate_carbon, excluded_oxygens, bonds, stereo = _validate_and_prepare_carboxylate(mol)
    return _name_acyclic_carboxylate(mol, carboxylate_carbon.GetIdx(), excluded_oxygens, bonds, stereo)
