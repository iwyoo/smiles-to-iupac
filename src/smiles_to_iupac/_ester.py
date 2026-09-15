"""Naming of esters (the '-oate' suffix, R-COO-R') on acyclic saturated or
unsaturated acyl chains, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-65.6.3, Table 3.3 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf
  for P-65.6.3; Table 3.3 lives in Chapter P-3,
  https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): an ester R-CO-O-R' is named
  as two words, "R'yl R-oate" — the alcohol part (R') cited first as a plain
  substituent-group name, then the acyl part (R) named the same way a
  carboxylic acid's R-COOH would be but with 'oate' instead of 'oic acid'.
  'oate' ranks junior to 'oic acid' and senior to 'amide' in Table 3.3.
- This module's scope (see also `_carboxylic_acid.py`/`_aldehyde.py`, the
  most similar existing modules): the acyl part (R) may be any acyclic
  saturated or unsaturated chain, with halogen substituents, the same as
  `_carboxylic_acid.py`'s R. The alcohol part (R') is restricted to a
  plain, unsubstituted, saturated, acyclic alkyl group (branched or
  unbranched) attached at its own chain terminus (e.g. 'methyl',
  'propan-2-yl', 'tert-butyl') — a substituted, unsaturated, or
  ring-bearing R' is deferred (see "Explicitly out of scope" below).
- The acyl carbon is always a chain terminus (its remaining two bonds, after
  the carbonyl and ester oxygens, allow at most one more substituent, which
  must be another chain carbon, or nothing for a formate ester), so it is
  never a genuine locant choice: it is always C1 of the acyl chain, the same
  way `_carboxylic_acid.py`'s -COOH carbon is, and that locant is never
  cited (P-14.3.3).
- Since this module only ever handles a single, isolated ester group (see
  scope below), the 'oate' word is never multiplied ('dioate' etc. is out
  of scope), unlike `_carboxylic_acid.py`'s 'dioic acid'.
- P-31.0 / P-31.1.1.1-.2: construction of the 'ene'/'yne' portion of the
  acyl part's name reuses the same mechanics as `_carboxylic_acid.py`.
- P-35.2.1: halogen substituents on the acyl chain are prefix-only and
  coexist freely with the 'oate' suffix, reusing `halogen_substituents`/
  `format_substituent_prefixes` unchanged.
- P-29.3.2.1: the alcohol part's name is a plain alkyl substituent-group
  name (P-13.2.1's "R'yl" role, not a locanted prefix), built with
  `name_branch` (P-29 PIN style, fixed project-wide by PR #237 --
  previously this module avoided `name_branch` for exactly this reason,
  but that blocker no longer applies; mirrors `_carbamate.py`'s
  identical fix, PR #328), e.g. 'propan-2-yl' for R' = isopropyl,
  matching the verified PIN 'propan-2-yl acetate' (PubChem CID 7915)
  and 'tert-butyl acetate' (CID 10908). R' is never parenthesized
  regardless of `name_branch`'s `is_compound` flag -- the "R'yl R-oate"
  two-word pattern has no nested-prefix ambiguity to guard against.
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters on the acyl chain --
  every one on the principal chain itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, cited
  immediately ahead of the acyl part specifically rather than the whole
  two-word name, e.g. 'ethyl (2R)-2-methylbutanoate' (PubChem CID
  7156991), same pattern as `_carboxylic_acid.py`/`_amide.py`. The
  alcohol part (R') is always a plain unbranched chain (above), so it can
  never itself hold a stereocenter -- this support only ever concerns the
  acyl side.

- P-65.6.3.2 (aryl esters): the alcohol part (R') may also be a single,
  otherwise-unsubstituted benzene ring attached directly to the ester
  oxygen (Ar-O-CO-R, e.g. 'phenyl ethanoate' for PhO-CO-CH3, confirmed via
  PubChem's own 'phenyl acetate'/'phenyl formate'/'phenyl propanoate' --
  this project's systematic-name convention carries over the same way it
  does for the acyl side). The acyl part (R) reuses the ordinary acyclic
  path unchanged (`_name_acyl_part`) since the ring and the acyl chain sit
  on opposite sides of the ester oxygen, never sharing a carbon-adjacency
  component. A ring-substituted phenol (more than one exocyclic
  attachment), a heteroaromatic or polycyclic alcohol part, or a genuine
  ring-embedded lactone (the carbonyl carbon itself inside the ring, a
  structurally different construction, P-65.6.3.3) are all still out of
  scope, deferred to a future pass.

- The alcohol part (R') may likewise be a single, otherwise-unsubstituted,
  saturated monocyclic ring attached directly to the ester oxygen (e.g.
  'cyclohexyl ethanoate', PubChem CID 12146's own 'cyclohexyl acetate') --
  the saturated counterpart of the aryl-ester case above, reusing
  `name_branch`'s existing plain-ring recognition (`_simple_ring_
  substituent`, already used by `_alcohol.py`/`_ketone.py`/`_carbamate.py`
  for the same ring-vs-chain shape) rather than new naming logic. The acyl
  part again reuses `_name_acyl_part` unchanged for the same
  opposite-sides-of-the-ester-oxygen reason. A ring bearing its own
  substituent/unsaturation, or a polycyclic shape, is still out of scope.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one ring, or a ring elsewhere in the molecule alongside the
  aryl-ester/cyclyl-ester alcohol-part ring above; a lone ring on the acyl
  side (see `_name_phenyl_acyl_ester`) is separately supported, but the
  ring paths don't combine, and a ring-embedded lactone (P-65.6.3.3) uses
  a different naming construction entirely.
- More than one ester group (a diester -- see `_diester_acyloxy.py` for a
  narrow, separately-scoped diester axis), or any oxygen that isn't part
  of the single ester's carbonyl/ester-oxygen pair (an ether, alcohol, or
  second carbonyl elsewhere).
- A substituted or unsaturated alcohol part (R'), or a cyclic one other
  than the plain-benzene aryl case or the plain-saturated-monocyclic
  cyclyl case above — only a plain saturated acyclic alkyl R' (branched or
  unbranched), a plain unsubstituted benzene ring, or a plain
  unsubstituted saturated monocyclic ring, is supported in this first
  pass.
- Any other heteroatom (N, S, ...).
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
    component_subgraph,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    longest_chains,
    lowest_locant_set,
    name_from_substituents,
    non_single_bonds,
    plain_saturated_ring_substituent_atoms,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._carboxylic_acid import _name_benzo_attached_carboxyl, _name_ring_attached_carboxyl
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
)

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def has_ester_shape(mol) -> bool:
    """True if some carbon carries both a doubly-bonded, monovalent oxygen
    and a singly-bonded oxygen that is itself bonded to a second carbon (a
    -C(=O)-O-C- pattern), regardless of whether the rest of the molecule is
    in scope. Used by `core.py` to route ahead of the carboxylic-acid/
    aldehyde/ketone dispatch, since an ester carbon would otherwise look
    carboxylic-acid- or carbonyl-shaped to those modules."""
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
        has_ester_oxygen = any(
            o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
            for o in oxygens
        )
        if has_carbonyl and has_ester_oxygen:
            return True
    return False


def _find_ester_group(mol):
    """Locate the molecule's single ester group and return
    (acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon) atoms, after
    checking the molecule has exactly one such group and no other oxygens
    (see module docstring)."""
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
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if carbonyls and ester_oxygens:
            matches.append((atom, carbonyls, ester_oxygens))

    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one ester group is required; zero or multiple ester "
            "groups (e.g. a diester) are not supported yet (P-65.6.3)"
        )
    acyl_carbon, carbonyls, ester_oxygens = matches[0]
    if len(carbonyls) != 1 or len(ester_oxygens) != 1:
        raise UnsupportedStructure(
            "an ester carbon with more than one carbonyl or ester oxygen "
            "does not match a simple ester group"
        )
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 2:
        raise UnsupportedStructure(
            "an oxygen outside the single ester group's carbonyl/ester pair "
            "(e.g. an ether or a hydroxyl) is out of scope for this module"
        )

    acyl_carbon_neighbors = [n for n in acyl_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(acyl_carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "an ester carbon with more than one carbon neighbor besides its "
            "two ester oxygens is not a valid ester group"
        )

    carbonyl_oxygen = carbonyls[0]
    ester_oxygen = ester_oxygens[0]
    alcohol_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != acyl_carbon.GetIdx())
    return acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon


def _alcohol_component(full_graph, alcohol_carbon_idx, ester_oxygen_idx):
    component = set()
    stack = [alcohol_carbon_idx]
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        for neighbor in full_graph[node]:
            if neighbor != ester_oxygen_idx and neighbor not in component:
                stack.append(neighbor)
    return component


def _name_alcohol_part(mol, alcohol_carbon, ester_oxygen_idx):
    full_graph = adjacency(mol)
    component = _alcohol_component(full_graph, alcohol_carbon.GetIdx(), ester_oxygen_idx)
    for idx in component:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "the alcohol part (R') must be a plain alkyl group; "
                "heteroatoms/halogens there are not supported yet (P-65.6.3)"
            )
    for a, b, _ in non_single_bonds(mol):
        if a in component and b in component:
            raise UnsupportedStructure(
                "unsaturation in the alcohol part (R') is not supported yet "
                "(P-65.6.3)"
            )

    name, _ = name_branch(full_graph, alcohol_carbon.GetIdx(), ester_oxygen_idx, {}, mol=mol)
    return name


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
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


def _substituents_for_chain(graph, chain, halogens, excluded_oxygens, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_oxygens]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_acyl_part(
    mol,
    acyl_carbon,
    carbonyl_oxygen_idx,
    ester_oxygen_idx,
    stereo=None,
    ring_atoms=frozenset(),
    extra_names=None,
    required_atoms=frozenset(),
):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylic_acid.py`/`_amide.py`'s identical treatment),
    and the winning candidate's own locants are used to format a
    "(<locant><R/S>,...)-" prefix onto the acyl name -- cited immediately
    ahead of the acyl part specifically (P-91.3), not the whole two-word
    ester name, confirmed against PubChem's own
    'ethyl (2R)-2-methylbutanoate' (CID 7156991). The alcohol part (R') is
    always a plain unbranched chain, or the single benzene ring
    `ring_atoms` names (`_name_phenol_ester`) -- either way it can never
    itself hold a stereocenter. `ring_atoms`: the alcohol-side aryl-ester
    ring's own atoms (empty by default), excluded here so this function's
    own unsaturation scan doesn't mistake the ring's internal aromatic
    bonds for acyl-chain unsaturation -- the acyl chain itself never
    shares a carbon-adjacency component with that ring regardless (they
    sit on opposite sides of the ester oxygen), so this exclusion is only
    needed for this scan, not the chain search below.

    `extra_names`/`required_atoms`: same coexisting-group injection point
    as `_amide.py`'s `_name_acyclic_amide`, reused by `_ester_amine.py` via
    `_coexisting_groups.py` -- both empty/None by default so existing
    callers are unaffected."""
    full_graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    halogens = {**halogen_substituents(mol), **(extra_names or {})}
    acyl_carbon_idx = acyl_carbon.GetIdx()
    excluded_oxygens = {carbonyl_oxygen_idx, ester_oxygen_idx}

    acyl_graph = component_subgraph(carbon_graph, acyl_carbon_idx)
    chains = longest_chains(acyl_graph)
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_oxygens
        and b[1] not in excluded_oxygens
        and not (b[0] in ring_atoms and b[1] in ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    eligible = []
    for chain in chains:
        if not required_atoms <= set(chain):
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            not bonds or bond_locants(c, bonds) is not None for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the acyl chain's unsaturation does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != acyl_carbon_idx:
                # The ester carbon must sit at C1 (see module docstring); a
                # direction that doesn't start there is never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(full_graph, candidate, halogens, excluded_oxygens, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of
    if best_name is None:
        raise UnsupportedStructure(
            "the ester's acyl carbon does not lie on a single longest "
            "carbon chain; a shorter principal chain is not supported yet"
        )

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _validate_ester_atoms(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule's atoms fit this module's scope (see module
    docstring): only C/O/halogen, no charges/isotopes, no aromatic carbon
    outside `aromatic_ring_atoms`, and at least one carbon.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_ester`'s benzene-ring-
    substituent path (see `_name_phenyl_acyl_ester`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the ester's own oxygens (P-65.6.3) "
                "and halogen substituents (P-35.2.1) are not supported yet"
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


def _name_phenyl_acyl_ester(mol, ring_atoms):
    """Name an ester whose acyl part (R-CO-) chain hangs off a single
    unbranched chain attached to an otherwise-plain, unsubstituted benzene
    ring -- e.g. 'methyl 3-phenylpropanoate'. Mirrors
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid` and
    `_alcohol.py`'s `_name_phenyl_chain_alcohol`, adapted for the ester's
    two-word (alcohol + acyl) name; `name_ester` only calls this when the
    ring sits on the acyl side (see `_name_phenol_ester` for the alcohol
    side instead) -- `_name_alcohol_part` is reused unchanged for the
    plain acyclic alcohol part here. Narrower than the general acyl path: no chain unsaturation and
    no specified stereocenter -- each a separate follow-up (see
    tasks/phenyl-substituent-on-ester-chain.md's scope note)."""
    _validate_ester_atoms(mol, aromatic_ring_atoms=ring_atoms)
    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "ester acyl chain is not supported yet"
        )
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_oxygens
        and b[1] not in excluded_oxygens
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent ester "
            "acyl chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain ester is not "
            "supported yet"
        )
    acyl_carbon_idx = acyl_carbon.GetIdx()
    chain, branches = longest_branched_chain(graph, acyl_carbon_idx, ring_atoms, excluded_oxygens, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an ester group directly attached to the benzene ring (no "
            "intervening chain carbon) is out of scope for this "
            "acyclic-chain-parent module"
        )

    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    acyl_name = _name_from_substituents(chain_length, [], [], grouped)
    return f"{alcohol_name} {acyl_name}"


def _name_benzoate_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon):
    """Name an ester whose acyl part (R-CO-) hangs directly off one carbon
    of an otherwise-plain (or substituted) benzene ring, with no
    intervening chain carbon (P-65.6.3.2's retained-name axis) -- e.g.
    'methyl benzoate' (PubChem CID 7150), '2-methylbenzoate' esters.
    Reuses `_carboxylic_acid.py`'s `_name_benzo_attached_carboxyl` kernel
    (the same ring-numbering search `_name_benzoic_acid` uses) with the
    'benzoate' suffix word in place of 'benzoic acid'; the alcohol part
    reuses `_name_alcohol_part` unchanged, same reasoning as
    `_name_phenol_ester`/`_name_cyclyl_ester`: the ring and the alcohol
    part sit on opposite sides of the ester oxygen."""
    _validate_ester_atoms(mol, aromatic_ring_atoms=ring_atoms)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzoate ester is not "
            "supported yet"
        )
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    all_non_single = non_single_bonds(mol)
    if any(
        a not in excluded_oxygens
        and b not in excluded_oxygens
        and (a not in ring_atoms or b not in ring_atoms)
        for a, b, _ in all_non_single
    ):
        raise UnsupportedStructure(
            "unsaturation outside the ring alongside a benzoate ester is "
            "not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_benzo_attached_carboxyl(graph, ring_atoms, acyl_carbon.GetIdx(), halogens, word="benzoate")
    return f"{alcohol_name} {acyl_name}"


def _name_ring_acyl_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon):
    """Saturated-ring counterpart of `_name_benzoate_ester`: the acyl part
    hangs directly off one atom of an otherwise-plain saturated
    monocyclic ring, with no intervening chain carbon (P-65.6.3.2's
    retained-name axis, extended to von Baeyer ring stems) -- e.g. 'methyl
    cyclohexanecarboxylate' (PubChem CID 20748). Reuses
    `_carboxylic_acid.py`'s `_name_ring_attached_carboxyl` kernel (the
    same ring-numbering search `_name_ring_carboxylic_acid` uses) with
    the 'carboxylate' suffix word in place of 'carboxylic acid'."""
    _validate_ester_atoms(mol)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a ring-attached ester "
            "acyl group is not supported yet"
        )
    graph = adjacency(mol)
    acyl_idx = acyl_carbon.GetIdx()
    (ring_atom,) = [n for n in graph[acyl_idx] if n in ring_atoms]
    other_ring_atom_branches = [n for n in graph[ring_atom] if n not in ring_atoms and n != acyl_idx]
    if other_ring_atom_branches:
        raise UnsupportedStructure(
            "a substituent on the same ring atom as the ester acyl group "
            "is not supported yet"
        )
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    all_non_single = non_single_bonds(mol)
    if any(a not in excluded_oxygens and b not in excluded_oxygens for a, b, _ in all_non_single):
        raise UnsupportedStructure(
            "an unsaturated ring alongside a ring-attached ester acyl "
            "group is not supported yet (see P-31.1.3)"
        )

    halogens = halogen_substituents(mol)
    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_ring_attached_carboxyl(graph, ring_atoms, acyl_idx, halogens, suffix="carboxylate")
    return f"{alcohol_name} {acyl_name}"


def _name_ring_acyl_chain_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon):
    """Saturated-ring counterpart of `_name_phenyl_acyl_ester`: the acyl
    part (R-CO-) chain hangs off a single unbranched chain attached to an
    otherwise-plain saturated monocyclic ring -- e.g. 'methyl
    2-cyclohexylacetate'. `name_branch`'s own plain-ring recognition
    (`_simple_ring_substituent`) names the ring as "cyclo..." with no
    `aromatic_atoms` needed, unlike the benzene case. Narrower than the
    benzene case: no ring alkyl/halogen substituents, matching this
    project's established saturated-ring-chain first-slice convention
    (see `_ketone.py`'s `_name_ring_substituent_chain_ketone`)."""
    _validate_ester_atoms(mol)
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a saturated-ring-"
            "substituent ester acyl chain is not supported yet"
        )
    excluded_oxygens = {carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx()}
    if any(a not in excluded_oxygens and b not in excluded_oxygens for a, b, _ in non_single_bonds(mol)):
        # Unlike the benzene case (whose ring bonds are inherently
        # aromatic and always correctly named "phenyl"), `name_branch`'s
        # plain-ring detector (`_simple_ring_substituent`) doesn't check
        # bond order at all, so an unsaturated saturated-ring shape
        # (cyclohexenyl) would otherwise be silently misnamed as the
        # plain "cyclohexyl" -- reject any non-single bond outside the
        # ester's own carbonyl, on the ring or the chain alike, rather
        # than risk that silent mismatch.
        raise UnsupportedStructure(
            "unsaturation on the ring or chain alongside a saturated-"
            "ring-substituent ester acyl chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch alongside a "
            "chain ester is not supported yet"
        )
    acyl_carbon_idx = acyl_carbon.GetIdx()
    chain, branches = longest_branched_chain(graph, acyl_carbon_idx, ring_atoms, excluded_oxygens, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an ester group directly attached to the ring (no intervening "
            "chain carbon) is out of scope for this acyclic-chain-parent "
            "module"
        )

    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())

    halogens = halogen_substituents(mol)
    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    acyl_name = _name_from_substituents(chain_length, [], [], grouped)
    return f"{alcohol_name} {acyl_name}"


def _name_phenol_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen):
    """Name an ester whose alcohol part (R') is a single, otherwise-plain
    benzene ring attached directly to the ester oxygen (Ar-O-CO-R,
    P-65.6.3.2) -- e.g. 'phenyl ethanoate' (PubChem CID 31229's own
    'phenyl acetate', this project's systematic-name convention applied).
    The acyl part reuses the ordinary acyclic path unchanged
    (`_name_acyl_part`): the ring and the acyl chain sit on opposite sides
    of the ester oxygen, so they never share a carbon-adjacency component
    and that function needs no change to handle this case."""
    _validate_ester_atoms(mol, aromatic_ring_atoms=ring_atoms)
    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None or attachment[1] != ester_oxygen.GetIdx():
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a phenol ester is not supported yet"
        )
    stereo = specified_stereocenters(mol)
    acyl_name = _name_acyl_part(
        mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx(), stereo, ring_atoms
    )
    return f"phenyl {acyl_name}"


def _name_cyclyl_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon):
    """Name an ester whose alcohol part (R') is a single, otherwise-plain,
    saturated monocyclic ring attached directly to the ester oxygen (e.g.
    'cyclohexyl ethanoate', PubChem CID 12146's own 'cyclohexyl acetate',
    this project's systematic-name convention applied) -- the saturated
    counterpart of `_name_phenol_ester`'s aromatic case. The acyl part
    reuses the ordinary acyclic path unchanged (`_name_acyl_part`), same
    reasoning as `_name_phenol_ester`: the ring and the acyl chain sit on
    opposite sides of the ester oxygen, so they never share a
    carbon-adjacency component."""
    _validate_ester_atoms(mol)
    graph = adjacency(mol)
    stereo = specified_stereocenters(mol)
    acyl_name = _name_acyl_part(
        mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx(), stereo, ring_atoms
    )
    alcohol_name, _ = name_branch(graph, alcohol_carbon.GetIdx(), ester_oxygen.GetIdx(), {}, mol=mol)
    return f"{alcohol_name} {acyl_name}"


def name_ester(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
            if alcohol_carbon.GetIdx() in ring_atoms:
                return _name_phenol_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen)
            graph = adjacency(mol)
            if any(n in ring_atoms for n in graph[acyl_carbon.GetIdx()]):
                return _name_benzoate_ester(
                    mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon
                )
            return _name_phenyl_acyl_ester(mol, ring_atoms)
        acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
        if acyl_carbon.GetIdx() in ring_atoms:
            # A ring-embedded lactone (the carbonyl carbon itself in the
            # ring, e.g. a macrolactone) is a structurally different
            # construction (P-65.6.3.3) from the ring-attached-acyl shapes
            # below, which all assume the acyl carbon sits outside the ring.
            raise UnsupportedStructure(
                "a ring-attached ester or lactone uses a different naming "
                "construction (P-65.6.3.2/P-65.6.3.3), out of scope for this "
                "acyclic-only module"
            )
        cyclyl_ring_atoms = plain_saturated_ring_substituent_atoms(
            mol, adjacency(mol), ester_oxygen.GetIdx(), alcohol_carbon.GetIdx()
        )
        if cyclyl_ring_atoms:
            return _name_cyclyl_ester(mol, cyclyl_ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon)
        graph = adjacency(mol)
        if any(n in ring_atoms for n in graph[acyl_carbon.GetIdx()]):
            return _name_ring_acyl_ester(mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon)
        acyl_component = _alcohol_component(graph, acyl_carbon.GetIdx(), ester_oxygen.GetIdx())
        if ring_atoms & acyl_component:
            return _name_ring_acyl_chain_ester(
                mol, ring_atoms, acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon
            )
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached ester or lactone uses a different naming "
            "construction (P-65.6.3.2/P-65.6.3.3), out of scope for this "
            "acyclic-only module"
        )
    _validate_ester_atoms(mol)

    acyl_carbon, carbonyl_oxygen, ester_oxygen, alcohol_carbon = _find_ester_group(mol)
    stereo = specified_stereocenters(mol)
    alcohol_name = _name_alcohol_part(mol, alcohol_carbon, ester_oxygen.GetIdx())
    acyl_name = _name_acyl_part(mol, acyl_carbon, carbonyl_oxygen.GetIdx(), ester_oxygen.GetIdx(), stereo)
    return f"{alcohol_name} {acyl_name}"
