"""Naming of a plain, saturated, monocyclic amine ring (piperidine/
pyrrolidine/azepane, and the 1,4-two-heteroatom morpholine/piperazine/
thiomorpholine siblings, P-22.2.1) whose nitrogen carries exactly one
substituent -- e.g. 1-methylpiperidine -- per the IUPAC 2013
Recommendations ("the Blue Book"):

- `_amine.py` explicitly defers "a secondary/tertiary amine nitrogen on
  or attached to a ring" (see that module's own docstring); this module
  fills in the most basic slice of that gap: the ring itself is the
  parent hydride (its retained/Hantzsch-Widman name, via
  `_hetero_monocyclic.saturated_ring_name`, already used by
  `_ketone.py`'s hetero-ring-ketone path and `_hidden_amide_ketone.py`'s
  pseudoketone path), and the nitrogen's one substituent is cited as an
  ordinary numbered prefix at locant '1' -- the ring's own heteroatom
  locant (P-22.2.1), not the special 'N-' prefix `_amine.py` uses for an
  acyclic secondary/tertiary amine. PubChem structure matches: 'CN1CCCCC1'
  -> '1-methylpiperidine' (CID 12291), 'CCN1CCCCC1' -> '1-ethylpiperidine'
  (CID 13007), 'CN1CCCC1' -> '1-methylpyrrolidine' (CID 8454),
  'ClCCN1CCCCC1' -> '1-(2-chloroethyl)piperidine' (CID 74827),
  'C1CC1N1CCCCC1' -> '1-cyclopropylpiperidine' (CID 10909573).
- The substituent's own name is built with `name_branch` unchanged
  (P-29 PIN style) -- a branched, halogenated, or plain-cyclic substituent
  is supported for free the same way `_carbamate.py`'s/`_ester.py`'s R
  parts are; a *compound* substituent (its own locant, e.g. 'propan-2-yl')
  is parenthesized, this project's usual convention over PubChem's own
  unparenthesized raw name ('1-propan-2-ylpiperidine' -> this module's
  '1-(propan-2-yl)piperidine').
- A second path handles the 1,4-two-heteroatom sibling rings
  (morpholine/piperazine/thiomorpholine, the same six-membered element
  pairs `_ketone.py`'s `_hetero_ring_two_heteroatoms` already recognizes)
  with one substituent on the ring's nitrogen -- e.g. 'CN1CCOCC1' ->
  '4-methylmorpholine' (PubChem CID 7972). The non-nitrogen heteroatom
  always wins locant 1 (P-22.2.1 element seniority O > S > N, mirroring
  `_ketone.py`'s own `_TWO_HETERO_PRIORITY` table -- not imported, since
  this project's convention accepts each module keeping its own small
  private copy of such a table, e.g. `_is_carbonyl_carbon` duplicated
  across `_amide.py`/`_hydrazide.py`), so the substituted nitrogen is
  always locant 4 there; for the symmetric N,N piperazine pair, the
  substituted nitrogen is always locant 1 (P-14.5.2 lowest locant to the
  cited substituent) and the plain nitrogen is 4. PubChem structure
  matches: 'CCN1CCOCC1' -> '4-ethylmorpholine' (CID 7525), 'CN1CCSCC1' ->
  '4-methylthiomorpholine' (CID 523249), 'CN1CCNCC1' ->
  '1-methylpiperazine' (CID 53167).

- A third path handles a sulfonyl group (P-65.3.1's -SO2- substituent
  prefix, not the '-sulfonic acid'/'-sulfonamide' suffix) as the ring
  nitrogen's substituent -- e.g. 'CS(=O)(=O)N1CCCCC1' ->
  '1-methylsulfonylpiperidine' (PubChem CID 273952). The prefix fuses the
  R group's own plain substituent name directly with 'sulfonyl' (not the
  '-ane'-stem acid form, i.e. 'methylsulfonyl' not 'methanesulfonyl'),
  inheriting R's own compound/parenthesization status unchanged -- e.g.
  'ClCS(=O)(=O)N1CCCCC1' -> '1-(chloromethylsulfonyl)piperidine' (CID
  1519936, since 'chloromethyl' is itself a compound substituent name,
  P-14.3.4.2(a)) but 'CCS(=O)(=O)N1CCCCC1' -> '1-ethylsulfonylpiperidine'
  (CID 3418537, no parens -- plain 'ethyl' isn't compound). Works
  uniformly with both the single- and two-heteroatom ring paths above:
  'CS(=O)(=O)N1CCOCC1' -> '4-methylsulfonylmorpholine' (CID 519344),
  'CS(=O)(=O)N1CCNCC1' -> '1-methylsulfonylpiperazine' (CID 709161).

`core.py` routes to this module from three places: the oxygen-free "any
nitrogen" branch (immediately before its final `name_amine` fallback --
piperidine/pyrrolidine/azepane and the N,N/N,S piperazine/thiomorpholine
siblings all have no oxygen when the substituent is plain alkyl); the
oxygen-gated branch ahead of `has_ether_shape`, for the N,O morpholine
sibling with a plain-alkyl substituent (morpholine's ring oxygen would
otherwise look like a perfectly ordinary ether oxygen to that check); and
-- via the narrower `has_ring_amine_sulfonyl_shape` -- ahead of
`has_sulfonamide_shape`, near the very start of `core.py`'s dispatch, for
*any* ring shape with a sulfonyl substituent (which brings its own
oxygen along regardless of ring type, and would otherwise be misclaimed
by that earlier, ring-unaware check). Never collides with
`_hidden_amide_ketone.py`'s acyl-on-ring-nitrogen path -- that module's
shape requires the substituent's own carbonyl oxygen and is routed
separately, and `has_ring_amine_sulfonyl_shape`'s extra sulfonyl check
(unlike the plain `has_ring_amine_shape`) never matches an acyl
substituent, so it's safe to route this one so early.

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring nitrogen with zero substituents (the plain unsubstituted ring
  itself, already named elsewhere) or more than one (structurally
  impossible for a neutral trivalent ring nitrogen with two ring bonds
  anyway).
- An acyl N-substituent -- handled by `_hidden_amide_ketone.py` instead
  (routed earlier).
- An aromatic-R or acyl-R sulfonyl N-substituent, or any other
  oxygen-bearing N-substituent (sulfinyl, phosphoryl, ...) -- a plain
  alkyl/halogenated/cyclic-R sulfonyl group is supported (see the
  P-65.3.1 note above), but that's the only oxygen-bearing shape this
  module recognizes so far.
- An aromatic N-substituent, a heteroatom other than N/O/S, a ring size
  other than 5/6/7 (single-heteroatom path) or other than 6
  (two-heteroatom path, `saturated_two_heteroatom_1_4_ring_name`'s own
  scope), a Se/Te two-heteroatom pair (no PubChem-registered N-substituted
  name to confirm the locant, same reasoning as `_ketone.py`'s identical
  exclusion), N,N'-disubstituted piperazine (both nitrogens substituted --
  a separate follow-up), any substituent on the ring itself, or any ring
  unsaturation.
- Any other heteroatom, charged/isotopically modified atom, or
  multi-fragment structure.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    non_single_bonds,
    ring_cycle,
)
from ._hetero_monocyclic import saturated_ring_name, saturated_two_heteroatom_1_4_ring_name
from ._substituents import name_branch

_RING_SIZES = (5, 6, 7)
_ALLOWED_SUBSTITUENT_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}
_TWO_HETERO_RING_SIZE = 6
_TWO_HETERO_ELEMENTS = {7: "N", 8: "O", 16: "S"}
_TWO_HETERO_ELEMENT_PAIRS = {
    frozenset(("N", "O")),
    frozenset(("N", "N")),
    frozenset(("N", "S")),
}
# P-22.2.1 element seniority for locant 1 (mirrors `_ketone.py`'s own
# `_TWO_HETERO_PRIORITY`, O > S > N) -- the non-nitrogen heteroatom always
# wins locant 1 for the N,O/N,S pairs, so the substituted nitrogen is
# always locant 4 there.
_TWO_HETERO_PRIORITY = {"O": 0, "S": 1, "N": 2}


def _ring_amine_shape(mol):
    """(n_idx, ring_atoms, substituent_root_idx) if `mol` has exactly one
    plain, saturated, monocyclic ring (5/6/7-membered) with one nitrogen
    heteroatom bearing exactly one exocyclic substituent -- else None. A
    *second*, unrelated ring may still exist inside that one substituent
    (e.g. a cyclopropyl N-substituent) -- `name_branch` recognizes that
    shape on its own; only the amine ring itself is required to be the
    sole N-heterocycle here."""
    ring_info = mol.GetRingInfo()
    candidates = [
        set(ring)
        for ring in ring_info.AtomRings()
        if len(ring) in _RING_SIZES
        and not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
        and sum(1 for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() != 6) == 1
    ]
    if len(candidates) != 1:
        return None
    (ring_atoms,) = candidates
    (n_idx,) = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    n_atom = mol.GetAtomWithIdx(n_idx)
    if n_atom.GetAtomicNum() != 7 or n_atom.GetFormalCharge() != 0 or n_atom.GetIsotope() != 0:
        return None
    exo = [n for n in n_atom.GetNeighbors() if n.GetIdx() not in ring_atoms]
    if len(exo) != 1:
        return None
    if mol.GetBondBetweenAtoms(n_idx, exo[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    for a in ring_atoms:
        if a == n_idx:
            continue
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if any(n.GetIdx() not in ring_atoms for n in atom.GetNeighbors()):
            return None
        if atom.GetTotalNumHs() != 2:
            return None
    if any(a in ring_atoms and b in ring_atoms for a, b, _ in non_single_bonds(mol)):
        return None
    return n_idx, ring_atoms, exo[0].GetIdx()


def _two_hetero_ring_amine_shape(mol):
    """(n_idx, ring_atoms, substituent_root_idx, n_locant, elements) for a
    plain, saturated, six-membered, 1,4-two-heteroatom ring (morpholine/
    piperazine/thiomorpholine's element pairs only) whose nitrogen carries
    exactly one exocyclic substituent -- else None. `elements` is the
    frozenset PubChem-verified pair for `saturated_two_heteroatom_1_4_ring_
    name`; `n_locant` is 4 for an N,O/N,S pair (the non-nitrogen heteroatom
    always wins locant 1) or 1 for the symmetric N,N piperazine pair (the
    substituted nitrogen wins the lowest locant)."""
    ring_info = mol.GetRingInfo()
    candidates = []
    for ring in ring_info.AtomRings():
        if len(ring) != _TWO_HETERO_RING_SIZE:
            continue
        if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            continue
        heteroatoms = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
        if len(heteroatoms) != 2:
            continue
        elements = [_TWO_HETERO_ELEMENTS.get(mol.GetAtomWithIdx(a).GetAtomicNum()) for a in heteroatoms]
        if None in elements or frozenset(elements) not in _TWO_HETERO_ELEMENT_PAIRS:
            continue
        candidates.append((set(ring), heteroatoms))
    if len(candidates) != 1:
        return None
    ring_atoms, heteroatoms = candidates[0]

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    het1, het2 = heteroatoms
    if abs(ring_order.index(het1) - ring_order.index(het2)) != _TWO_HETERO_RING_SIZE // 2:
        return None

    nitrogens = [a for a in heteroatoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    substituted = []
    for n in nitrogens:
        n_atom = mol.GetAtomWithIdx(n)
        if n_atom.GetFormalCharge() != 0 or n_atom.GetIsotope() != 0:
            return None
        exo = [x for x in n_atom.GetNeighbors() if x.GetIdx() not in ring_atoms]
        if len(exo) > 1:
            return None
        if exo:
            if mol.GetBondBetweenAtoms(n, exo[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
                return None
            substituted.append((n, exo[0].GetIdx()))
        elif n_atom.GetTotalNumHs() != 1:
            return None
    if len(substituted) != 1:
        return None
    (n_idx, root) = substituted[0]

    other_heteroatom = het2 if n_idx == het1 else het1
    other_atom = mol.GetAtomWithIdx(other_heteroatom)
    if other_atom.GetAtomicNum() != 7:
        if other_atom.GetFormalCharge() != 0 or other_atom.GetIsotope() != 0:
            return None
        if any(x.GetIdx() not in ring_atoms for x in other_atom.GetNeighbors()):
            return None

    for a in ring_atoms:
        if a in heteroatoms:
            continue
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if any(x.GetIdx() not in ring_atoms for x in atom.GetNeighbors()):
            return None
        if atom.GetTotalNumHs() != 2:
            return None
    if any(a in ring_atoms and b in ring_atoms for a, b, _ in non_single_bonds(mol)):
        return None

    other_element = _TWO_HETERO_ELEMENTS[other_atom.GetAtomicNum()]
    n_locant = 1 if other_element == "N" else 4
    elements = frozenset(("N", other_element))
    return n_idx, ring_atoms, root, n_locant, elements


def has_ring_amine_shape(mol) -> bool:
    return _ring_amine_shape(mol) is not None or _two_hetero_ring_amine_shape(mol) is not None


def has_ring_amine_sulfonyl_shape(mol) -> bool:
    """True if `mol` fits `_ring_amine_shape`/`_two_hetero_ring_amine_shape`
    *and* the ring nitrogen's one substituent is specifically a sulfonyl
    group. Exposed separately (narrower than `has_ring_amine_shape`) so
    `core.py` can route it ahead of `has_sulfonamide_shape` -- which would
    otherwise misclaim this exact shape (a sulfonyl group on a plain ring
    nitrogen looks like an ordinary, but ring-unaware, sulfonamide to that
    check) -- without also preempting `_hidden_amide_ketone.py`'s
    acyl-on-ring-nitrogen shape, which `has_ring_amine_shape` alone can't
    distinguish from any other single N-substituent (it only checks that
    exactly one exists, not what it is)."""
    shape = _ring_amine_shape(mol)
    if shape is not None:
        n_idx, _, root = shape
        return _sulfonyl_root_shape(mol, root, n_idx) is not None
    shape2 = _two_hetero_ring_amine_shape(mol)
    if shape2 is not None:
        n_idx, _, root, _, _ = shape2
        return _sulfonyl_root_shape(mol, root, n_idx) is not None
    return False


def _sulfonyl_root_shape(mol, root, n_idx):
    """(r_carbon_idx, {oxygen_idx, oxygen_idx}) if `root` is a sulfonyl
    sulfur (P-65.3.1's -SO2- group: two double-bonded, monovalent oxygens
    plus exactly one carbon neighbor besides the ring nitrogen) -- else
    None."""
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() != 16 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return None
    if mol.GetBondBetweenAtoms(root, n_idx).GetBondTypeAsDouble() != 1.0:
        return None
    neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() != n_idx]
    oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
    carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
    if len(oxygens) != 2 or len(carbons) != 1:
        return None
    if any(
        o.GetDegree() != 1 or mol.GetBondBetweenAtoms(root, o.GetIdx()).GetBondTypeAsDouble() != 2.0
        for o in oxygens
    ):
        return None
    (r_carbon,) = carbons
    if mol.GetBondBetweenAtoms(root, r_carbon.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    return r_carbon.GetIdx(), {o.GetIdx() for o in oxygens}


def _validate_and_name_substituent(mol, ring_atoms, n_idx, root):
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    sulfonyl = _sulfonyl_root_shape(mol, root, n_idx)
    extra_allowed = {root, *sulfonyl[1]} if sulfonyl is not None else set()

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms or idx in extra_allowed:
            continue
        if atom.GetAtomicNum() not in _ALLOWED_SUBSTITUENT_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "an N-substituent containing anything other than carbon, "
                "halogens, and a sulfonyl group is not supported yet for "
                "this ring-amine path"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic():
            raise UnsupportedStructure(
                "an aromatic N-substituent is not supported yet for this "
                "ring-amine path"
            )
        if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1)"
            )

    if any(
        a not in ring_atoms and b not in ring_atoms and not ({a, b} <= extra_allowed)
        for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure("unsaturation in the N-substituent is not supported yet")

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    if sulfonyl is not None:
        r_carbon, _ = sulfonyl
        r_name, r_is_compound = name_branch(graph, r_carbon, root, halogens)
        return f"{r_name}sulfonyl", r_is_compound
    return name_branch(graph, root, n_idx, halogens)


def name_ring_amine(mol) -> str:
    shape = _ring_amine_shape(mol)
    if shape is not None:
        n_idx, ring_atoms, root = shape
        name, is_compound = _validate_and_name_substituent(mol, ring_atoms, n_idx, root)
        stem = saturated_ring_name("N", len(ring_atoms))
        if stem is None:
            raise UnsupportedStructure(
                f"no retained/Hantzsch-Widman name for a {len(ring_atoms)}-membered "
                f"N-heteroatom saturated ring (P-22.2.1)"
            )
        sub_name = f"({name})" if is_compound else name
        return f"1-{sub_name}{stem}"

    shape2 = _two_hetero_ring_amine_shape(mol)
    if shape2 is not None:
        n_idx, ring_atoms, root, n_locant, elements = shape2
        name, is_compound = _validate_and_name_substituent(mol, ring_atoms, n_idx, root)
        stem = saturated_two_heteroatom_1_4_ring_name(elements)
        if stem is None:
            raise UnsupportedStructure(
                "no retained name for this 1,4-two-heteroatom saturated "
                "ring element pair (P-22.2.1)"
            )
        sub_name = f"({name})" if is_compound else name
        return f"{n_locant}-{sub_name}{stem}"

    raise UnsupportedStructure(
        "no plain saturated monocyclic ring with a single N-substituted "
        "nitrogen shape found"
    )
