"""Naming of the single-atom-bridged fused-aromatic case (1,4-type methano/
epoxy bridges on a terminal ring of any plain all-carbon ortho-fused
mancude polycyclic aromatic parent this project can already name --
naphthalene ("benzonorbornadiene"), phenanthrene, tetracene, etc. -- plus
the anthracene analogue bridging the 9,10 meso positions,
"9,10-dihydro-9,10-methanoanthracene"), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-25.4 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): a
  one-atom bridge across two nonadjacent ('1,4'-type) positions of a
  mancude ring system is named by citing the bridge prefix and its own
  attachment locants after the ring system's name.
- P-25.4.2.1.4: the preselected bridge prefix for a divalent -O- bridge is
  'epoxy', confirmed directly in the primary text ("epoxy (preselected
  prefix) -O- (not epoxidano)") -- not a guess. A single-carbon -CH2-
  bridge is 'methano' (P-25.4.2.1.1), already supported below.
- Structural necessity for 'dihydro' (P-31.1.4.2): bridging naphthalene's
  1,4-positions forces both bridgehead carbons to sp3 (a bridged atom can't
  stay part of a mancude/aromatic ring), so the correct name is
  '1,4-dihydro-1,4-methanonaphthalene'/'1,4-dihydro-1,4-epoxynaphthalene',
  not the bare bridge name -- confirmed against CAS 4453-90-1 (methano)
  and CAS 573-57-9 (epoxy). Both locant sets ('...-dihydro-' and the bridge's own) are always
  numerically identical for this exact shape (the bridge only ever spans
  the two bridgehead positions the hydro prefix already names), so this
  module computes the pair once and reuses it for both.
- P-25.3.3.1.1: the numbering of the ring system containing the bridge
  reuses `_aromatic.py`'s own shape classifier and candidate-numbering
  machinery (`_classify_shape`/`_retained_chain_name` picks the parent
  name; `_straight_chain_candidates`/`_anthracene_candidates`/
  `_phenanthrene_candidates` generate its numbering candidates)
  unchanged, applied to the bridged ring system with the bridge atom
  excluded -- exactly as `_dihydro_aromatic.py` does for its own
  (unbridged) partially-saturated case, since the bridge doesn't change
  which atoms are "1" through "8a" (or however many positions the parent
  has).
- Deliberately not built on RDKit's own SSSR for the bridged ring itself:
  for this bridged shape RDKit picks two 5-membered rings (via the bridge)
  instead of the bridged ring's natural 6-membered shape (confirmed by
  direct inspection), so this module identifies that one ring itself from
  atom roles (bridgeheads, bridge, fusion atoms, ene atoms) rather than
  trusting `mol.GetRingInfo()`'s ring choice there, then combines it with
  whichever *other* rings RDKit's SSSR reports correctly (every other ring
  in the parent is untouched by the bridge and perceived normally).
  PubChem's own computed "IUPACName" for both the methano and epoxy
  compounds is a von Baeyer bridged-ring name (e.g.
  '11-oxatricyclo[6.2.1.0^2,7]undeca-2,4,6,9-tetraene' for the naphthalene
  epoxy case), not a fusion+bridge name, so it isn't usable for verifying
  either name this module returns -- CAS/NIST WebBook/reagent-catalog
  naming (which actually use fusion+bridge nomenclature for this shape) is
  the cross-check instead.
- Anthracene 9,10-bridge: a one-atom bridge across anthracene's own meso positions
  (C9/C10, the middle ring's two non-fusion atoms) is structurally the
  same P-25.4 shape, just on a 3-ring parent instead of 2-ring
  naphthalene. Unlike the naphthalene case, no locant search/tie-break is
  actually needed here: `_aromatic.py`'s `_anthracene_candidates` always
  assigns anthracene's meso positions locants 9 and 10 in every candidate
  (that's what "meso" numbering means), so the bridgeheads are always
  '9,10' regardless of which of the 4 symmetric candidates is picked --
  the search below is kept anyway only for structural symmetry with the
  naphthalene function and as a defensive check. Confirmed against the
  literature name "9,10-Dihydro-9,10-methanoanthracene" (J. Org. Chem.)
  and PubChem CID 12651785 (structure only, cross-checked by
  ConnectivitySMILES; PubChem's own computed IUPACName is von Baeyer-style
  and unusable for verifying the name itself, same limitation as the
  naphthalene case above).
- Halogen substituents anywhere on the aromatic portion (not the bridge
  atom or its two bridgehead/ene neighbors): `find_bridged_aromatic_core`/
  `name_bridged_aromatic` already compute their own ring atom sets and
  locants directly (unlike a hardcoded whole-molecule-match module), so no
  new numbering mechanism was needed -- a halogen substituent is just
  another `graph[atom]` neighbor outside the ring-atom set, exactly as
  `_aromatic.py`'s plain fused-ring module already handles it. Confirmed
  against a from-scratch-built chlorobenzonorbornadiene (structure only;
  PubChem CID 12473502's own computed IUPACName is von-Baeyer-style, same
  limitation as the unsubstituted cases above) and its epoxy analogue
  (CID 14208771).

Explicitly out of scope (raise `UnsupportedStructure` via the generic
fallback in `core.py`, since `find_bridged_aromatic_core`/
`find_bridged_anthracene_core` below simply return None for any of these):
- More than one bridge, a bridge longer/shorter than one atom, or any
  bridge atom other than carbon ('methano') or oxygen ('epoxy') -- e.g.
  'imino' (-NH-), 'etheno' (-CH=CH-), 'epidioxy' (-O-O-, two atoms). On
  anthracene specifically, only 'methano' is supported (see
  `find_bridged_anthracene_core`'s own docstring for why 'epoxy' is left
  for a follow-up task).
- Any fused-ring chain shape `_aromatic.py` itself doesn't have a retained
  name for yet (a branched ring-fusion, a longer angular/polyphene chain,
  a peri-fused system such as pyrene, etc. -- see
  `find_aromatic_fused_core`/`_retained_chain_name`), or a bridge
  positioned anywhere on anthracene other than the 9,10 meso positions
  (e.g. a 1,4-type bridge within one terminal ring, which this module
  doesn't attempt to distinguish from the same 1,4-type bridge on any
  other parent's terminal ring -- handled by `find_bridged_aromatic_core`
  instead, not `find_bridged_anthracene_core`).
- A bridge spanning adjacent ('1,2'-type) ring positions (a structurally
  different, cyclopropa-fused system, not a P-25.4 bridge at all).
- Any substituent other than a halogen anywhere on the aromatic portion
  (see above); any substituent at all on the bridge atom itself, the
  reduced-ring bridgeheads/alkene carbons, or anywhere on the anthracene
  case (still fully unsubstituted-only).
"""

from ._aromatic import (
    _anthracene_candidates,
    _classify_shape,
    _phenanthrene_candidates,
    _retained_chain_name,
    _ring_adjacency,
    _ring_path_order,
    _straight_chain_candidates,
)
from ._common import UnsupportedStructure, adjacency, group_substituents, halogen_substituents, lowest_locant_set, ring_cycle
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_BRIDGE_PREFIXES = {6: "methano", 8: "epoxy"}


def _bridge_candidates(mol):
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() or atom.GetDegree() != 2:
            continue
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6 and atom.GetTotalNumHs() == 2 and atom.GetHybridization().name == "SP3":
            candidates.append(atom)
        elif atomic_num == 8 and atom.GetTotalNumHs() == 0:
            candidates.append(atom)
    return candidates


def find_bridged_aromatic_core(mol):
    """Return (atom_rings, ring_atom_sets, fusion_bond_idxs, bridgeheads,
    bridge_prefix, bridge_idx) if `mol` is a plain all-carbon ortho-fused
    mancude polycyclic aromatic parent -- any parent `_aromatic.py`'s
    `find_aromatic_fused_core`/`name_aromatic_fused` can already name
    (naphthalene, anthracene, phenanthrene, tetracene..nonacene) -- plus
    exactly one -CH2- or -O- bridge across one ring's 1,4-type positions
    (see module docstring), optionally with halogen substituents on the
    rest of the aromatic system, else None.

    RDKit's own ring perception can't be trusted for the bridged ring
    itself (see module docstring: it finds two 5-membered artifact rings
    running through the bridge atom instead of the bridged ring's natural
    6-membered shape), so that one ring is reconstructed manually from the
    bridgehead/fusion/ene atom roles below and combined with whichever
    *other* rings RDKit's SSSR reports correctly (any 6-membered
    all-aromatic ring elsewhere in the molecule, untouched by the bridge)
    before handing the combined ring list to `_aromatic.py`'s shape
    classifier and numbering machinery. Those are purely
    connectivity-based (ring-membership sets and fusion-bond indices, not
    RDKit's own aromaticity perception), so they work unchanged on this
    hybrid ring list -- only the *finding* of the bridged ring itself ever
    needed naphthalene-specific code; the numbering was always
    parent-agnostic."""
    halogens = halogen_substituents(mol)
    for idx in halogens:
        if mol.GetAtomWithIdx(idx).GetDegree() != 1:
            return None

    bridge_candidates = _bridge_candidates(mol)
    if len(bridge_candidates) != 1:
        return None
    bridge_atom = bridge_candidates[0]
    if bridge_atom.GetFormalCharge() != 0 or bridge_atom.GetIsotope() != 0:
        return None
    bridge_prefix = _BRIDGE_PREFIXES[bridge_atom.GetAtomicNum()]

    for atom in mol.GetAtoms():
        if atom.GetIdx() == bridge_atom.GetIdx() or atom.GetIdx() in halogens:
            continue
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None

    bridgeheads = list(bridge_atom.GetNeighbors())
    if len(bridgeheads) != 2:
        return None
    for bh in bridgeheads:
        if bh.GetIsAromatic() or bh.GetDegree() != 3 or bh.GetTotalNumHs() != 1:
            return None
        if bh.GetHybridization().name != "SP3":
            return None

    fusion_atoms, ene_atoms = [], []
    for bh in bridgeheads:
        others = [n for n in bh.GetNeighbors() if n.GetIdx() != bridge_atom.GetIdx()]
        if len(others) != 2:
            return None
        found_fusion = [n for n in others if n.GetIsAromatic() and n.GetDegree() == 3 and n.GetTotalNumHs() == 0]
        found_ene = [
            n
            for n in others
            if not n.GetIsAromatic() and n.GetDegree() == 2 and n.GetTotalNumHs() == 1
        ]
        if len(found_fusion) != 1 or len(found_ene) != 1:
            return None
        fusion_atoms.append(found_fusion[0])
        ene_atoms.append(found_ene[0])

    if fusion_atoms[0].GetIdx() == fusion_atoms[1].GetIdx() or ene_atoms[0].GetIdx() == ene_atoms[1].GetIdx():
        return None
    ene_bond = mol.GetBondBetweenAtoms(ene_atoms[0].GetIdx(), ene_atoms[1].GetIdx())
    if ene_bond is None or ene_bond.GetIsAromatic() or ene_bond.GetBondTypeAsDouble() != 2.0:
        return None

    fusion_bond = mol.GetBondBetweenAtoms(fusion_atoms[0].GetIdx(), fusion_atoms[1].GetIdx())
    if fusion_bond is None or not fusion_bond.GetIsAromatic():
        return None

    bridged_ring_atoms = {
        bridgeheads[0].GetIdx(),
        bridgeheads[1].GetIdx(),
        ene_atoms[0].GetIdx(),
        ene_atoms[1].GetIdx(),
        fusion_atoms[0].GetIdx(),
        fusion_atoms[1].GetIdx(),
    }
    other_rings = [
        set(ring)
        for ring in mol.GetRingInfo().AtomRings()
        if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
    ]
    ring_atom_sets = other_rings + [bridged_ring_atoms]

    full_ring_atoms = set()
    for s in ring_atom_sets:
        full_ring_atoms |= s
    other_atoms = {
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() != bridge_atom.GetIdx() and atom.GetIdx() not in halogens
    }
    if other_atoms - full_ring_atoms:
        # A substituent other than a halogen (see module docstring) --
        # out of scope, not a different aromatic parent this function
        # should try to recognize.
        return None

    graph = adjacency(mol)
    try:
        atom_rings = [tuple(ring_cycle(graph, list(s))) for s in ring_atom_sets]
    except StopIteration:
        return None

    bond_ring_count = {}
    for atom_set in ring_atom_sets:
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in atom_set and b in atom_set:
                bond_ring_count[bond.GetIdx()] = bond_ring_count.get(bond.GetIdx(), 0) + 1
    fusion_bond_idxs = {b for b, c in bond_ring_count.items() if c >= 2}

    bridgehead_idxs = (bridgeheads[0].GetIdx(), bridgeheads[1].GetIdx())
    return atom_rings, ring_atom_sets, fusion_bond_idxs, bridgehead_idxs, bridge_prefix, bridge_atom.GetIdx()


def name_bridged_aromatic(mol, core) -> str:
    atom_rings, ring_atom_sets, fusion_bond_idxs, bridgeheads, bridge_prefix, bridge_idx = core
    n = len(atom_rings)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    full_ring_atoms = set()
    for s in ring_atom_sets:
        full_ring_atoms |= s
    excluded = full_ring_atoms | {bridge_idx}

    adj, fusion_bonds_by_pair = _ring_adjacency(atom_rings, ring_atom_sets, fusion_bond_idxs, mol)
    ring_order = _ring_path_order(adj, n)
    shape = _classify_shape(graph, atom_rings, ring_order, fusion_bonds_by_pair)
    parent = _retained_chain_name(n, shape)
    if parent is None:
        raise UnsupportedStructure(
            "this bridged aromatic ring chain does not have one of the "
            "retained names supported so far (naphthalene, anthracene, "
            "phenanthrene, or a straight tetracene..nonacene polyacene); "
            "see P-25.3.1-P-25.3.3"
        )

    if parent == "anthracene":
        candidates = _anthracene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order)
    elif parent == "phenanthrene":
        candidates = _phenanthrene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, ring_order)
    else:
        candidates = _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)

    best_key = None
    for locants in candidates:
        pair = tuple(sorted(locants[a] for a in bridgeheads))
        substituents = {}
        for atom, position in locants.items():
            branch_roots = [n for n in graph[atom] if n not in excluded]
            if branch_roots:
                substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
        grouped = group_substituents(substituents)
        locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
        citation_locants = tuple(
            loc
            for name in sorted(grouped, key=alpha_sort_key)
            for loc in sorted(grouped[name]["locants"])
        )
        a, b = pair
        base = f"{a},{b}-dihydro-{a},{b}-{bridge_prefix}{parent}"
        name = base if not grouped else format_substituent_prefixes(grouped) + "-" + base
        key = (pair, locant_set, citation_locants, name)
        if best_key is None or key < best_key:
            best_key = key

    return best_key[-1]


def find_bridged_anthracene_core(mol):
    """Return (ring_atom_sets, fusion_bonds_by_pair, bridgeheads,
    bridge_prefix) if `mol` is anthracene's carbon skeleton plus exactly
    one -CH2- bridge ('methano') across the middle ring's two meso
    (9,10-type) positions, else None. Unlike naphthalene's bridgeheads
    (which sit between one fusion atom and one plain -CH= atom),
    anthracene's meso bridgeheads sit between two fusion atoms, one shared
    with each terminal ring -- so this function's atom-role matching is
    deliberately different from `find_bridged_aromatic_core` above, not
    a parameterized reuse of it.

    Deliberately methano-only (unlike `find_bridged_aromatic_core`, which
    also accepts an 'epoxy' -O- bridge): no worked example or independently
    verifiable name for the anthracene 9,10-epoxy analogue was found --
    only the structure was confirmed (PubChem, matching connectivity),
    which isn't enough to assert the name under this project's
    test-writing policy. Left for a follow-up task if a real worked
    example turns up."""
    if mol.GetNumAtoms() != 15:
        return None

    bridge_candidates = [a for a in _bridge_candidates(mol) if a.GetAtomicNum() == 6]
    if len(bridge_candidates) != 1:
        return None
    bridge_atom = bridge_candidates[0]
    if bridge_atom.GetFormalCharge() != 0 or bridge_atom.GetIsotope() != 0:
        return None
    bridge_prefix = _BRIDGE_PREFIXES[bridge_atom.GetAtomicNum()]

    for atom in mol.GetAtoms():
        if atom.GetIdx() == bridge_atom.GetIdx():
            continue
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None

    bridgeheads = list(bridge_atom.GetNeighbors())
    if len(bridgeheads) != 2:
        return None
    for bh in bridgeheads:
        if bh.GetIsAromatic() or bh.GetDegree() != 3 or bh.GetTotalNumHs() != 1:
            return None
        if bh.GetHybridization().name != "SP3":
            return None

    neighbor_sets = []
    for bh in bridgeheads:
        others = [n for n in bh.GetNeighbors() if n.GetIdx() != bridge_atom.GetIdx()]
        if len(others) != 2:
            return None
        for n in others:
            if not n.GetIsAromatic() or n.GetDegree() != 3 or n.GetTotalNumHs() != 0:
                return None
        neighbor_sets.append({n.GetIdx() for n in others})

    fusion_pairs = [(a, b) for a in neighbor_sets[0] for b in neighbor_sets[1] if mol.GetBondBetweenAtoms(a, b)]
    if len(fusion_pairs) != 2:
        return None
    used = {atom for pair in fusion_pairs for atom in pair}
    if len(used) != 4:
        return None

    terminal_rings = []
    for a, b in fusion_pairs:
        bond = mol.GetBondBetweenAtoms(a, b)
        if not bond.GetIsAromatic():
            return None
        found = None
        for ring in mol.GetRingInfo().AtomRings():
            idx_set = set(ring)
            if len(ring) == 6 and a in idx_set and b in idx_set and all(
                mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring
            ):
                found = idx_set
                break
        if found is None:
            return None
        terminal_rings.append(found)
    if terminal_rings[0] == terminal_rings[1]:
        return None

    middle_ring_atoms = {bridgeheads[0].GetIdx(), bridgeheads[1].GetIdx()} | used
    ring_atom_sets = [terminal_rings[0], middle_ring_atoms, terminal_rings[1]]
    fusion_bonds_by_pair = {
        frozenset((0, 1)): fusion_pairs[0],
        frozenset((1, 2)): fusion_pairs[1],
    }
    bridgehead_idxs = (bridgeheads[0].GetIdx(), bridgeheads[1].GetIdx())
    return ring_atom_sets, fusion_bonds_by_pair, bridgehead_idxs, bridge_prefix


def name_bridged_anthracene(mol, core) -> str:
    ring_atom_sets, fusion_bonds_by_pair, bridgeheads, bridge_prefix = core
    graph = adjacency(mol)

    candidates = _anthracene_candidates(graph, ring_atom_sets, fusion_bonds_by_pair, [0, 1, 2])
    best_locants = None
    for locants in candidates:
        pair = tuple(sorted(locants[a] for a in bridgeheads))
        if best_locants is None or pair < best_locants:
            best_locants = pair

    a, b = best_locants
    return f"{a},{b}-dihydro-{a},{b}-{bridge_prefix}anthracene"
