"""Naming of the minimal single-atom-bridged fused-aromatic case
(1,4-dihydro-1,4-methanonaphthalene "benzonorbornadiene" and its 'epoxy'
oxygen-bridged analogue), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-25.4 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): a
  one-atom bridge across two nonadjacent ('1,4'-type) positions of a
  mancude ring system is named by citing the bridge prefix and its own
  attachment locants (here, always the naphthalene locants unchanged --
  the bridge atom's own locant, '9', doesn't appear in this module's
  unsubstituted scope) after the ring system's name.
- P-25.4.2.1.4 (as of `tasks/bridged-naphthalene-oxa-bridge-naming.md`,
  2026-08-25): the preselected bridge prefix for a divalent -O- bridge is
  'epoxy', confirmed directly in the primary text ("epoxy (preselected
  prefix) -O- (not epoxidano)") -- not a guess. A single-carbon -CH2-
  bridge is 'methano' (P-25.4.2.1.1), already supported below.
- Structural necessity for 'dihydro' (P-31.1.4.2): bridging naphthalene's
  1,4-positions forces both bridgehead carbons to sp3 (a bridged atom can't
  stay part of a mancude/aromatic ring), so the correct name is
  '1,4-dihydro-1,4-methanonaphthalene'/'1,4-dihydro-1,4-epoxynaphthalene',
  not the bare bridge name -- confirmed against CAS 4453-90-1 (methano,
  `tasks/bridged-naphthalene-naming.md`) and CAS 573-57-9 (epoxy, this
  task). Both locant sets ('...-dihydro-' and the bridge's own) are always
  numerically identical for this exact shape (the bridge only ever spans
  the two bridgehead positions the hydro prefix already names), so this
  module computes the pair once and reuses it for both.
- P-25.3.3.1.1: the underlying naphthalene numbering is reused unchanged
  from `_aromatic.py`'s naphthalene machinery (`_periphery_cycle`/
  `_straight_chain_candidates`), applied to the 6+6 ring skeleton with the
  bridge atom excluded -- exactly as `_dihydro_aromatic.py` does, since the
  bridge doesn't change which atoms are "1" through "8a".
- Deliberately not built on RDKit's own SSSR: for this bridged shape RDKit
  picks two 5-membered rings (via the bridge) instead of naphthalene's
  natural two 6-membered rings (confirmed by direct inspection), so this
  module identifies the two "naphthalene-shape" 6-rings itself from atom
  roles (bridgeheads, bridge, fusion atoms, intact-ring atoms) rather than
  trusting `mol.GetRingInfo()`'s ring choice. PubChem's own computed
  "IUPACName" for both the methano and epoxy compounds is a von Baeyer
  bridged-ring name (e.g. '11-oxatricyclo[6.2.1.0^2,7]undeca-2,4,6,9-
  tetraene' for the epoxy case), not a fusion+bridge name, so it isn't
  usable for verifying either name this module returns -- CAS/NIST WebBook/
  reagent-catalog naming (which actually use fusion+bridge nomenclature for
  this shape) is the cross-check instead.

Explicitly out of scope (raise `UnsupportedStructure` via the generic
fallback in `core.py`, since `find_bridged_naphthalene_core` below simply
returns None for any of these):
- More than one bridge, a bridge longer/shorter than one atom, or any
  bridge atom other than carbon ('methano') or oxygen ('epoxy') -- e.g.
  'imino' (-NH-), 'etheno' (-CH=CH-), 'epidioxy' (-O-O-, two atoms).
- Any fused ring system other than plain 2-ring naphthalene.
- A bridge spanning adjacent ('1,2'-type) ring positions (a structurally
  different, cyclopropa-fused system, not a P-25.4 bridge at all).
- Substituents of any kind, including on the bridge atom itself (which
  would need its own locant, '9', not used in this unsubstituted scope).
"""

from ._aromatic import _straight_chain_candidates

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


def find_bridged_naphthalene_core(mol):
    """Return (ring_atom_sets, fusion_bond_idx, bridgeheads, bridge_prefix)
    if `mol` is naphthalene's carbon skeleton plus exactly one -CH2- or -O-
    bridge across one ring's 1,4-type positions (see module docstring),
    else None."""
    if mol.GetNumAtoms() != 11:
        return None

    bridge_candidates = _bridge_candidates(mol)
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

    intact_ring_atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        idx_set = set(ring)
        if len(ring) == 6 and fusion_atoms[0].GetIdx() in idx_set and fusion_atoms[1].GetIdx() in idx_set:
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                intact_ring_atoms = idx_set
                break
    if len(intact_ring_atoms) != 6:
        return None

    reduced_ring_atoms = {
        bridgeheads[0].GetIdx(),
        bridgeheads[1].GetIdx(),
        ene_atoms[0].GetIdx(),
        ene_atoms[1].GetIdx(),
        fusion_atoms[0].GetIdx(),
        fusion_atoms[1].GetIdx(),
    }
    ring_atom_sets = [reduced_ring_atoms, intact_ring_atoms]
    fusion_bond_idxs = {fusion_bond.GetIdx()}
    bridgehead_idxs = (bridgeheads[0].GetIdx(), bridgeheads[1].GetIdx())
    return ring_atom_sets, fusion_bond_idxs, bridgehead_idxs, bridge_prefix


def name_bridged_naphthalene(mol, core) -> str:
    ring_atom_sets, fusion_bond_idxs, bridgeheads, bridge_prefix = core
    atom_rings = [tuple(s) for s in ring_atom_sets]

    candidates = _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, [0, 1])
    best_locants = None
    for locants in candidates:
        pair = tuple(sorted(locants[a] for a in bridgeheads))
        if best_locants is None or pair < best_locants:
            best_locants = pair

    a, b = best_locants
    return f"{a},{b}-dihydro-{a},{b}-{bridge_prefix}naphthalene"
