"""Naming of the minimal partially-hydrogenated fused-aromatic case
(dihydronaphthalene), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-31.1.4.2 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf):
  the nondetachable 'hydro' prefix ('dihydro-', 'tetrahydro-', etc.) reduces
  the number of noncumulative double bonds of a mancude-ring-system parent
  hydride (here, naphthalene, P-25.1.2.1) by saturating a pair of ring
  atoms; lowest locants go to the hydro prefix's own locant set first (this
  module's whole scope, since there is nothing else -- no indicated
  hydrogen, no substituents -- to compete with it here). This is a distinct
  mechanism from 'indicated hydrogen' (an italic "1H"-style locant
  disambiguating a tautomer a mancude ring system can already have on its
  own, e.g. `_heteroaromatic_fused.py`'s "1H-indole") and from 'added
  hydrogen' (a suffix-triggered parenthetical locant, e.g.
  "naphthalen-1(2H)-one"): naphthalene has no indicated-hydrogen tautomer
  of its own, and no suffix is involved here, so 'hydro' is the only
  mechanism in play (confirmed against secondary sources cross-referencing
  FR-9.3's indicated/added/hydro-prefix distinction).
- P-25.3.3.1.1: dihydronaphthalene reuses naphthalene's own peripheral
  numbering (1,2,3,4,4a,5,6,7,8,8a) unchanged -- the saturation doesn't
  renumber anything, it just adds a locant pair citation. Numbering
  candidates (which physical ring is "1", which walk direction) are
  reused unchanged from `_aromatic.py`'s naphthalene machinery
  (`_periphery_cycle`/`_straight_chain_candidates`), which is purely
  connectivity-based and doesn't care that one ring is no longer aromatic.
- Structural shape: naphthalene's two fused benzo rings share exactly one
  bond (two atoms, always aromatic since they stay part of the untouched
  ring). In the reduced ring, saturating any two *adjacent* skeletal
  positions leaves a well-formed structure (the remaining two nonfusion
  atoms keep a single, isolated C=C between them) -- giving either
  '1,2-dihydronaphthalene' (CAS 447-53-0) or '1,4-dihydronaphthalene'
  (CAS 612-17-9) depending on which adjacent pair, both real, independently
  named compounds (NIST WebBook, Sigma-Aldrich). Saturating a
  *non-adjacent* pair (e.g. the two "middle" nonfusion positions) is not
  reachable this way: the remaining nonfusion atom would need a double
  bond to a fusion atom, which would break that fusion atom's own
  aromaticity and therefore the untouched ring's aromaticity too -- so it
  never arises for a molecule that (by this module's own detection) already
  has one fully intact aromatic ring. This module doesn't special-case
  '1,2-' vs '1,4-': it structurally detects "exactly two sp3 ring carbons,
  the other two nonfusion ring carbons doubly bonded to each other", which
  covers both shapes uniformly and rejects everything else (including any
  shape requiring a fusion atom to lose aromaticity) by construction.

Explicitly out of scope (raise `UnsupportedStructure` via the generic
P-23/P-24/P-25 fallback in `core.py`, since `find_dihydronaphthalene_core`
below simply returns None for any of these):
- More than one 'hydro' pair short of full saturation (tetrahydro-/
  hexahydro-/octahydronaphthalene -- a genuinely bigger generalization,
  needing locant computation and lowest-locant tie-breaking this module's
  two-function pair doesn't attempt).
- Any fused ring system other than the plain 2-ring naphthalene shape
  (anthracene, indole, heteroaromatics, ...).
- Any indicated-hydrogen-requiring parent (not applicable to naphthalene
  itself, but future modules built on this one will need to handle it).
- Substituents of any kind, including halogens (P-35.2.1) -- every ring
  atom must have exactly its "bare" degree (2 for a CH/CH2 position, 3 for
  a fusion carbon), so any substituent is rejected by construction.
- Bridged derivatives (methanonaphthalene, etc.), a separate, larger
  prerequisite task.

- **Full saturation (`decahydronaphthalene`)**: every ring double bond
  removed, not just one pair. Unlike the partial-hydro case above, this
  needs no locant computation at all -- P-31.2.3.3.2's own worked example
  cites total hydrogenation with no locants (P-14.3.4.5: a fully-cited
  hydro count is redundant once every ring position is saturated), so
  `find_decahydronaphthalene_core`/`name_decahydronaphthalene` below just
  detect the bare skeleton and return the literal string unconditionally.
  This landed before the general partial-hydrogenation case above because
  it's structurally simpler, not because it's a subset of it -- confirmed
  via PubChem PUG REST (CID 7044: `1,2,3,4,4a,5,6,7,8,8a-
  decahydronaphthalene`, though the Blue Book's own plain `decahydro-
  naphthalene (PIN)` omits the redundant locants PubChem always cites).
  Reuses `_bicyclic.py`'s `find_bicyclic_core` for the actual skeleton
  detection (a 0-bridge 6,6 bicyclic is exactly this shape) rather than
  reimplementing ring-walking here; this is deliberately routed ahead of
  `_bicyclic.py`'s own von Baeyer naming in `core.py`; a von Baeyer name
  for this skeleton has been wrong all along (P-23's own naming method is
  reserved for skeletons with no competing mancude-ring-system name).
"""

from rdkit import Chem

from ._aromatic import _periphery_cycle, _straight_chain_candidates
from ._bicyclic import find_bicyclic_core
from ._common import multiplied_word


def find_dihydronaphthalene_core(mol):
    """Return (ring_atom_sets, fusion_bond_idxs, ring_order, sp3_atoms) if
    `mol` is naphthalene's carbon skeleton with exactly one adjacent pair of
    ring atoms saturated (see module docstring), else None.
    """
    if mol.GetNumAtoms() != 10:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2 or any(len(r) != 6 for r in atom_rings):
        return None
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None

    bond_ring_count = {}
    for bonds in ring_info.BondRings():
        for b in bonds:
            bond_ring_count[b] = bond_ring_count.get(b, 0) + 1
    fusion_bond_idxs = {b for b, c in bond_ring_count.items() if c >= 2}
    if len(fusion_bond_idxs) != 1:
        return None

    ring_atom_sets = [set(r) for r in atom_rings]
    fusion_atoms = ring_atom_sets[0] & ring_atom_sets[1]
    if len(fusion_atoms) != 2:
        return None
    for a in fusion_atoms:
        atom = mol.GetAtomWithIdx(a)
        if not atom.GetIsAromatic() or atom.GetDegree() != 3:
            return None

    # The reduced ring's only aromatic atoms are the two it shares with the
    # intact ring (a fusion atom is aromatic whenever it belongs to any
    # aromatic ring, regardless of the other ring's own status), so it has
    # an aromatic count of 2 (not 4); the intact ring has all 6.
    aromatic_counts = [
        sum(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring) for ring in atom_rings
    ]
    if sorted(aromatic_counts) != [2, 6]:
        return None
    reduced_idx = aromatic_counts.index(2)
    intact_idx = 1 - reduced_idx

    intact_non_fusion = ring_atom_sets[intact_idx] - fusion_atoms
    for a in intact_non_fusion:
        atom = mol.GetAtomWithIdx(a)
        if not atom.GetIsAromatic() or atom.GetDegree() != 2:
            return None

    reduced_non_fusion = ring_atom_sets[reduced_idx] - fusion_atoms
    if len(reduced_non_fusion) != 4:
        return None

    sp3_atoms, ene_atoms = [], []
    for a in reduced_non_fusion:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetIsAromatic() or atom.GetDegree() != 2:
            return None
        if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetTotalNumHs() == 2:
            sp3_atoms.append(a)
        else:
            ene_atoms.append(a)
    if len(sp3_atoms) != 2 or len(ene_atoms) != 2:
        return None

    ene_bond = mol.GetBondBetweenAtoms(*ene_atoms)
    if ene_bond is None or ene_bond.GetBondTypeAsDouble() != 2.0:
        return None

    return ring_atom_sets, fusion_bond_idxs, [0, 1], tuple(sp3_atoms)


def name_dihydronaphthalene(mol, core) -> str:
    ring_atom_sets, fusion_bond_idxs, ring_order, sp3_atoms = core
    atom_rings = [tuple(s) for s in ring_atom_sets]

    candidates = _straight_chain_candidates(mol, atom_rings, ring_atom_sets, fusion_bond_idxs, ring_order)
    best_locants = None
    for locants in candidates:
        pair = tuple(sorted(locants[a] for a in sp3_atoms))
        if best_locants is None or pair < best_locants:
            best_locants = pair

    return f"{best_locants[0]},{best_locants[1]}-dihydronaphthalene"


def find_decahydronaphthalene_core(mol):
    """Return `_bicyclic.py`'s bicyclic core if `mol` is naphthalene's
    carbon skeleton fully saturated (a bare, unsubstituted 0-bridge 6,6
    bicyclic all-carbon hydrocarbon), else None."""
    if mol.GetNumAtoms() != 10:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            return None
    core = find_bicyclic_core(mol)
    if core is None:
        return None
    _, _, bridges = core
    if sorted(len(bridge) for bridge in bridges) != [0, 4, 4]:
        return None
    return core


def name_decahydronaphthalene(mol, core) -> str:
    return "decahydronaphthalene"


def _hydro_locant_sort_key(label):
    digits = "".join(c for c in label if c.isdigit())
    letters = "".join(c for c in label if c.isalpha())
    return int(digits), letters


# Naphthalene's own reference Kekule structure (P-25.3.3.1.1's fixed
# numbering, the symmetric form with the 4a-8a fusion bond double) -- the
# 5-edge perfect matching every valid partial-hydrogenation state below
# (M2 step 1, #828) is checked as a subset of, per the issue's own scope
# ("removing double bonds from naphthalene's own Kekule structure one at a
# time" -- this project's practical scope, not every one of naphthalene's
# 3 possible Kekule structures).
_REFERENCE_KEKULE_PAIRS = ((1, 2), (3, 4), ("4a", "8a"), (5, 6), (7, 8))


def find_partially_unsaturated_naphthalene_core(mol):
    """Return a sorted hydro-locant list (e.g. `["1", "2", "3", "4", "4a",
    "8a"]`) if `mol` is naphthalene's carbon skeleton with 1-4 ring double
    bonds forming a valid partial-hydrogenation state of
    `_REFERENCE_KEKULE_PAIRS` (see module docstring), else None. Unlike
    `find_dihydronaphthalene_core` above, no ring atom may be RDKit
    aromatic-flagged -- a literally-aromatic-ring case is that function's
    own, narrower shape."""
    if mol.GetNumAtoms() != 10:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2 or any(len(r) != 6 for r in atom_rings):
        return None
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
                return None
            if atom.GetIsAromatic():
                return None

    bond_ring_count = {}
    for bonds in ring_info.BondRings():
        for b in bonds:
            bond_ring_count[b] = bond_ring_count.get(b, 0) + 1
    fusion_bond_idxs = {b for b, c in bond_ring_count.items() if c >= 2}
    if len(fusion_bond_idxs) != 1:
        return None

    ring_atom_sets = [set(r) for r in atom_rings]
    fusion_atoms = ring_atom_sets[0] & ring_atom_sets[1]
    if len(fusion_atoms) != 2:
        return None
    for a in fusion_atoms:
        if mol.GetAtomWithIdx(a).GetDegree() != 3:
            return None
    non_fusion_atoms = (ring_atom_sets[0] | ring_atom_sets[1]) - fusion_atoms
    for a in non_fusion_atoms:
        if mol.GetAtomWithIdx(a).GetDegree() != 2:
            return None

    double_bond_pairs = set()
    for bond in mol.GetBonds():
        order = bond.GetBondTypeAsDouble()
        if order == 1.0:
            continue
        if order != 2.0:
            return None
        double_bond_pairs.add(frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    if not 1 <= len(double_bond_pairs) <= 4:
        return None

    atom_rings_tuples = [tuple(s) for s in ring_atom_sets]
    candidates = _straight_chain_candidates(mol, atom_rings_tuples, ring_atom_sets, fusion_bond_idxs, [0, 1])

    best_labels = None
    for locants in candidates:
        atom_by_locant = {v: k for k, v in locants.items()}
        atom4a = next(
            a for a in fusion_atoms if a in {n.GetIdx() for n in mol.GetAtomWithIdx(atom_by_locant[4]).GetNeighbors()}
        )
        atom8a = next(a for a in fusion_atoms if a != atom4a)
        label_of_atom = {atom: str(locant) for atom, locant in locants.items()}
        label_of_atom[atom4a] = "4a"
        label_of_atom[atom8a] = "8a"

        reference_edges = {
            frozenset(
                (
                    atom8a if pos == "8a" else atom4a if pos == "4a" else atom_by_locant[pos]
                    for pos in pair
                )
            )
            for pair in _REFERENCE_KEKULE_PAIRS
        }
        if not double_bond_pairs <= reference_edges:
            continue

        removed_pairs = reference_edges - double_bond_pairs
        sp3_atoms = {atom for pair in removed_pairs for atom in pair}
        labels = sorted((label_of_atom[atom] for atom in sp3_atoms), key=_hydro_locant_sort_key)
        key = tuple(_hydro_locant_sort_key(label) for label in labels)
        if best_labels is None or key < tuple(_hydro_locant_sort_key(label) for label in best_labels):
            best_labels = labels

    return best_labels


def name_partially_unsaturated_naphthalene(mol, core) -> str:
    hydro_word = multiplied_word(len(core), "hydro")
    return f"{','.join(core)}-{hydro_word}naphthalene"
