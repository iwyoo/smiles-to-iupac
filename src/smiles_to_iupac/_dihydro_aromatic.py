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
- More than one 'hydro' pair (tetrahydronaphthalene/"tetralin" and beyond).
- Any fused ring system other than the plain 2-ring naphthalene shape
  (anthracene, indole, heteroaromatics, ...).
- Any indicated-hydrogen-requiring parent (not applicable to naphthalene
  itself, but future modules built on this one will need to handle it).
- Substituents of any kind, including halogens (P-35.2.1) -- every ring
  atom must have exactly its "bare" degree (2 for a CH/CH2 position, 3 for
  a fusion carbon), so any substituent is rejected by construction.
- Bridged derivatives (methanonaphthalene, etc.) -- `tasks/
  bridged-fused-ring-naming.md`'s territory, a separate, larger task that
  this one is a prerequisite for.
"""

from rdkit import Chem

from ._aromatic import _periphery_cycle, _straight_chain_candidates


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
