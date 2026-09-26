"""Naming of a common alpha-amino carboxylic acid using its P-103 retained
name and L/D stereodescriptor, for the side chains below, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained names, keyed by a canonical
  fragment SMILES of the side chain itself (see `_SIDE_CHAIN_TABLE`)
  rather than a hand-coded graph walk per shape -- a growing per-shape
  `if`/`elif` chain (this module's own earlier form, #1057/#1068/#1070)
  is a narrow-enumeration smell even when each addition individually
  looks justified; a canonical-fragment table lookup is the general
  mechanism that actually covers the whole family in one place. New side
  chains are added by extending the table, not by writing a new
  structural walk.
- P-103.1.3.1: 'L' corresponds to the CIP 'S' configuration at the
  alpha-carbon for every common amino acid *except* cysteine, whose
  side-chain sulfur outranks the ring-ward carbon in CIP priority and so
  flips the correspondence to L='R'/D='S' -- confirmed against
  cysteine's real PubChem L-/D- structures (CIDs 92851, 5862). 'D'
  corresponds to 'R' (or 'S' for cysteine). Confirmed against real
  PubChem L-/D- structures for alanine/valine/leucine/serine/aspartic
  acid/glutamic acid during scoping (#1056 M1 steps 1-3) -- no exception
  for aspartic/glutamic acid, the ordinary S=L/R=D rule holds. Glycine's
  alpha-carbon bears two hydrogens, so it's never a stereocenter and gets
  no L-/D- prefix at all.

Scope, deliberately narrow (table entries only): exactly the side chains
in `_SIDE_CHAIN_TABLE`, each already correctly recognized as a plain
alpha-amino-acid backbone by `_carboxylic_acid_amine.py`'s shared
detection -- this module only intercepts those specific side chains
ahead of that module's own generic (CIP-only, no-retained-name) fallback
in `core.py`'s dispatch order, unchanged for every other amino acid or
amine/acid combination (including isoleucine/threonine, which have a
second side-chain stereocenter and their own 'allo' complexity, and every
ring-containing side chain -- the top-level `NumRings() > 0` guard below
rejects those regardless of the table, deferred to a later step once
this mechanism is extended to handle a ring-bearing side chain).
"""

from rdkit import Chem

from ._common import adjacency, find_primary_amines, specified_stereocenters

_ALPHA_TO_LD = {"S": "L", "R": "D"}
_ALPHA_TO_LD_CYSTEINE = {"R": "L", "S": "D"}  # P-103.1.3.1's stated exception

_SIDE_CHAIN_SMILES = {
    "alanine": "*C",
    "valine": "*C(C)C",
    "leucine": "*CC(C)C",
    "serine": "*CO",
    "cysteine": "*CS",
    "aspartic acid": "*CC(=O)O",
    "glutamic acid": "*CCC(=O)O",
}
# Canonicalized at import time (rather than hardcoding the already-
# canonical strings above) so a future RDKit version's canonicalization
# doesn't silently desync the table from what `_side_chain_fragment`
# actually produces.
_SIDE_CHAIN_TABLE = {Chem.CanonSmiles(smiles): name for name, smiles in _SIDE_CHAIN_SMILES.items()}


def _is_plain_carbon(atom):
    return atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetIsotope() == 0


def _is_terminal_carboxylic_acid_carbon(mol, graph, atom_idx, coming_from):
    """True for a terminal -COOH carbon (bonded only to `coming_from` and
    its own carbonyl/hydroxyl oxygens) hanging off `coming_from` -- used
    to find the *main-chain* acid carbon among the alpha-carbon's own
    neighbors, not for side-chain recognition (aspartic/glutamic acid's
    second, side-chain-terminal carboxylic acid is matched via
    `_SIDE_CHAIN_TABLE` instead, like every other side chain)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if not _is_plain_carbon(atom):
        return False
    neighbors = [n for n in graph[atom_idx] if n != coming_from]
    if len(neighbors) != 2:
        return False
    oxygens = [n for n in neighbors if mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return False
    carbonyls = [
        o
        for o in oxygens
        if mol.GetAtomWithIdx(o).GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom_idx, o).GetBondTypeAsDouble() == 2.0
    ]
    hydroxyls = [
        o
        for o in oxygens
        if mol.GetAtomWithIdx(o).GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom_idx, o).GetBondTypeAsDouble() == 1.0
        and mol.GetAtomWithIdx(o).GetTotalNumHs() == 1
    ]
    return len(carbonyls) == 1 and len(hydroxyls) == 1


def _side_chain_fragment(mol, alpha_carbon, other_alpha_neighbors):
    """The side-chain fragment hanging off `alpha_carbon`, as an RDKit Mol
    plus its canonical SMILES -- the cut bond to `alpha_carbon` is marked
    by a dummy atom (`*`) so the fragment's canonical form is stable
    regardless of what's on the other side of that bond (the amine/acid
    this module already excludes). `(None, None)` if `alpha_carbon` isn't
    itself a plain, unbranched attachment point (shouldn't happen given
    how callers use this, but keeps this function total)."""
    rw = Chem.RWMol(mol)
    for neighbor in other_alpha_neighbors:
        rw.RemoveBond(alpha_carbon, neighbor)
    dummy = rw.GetAtomWithIdx(alpha_carbon)
    dummy.SetAtomicNum(0)
    dummy.SetFormalCharge(0)
    dummy.SetNoImplicit(True)
    dummy.SetNumExplicitHs(0)
    frag_mol = rw.GetMol()
    Chem.SanitizeMol(frag_mol)
    for frag in Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True):
        if any(a.GetAtomicNum() == 0 for a in frag.GetAtoms()):
            return frag, Chem.MolToSmiles(frag)
    return None, None


def _match(mol):
    """(retained_name, alpha_carbon_idx_or_None) if `mol` is a plain
    acyclic alpha-amino acid whose side chain matches a `_SIDE_CHAIN_TABLE`
    entry, else None. `alpha_carbon_idx` is None only for glycine (no
    stereocenter to look up an L/D descriptor for).

    Finds the *main-chain* acid carbon by looking directly at the
    alpha-carbon's own neighbors, rather than `_find_carboxylic_acid_carbon`'s
    whole-molecule first-match search -- aspartic/glutamic acid also carry
    a second, side-chain-terminal carboxylic-acid-shaped carbon, and that
    whole-molecule search can find the wrong one first (confirmed:
    for 'OC(=O)C[C@H](N)C(=O)O', it returned the side-chain acid carbon,
    not the one bonded to the alpha carbon)."""
    if mol.GetRingInfo().NumRings() > 0:
        return None
    amines = find_primary_amines(mol)
    if len(amines) != 1:
        return None
    (amine_n,) = amines

    graph = adjacency(mol)
    alpha_carbon = graph[amine_n][0]
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if not _is_plain_carbon(alpha_atom):
        return None

    acid_candidates = [
        n
        for n in graph[alpha_carbon]
        if n != amine_n and _is_terminal_carboxylic_acid_carbon(mol, graph, n, alpha_carbon)
    ]
    if len(acid_candidates) != 1:
        return None
    (acid_carbon_idx,) = acid_candidates

    side_neighbors = [n for n in graph[alpha_carbon] if n not in (amine_n, acid_carbon_idx)]

    if not side_neighbors:
        if alpha_atom.GetTotalNumHs() != 2:
            return None
        # amine N + acid C + its 2 oxygens + alpha C, nothing else anywhere
        # in the molecule (rejects e.g. a disconnected salt fragment).
        if mol.GetNumAtoms() != 5:
            return None
        return "glycine", None

    if len(side_neighbors) != 1 or alpha_atom.GetTotalNumHs() != 1:
        return None

    frag, frag_smiles = _side_chain_fragment(mol, alpha_carbon, (amine_n, acid_carbon_idx))
    if frag_smiles is None:
        return None
    name = _SIDE_CHAIN_TABLE.get(frag_smiles)
    if name is None:
        return None
    # amine N + acid C + its 2 oxygens + alpha C + the side chain's own
    # atoms (frag includes a dummy atom standing in for alpha_carbon, so
    # its real atom count is frag.GetNumAtoms() - 1) - same disconnected-
    # fragment guard as glycine's case above, generalized via the
    # fragment's own atom count instead of a per-name lookup table.
    if mol.GetNumAtoms() != 4 + frag.GetNumAtoms():
        return None
    return name, alpha_carbon


def has_amino_acid_shape(mol) -> bool:
    return _match(mol) is not None


def name_amino_acid(mol) -> str:
    name, alpha_carbon = _match(mol)
    if alpha_carbon is None:
        return name
    stereo = specified_stereocenters(mol)
    if not stereo:
        return name
    labels = {atom_idx: label for atom_idx, label in stereo}
    if alpha_carbon not in labels:
        return name
    mapping = _ALPHA_TO_LD_CYSTEINE if name == "cysteine" else _ALPHA_TO_LD
    return f"{mapping[labels[alpha_carbon]]}-{name}"
