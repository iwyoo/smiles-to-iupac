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
amine/acid combination. The backbone-anchor check in `_match` (see its
own docstring) tolerates a side chain with its own extra nitrogen
(lysine, arginine) or its own ring (phenylalanine, tyrosine, tryptophan)
as long as the table recognizes it -- an *unrecognized* multi-amine or
ring-bearing shape still falls through unmatched, same as any other
unrecognized side chain. Isoleucine/threonine (a second side-chain
stereocenter and 'allo' complexity) remain out of scope regardless, since
this whole mechanism only ever looks up *one* alpha-stereocenter.
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
    "asparagine": "*CC(=O)N",
    "glutamine": "*CCC(=O)N",
    "methionine": "*CCSC",
    "lysine": "*CCCCN",
    "arginine": "*CCCNC(N)=N",
    "phenylalanine": "*Cc1ccccc1",
    "tyrosine": "*Cc1ccc(O)cc1",
    "tryptophan": "*Cc1c[nH]c2ccccc12",
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


def _backbone_candidate(mol, graph, amine_n):
    """(alpha_carbon, acid_carbon_idx) if `amine_n` directly anchors a
    plain alpha-amino-acid backbone (bonded to a plain carbon that itself
    has exactly one terminal-carboxylic-acid-carbon neighbor), else None.
    Doesn't care whether `mol` has other amines or rings elsewhere -- a
    side chain with its own extra nitrogen (lysine, arginine) or its own
    ring (phenylalanine, tyrosine, tryptophan) is fine, as long as
    exactly one amine in the whole molecule anchors a backbone this way;
    that side chain is handled entirely by `_side_chain_fragment` +
    `_SIDE_CHAIN_TABLE` below, not by this function."""
    neighbors = graph[amine_n]
    if len(neighbors) != 1:
        return None
    (alpha_carbon,) = neighbors
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
    return alpha_carbon, acid_candidates[0]


def _match(mol):
    """(retained_name, alpha_carbon_idx_or_None) if `mol` is a plain
    alpha-amino acid whose side chain matches a `_SIDE_CHAIN_TABLE` entry,
    else None. `alpha_carbon_idx` is None only for glycine (no
    stereocenter to look up an L/D descriptor for).

    Anchors the backbone by finding which amine (of possibly several --
    lysine/arginine's side chains carry their own extra nitrogen) is
    directly bonded to a plain carbon that itself neighbors a terminal
    carboxylic acid carbon, rather than requiring the *whole molecule*
    have exactly one amine and no ring anywhere (confirmed this session:
    that whole-molecule check rejected lysine/arginine/asparagine/
    glutamine outright over their side chain's own nitrogen, and
    phenylalanine/tyrosine/tryptophan over their side chain's own ring,
    long before side-chain matching was ever reached)."""
    graph = adjacency(mol)
    amines = find_primary_amines(mol)
    backbones = []
    for amine_n in amines:
        found = _backbone_candidate(mol, graph, amine_n)
        if found is not None:
            backbones.append((amine_n, *found))
    if len(backbones) != 1:
        return None
    amine_n, alpha_carbon, acid_carbon_idx = backbones[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)

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
