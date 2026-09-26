"""Naming of a common alpha-amino carboxylic acid using its P-103 retained
name and L/D stereodescriptor, for the plain-hydrocarbon side chains
(glycine, alanine, valine, leucine), the two single-heteroatom-terminated
side chains (serine, cysteine), and the two second-carboxylic-acid side
chains (aspartic acid, glutamic acid), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained names for these eight.
- P-103.1.3.1: 'L' corresponds to the CIP 'S' configuration at the
  alpha-carbon for every common amino acid *except* cysteine, whose
  side-chain sulfur outranks the ring-ward carbon in CIP priority and so
  flips the correspondence to L='R'/D='S' -- confirmed this session
  against cysteine's real PubChem L-/D- structures (CIDs 92851, 5862).
  'D' corresponds to 'R' (or 'S' for cysteine). Confirmed against real
  PubChem L-/D- structures for alanine/valine/leucine/serine/aspartic
  acid/glutamic acid during scoping (#1056 M1 steps 1-3) -- no exception
  for aspartic/glutamic acid, the ordinary S=L/R=D rule holds. Glycine's
  alpha-carbon bears two hydrogens, so it's never a stereocenter and gets
  no L-/D- prefix at all.

Scope, deliberately narrow (#1057, extended by #1068 and #1070): exactly
the eight side-chain shapes below, each already correctly recognized as a
plain alpha-amino-acid backbone by `_carboxylic_acid_amine.py`'s shared
detection -- this module only intercepts those specific side chains
ahead of that module's own generic (CIP-only, no-retained-name) fallback
in `core.py`'s dispatch order, unchanged for every other amino acid or
amine/acid combination (including isoleucine/threonine, which have a
second side-chain stereocenter and their own 'allo' complexity, and every
other heteroatom-bearing or aromatic side chain, all deferred to a later
batch-rollout step).
"""

from ._common import adjacency, find_primary_amines, specified_stereocenters

_UNRECOGNIZED = None

_ALPHA_TO_LD = {"S": "L", "R": "D"}
_ALPHA_TO_LD_CYSTEINE = {"R": "L", "S": "D"}  # P-103.1.3.1's stated exception


def _is_plain_carbon(atom):
    return atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetIsotope() == 0


def _is_terminal_methyl(mol, graph, atom_idx, coming_from):
    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = [n for n in graph[atom_idx] if n != coming_from]
    return _is_plain_carbon(atom) and not neighbors and atom.GetTotalNumHs() == 3


def _is_terminal_heteroatom(mol, graph, atom_idx, coming_from):
    """True for a terminal -OH/-SH heteroatom (degree 1, one hydrogen, no
    charge/isotope) hanging off `coming_from` -- serine/cysteine's own
    side-chain terminus, distinguished from an ether/thioether or a
    charged/isotopically modified variant."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() not in (8, 16):
        return False
    if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return False
    neighbors = [n for n in graph[atom_idx] if n != coming_from]
    return not neighbors and atom.GetTotalNumHs() == 1


def _is_terminal_carboxylic_acid_carbon(mol, graph, atom_idx, coming_from):
    """True for a terminal -COOH carbon (bonded only to `coming_from` and
    its own carbonyl/hydroxyl oxygens) hanging off `coming_from` --
    aspartic/glutamic acid's own side-chain terminus."""
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


def _side_chain_name(mol, graph, root, coming_from):
    """'alanine', 'valine', 'leucine', 'serine', 'cysteine', 'aspartic
    acid', 'glutamic acid', or `_UNRECOGNIZED` for the side-chain branch
    hanging off `root` (the alpha-carbon's own side-chain neighbor,
    already confirmed to exist -- glycine's no-side-chain case is handled
    by the caller before this is reached)."""
    atom = mol.GetAtomWithIdx(root)
    if not _is_plain_carbon(atom):
        return _UNRECOGNIZED
    neighbors = [n for n in graph[root] if n != coming_from]
    if not neighbors:
        return "alanine" if atom.GetTotalNumHs() == 3 else _UNRECOGNIZED
    if len(neighbors) == 2 and atom.GetTotalNumHs() == 1:
        if all(_is_terminal_methyl(mol, graph, n, root) for n in neighbors):
            return "valine"
        return _UNRECOGNIZED
    if len(neighbors) == 1 and atom.GetTotalNumHs() == 2:
        (next_atom,) = neighbors
        next_obj = mol.GetAtomWithIdx(next_atom)
        if _is_terminal_heteroatom(mol, graph, next_atom, root):
            return "serine" if next_obj.GetAtomicNum() == 8 else "cysteine"
        if _is_terminal_carboxylic_acid_carbon(mol, graph, next_atom, root):
            return "aspartic acid"
        next_neighbors = [n for n in graph[next_atom] if n != root]
        if (
            _is_plain_carbon(next_obj)
            and next_obj.GetTotalNumHs() == 1
            and len(next_neighbors) == 2
            and all(_is_terminal_methyl(mol, graph, n, next_atom) for n in next_neighbors)
        ):
            return "leucine"
        if _is_plain_carbon(next_obj) and next_obj.GetTotalNumHs() == 2 and len(next_neighbors) == 1:
            (next_next_atom,) = next_neighbors
            if _is_terminal_carboxylic_acid_carbon(mol, graph, next_next_atom, next_atom):
                return "glutamic acid"
        return _UNRECOGNIZED
    return _UNRECOGNIZED


def _match(mol):
    """(retained_name, alpha_carbon_idx_or_None) if `mol` is a plain
    acyclic alpha-amino acid whose side chain is glycine/alanine/valine/
    leucine/serine/cysteine/aspartic acid/glutamic acid's shape, else
    None. `alpha_carbon_idx` is None only for glycine (no stereocenter to
    look up an L/D descriptor for).

    Finds the *main-chain* acid carbon by looking directly at the
    alpha-carbon's own neighbors, rather than `_find_carboxylic_acid_carbon`'s
    whole-molecule first-match search -- aspartic/glutamic acid also carry
    a second, side-chain-terminal carboxylic-acid-shaped carbon, and that
    whole-molecule search can find the wrong one first (confirmed this
    session: for 'OC(=O)C[C@H](N)C(=O)O', it returned the side-chain acid
    carbon, not the one bonded to the alpha carbon)."""
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
    acid_oxygens = [n for n in graph[acid_carbon_idx] if n != alpha_carbon]
    carbonyl_oxygen = next(
        o for o in acid_oxygens if mol.GetBondBetweenAtoms(acid_carbon_idx, o).GetBondTypeAsDouble() == 2.0
    )
    hydroxyl_oxygen = next(o for o in acid_oxygens if o != carbonyl_oxygen)

    exclude = {amine_n, acid_carbon_idx}
    side_neighbors = [n for n in graph[alpha_carbon] if n not in exclude]

    excluded_acid_oxygens = {carbonyl_oxygen, hydroxyl_oxygen}
    accounted = exclude | excluded_acid_oxygens

    if not side_neighbors:
        if alpha_atom.GetTotalNumHs() != 2:
            return None
        if mol.GetNumAtoms() != len(accounted) + 1:
            return None
        return "glycine", None

    if len(side_neighbors) != 1 or alpha_atom.GetTotalNumHs() != 1:
        return None
    (side_root,) = side_neighbors
    name = _side_chain_name(mol, graph, side_root, alpha_carbon)
    if name is None:
        return None

    side_chain_atoms = {
        "alanine": 1,
        "valine": 3,
        "leucine": 4,
        "serine": 2,
        "cysteine": 2,
        "aspartic acid": 4,
        "glutamic acid": 5,
    }[name]
    if mol.GetNumAtoms() != len(accounted) + 1 + side_chain_atoms:
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
