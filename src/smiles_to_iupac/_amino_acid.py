"""Naming of a common alpha-amino carboxylic acid using its P-103 retained
name and L/D stereodescriptor, for the plain-hydrocarbon-side-chain subset
(glycine, alanine, valine, leucine), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained names for these four.
- P-103.1.3.1: 'L' corresponds to the CIP 'S' configuration at the
  alpha-carbon for every common amino acid except cysteine (out of scope
  here -- cysteine's side-chain sulfur flips the correspondence to L='R',
  #1057's own explicit exclusion); 'D' corresponds to 'R'. Confirmed
  against real PubChem L-/D- structures for alanine/valine/leucine during
  scoping (#1056 M1 step 1). Glycine's alpha-carbon bears two hydrogens,
  so it's never a stereocenter and gets no L-/D- prefix at all.

Scope, deliberately narrow (#1057): exactly the four side-chain shapes
below, each already correctly recognized as a plain alpha-amino-acid
backbone by `_carboxylic_acid_amine.py`'s shared detection -- this module
only intercepts those four specific side chains ahead of that module's
own generic (CIP-only, no-retained-name) fallback in `core.py`'s dispatch
order, unchanged for every other amino acid or amine/acid combination
(including isoleucine, which has a second side-chain stereocenter and its
own 'allo' complexity, and every heteroatom-bearing or aromatic side
chain, all deferred to a later batch-rollout step).
"""

from ._carboxylic_acid_amine import _amine_on_a_different_carbon, _find_carboxylic_acid_carbon
from ._common import adjacency, find_primary_amines, specified_stereocenters

_UNRECOGNIZED = None

_ALPHA_TO_LD = {"S": "L", "R": "D"}


def _is_plain_carbon(atom):
    return atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetIsotope() == 0


def _is_terminal_methyl(mol, graph, atom_idx, coming_from):
    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = [n for n in graph[atom_idx] if n != coming_from]
    return _is_plain_carbon(atom) and not neighbors and atom.GetTotalNumHs() == 3


def _side_chain_name(mol, graph, root, coming_from):
    """'alanine', 'valine', 'leucine', or `_UNRECOGNIZED` for the side-
    chain branch hanging off `root` (the alpha-carbon's own side-chain
    neighbor, already confirmed to exist -- glycine's no-side-chain case
    is handled by the caller before this is reached)."""
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
        next_neighbors = [n for n in graph[next_atom] if n != root]
        if (
            _is_plain_carbon(next_obj)
            and next_obj.GetTotalNumHs() == 1
            and len(next_neighbors) == 2
            and all(_is_terminal_methyl(mol, graph, n, next_atom) for n in next_neighbors)
        ):
            return "leucine"
        return _UNRECOGNIZED
    return _UNRECOGNIZED


def _match(mol):
    """(retained_name, alpha_carbon_idx_or_None) if `mol` is a plain
    acyclic alpha-amino acid whose side chain is glycine/alanine/valine/
    leucine's shape, else None. `alpha_carbon_idx` is None only for
    glycine (no stereocenter to look up an L/D descriptor for)."""
    if mol.GetRingInfo().NumRings() > 0:
        return None
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return None
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    amines = find_primary_amines(mol)
    if len(amines) != 1:
        return None
    if not _amine_on_a_different_carbon(mol, acid_carbon.GetIdx(), amines):
        return None
    (amine_n,) = amines

    graph = adjacency(mol)
    acid_carbon_idx = acid_carbon.GetIdx()
    alpha_carbon = graph[amine_n][0]
    if acid_carbon_idx not in graph[alpha_carbon]:
        return None

    exclude = {amine_n, acid_carbon_idx}
    side_neighbors = [n for n in graph[alpha_carbon] if n not in exclude]
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if not _is_plain_carbon(alpha_atom):
        return None

    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
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

    side_chain_atoms = {"alanine": 1, "valine": 3, "leucine": 4}[name]
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
    return f"{_ALPHA_TO_LD[labels[alpha_carbon]]}-{name}"
