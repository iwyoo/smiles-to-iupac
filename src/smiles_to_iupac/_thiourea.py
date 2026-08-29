"""Naming of thiourea (H2N-C(=S)-NH2), the sulfur analogue of urea, and its
N-substituted derivatives, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'thiourea' is
  a retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- mirroring `_urea.py`'s own 'urea' exactly,
  just with sulfur instead of oxygen. Confirmed via PubChem structure
  match: `NC(=S)N` -> "thiourea" (`CNC(=S)N` -> "methylthiourea" -- a
  structure-only match; PubChem's own generated name doesn't carry a
  locant at all here, an ambiguity that doesn't arise once a second
  substituent forces one). Like urea, thiourea has no parent-hydride chain
  to number and no chain-length logic of its own.
- Each nitrogen may carry 0, 1, or 2 plain, unbranched, saturated alkyl
  substituents, cited exactly the way `_urea.py` already established for
  its own N-/N,N-/N,N'- letter-locant convention (confirmed from the Blue
  Book's own text for urea itself, `tmp/bluebook/P6.txt` lines 630,
  1373-1378) -- thiourea and urea are named by the same P-6 retained-name
  mechanism (differing only in the chalcogen), so the identical
  N-substitution citation logic is reused unchanged rather than
  re-verified from a thiourea-specific worked example.

Scope, deliberately narrow, identical to `_urea.py`'s own: substituents
landing on a single nitrogen (one or two, using the same 'N-'/'N,N-di'
citation), or an identical single substituent on each of the two
different nitrogens (symmetric 'N,N'-di...' citation). Explicitly out of
scope (raise `UnsupportedStructure`): two DIFFERENT substituents split
across the two different nitrogens (no confirmed worked example settles
which physical nitrogen becomes 'N' vs 'N'' in that case), a
branched/unsaturated/ring-bearing N-substituent, a ring-fused thiourea,
and the selenium/tellurium analogues (selenourea/tellurourea).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, carbon_adjacency, linear_branch, non_single_bonds
from ._numerals import alkyl_name


def _thiourea_core(mol):
    """(carbon_idx, (nitrogen1_idx, nitrogen2_idx)) for the thiourea carbonyl
    carbon and its two nitrogens, or None if the molecule isn't shaped like
    a thiourea core at all (a carbon with exactly one double-bonded, terminal
    sulfur and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        neighbors = atom.GetNeighbors()
        sulfurs = [n for n in neighbors if n.GetAtomicNum() == 16]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(sulfurs) != 1 or len(nitrogens) != 2:
            continue
        (sulfur,) = sulfurs
        if sulfur.GetDegree() != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), sulfur.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
            continue
        if any(nn.GetAtomicNum() != 6 for n in nitrogens for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx()):
            continue
        return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    return None


def has_thiourea_shape(mol) -> bool:
    return _thiourea_core(mol) is not None


def _n_substituent_carbons(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(
        n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx
    )


def _substituent_names(carbon_graph, substituent_carbons):
    names = []
    for c in substituent_carbons:
        length = linear_branch(carbon_graph, c, None)
        if length is None:
            raise UnsupportedStructure("a branched N-substituent is not supported yet")
        names.append(alkyl_name(length))
    return names


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        previous, current = None, root
        while current is not None:
            atoms.add(current)
            neighbors = [n for n in carbon_graph[current] if n != previous]
            previous, current = current, (neighbors[0] if neighbors else None)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def _n_prefix(letter, names):
    if not names:
        return ""
    if len(names) == 1:
        return f"{letter}-{names[0]}"
    if names[0] == names[1]:
        return f"{letter},{letter}-di{names[0]}"
    a, b = sorted(names)
    return f"{letter}-{a}-{letter}-{b}"


def name_thiourea(mol) -> str:
    core = _thiourea_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "no thiourea (H2N-C(=S)-NH2 or an N-substituted derivative) shape "
            "found; this module only handles thiourea and simple N-substituted "
            "thioureas"
        )
    carbon_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-fused thiourea (e.g. hydantoin) is out of scope for this module"
        )

    (sulfur_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 16)
    n1_carbons = _n_substituent_carbons(mol, n1_idx, carbon_idx)
    n2_carbons = _n_substituent_carbons(mol, n2_idx, carbon_idx)

    carbon_graph = carbon_adjacency(mol)
    n1_chain_atoms = _substituent_chain_atoms(carbon_graph, n1_carbons)
    n2_chain_atoms = _substituent_chain_atoms(carbon_graph, n2_carbons)
    known_atoms = {carbon_idx, sulfur_idx, n1_idx, n2_idx} | n1_chain_atoms | n2_chain_atoms
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in known_atoms:
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the "
                "thiourea core and its plain N-alkyl substituents is not "
                "supported yet"
            )

    _reject_unsaturated_substituents(mol, n1_chain_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms)

    n1_names = _substituent_names(carbon_graph, n1_carbons)
    n2_names = _substituent_names(carbon_graph, n2_carbons)

    if not n1_names and not n2_names:
        return "thiourea"

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1 or n1_names[0] != n2_names[0]:
            raise UnsupportedStructure(
                "different substituents split across thiourea's two nitrogens "
                "is not supported yet (no confirmed worked example settles "
                "which nitrogen becomes N vs N' in that case)"
            )
        return f"N,N'-di{n1_names[0]}thiourea"

    names = n1_names or n2_names
    return f"{_n_prefix('N', names)}thiourea"
