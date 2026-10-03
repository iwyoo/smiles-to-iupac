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
citation), an identical single substituent on each of the two different
nitrogens (symmetric 'N,N'-di...' citation), or one DIFFERENT substituent
on each of the two nitrogens -- the alphanumerical order (P-14.5.2,
locants and italicized prefixes like 'tert-' ignored) decides which
substituent becomes 'N-' and which becomes 'N''-', same rule as `_urea.py`
(PubChem structure match: `CCNC(=S)NC` -> '1-ethyl-3-methylthiourea', CID
15568242). Each N-substituent's own name is built with `name_branch`
(P-29 PIN style, fixed project-wide by PR #237; mirrors `_urea.py`'s
identical fix, PR #332) -- a branched N-substituent is supported (e.g.
'N-tert-butylthiourea', CID 737374, a retained non-compound name), and a
*compound* one is always parenthesized -- 'N-(propan-2-yl)thiourea', not
PubChem's own raw 'N-propan-2-ylthiourea' (CID 1711921), same correction
as `_urea.py` (see that module's docstring for the Blue Book citations).
Explicitly out of scope (raise `UnsupportedStructure`): a different
substituent *count* on each nitrogen (no confirmed worked example settles
that locant tie-break), an unsaturated N-substituent, a ring-bearing
N-substituent other than a single plain (unsubstituted) benzene ring, a
ring-fused thiourea, and the selenium/tellurium analogues (selenourea/
tellurourea).

- A plain benzene ring bonded directly to one nitrogen is cited as
  'phenyl', mirroring `_urea.py`'s identical extension (PR #382):
  PubChem structure match `NC(=S)Nc1ccccc1` -> "phenylthiourea", CID
  676454. Combines with the existing symmetric/asymmetric machinery
  unchanged (`N-methyl-N'-phenylthiourea`, CID 698294). A *substituted*
  phenyl ring or a second substituent sharing that same nitrogen is out
  of scope.
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    plain_phenyl_substituent_atoms,
    reject_unsaturated_substituents,
)
from ._substituents import alpha_sort_key, name_branch


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


def _substituent_names(full_graph, nitrogen_idx, substituent_carbons, aromatic_atoms=frozenset(), mol=None):
    return [name_branch(full_graph, c, nitrogen_idx, {}, aromatic_atoms, mol=mol) for c in substituent_carbons]


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        reached, _ = bfs(carbon_graph, root)
        atoms.update(reached)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    reject_unsaturated_substituents(mol, atoms)


def _di_name(name, is_compound):
    return enclose(name) if is_compound else name


def _n_letter_entry(letter, name, is_compound):
    return f"{letter}-({name})" if is_compound else f"{letter}-{name}"


def _n_prefix(letter, entries):
    if not entries:
        return ""
    if len(entries) == 1:
        (name, is_compound), = entries
        return _n_letter_entry(letter, name, is_compound)
    (name_a, compound_a), (name_b, compound_b) = entries
    if name_a == name_b:
        return f"{letter},{letter}-di{_di_name(name_a, compound_a)}"
    (a, ca), (b, cb) = sorted(entries, key=lambda e: alpha_sort_key(e[0]))
    return f"{_n_letter_entry(letter, a, ca)}-{_n_letter_entry(letter, b, cb)}"


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

    (sulfur_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 16)
    n1_carbons = _n_substituent_carbons(mol, n1_idx, carbon_idx)
    n2_carbons = _n_substituent_carbons(mol, n2_idx, carbon_idx)

    full_graph = adjacency(mol)
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, full_graph, n1_carbons + n2_carbons)
    if phenyl_atoms and (
        (any(c in phenyl_atoms for c in n1_carbons) and len(n1_carbons) > 1)
        or (any(c in phenyl_atoms for c in n2_carbons) and len(n2_carbons) > 1)
    ):
        raise UnsupportedStructure(
            "a phenyl N-substituent alongside another substituent on the "
            "same nitrogen is not supported yet"
        )

    if mol.GetRingInfo().NumRings() > 0:
        all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
        if all_ring_atoms - phenyl_atoms:
            raise UnsupportedStructure(
                "a ring-fused thiourea (e.g. hydantoin) or a ring "
                "N-substituent other than a plain, unsubstituted benzene "
                "ring is out of scope for this module"
            )

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

    _reject_unsaturated_substituents(mol, n1_chain_atoms - phenyl_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms - phenyl_atoms)

    n1_names = _substituent_names(full_graph, n1_idx, n1_carbons, frozenset(phenyl_atoms), mol=mol)
    n2_names = _substituent_names(full_graph, n2_idx, n2_carbons, frozenset(phenyl_atoms), mol=mol)

    if not n1_names and not n2_names:
        return "thiourea"

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1:
            raise UnsupportedStructure(
                "a different substituent count on each of thiourea's two "
                "nitrogens is not supported yet (no confirmed worked "
                "example settles the locant tie-break for that case)"
            )
        (name_a, compound_a), (name_b, compound_b) = n1_names[0], n2_names[0]
        if name_a == name_b:
            return f"N,N'-di{_di_name(name_a, compound_a)}thiourea"
        (first, first_compound), (second, second_compound) = sorted(
            (n1_names[0], n2_names[0]), key=lambda e: alpha_sort_key(e[0])
        )
        first_entry = _n_letter_entry("N", first, first_compound)
        second_entry = _n_letter_entry("N'", second, second_compound)
        return f"{first_entry}-{second_entry}thiourea"

    names = n1_names or n2_names
    return f"{_n_prefix('N', names)}thiourea"
