"""Naming of selenourea (H2N-C(=Se)-NH2), the selenium analogue of urea,
and its N-substituted derivatives, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-66.1.6.1.3.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  "Chalcogen analogues of urea are named by functional replacement
  nomenclature using the prefixes 'thio', 'seleno', and 'telluro'." --
  'selenourea' is a retained name that is itself the preferred IUPAC name,
  mirroring `_thiourea.py`'s own 'thiourea' exactly, just with selenium
  instead of sulfur. Confirmed via the Blue Book's own worked example
  (`tmp/bluebook/P6a.txt` line 1143): "N-(butan-2-yl)selenourea (PIN)".
  Confirmed via PubChem structure match: `NC(=[Se])N` is CID 6327594,
  whose synonym list includes "Selenourea" (PubChem's own computed
  IUPACName property is unavailable for this structure, the same PubChem
  engine limitation already seen for other selenium/tellurium compounds
  this project has covered -- structure match plus the Blue Book's own
  explicit rule text and worked example are the evidence here).
- Each nitrogen may carry 0, 1, or 2 plain, unsubstituted, saturated,
  acyclic alkyl substituents (branched or unbranched), cited exactly the
  way `_thiourea.py`/`_urea.py` already established for their own
  N-/N,N-/N,N'- letter-locant convention. Each N-substituent's own name
  is built with `name_branch` (P-29 PIN style, fixed project-wide by PR
  #237) -- a branched N-substituent is directly confirmed by this
  chapter's own worked example: 'N-(butan-2-yl)selenourea (PIN)'
  (`tmp/bluebook/P6a.txt` line 1143) -- a compound (has its own locant)
  N-substituent is always parenthesized, even alone, matching
  `_urea.py`'s identical correction (PR #341).

Scope, deliberately narrow, identical to `_thiourea.py`'s own:
substituents landing on a single nitrogen (one or two, using the same
'N-'/'N,N-di' citation), or an identical single substituent on each of the
two different nitrogens (symmetric 'N,N'-di...' citation). Explicitly out
of scope (raise `UnsupportedStructure`): two DIFFERENT substituents split
across the two different nitrogens, an unsaturated N-substituent, a
ring-bearing N-substituent other than a single plain (unsubstituted)
benzene ring, and a ring-fused selenourea.

- A plain benzene ring bonded directly to a nitrogen is cited as
  'phenyl', mirroring `_urea.py`'s/`_thiourea.py`'s/`_guanidine.py`'s
  identical fix (PR #382/#383/#384). PubChem's IUPACName property is
  unavailable for this selenium compound (the same engine limitation
  noted above for plain selenourea), so this is a structure match plus
  synonym confirmation instead: `NC(=[Se])Nc1ccccc1` is CID 6329315,
  synonym "Phenylselenourea"; `c1ccc(NC(=[Se])Nc2ccccc2)cc1` is CID
  6327858, synonym "N,N'-Diphenylselenourea". A *substituted* phenyl
  ring or a second substituent sharing that same nitrogen is out of
  scope.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
from ._substituents import alpha_sort_key, name_branch


def _selenourea_core(mol):
    """(carbon_idx, (nitrogen1_idx, nitrogen2_idx)) for the selenourea
    carbonyl carbon and its two nitrogens, or None if the molecule isn't
    shaped like a selenourea core at all (a carbon with exactly one
    double-bonded, terminal selenium and two singly-bonded nitrogens, each
    nitrogen bonded only to that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        neighbors = atom.GetNeighbors()
        seleniums = [n for n in neighbors if n.GetAtomicNum() == 34]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(seleniums) != 1 or len(nitrogens) != 2:
            continue
        (selenium,) = seleniums
        if selenium.GetDegree() != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), selenium.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
            continue
        if any(nn.GetAtomicNum() != 6 for n in nitrogens for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx()):
            continue
        return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    return None


def has_selenourea_shape(mol) -> bool:
    return _selenourea_core(mol) is not None


def _n_substituent_carbons(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(
        n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx
    )


def _substituent_names(full_graph, nitrogen_idx, substituent_carbons, aromatic_atoms=frozenset(), mol=None):
    return [
        name_branch(full_graph, c, nitrogen_idx, {}, aromatic_atoms, mol=mol)
        for c in substituent_carbons
    ]


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        reached, _ = bfs(carbon_graph, root)
        atoms.update(reached)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def _di_name(name, is_compound):
    return f"({name})" if is_compound else name


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


def name_selenourea(mol) -> str:
    core = _selenourea_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "no selenourea (H2N-C(=Se)-NH2 or an N-substituted derivative) "
            "shape found; this module only handles selenourea and simple "
            "N-substituted selenoureas"
        )
    carbon_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (selenium_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 34)
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
                "a ring-fused selenourea or a ring N-substituent other "
                "than a plain, unsubstituted benzene ring is out of scope "
                "for this module"
            )

    carbon_graph = carbon_adjacency(mol)
    n1_chain_atoms = _substituent_chain_atoms(carbon_graph, n1_carbons)
    n2_chain_atoms = _substituent_chain_atoms(carbon_graph, n2_carbons)
    known_atoms = {carbon_idx, selenium_idx, n1_idx, n2_idx} | n1_chain_atoms | n2_chain_atoms
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in known_atoms:
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the "
                "selenourea core and its plain N-alkyl substituents is not "
                "supported yet"
            )

    _reject_unsaturated_substituents(mol, n1_chain_atoms - phenyl_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms - phenyl_atoms)

    aromatic_atoms = frozenset(phenyl_atoms)
    n1_names = _substituent_names(full_graph, n1_idx, n1_carbons, aromatic_atoms, mol=mol)
    n2_names = _substituent_names(full_graph, n2_idx, n2_carbons, aromatic_atoms, mol=mol)

    if not n1_names and not n2_names:
        return "selenourea"

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1 or n1_names[0][0] != n2_names[0][0]:
            raise UnsupportedStructure(
                "different substituents split across selenourea's two "
                "nitrogens is not supported yet (no confirmed worked "
                "example settles which nitrogen becomes N vs N' in that "
                "case)"
            )
        name, is_compound = n1_names[0]
        return f"N,N'-di{_di_name(name, is_compound)}selenourea"

    names = n1_names or n2_names
    return f"{_n_prefix('N', names)}selenourea"
