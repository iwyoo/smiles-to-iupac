"""Naming of a plain, saturated, monocyclic amine ring (piperidine/
pyrrolidine/azepane, P-22.2.1) whose sole nitrogen carries exactly one
substituent -- e.g. 1-methylpiperidine -- per the IUPAC 2013
Recommendations ("the Blue Book"):

- `_amine.py` explicitly defers "a secondary/tertiary amine nitrogen on
  or attached to a ring" (see that module's own docstring); this module
  fills in the most basic slice of that gap: the ring itself is the
  parent hydride (its retained/Hantzsch-Widman name, via
  `_hetero_monocyclic.saturated_ring_name`, already used by
  `_ketone.py`'s hetero-ring-ketone path and `_hidden_amide_ketone.py`'s
  pseudoketone path), and the nitrogen's one substituent is cited as an
  ordinary numbered prefix at locant '1' -- the ring's own heteroatom
  locant (P-22.2.1), not the special 'N-' prefix `_amine.py` uses for an
  acyclic secondary/tertiary amine. PubChem structure matches: 'CN1CCCCC1'
  -> '1-methylpiperidine' (CID 12291), 'CCN1CCCCC1' -> '1-ethylpiperidine'
  (CID 13007), 'CN1CCCC1' -> '1-methylpyrrolidine' (CID 8454),
  'ClCCN1CCCCC1' -> '1-(2-chloroethyl)piperidine' (CID 74827),
  'C1CC1N1CCCCC1' -> '1-cyclopropylpiperidine' (CID 10909573).
- The substituent's own name is built with `name_branch` unchanged
  (P-29 PIN style) -- a branched, halogenated, or plain-cyclic substituent
  is supported for free the same way `_carbamate.py`'s/`_ester.py`'s R
  parts are; a *compound* substituent (its own locant, e.g. 'propan-2-yl')
  is parenthesized, this project's usual convention over PubChem's own
  unparenthesized raw name ('1-propan-2-ylpiperidine' -> this module's
  '1-(propan-2-yl)piperidine').

`core.py` must route to this module in the oxygen-free "any nitrogen"
branch, immediately before its final `name_amine` fallback -- a plain
alkyl/cyclic N-substituent has no oxygen at all, so it never reaches any
oxygen-gated check (including `_hidden_amide_ketone.py`'s acyl-on-ring-
nitrogen path, which requires the substituent's own carbonyl oxygen and
is therefore never a routing collision with this module).

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring nitrogen with zero substituents (the plain unsubstituted ring
  itself, already named elsewhere) or more than one (structurally
  impossible for a neutral trivalent ring nitrogen with two ring bonds
  anyway).
- A sulfonyl (or any other oxygen-bearing) N-substituent (e.g.
  'CS(=O)(=O)N1CCCCC1', PubChem's own '1-methylsulfonylpiperidine') --
  `name_branch` doesn't recognize a sulfonyl-rooted substituent at all
  yet; a separate follow-up (overlaps the "sulfone combined with other
  groups" real-data domain).
- An acyl N-substituent -- handled by `_hidden_amide_ketone.py` instead
  (routed earlier, in the oxygen-gated branch).
- An aromatic N-substituent, a ring with more than one heteroatom (e.g.
  morpholine/piperazine -- different locant numbering), a heteroatom
  other than nitrogen, a ring size other than 5/6/7 (`saturated_ring_name`'s
  own P-22.2.1 scope), any substituent on the ring itself, or any ring
  unsaturation.
- Any other heteroatom, charged/isotopically modified atom, or
  multi-fragment structure.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    non_single_bonds,
)
from ._hetero_monocyclic import saturated_ring_name
from ._substituents import name_branch

_RING_SIZES = (5, 6, 7)
_ALLOWED_SUBSTITUENT_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}


def _ring_amine_shape(mol):
    """(n_idx, ring_atoms, substituent_root_idx) if `mol` has exactly one
    plain, saturated, monocyclic ring (5/6/7-membered) with one nitrogen
    heteroatom bearing exactly one exocyclic substituent -- else None. A
    *second*, unrelated ring may still exist inside that one substituent
    (e.g. a cyclopropyl N-substituent) -- `name_branch` recognizes that
    shape on its own; only the amine ring itself is required to be the
    sole N-heterocycle here."""
    ring_info = mol.GetRingInfo()
    candidates = [
        set(ring)
        for ring in ring_info.AtomRings()
        if len(ring) in _RING_SIZES
        and not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
        and sum(1 for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() != 6) == 1
    ]
    if len(candidates) != 1:
        return None
    (ring_atoms,) = candidates
    (n_idx,) = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    n_atom = mol.GetAtomWithIdx(n_idx)
    if n_atom.GetAtomicNum() != 7 or n_atom.GetFormalCharge() != 0 or n_atom.GetIsotope() != 0:
        return None
    exo = [n for n in n_atom.GetNeighbors() if n.GetIdx() not in ring_atoms]
    if len(exo) != 1:
        return None
    if mol.GetBondBetweenAtoms(n_idx, exo[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    for a in ring_atoms:
        if a == n_idx:
            continue
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if any(n.GetIdx() not in ring_atoms for n in atom.GetNeighbors()):
            return None
        if atom.GetTotalNumHs() != 2:
            return None
    if any(a in ring_atoms and b in ring_atoms for a, b, _ in non_single_bonds(mol)):
        return None
    return n_idx, ring_atoms, exo[0].GetIdx()


def has_ring_amine_shape(mol) -> bool:
    return _ring_amine_shape(mol) is not None


def name_ring_amine(mol) -> str:
    shape = _ring_amine_shape(mol)
    if shape is None:
        raise UnsupportedStructure(
            "no plain saturated monocyclic ring with a single N-substituted "
            "nitrogen shape found"
        )
    n_idx, ring_atoms, root = shape
    ring_size = len(ring_atoms)

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            continue
        if atom.GetAtomicNum() not in _ALLOWED_SUBSTITUENT_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "an N-substituent containing anything other than carbon "
                "and halogens is not supported yet for this ring-amine "
                "path"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic():
            raise UnsupportedStructure(
                "an aromatic N-substituent is not supported yet for this "
                "ring-amine path"
            )
        if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1)"
            )

    if any(a not in ring_atoms and b not in ring_atoms for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure("unsaturation in the N-substituent is not supported yet")

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    name, is_compound = name_branch(graph, root, n_idx, halogens)

    stem = saturated_ring_name("N", ring_size)
    if stem is None:
        raise UnsupportedStructure(
            f"no retained/Hantzsch-Widman name for a {ring_size}-membered "
            f"N-heteroatom saturated ring (P-22.2.1)"
        )
    sub_name = f"({name})" if is_compound else name
    return f"1-{sub_name}{stem}"
