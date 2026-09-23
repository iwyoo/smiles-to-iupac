"""Naming of simple carbenium cations ('methylium', 'propylium',
'cyclobutylium', ...), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-73.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf),
  the "specific method": "A cation formally derived by the removal of a
  hydride ion, H-, ... from a terminal atom of a saturated unbranched
  acyclic hydrocarbon [or] a saturated monocyclic hydrocarbon ... is named
  by replacing the 'ane' ending in the name of the parent hydride by the
  suffix 'ylium'." Confirmed worked examples (`tmp/bluebook/P7.txt`):
  `[CH3+]` -> 'methylium (PIN)'; a terminal cation on propane ->
  'propylium (PIN)'; a cyclobutane ring cation -> 'cyclobutylium (PIN)'.
  This is structurally identical to `_radical.py`'s own "ane"->"yl"
  mechanism (P-71.2.1.1), just with 'ylium' instead of 'yl' -- this
  module reuses `_numerals.py`'s `alkyl_name` the same way, appending
  'ium' to its own output.

- P-73.2.2.1.2, the "general method": a cation carbon that is itself a
  branch point (not a chain terminus) is named by adding '-ylium' to the
  PIN of the parent hydride (P-2/P-5 parent-hydride selection), with
  elision of the final 'e' -- confirmed worked examples in
  `tmp/bluebook/P7.txt` (1948-1995) are silane/furan/spiro-based (e.g.
  'heptamethyltrisilan-2-ylium (PIN)'), with no plain-carbon-chain
  worked example, but the mechanism ("parent hydride PIN + '-ylium'")
  parallels `_radical.py`'s own P-29.3.2.2 branch-point mechanism
  ('yl' instead of 'ylium') closely enough for the two-branch case: both
  reduce to picking the longest chain through the branch-point atom, that
  atom getting the lowest possible locant. Scoped conservatively to
  exactly two branches off the cation (e.g. 'propan-2-ylium',
  'butan-2-ylium') -- a three-branch shape (needing an extra substituent
  prefix, or the `_radical.py`-style 'tert-butyl' retained-name
  exception) is deferred: neither has a confirmed carbon-only PIN worked
  example here, and PubChem doesn't reliably register/verify these
  cationic structures (tested directly: it silently returns the neutral
  isomer's name instead of erroring), so this project's usual
  verification bar isn't met for that wider case yet.

- P-23.2.1/P-24.2.1: a single cation on a von Baeyer bicyclic/polycyclic
  or monospiro ring system is named the same way as `_alcohol.py`'s/
  `_amine.py`'s own polycyclic/spiro suffix support -- the parent's
  already-established skeleton numbering picks the cation's locant, and
  `-ylium` is appended per P-73.2.2.1.1/.2 above (`_polycyclic_suffix.py`,
  shared mechanism). PubChem doesn't reliably register/verify these
  charged structures either (see above), so verification for this case
  cross-checks the neutral parent-hydride skeleton's own name instead of
  the charged species' name directly.

Explicitly out of scope (raise `UnsupportedStructure`):
- A cation carbon that is a branch point with other than exactly two
  branches (three-branch shapes, and the associated 'tert-butyl'-style
  retained-name question, are deferred -- see above), or where either
  branch is itself further branched (P-29.5-style "complex substituent
  groups").
- More than one cationic center, a cation on an aromatic ring, a
  specified stereocenter alongside a polycyclic/spiro cation, or
  coexisting with any heteroatom, halogen, unsaturation, or isotopic
  modification.
"""

from rdkit import Chem

from ._bicyclic import find_bicyclic_core
from ._common import UnsupportedStructure, adjacency, linear_branch, specified_stereocenters
from ._numerals import alkane_name, alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom


def has_carbenium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged carbon shaped
    like a genuine carbenium (a carbon with total substituent+hydrogen
    count of 3, i.e. one fewer bond than neutral carbon's four). Used by
    `core.py` to route here before any other branch, none of which
    recognize a charged atom."""
    charged_carbons = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 1
    ]
    if len(charged_carbons) != 1:
        return False
    carbon = charged_carbons[0]
    if carbon.GetIsotope() != 0:
        return False
    degree = carbon.GetDegree()
    if degree > 3 or carbon.GetTotalNumHs() + degree != 3:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in carbon.GetNeighbors()
    )


def name_carbenium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_carbons = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() != 0
    ]
    (cation,) = charged_carbons
    if cation.GetFormalCharge() != 1 or cation.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "carbenium carbon is supported (P-73.2.2.1.1)"
        )

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "only an all-carbon skeleton is supported yet (P-73.2.2.1.1)"
            )
        if atom.GetIdx() != cation.GetIdx() and atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("more than one charged atom is not supported yet")
        if atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("unsaturated skeletons are not supported yet")

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_chain_carbenium(mol, cation)
    if num_rings == 1:
        return _name_ring_carbenium(mol, ring_info)
    return _name_von_baeyer_or_spiro_carbenium(mol, cation)


def _name_von_baeyer_or_spiro_carbenium(mol, cation):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single carbenium suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`
    (the shared mechanism `_alcohol.py`/`_amine.py` already use). The
    cation carries no separate heteroatom to exclude from substituent
    enumeration -- the charge sits directly on a ring carbon, unlike
    alcohol's oxygen or amine's nitrogen -- so `excluded` is empty.
    Restricted to no specified stereocenter, matching every other WS2
    step's initial scope (unsaturation is already rejected earlier in
    `name_carbenium`, for every ring count)."""
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a von Baeyer bicyclic/"
            "polycyclic or monospiro carbenium cation is not supported "
            "yet (see P-92)"
        )

    bicyclic_core = find_bicyclic_core(mol)
    polycyclic_core = None
    von_baeyer_ring_count = None
    if bicyclic_core is None:
        for candidate_ring_count in (3, 4, 5, 6):
            polycyclic_core = find_polycyclic_core(mol, candidate_ring_count)
            if polycyclic_core is not None:
                von_baeyer_ring_count = candidate_ring_count
                break
    if bicyclic_core is not None or polycyclic_core is not None:
        return name_von_baeyer_suffix(
            mol, cation.GetIdx(), set(), "ylium", "carbenium", bicyclic_core, polycyclic_core, von_baeyer_ring_count
        )

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro_suffix(mol, cation.GetIdx(), set(), "ylium", "carbenium", spiro_atom)

    raise UnsupportedStructure(
        "polycyclic and fused-ring carbenium cations are not supported "
        "yet (P-23/P-25 numbering integration with a suffix group is "
        "future work)"
    )


def _name_chain_carbenium(mol, cation) -> str:
    if mol.GetNumAtoms() == 1:
        # P-73.2.2.1.1's own wording covers this directly: a mononuclear
        # parent hydride (methane) -- methylium (CH3+) has no bond to
        # another atom, so degree 0 rather than 1.
        return alkyl_name(1) + "ium"
    if cation.GetDegree() != 1:
        return _name_branch_point_carbenium(mol, cation)
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            raise UnsupportedStructure(
                "a branched chain is out of scope for this module "
                "(P-73.2.2.1.2, the 'general method')"
            )
    return alkyl_name(mol.GetNumAtoms()) + "ium"


def _name_branch_point_carbenium(mol, cation) -> str:
    """P-73.2.2.1.2: the cation carbon itself is a branch point (not a
    chain terminus) -- see module docstring for the full derivation.
    Scoped to exactly two branches (mirrors `_radical.py`'s
    `_name_branch_point_radical`'s two-longest-branches mechanism, 'yl'
    replaced by 'ylium'); a third branch or a further-branched branch is
    out of scope."""
    graph = adjacency(mol)
    root_idx = cation.GetIdx()
    branch_roots = list(graph[root_idx])
    if len(branch_roots) != 2:
        raise UnsupportedStructure(
            "a carbenium branch point with other than exactly two "
            "branches is not supported yet (P-73.2.2.1.2, the 'general "
            "method')"
        )
    lengths = []
    for branch_root in branch_roots:
        length = linear_branch(graph, branch_root, root_idx)
        if length is None:
            raise UnsupportedStructure(
                "a substituent branch with its own branch point is out "
                "of scope for this module (P-29.5, complex substituent "
                "groups)"
            )
        lengths.append(length)

    chain_length = lengths[0] + lengths[1] + 1
    root_locant = min(lengths[0] + 1, lengths[1] + 1)
    stem = alkane_name(chain_length)[:-1]
    return f"{stem}-{root_locant}-ylium"


def _name_ring_carbenium(mol, ring_info) -> str:
    (ring_atoms,) = ring_info.AtomRings()
    if len(ring_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a substituent hanging off the ring (other than the "
            "carbenium cation itself) is out of scope for this module"
        )
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 2:
            raise UnsupportedStructure(
                "a monocyclic carbenium ring must otherwise be "
                "unsubstituted (P-73.2.2.1.1)"
            )
    return "cyclo" + alkyl_name(len(ring_atoms)) + "ium"
