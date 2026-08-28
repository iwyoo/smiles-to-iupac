"""Naming of simple carbon-centered radicals ('methyl', 'propyl',
'cyclobutyl', ...), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-71.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf),
  the "specific method": "A radical formally derived by the removal of one
  hydrogen atom from a mononuclear parent hydride of an element of Group
  14, from a terminal atom of an unbranched acyclic hydrocarbon, or from
  any position of a monocyclic saturated hydrocarbon ring is named by
  replacing the 'ane' ending of the systematic name of the parent hydride
  by 'yl'." Confirmed worked examples: *CH3 -> 'methyl (PIN)'; a terminal
  radical on propane -> 'propyl (PIN)'; a cyclobutane ring radical ->
  'cyclobutyl (PIN)'.
- This reuses `_numerals.py`'s existing `alkyl_name` (already used by
  `_substituents.py` for exactly this "ane"->"yl" replacement, including
  the retained one-to-four-carbon names 'methyl'/'ethyl'/'propyl'/'butyl'),
  applied to a bare, otherwise-unsubstituted unbranched chain or
  monocyclic ring -- no locant is ever cited, since a terminal chain
  position is always locant 1, and a monocyclic ring's radical position is
  symmetric ("any position").

- P-71.2.1.2 (the "general method"), as of
  `tasks/radical-branch-point-naming.md` (2026-08-28): when the radical
  carbon is itself a branch point (not a chain terminus), reusing
  `_substituents.py`'s `name_branch` for this was tried first and found to
  return the pre-2013 substituent name ("1-methylethyl") rather than the
  Blue Book PIN ("propan-2-yl") -- fixing `name_branch` itself is a much
  larger, riskier axis (`tasks/parent-derived-substituent-prefixes.md`,
  parked: at least 8 existing test files assert the old-style name as a
  *substituent* prefix, which is still correct general nomenclature there,
  just not radical PIN naming). Instead, this module implements P-29.3.2.2
  directly and independently: number the two longest branches plus the
  root as one parent chain, citing the free valence's own locant (e.g.
  'propan-2-yl', 'butan-2-yl', never the elided-locant 'prop-2-yl'); any
  third branch becomes an ordinary substituent prefix at that same locant
  (e.g. 'sec-butyl'/'tert-pentyl' are general-nomenclature-only names for
  radicals this module instead renders as their PINs, 'butan-2-yl'/
  '2-methylbutan-2-yl' -- P-29.6.2.2/P-29.6.3 confirm neither retained
  name is a PIN). The sole exception is P-29.6.1: unsubstituted (CH3)3C-
  keeps its retained PIN 'tert-butyl' rather than the rule's own
  '2-methylpropan-2-yl'.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any branch off the radical carbon that is itself further branched
  (P-29.5, "complex substituent groups") -- only the radical carbon itself
  may be a branch point.
- Divalent/trivalent radicals ('-ylidene'/'-ylidyne', P-71.2.2), more than
  one radical center (P-71.2.3), a radical on a functional group (P-71.3),
  on an aromatic ring, on a polycyclic/spiro skeleton, or coexisting with
  any heteroatom, halogen, charge, or isotopic modification.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, linear_branch
from ._numerals import alkane_name, alkyl_name
from ._substituents import format_substituent_prefixes


def has_radical_shape(mol) -> bool:
    """True if the molecule contains any atom with a nonzero radical
    electron count, regardless of whether the rest of the molecule is in
    scope. Used by `core.py` to route here before every other branch, none
    of which recognize a radical center at all."""
    return any(atom.GetNumRadicalElectrons() != 0 for atom in mol.GetAtoms())


def name_radical(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    radicals = [atom for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() != 1:
        raise UnsupportedStructure(
            "only a single, monovalent radical center is supported yet "
            "(P-71.2.1.1); zero, multiple, or higher-valence radical "
            "centers are not supported"
        )
    (radical,) = radicals

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "only an all-carbon skeleton is supported yet (P-71.2.1.1)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("unsaturated skeletons are not supported yet")

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_chain_radical(mol, radical)
    if num_rings == 1:
        return _name_ring_radical(mol, ring_info)
    raise UnsupportedStructure("polycyclic and spiro radicals are not supported yet")


def _name_chain_radical(mol, radical) -> str:
    if mol.GetNumAtoms() == 1:
        # P-71.2.1.1's own wording covers this directly: "a mononuclear
        # parent hydride of an element of Group 14" -- methyl (*CH3) has
        # no bond to another atom, so degree 0 rather than 1.
        return alkyl_name(1)
    if radical.GetDegree() != 1:
        return _name_branch_point_radical(mol, radical)
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            raise UnsupportedStructure(
                "a branched chain is out of scope for this module "
                "(P-71.2.1.2, the 'general method')"
            )
    return alkyl_name(mol.GetNumAtoms())


def _name_branch_point_radical(mol, radical) -> str:
    """P-29.3.2.2: the radical carbon itself is a branch point (not a
    chain terminus) -- see module docstring for the full derivation."""
    graph = adjacency(mol)
    root_idx = radical.GetIdx()
    branch_roots = list(graph[root_idx])
    lengths = []
    for branch_root in branch_roots:
        length = linear_branch(graph, branch_root, root_idx)
        if length is None:
            raise UnsupportedStructure(
                "a substituent branch with its own branch point is out of "
                "scope for this module (P-29.5, complex substituent groups)"
            )
        lengths.append(length)

    if len(branch_roots) == 3 and all(length == 1 for length in lengths):
        # P-29.6.1: the retained name 'tert-butyl' is the PIN for the
        # unsubstituted (CH3)3C- radical, never the general rule's own
        # '2-methylpropan-2-yl'.
        return "tert-butyl"

    order = sorted(range(len(branch_roots)), key=lambda i: -lengths[i])
    chain_indices = order[:2]
    extra_indices = order[2:]
    side_lengths = [lengths[i] for i in chain_indices]
    chain_length = side_lengths[0] + side_lengths[1] + 1
    root_locant = min(side_lengths[0] + 1, side_lengths[1] + 1)
    stem = alkane_name(chain_length)[:-1]

    if not extra_indices:
        return f"{stem}-{root_locant}-yl"

    (extra_idx,) = extra_indices
    extra_name = alkyl_name(lengths[extra_idx])
    grouped = {extra_name: {"locants": [root_locant], "compound": False}}
    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{stem}-{root_locant}-yl"


def _name_ring_radical(mol, ring_info) -> str:
    (ring_atoms,) = ring_info.AtomRings()
    if len(ring_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a substituent hanging off the ring (other than the radical "
            "itself) is out of scope for this module"
        )
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 2:
            raise UnsupportedStructure(
                "a monocyclic radical ring must otherwise be unsubstituted "
                "(P-71.2.1.1)"
            )
    return "cyclo" + alkyl_name(len(ring_atoms))
