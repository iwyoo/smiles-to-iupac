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

- P-71.2.1.2 (the "general method"): when the radical
  carbon is itself a branch point (not a chain terminus), reusing
  `_substituents.py`'s `name_branch` for this was tried first and found to
  return the pre-2013 substituent name ("1-methylethyl") rather than the
  Blue Book PIN ("propan-2-yl") -- fixing `name_branch` itself is a much
  larger, riskier axis (parked: at least 8 existing test files assert the
  old-style name as a *substituent* prefix, which is still correct general
  nomenclature there, just not radical PIN naming). Instead, this module
  implements P-29.3.2.2
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

- P-71.2.2.1: a divalent or trivalent radical center (`=CH2` methylidene,
  `#CH` methylidyne, ...) on the same unbranched-chain-terminus or
  monocyclic-ring shape is named the identical way, just appending
  'idene'/'idyne' after the '-yl' name instead of using it bare -- e.g.
  'methyl' + 'idene' -> 'methylidene', 'cyclohexyl' + 'idene' ->
  'cyclohexylidene' (worked examples confirmed against `tmp/bluebook/
  P7.txt` lines ~284-315). RDKit's `GetNumRadicalElectrons()` reports this
  free valence directly (2 or 3), with no bond-order difference from the
  monovalent case -- `[CH]C` (ethylidene) and `[C]C` (ethylidyne) both
  have a degree-1 radical carbon, same as a monovalent chain terminus.

- P-71.2.3: two separate monovalent radical centers on different atoms of
  the same unbranched acyclic chain or monocyclic ring is a distinct
  '-diyl' mechanism (`_name_chain_diradical`/`_name_ring_diradical`), not
  a generalization of the single-center '-yl'/'-ylidene'/'-ylidyne' cases
  above: the parent's numbering is chosen to give the lowest *combined*
  locant set to both radical positions together (mirroring how
  `_isotope.py`'s multi-position deuterium and `_alcohol.py`'s multi-
  hydroxyl chains already pick numbering direction), and both locants are
  always cited explicitly (e.g. 'ethane-1,2-diyl', 'propane-1,3-diyl',
  'butane-1,4-diyl'), never elided the way a single terminal '-yl' is --
  P-14.3.3's single-position omission never applies once there are two
  positions to distinguish. No elision of the parent stem's trailing 'e'
  before '-diyl' either, mirroring the analogous '-ylidene'/'-ylidyne'
  non-elided pattern.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any branch off the radical carbon that is itself further branched
  (P-29.5, "complex substituent groups") -- only the radical carbon itself
  may be a branch point.
- A divalent/trivalent radical carbon that is itself a branch point
  (P-29.3.2.2's general method is monovalent-only here; the branched case
  for '-ylidene'/'-ylidyne' needs its own locant-citation research, not
  done yet).
- A '-diyl' radical carbon that is itself a branch point (mirrors the
  single-center exclusion above), three or more radical centers, a mixed
  monovalent+divalent/trivalent combination on the same molecule (P-71.6's
  'ethan-1-yl-2-ylidene'-shaped case), a radical on a functional group
  (P-71.3), on an aromatic ring, on a polycyclic/spiro skeleton, or
  coexisting with any heteroatom, halogen, charge, or isotopic
  modification.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, linear_branch, ring_cycle
from ._numerals import alkane_name, alkyl_name
from ._substituents import format_substituent_prefixes


def has_radical_shape(mol) -> bool:
    """True if the molecule contains any atom with a nonzero radical
    electron count, regardless of whether the rest of the molecule is in
    scope. Used by `core.py` to route here before every other branch, none
    of which recognize a radical center at all."""
    return any(atom.GetNumRadicalElectrons() != 0 for atom in mol.GetAtoms())


def _validate_carbon_skeleton(mol):
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


def name_radical(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    radicals = [atom for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() != 0]
    if len(radicals) == 2 and all(r.GetNumRadicalElectrons() == 1 for r in radicals):
        _validate_carbon_skeleton(mol)
        ring_info = mol.GetRingInfo()
        num_rings = ring_info.NumRings()
        if num_rings == 0:
            return _name_chain_diradical(mol, radicals)
        if num_rings == 1:
            return _name_ring_diradical(mol, ring_info, radicals)
        raise UnsupportedStructure("polycyclic and spiro radicals are not supported yet")

    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() not in (1, 2, 3):
        raise UnsupportedStructure(
            "only a single radical center of valence 1, 2, or 3, or two "
            "monovalent radical centers (P-71.2.3), is supported yet; "
            "zero, three or more, or a mixed-valence combination of "
            "radical centers is not supported"
        )
    (radical,) = radicals
    valence = radical.GetNumRadicalElectrons()

    _validate_carbon_skeleton(mol)

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_chain_radical(mol, radical, valence)
    if num_rings == 1:
        return _name_ring_radical(mol, ring_info, valence)
    raise UnsupportedStructure("polycyclic and spiro radicals are not supported yet")


def _radical_suffix(yl_name: str, valence: int) -> str:
    """P-71.2.2.1: the divalent/trivalent suffix is formed by appending
    'idene'/'idyne' after the '-yl' name (not replacing it), e.g.
    'methyl' -> 'methylidene'/'methylidyne'."""
    if valence == 1:
        return yl_name
    return yl_name + ("idene" if valence == 2 else "idyne")


def _name_chain_radical(mol, radical, valence) -> str:
    if mol.GetNumAtoms() == 1:
        # P-71.2.1.1's own wording covers this directly: "a mononuclear
        # parent hydride of an element of Group 14" -- methyl (*CH3) has
        # no bond to another atom, so degree 0 rather than 1.
        return _radical_suffix(alkyl_name(1), valence)
    if radical.GetDegree() != 1:
        if valence != 1:
            raise UnsupportedStructure(
                "a divalent/trivalent radical carbon that is itself a "
                "branch point is out of scope for this module (see module "
                "docstring)"
            )
        return _name_branch_point_radical(mol, radical)
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            raise UnsupportedStructure(
                "a branched chain is out of scope for this module "
                "(P-71.2.1.2, the 'general method')"
            )
    return _radical_suffix(alkyl_name(mol.GetNumAtoms()), valence)


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


def _name_ring_radical(mol, ring_info, valence) -> str:
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
    return _radical_suffix("cyclo" + alkyl_name(len(ring_atoms)), valence)


def _name_chain_diradical(mol, radicals) -> str:
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            raise UnsupportedStructure(
                "a branched chain is out of scope for this module "
                "(P-71.2.1.2, the 'general method')"
            )

    n = mol.GetNumAtoms()
    graph = adjacency(mol)
    termini = [idx for idx, neighbors in graph.items() if len(neighbors) <= 1]
    if len(termini) != 2:
        raise UnsupportedStructure(
            "not a single unbranched chain (P-71.2.1.1, unbranched chains only)"
        )
    (start, _) = termini
    order = [start]
    previous, current = None, start
    while len(order) < n:
        next_atoms = [a for a in graph[current] if a != previous]
        if not next_atoms:
            break
        previous, current = current, next_atoms[0]
        order.append(current)
    if len(order) != n:
        raise UnsupportedStructure(
            "not a single unbranched chain (P-71.2.1.1, unbranched chains only)"
        )

    r1, r2 = radicals[0].GetIdx(), radicals[1].GetIdx()
    positions_fwd = {atom: i + 1 for i, atom in enumerate(order)}
    positions_rev = {atom: n - i for i, atom in enumerate(order)}
    locants_fwd = tuple(sorted((positions_fwd[r1], positions_fwd[r2])))
    locants_rev = tuple(sorted((positions_rev[r1], positions_rev[r2])))
    lo, hi = min(locants_fwd, locants_rev)
    return f"{alkane_name(n)}-{lo},{hi}-diyl"


def _name_ring_diradical(mol, ring_info, radicals) -> str:
    (ring_atoms,) = ring_info.AtomRings()
    if len(ring_atoms) != mol.GetNumAtoms():
        raise UnsupportedStructure(
            "a substituent hanging off the ring (other than the radical "
            "centers) is out of scope for this module"
        )
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 2:
            raise UnsupportedStructure(
                "a monocyclic radical ring must otherwise be unsubstituted "
                "(P-71.2.1.1)"
            )

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    r1, r2 = radicals[0].GetIdx(), radicals[1].GetIdx()

    best = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locants = tuple(sorted((position_of[r1], position_of[r2])))
            if best is None or locants < best:
                best = locants
    lo, hi = best
    return f"cyclo{alkane_name(ring_size)}-{lo},{hi}-diyl"
