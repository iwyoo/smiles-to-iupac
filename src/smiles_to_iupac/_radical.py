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
  'cyclohexylidene' (worked examples confirmed against the Blue Book). RDKit's `GetNumRadicalElectrons()` reports this
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

- P-71.3.1: an acid-derived radical ("acyl" radical -- a carbon bearing a
  formal C=O double bond and one radical electron) is named by replacing
  the parent acid's '-oic acid'/'carboxylic acid' ending with '-oyl'/
  '-carbonyl' -- `has_radical_shape` (`core.py`) already sends any
  radical-bearing molecule here first, ahead of every other module.
  One- and two-carbon chains take the retained 'formyl'/'acetyl'
  (`_retained_acids.py`, P-65.1.7.2.1). Two shapes, both reusing
  `_carboxylic_acid.py`'s own chain-/ring-numbering machinery directly
  rather than duplicating it
  (that module's ring-attached kernels are already suffix/word-
  parameterized and reused the same way by `_ester.py`):
  - An acyclic chain (unbranched or branched), radical carbon fixed at
    C1 like `-COOH`'s own carbon: `_common.py`'s generic
    `name_from_substituents`/`longest_chains`/`substituents_for_chain`
    assemble the name, own_word `"oyl"`. Confirmed worked example
    `hexanoyl (PIN)`, the Blue Book.
  - A ring-attached acyl carbon (P-65.1.7's '-carbonyl'/'benzoyl'
    construction): `_carboxylic_acid.py`'s `_name_ring_attached_carboxyl`/
    `_name_benzo_attached_carboxyl` reused with `suffix="carbonyl"`/
    `word="benzoyl"`. Confirmed worked examples `benzoyl (PIN)`,
    `cyclohexanecarbonyl (PIN)`, the Blue Book.
  The hydroxy-derived radical family (P-71.3.4) is a separate follow-up
  step, not this one.

- P-71.3.2: a radical derived from an amine, imine, or amide
  characteristic group (a single nitrogen bearing one radical electron)
  is named by taking the neutral parent's own name and eliding its final
  'e', then appending 'yl' (`methanamine` -> `methanaminyl`) --
  `_characteristic_group_radical_name` reconstructs that neutral parent
  (one additional explicit hydrogen replacing the radical electron,
  mirroring `_dipole_oxide.py`'s own "strip the dipole atom's charge/
  electron, sanitize, delegate to the plain neutral namer" pattern) and
  delegates to whichever of `_imine.py`/`_amide.py`/`_amine.py` claims
  the reconstructed shape. The same string-transformation mechanism
  `_ylide.py`'s own `amine_name[:-1] + "ium"` already proves works for an
  arbitrary substituted amine name (confirmed worked example
  `methanaminyl (PIN)`, the Blue Book). P-71.3.3's own
  polyamine/polyimine/polyamide multiplicative radicals (two or more
  radical centers on separate characteristic groups) and a divalent
  '-ylidene' version of this same suffix family are each a separate
  follow-up step, not this one.

- P-74.2.2.3.4 (vinyl carbenes): a divalent/trivalent radical at one
  terminus of an otherwise-plain unbranched all-carbon chain that also
  carries exactly one C=C/C#C multiple bond elsewhere on the chain is
  named via the ordinary '-ylidene'/'-ylidyne' suffix combined with the
  chain's own ene/yne locant -- `_vinyl_carbene_name` reuses `_common.py`'s
  generic `name_from_substituents(chain_length, ene_locants, yne_locants,
  own_word, own_locants=[1])` directly (the radical terminus is always
  numbered first, P-74.2.2.3.4's own "low locants... to the suffix"
  rule), the same shared ene/yne-vs-suffix priority logic every other
  suffix module in this project already uses -- no new locant research
  needed. Confirmed worked example `prop-2-en-1-ylidene (PIN)`,
  the Blue Book. The other three "carbene type"
  subtypes in this same Blue Book section (acyl carbenes, imidoyl
  carbenes, imidoyl nitrenes, P-74.2.2.3.1-.3) each need a coexisting
  characteristic-group substituent this project's general "coexisting
  groups" machinery doesn't reach from `_radical.py` yet -- separate,
  larger follow-up steps, not this one.

- P-71.3.4: a radical derived by removing the hydrogen from a hydroxy
  group or its chalcogen analogue (a single terminal O/S/Se radical
  bonded to one plain hydrocarbon substituent) is named additively --
  oxygen uses one of seven retained short names (`_OXYL_RETAINED`,
  confirmed worked examples `methoxyl (PIN)`, `phenoxyl (PIN)`),
  sulfur/selenium have no retained contractions and are always
  systematic (`<R>sulfanyl`/`<R>selanyl`, confirmed worked examples
  `phenylsulfanyl (PIN)`, `methylselanyl (PIN)`). `_chalcogen_radical_
  name` reuses `_substituents.py`'s own `name_branch` directly for R
  (the same plain alkyl/branched/halogenated/phenyl scope `_nitrate_
  ester.py`/`_sulfate.py`'s own R already covers) -- no new substituent-
  naming logic needed. The acid-derived acyloxy radical (R-CO-O•, needs
  a coexisting carbonyl), peroxyl/chalcogen-chain radicals (R-O-O•/
  R-S-S•, a 2-chalcogen chain terminating in a radical instead of
  `_disulfide.py`'s own -SH/-R' termini), and aminoxyl (an amine
  substituent on the radical oxygen) are each a separate follow-up step,
  not this one.

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
  'ethan-1-yl-2-ylidene'-shaped case), a radical center that itself sits
  on an aromatic ring or a polycyclic/spiro skeleton (the acyl case
  above's radical carbon is always exocyclic to its own ring, never on
  it), a radical on a functional group other than the acyl case above
  (P-71.3.1), or coexisting with any heteroatom other than that one acyl
  oxygen, any halogen, charge, or isotopic modification.
"""

from rdkit import Chem

from ._amide import has_amide_shape, name_amide
from ._amine import name_amine
from ._bicyclic import find_bicyclic_core
from ._carboxylic_acid import _name_benzo_attached_carboxyl, _name_ring_attached_carboxyl
from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    linear_branch,
    longest_chains,
    name_from_substituents,
    non_single_bonds,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._imine import has_simple_imine_shape, name_imine
from ._numerals import alkane_name, alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._retained_acids import retained_chain_acid
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain


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


def _acyl_radical_core(mol):
    """(radical_atom, oxygen_atom) if `mol` has P-71.3.1's basic acyl-
    radical charge/bond pattern -- a single monovalent carbon radical
    double-bonded to one terminal oxygen -- else None. Doesn't itself
    constrain the rest of the skeleton (chain branching, ring membership)
    -- `_acyl_radical_acyclic_name`/`_acyl_radical_ring_name` do that."""
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() != 1:
        return None
    (radical,) = radicals
    if (
        radical.GetAtomicNum() != 6
        or radical.GetFormalCharge() != 0
        or radical.GetIsotope() != 0
        or radical.GetIsAromatic()
    ):
        return None

    oxygens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 8]
    if len(oxygens) != 1:
        return None
    (oxygen,) = oxygens
    carbonyl_bond = mol.GetBondBetweenAtoms(radical.GetIdx(), oxygen.GetIdx())
    if (
        carbonyl_bond is None
        or carbonyl_bond.GetBondTypeAsDouble() != 2.0
        or oxygen.GetDegree() != 1
        or oxygen.GetFormalCharge() != 0
        or oxygen.GetIsotope() != 0
    ):
        return None
    return radical, oxygen


def _acyl_radical_acyclic_name(mol):
    """Name if `mol` is P-71.3.1's acyclic (unbranched or branched) acid-
    derived acyl radical shape -- the radical carbon fixed at C1, mirroring
    `_carboxylic_acid.py`'s own -COOH-carbon-at-C1 convention -- else None.
    Reuses `_common.py`'s generic `name_from_substituents`/`longest_chains`/
    `substituents_for_chain` (the same building blocks `_carboxylic_acid.py`'s
    own acyclic path wraps) directly, since P-71.3.1's 'oyl' suffix needs
    the identical chain-numbering/substituent-citation logic, just a
    different trailing suffix word than '-oic acid'."""
    core = _acyl_radical_core(mol)
    if core is None:
        return None
    radical, oxygen = core
    if mol.GetRingInfo().NumRings() != 0:
        return None

    chain_atoms = [a for a in mol.GetAtoms() if a.GetIdx() != oxygen.GetIdx()]
    for atom in chain_atoms:
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            return None
    carbonyl_bond_idx = mol.GetBondBetweenAtoms(radical.GetIdx(), oxygen.GetIdx()).GetIdx()
    for bond in mol.GetBonds():
        if bond.GetIdx() != carbonyl_bond_idx and bond.GetBondTypeAsDouble() != 1.0:
            return None

    graph = carbon_adjacency(mol)
    chains = [chain for chain in longest_chains(graph) if radical.GetIdx() in chain]
    if not chains:
        return None

    best_key = None
    best_name = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != radical.GetIdx():
                # The radical carbon must sit at C1, same as `-COOH`'s own
                # carbon in `_carboxylic_acid.py` -- a direction that
                # doesn't start there is never valid.
                continue
            chain_length = len(candidate)
            substituents = substituents_for_chain(graph, candidate, {}, mol=mol)
            grouped = group_substituents(substituents)
            locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
            name = retained_chain_acid(grouped, chain_length, [], [], 1, "acyl")
            if name is None:
                name = format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "oyl")
            key = (locant_set, citation_locants, name)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def _acyl_radical_ring_name(mol):
    """Name if `mol` is P-71.3.1's ring-attached acyl radical shape (the
    '-carbonyl'/'benzoyl' construction) -- the radical carbon hanging off
    exactly one atom of an otherwise-plain carbocyclic ring -- else None.
    Reuses `_carboxylic_acid.py`'s own ring-numbering kernels directly
    (`_name_ring_attached_carboxyl`/`_name_benzo_attached_carboxyl`, both
    already suffix/word-parameterized and reused the same way by
    `_ester.py`), only the trailing suffix word ('carbonyl'/'benzoyl'
    instead of 'carboxylic acid'/'benzoic acid') differs."""
    core = _acyl_radical_core(mol)
    if core is None:
        return None
    radical, oxygen = core
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = set(ring_info.AtomRings()[0])
    if radical.GetIdx() in ring_atoms:
        return None

    graph = adjacency(mol)
    other_neighbors = [n for n in graph[radical.GetIdx()] if n != oxygen.GetIdx()]
    if other_neighbors == [] or len(other_neighbors) != 1 or other_neighbors[0] not in ring_atoms:
        return None
    ring_atom = other_neighbors[0]
    if any(n for n in graph[ring_atom] if n not in ring_atoms and n != radical.GetIdx()):
        return None
    if any(mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 for atom in ring_atoms):
        return None

    carbonyl_bond_idx = mol.GetBondBetweenAtoms(radical.GetIdx(), oxygen.GetIdx()).GetIdx()
    ring_is_aromatic = all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms)
    if ring_is_aromatic:
        if len(ring_atoms) != 6:
            return None
    else:
        if any(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms):
            return None
        if any(
            bond.GetIdx() != carbonyl_bond_idx and (bond.GetBeginAtomIdx() in ring_atoms or bond.GetEndAtomIdx() in ring_atoms)
            for bond in mol.GetBonds()
            if bond.GetBondTypeAsDouble() != 1.0
        ):
            return None

    halogens = halogen_substituents(mol)
    if ring_is_aromatic:
        return _name_benzo_attached_carboxyl(graph, ring_atoms, radical.GetIdx(), halogens, word="benzoyl", mol=mol)
    return _name_ring_attached_carboxyl(graph, ring_atoms, radical.GetIdx(), halogens, suffix="carbonyl", mol=mol)


def _characteristic_group_radical_name(mol):
    """Name if `mol` is P-71.3.2's amine/imine/amide radical shape -- a
    single nitrogen bearing exactly one radical electron, formal charge
    0, not aromatic -- else None. Reconstructs the neutral parent
    (replacing the radical electron with one additional explicit
    hydrogen, mirroring `_dipole_oxide.py`'s own "strip the dipole atom's
    charge/electron, sanitize, delegate to the plain neutral namer"
    pattern) and delegates to whichever of `_imine.py`/`_amide.py`/
    `_amine.py` claims the reconstructed shape, then transforms that
    name: strip the trailing 'e' and append 'yl' (the same string
    transformation `_ylide.py`'s own `amine_name[:-1] + "ium"` already
    proves works for an arbitrary substituted amine name)."""
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() != 1:
        return None
    (radical,) = radicals
    if radical.GetAtomicNum() != 7 or radical.GetFormalCharge() != 0 or radical.GetIsotope() != 0 or radical.GetIsAromatic():
        return None

    rw = Chem.RWMol(mol)
    n_idx = radical.GetIdx()
    atom = rw.GetAtomWithIdx(n_idx)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(mol.GetAtomWithIdx(n_idx).GetTotalNumHs() + 1)
    atom.SetNumRadicalElectrons(0)
    neutral_mol = rw.GetMol()
    try:
        Chem.SanitizeMol(neutral_mol)
    except (Chem.rdchem.AtomValenceException, Chem.rdchem.KekulizeException):
        return None

    if has_simple_imine_shape(neutral_mol):
        neutral_name = name_imine(neutral_mol)
    elif has_amide_shape(neutral_mol):
        neutral_name = name_amide(neutral_mol)
    else:
        try:
            neutral_name = name_amine(neutral_mol)
        except UnsupportedStructure:
            return None
    return neutral_name[:-1] + "yl"


def _vinyl_carbene_name(mol, radical, valence):
    """Name if `mol` is P-74.2.2.3.4's vinyl-carbene shape -- a divalent
    or trivalent radical at one terminus of an otherwise-plain unbranched
    all-carbon chain that also carries exactly one C=C/C#C multiple bond
    -- else None. Reuses `_common.py`'s generic `name_from_substituents`
    (own_word 'ylidene'/'ylidyne', own_locants=[1] since the radical
    terminus is always numbered first per P-74.2.2.3.4's own "low locants
    ... to the suffix" rule) rather than any new locant-priority logic."""
    if valence not in (2, 3) or radical.GetDegree() != 1:
        return None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            return None
        if atom.GetDegree() > 2:
            return None
    if mol.GetRingInfo().NumRings() != 0:
        return None

    bonds = non_single_bonds(mol)
    if len(bonds) != 1 or bonds[0][2] not in (2.0, 3.0):
        return None
    a, b, order = bonds[0]
    (neighbor,) = [n.GetIdx() for n in radical.GetNeighbors()]
    if mol.GetBondBetweenAtoms(radical.GetIdx(), neighbor).GetBondTypeAsDouble() != 1.0:
        return None

    graph = adjacency(mol)
    n = mol.GetNumAtoms()
    order_list = [radical.GetIdx()]
    previous, current = None, radical.GetIdx()
    while len(order_list) < n:
        next_atoms = [atom for atom in graph[current] if atom != previous]
        if not next_atoms:
            return None
        previous, current = current, next_atoms[0]
        order_list.append(current)
    position_of = {atom: i + 1 for i, atom in enumerate(order_list)}
    bond_locant = min(position_of[a], position_of[b])
    own_word = "ylidene" if valence == 2 else "ylidyne"
    ene_locants = [bond_locant] if order == 2.0 else []
    yne_locants = [bond_locant] if order == 3.0 else []
    return name_from_substituents(n, ene_locants, yne_locants, own_word, own_locants=[1])


# P-71.3.4's own seven retained short 'oxyl' names, minus 'aminoxyl'
# (out of scope, see module docstring) -- only the oxygen case has
# retained contractions; sulfur/selenium are always systematic.
_OXYL_RETAINED = {
    "methyl": "methoxyl",
    "ethyl": "ethoxyl",
    "propyl": "propoxyl",
    "butyl": "butoxyl",
    "tert-butyl": "tert-butoxyl",
    "phenyl": "phenoxyl",
    "amino": "aminoxyl",
}
_CHALCOGEN_RADICAL_SUFFIXES = {8: "oxyl", 16: "sulfanyl", 34: "selanyl"}


def _chalcogen_radical_name(mol):
    """Name if `mol` is P-71.3.4's single-chalcogen terminal radical shape
    -- one O/S/Se atom bearing exactly one radical electron, bonded to
    exactly one carbon root, the rest of the molecule matching
    `name_branch`'s own plain-hydrocarbon scope -- else None."""
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons() != 0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() != 1:
        return None
    (radical,) = radicals
    suffix = _CHALCOGEN_RADICAL_SUFFIXES.get(radical.GetAtomicNum())
    if suffix is None or radical.GetFormalCharge() != 0 or radical.GetIsotope() != 0 or radical.GetDegree() != 1:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    graph = adjacency(mol)
    (root,) = graph[radical.GetIdx()]
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    try:
        name, compound = name_branch(graph, root, radical.GetIdx(), halogens, aromatic_atoms, mol=mol)
    except UnsupportedStructure:
        return None
    if compound:
        return None

    if radical.GetAtomicNum() == 8 and name in _OXYL_RETAINED:
        return _OXYL_RETAINED[name]
    return name + suffix


def name_radical(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    from ._radical_hetero import hetero_radical_name

    hetero_name = hetero_radical_name(mol)
    if hetero_name is not None:
        return hetero_name

    characteristic_group_name = _characteristic_group_radical_name(mol)
    if characteristic_group_name is not None:
        return characteristic_group_name

    chalcogen_radical_name = _chalcogen_radical_name(mol)
    if chalcogen_radical_name is not None:
        return chalcogen_radical_name

    acyl_name = _acyl_radical_acyclic_name(mol)
    if acyl_name is not None:
        return acyl_name
    acyl_name = _acyl_radical_ring_name(mol)
    if acyl_name is not None:
        return acyl_name

    radicals = [atom for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() != 0]
    if len(radicals) == 1 and radicals[0].GetNumRadicalElectrons() in (2, 3):
        vinyl_carbene_name = _vinyl_carbene_name(mol, radicals[0], radicals[0].GetNumRadicalElectrons())
        if vinyl_carbene_name is not None:
            return vinyl_carbene_name
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
    return _name_von_baeyer_or_spiro_radical(mol, radical, valence)


def _radical_suffix(yl_name: str, valence: int) -> str:
    """P-71.2.2.1: the divalent/trivalent suffix is formed by appending
    'idene'/'idyne' after the '-yl' name (not replacing it), e.g.
    'methyl' -> 'methylidene'/'methylidyne'."""
    if valence == 1:
        return yl_name
    return yl_name + ("idene" if valence == 2 else "idyne")


def _name_von_baeyer_or_spiro_radical(mol, radical, valence) -> str:
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single radical suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`
    (the shared mechanism `_alcohol.py`/`_carbenium.py` already use, the
    latter structurally parallel to this one -- both cite a single
    free-valence/charge locant on the ring skeleton with no suffix
    elision needed). The '-yl' name is computed first, then P-71.2.2.1's
    'idene'/'idyne' suffix is appended on top for a divalent/trivalent
    radical center, exactly as `_radical_suffix` already does for the
    chain/monocyclic-ring cases."""
    stereo = specified_stereocenters(mol)

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
        yl_name = name_von_baeyer_suffix(
            mol, radical.GetIdx(), set(), "yl", "radical", bicyclic_core, polycyclic_core, von_baeyer_ring_count,
            stereo=stereo,
        )
        return _radical_suffix(yl_name, valence)

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        yl_name = name_monospiro_suffix(mol, radical.GetIdx(), set(), "yl", "radical", spiro_atom, stereo=stereo)
        return _radical_suffix(yl_name, valence)

    raise UnsupportedStructure(
        "polycyclic and fused-ring radicals are not supported yet (P-23/"
        "P-25 numbering integration with a suffix group is future work)"
    )


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
