"""Naming of simple Group 1/2 organometallics (a single lithium, sodium,
potassium, magnesium, or calcium atom bearing 1 (Group 1) or 1-2 (Group 2)
alkyl/phenyl substituents), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-69.3 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  organometallic compounds of Groups 1 and 2 are named additively -- the
  ligand name(s) are cited directly before the bare metal element name,
  with no parent-hydride ('-ane') construction at all (unlike P-69.1's
  Group 13-16 mechanism, `_group13_hydride.py`). Confirmed worked example
  (`tmp/bluebook/P6a.txt` line 8843-8845): `[LiMe]` -> 'methyllithium'.
  Doubly-substituted Group 2 case (unconfirmed by a worked example in the
  cached excerpt) follows the same ordinary P-14.2.1 multiplying-prefix
  convention already used throughout this project for an otherwise
  unremarkable identical-substituent citation, e.g. 'diethylmagnesium'.
- PubChem is not usable to verify this shape: its own registered entries
  for these compounds use a fully-dissociated ionic representation
  (`[Li+].[CH3-]`) rather than the covalent single-molecule structure the
  Blue Book itself names, and its auto-generated names for the covalent
  form are simply wrong ('lithium carbanide', 'magnesium ethane' -- not
  real nomenclature, confirmed via PubChem PUG REST structure search,
  CIDs 2724049/11183/13351308). Verification here is therefore against
  the primary source's own worked example plus mechanical, self-
  consistent structural reasoning, the same situation `_radical.py`'s
  polycyclic/spiro step was already in.
- This module is structurally similar to `_group13_hydride.py`
  (mononuclear metal parent + substituent citation) but is NOT the same
  mechanism: no parent-hydride stem is constructed here (Group 1/2
  elements have no borane/alumane-style '-ane' hydride name in this
  additive scheme), and the substituent-name-then-plain-multiplying-
  prefix pattern is the ester-word style already established in
  `_phosphate.py` (plain 'di'/'tri' concatenated directly, enclosing
  marks only when the substituent's own name starts with a locant digit)
  -- not `format_mononuclear_prefixes`'s substituent-prefix-on-a-parent-
  hydride style.

Scope: a plain/branched alkyl substituent or a plain phenyl substituent,
1 substituent for Group 1 (Li/Na/K), 1-2 identical substituents for
Group 2 (Mg/Ca). Also: a Group 2 metal bearing exactly one organic
substituent and one halogen (the Grignard-reagent-shaped R-M-X case) --
P-69.3's own worked example (`tmp/bluebook/P6a.txt` lines 8859-8863):
`[MgMe]I` -> 'methylmagnesium iodide (compositional name; the formally
electropositive component named by additive nomenclature)'. This is a
*different* citation style from the plain R-M/R2-M case above: the
organic-group-name-fused-to-metal word (identical to the single-
substituent case), then a separate, space-delimited word for the halide's
anion name (`HALIDE_WORDS`, the same `fluoride`/`chloride`/`bromide`/
`iodide` dict `_acyl_halide.py` already uses), not the `fluoro`/`chloro`/
... substituent-prefix form. Group 1 metals are monovalent, so an R-M-X
shape isn't chemically possible for them (would need valence 2).
PubChem's own registered entries for this shape are the same
fully-dissociated-ionic, auto-generated-name-not-real-nomenclature
situation already documented above for the plain R-M case -- their
existence (e.g. methylmagnesium iodide, ethylmagnesium bromide,
phenylmagnesium chloride are all well-known, commonly-registered Grignard
reagents) confirms these are real compounds, but verification is against
the primary source's own worked example plus structural consistency, not
an exact PubChem name match.

Explicitly out of scope (raise `UnsupportedStructure`): a halogen on a
Group 1 metal, two organic substituents alongside a halogen, two
different organic substituents on a Group 2 metal (with or without a
halogen), more than one metal atom, an unsaturated substituent, an
aromatic substituent other than plain phenyl, charged/isotopically
modified atoms, and coexistence with any other heteroatom.
"""

from rdkit import Chem

from ._common import (
    HALIDE_WORDS,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
from ._numerals import multiplying_prefix
from ._substituents import name_branch

_GROUP1_ELEMENTS = {3: "lithium", 11: "sodium", 19: "potassium"}
_GROUP2_ELEMENTS = {12: "magnesium", 20: "calcium"}
_ELEMENT_NAMES = {**_GROUP1_ELEMENTS, **_GROUP2_ELEMENTS}


def has_group1_2_organometallic_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in _ELEMENT_NAMES for atom in mol.GetAtoms())


def name_group1_2_organometallic(mol) -> str:
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in _ELEMENT_NAMES]
    if len(metal_atoms) != 1:
        raise UnsupportedStructure("more than one Group 1/2 metal atom is not supported yet")
    (metal,) = metal_atoms
    atomic_num = metal.GetAtomicNum()
    element_name = _ELEMENT_NAMES[atomic_num]

    if metal.GetFormalCharge() != 0 or metal.GetIsotope() != 0:
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")

    graph = adjacency(mol)
    metal_neighbors = list(metal.GetNeighbors())
    halogen_neighbors = [n for n in metal_neighbors if n.GetAtomicNum() in HALOGEN_PREFIXES]
    has_halide = bool(halogen_neighbors)

    if has_halide and atomic_num in _GROUP1_ELEMENTS:
        raise UnsupportedStructure(
            "a halogen on a Group 1 metal is not a valid organometallic "
            "shape for this module (a monovalent metal has no room for "
            "both an organic group and a halide)"
        )
    if len(halogen_neighbors) > 1:
        raise UnsupportedStructure("more than one halogen on the same metal atom is not supported yet")

    # A neutral Group 1 metal is monovalent (exactly 1 organic
    # substituent) and a neutral Group 2 metal is divalent -- either a
    # closed-shell R2M species (2 organic substituents) or, if one
    # neighbor is a halogen, the Grignard-shaped R-M-X compositional case
    # (1 organic + 1 halide). A Group 2 metal with only 1 organic
    # substituent and no halide would be an open-shell radical, not a
    # real compound.
    required_substituents = 1 if atomic_num in _GROUP1_ELEMENTS else 2
    if metal.GetDegree() != required_substituents:
        raise UnsupportedStructure(
            f"a {element_name} atom must have exactly {required_substituents} "
            "organic/halide substituent(s) to be supported here"
        )

    roots = {n.GetIdx() for n in metal_neighbors} - {n.GetIdx() for n in halogen_neighbors}
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots)

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        is_metal_halogen = has_halide and idx == halogen_neighbors[0].GetIdx()
        if atom.GetAtomicNum() not in (atomic_num, 6) and not is_metal_halogen:
            raise UnsupportedStructure(
                f"heteroatoms other than the {element_name} atom itself are "
                "not supported yet (see P-69.3)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and idx not in phenyl_atoms:
            raise UnsupportedStructure(
                "an aromatic substituent other than a plain phenyl group is "
                "out of scope for this module"
            )
    all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if all_ring_atoms - phenyl_atoms:
        raise UnsupportedStructure(
            "a ring other than a plain phenyl substituent directly on the "
            "metal is out of scope for this module"
        )
    non_ring_unsaturation = [b for b in non_single_bonds(mol) if b[0] not in phenyl_atoms and b[1] not in phenyl_atoms]
    if non_ring_unsaturation:
        raise UnsupportedStructure("an unsaturated substituent is out of scope for this module (see P-69.3)")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    names = []
    for root in roots:
        if root in phenyl_atoms:
            names.append("phenyl")
            continue
        names.append(name_branch(graph, root, metal.GetIdx(), {}, mol=mol)[0])

    (first_name, *rest) = names
    if any(name != first_name for name in rest):
        raise UnsupportedStructure(
            "two different substituents on the same Group 2 metal is not "
            "supported yet"
        )

    if has_halide:
        # P-69.3's compositional-name worked example (`[MgMe]I` ->
        # 'methylmagnesium iodide') cites the organic group fused to the
        # metal exactly like the plain single-substituent case, then a
        # separate, space-delimited word for the halide's own anion name
        # -- not the `fluoro`/`chloro`/... substituent-prefix form.
        halide_word = HALIDE_WORDS[halogen_neighbors[0].GetAtomicNum()]
        return f"{first_name}{element_name} {halide_word}"

    if len(names) == 1:
        return f"{first_name}{element_name}"

    # P-69.3's additive naming cites the substituent(s) directly before
    # the bare metal name, the same word-citation style `_phosphate.py`
    # already established for P-67.1.3.2's ester names: a plain
    # multiplying prefix concatenated directly, enclosing marks only when
    # the substituent's own name starts with a locant digit.
    needs_enclosure = first_name[0].isdigit()
    prefix = multiplying_prefix(len(names), compound=needs_enclosure)
    group = f"({first_name})" if needs_enclosure else first_name
    return f"{prefix}{group}{element_name}"
