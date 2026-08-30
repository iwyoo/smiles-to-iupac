"""Naming of simple ammonium cations (the '-aminium'/'azanium' suffix
family, P-73.1.1.2), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium'.
  Applied to a primary amine (`_amine.py`'s '-amine' names, P-33.1), this
  gives e.g. 'methanamine' -> 'methanaminium', 'propan-1-amine' ->
  'propan-1-aminium'. Unsubstituted NH4+ is the hydron-added form of
  'azane' (NH3) itself, giving 'azanium' -- both the Blue Book PIN and
  PubChem's own auto-generated name agree on this one (an unusually
  complete match for this project).
- This 'aminium' derivation -- not the alternative 'azanium' substitutive
  style (a multiplying prefix directly on the 'azanium' parent, e.g.
  'tetramethylazanium') -- is the PIN all the way up through a
  *quaternary* ammonium cation too: the primary source's own Table 7.3
  worked example for (CH3)4N+ lists 'tetramethylazanium' right next to
  the PIN 'N,N,N-trimethylmethanaminium', confirming the 'azanium' form
  is explicitly non-preferred (unlike `_phosphonium.py`/`_sulfonium.py`,
  where the analogous direct-substituent style *is* the PIN -- nitrogen
  is the odd one out here, confirmed directly from the primary source
  rather than assumed by analogy).
- Secondary/tertiary ammonium (nitrogen bonded to 2-3 carbons) reuses
  `_amine.py`'s existing secondary/tertiary amine naming wholesale: the
  cation is neutralized (formal charge 0, one more hydrogen than carbon
  neighbors leaves), named via `name_amine`, and the resulting name's
  terminal 'e' is replaced with 'ium' -- e.g. diethylamine's own
  'N-ethylethanamine' -> 'N-ethylethanaminium' (diethylammonium), and
  triethylamine's 'N,N-diethylethanamine' -> 'N,N-diethylethanaminium'
  (confirmed against the primary source's own
  'N,N-diethylethanaminium hydrogen sulfate (PIN)' salt example).
- Quaternary ammonium (nitrogen bonded to 4 carbons) has no neutral
  counterpart to derive an '-ium' name from at all (a neutral nitrogen
  can carry at most 3 substituents) -- instead this module calls
  `_amine.py`'s own `_name_acyclic_secondary_tertiary_amine` directly on
  the still-charged molecule (that helper's internals are agnostic to
  the nitrogen's formal charge; see its own docstring), producing the
  same 'parent chain + N,N,N-prefix' shape as an '-amine' name would,
  then applies the same terminal 'e' -> 'ium' swap.
- P-91.3/P-92 stereocenters (`tasks/ammonium-stereocenter-naming.md`): the
  1-3-substituted paths delegate to `name_amine` on the neutralized
  molecule, so they inherited full R/S support automatically once
  `_amine.py` gained it. The quaternary (4-substituted) path calls
  `_name_acyclic_secondary_tertiary_amine` directly and now passes
  `specified_stereocenters(mol)` through the same way, e.g.
  '(2R)-N,N,N-trimethylbutan-2-aminium'.

Explicitly out of scope (raise `UnsupportedStructure`):
- Anything `_amine.py` itself would reject for the neutralized (or,
  for the quaternary case, still-charged) molecule -- a ring, a
  branched/unsaturated N-substituent, a halogen substituent coexisting
  with a secondary/tertiary/quaternary nitrogen, etc.
- Any ammonium nitrogen not shaped like a simple, singly-charged R-NH3+,
  NH4+, or a nitrogen bonded to 2-4 carbons (e.g. formal charge other
  than +1, more than one charged atom, isotopic modification, a nitrogen
  double/triple-bonded to carbon).
- Any other cation-forming parent (oxonium R3O+, sulfonium R3S+, ...) --
  each is a separate P-73 subsection with its own derivation rule.
"""

from rdkit import Chem

from ._amine import _name_acyclic_secondary_tertiary_amine, name_amine
from ._common import UnsupportedStructure, non_single_bonds, specified_stereocenters


def has_ammonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged nitrogen shaped
    like a genuine ammonium: NH4+, or a nitrogen singly bonded to 1-4
    carbons with the rest of its valence as hydrogens. Used by `core.py` to
    route here before `_amine.py`, which rejects any charged atom outright.

    Deliberately narrower than "any charged nitrogen exists": a nitro
    group's canonical Lewis structure (-[N+](=O)[O-]) and an isocyanide's
    (-[N+]#[C-]) both also carry a formally charged nitrogen, but neither
    is singly-bonded-to-carbon/all-remaining-valence-as-H shaped, so this
    predicate correctly leaves them to their own modules."""
    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1
    ]
    if len(charged_nitrogens) != 1:
        return False
    nitrogen = charged_nitrogens[0]
    if nitrogen.GetIsotope() != 0:
        return False
    degree = nitrogen.GetDegree()
    if degree == 0:
        return nitrogen.GetTotalNumHs() == 4
    if degree > 4:
        return False
    if nitrogen.GetTotalNumHs() != 4 - degree:
        return False
    return all(
        neighbor.GetAtomicNum() == 6
        and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 1.0
        for neighbor in nitrogen.GetNeighbors()
    )


def name_ammonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() != 0
    ]
    (nitrogen,) = charged_nitrogens
    if nitrogen.GetFormalCharge() != 1 or nitrogen.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "ammonium nitrogen is supported (P-73.1.1.2)"
        )

    degree = nitrogen.GetDegree()
    if degree == 0:
        if nitrogen.GetTotalNumHs() != 4:
            raise UnsupportedStructure("an unsupported unsubstituted ammonium shape")
        return "azanium"

    neighbors = list(nitrogen.GetNeighbors())
    if degree > 4 or nitrogen.GetTotalNumHs() != 4 - degree:
        raise UnsupportedStructure(
            "only NH4+, or a nitrogen bonded to 1-4 carbons with the rest "
            "of its valence as hydrogens, is supported (P-73.1.1.2)"
        )
    if any(neighbor.GetAtomicNum() != 6 for neighbor in neighbors):
        raise UnsupportedStructure("an ammonium nitrogen must be attached only to carbon atoms")
    if any(
        mol.GetBondBetweenAtoms(nitrogen.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() != 1.0
        for neighbor in neighbors
    ):
        raise UnsupportedStructure("the ammonium nitrogen must be singly bonded to carbon")

    if degree == 4:
        # No neutral nitrogen can carry 4 substituents, so there's no
        # '-amine' name to derive 'ium' from by neutralizing -- instead
        # this reuses _amine.py's own secondary/tertiary machinery
        # directly on the still-charged molecule (see module docstring).
        all_non_single = non_single_bonds(mol)
        bonds = [b for b in all_non_single if b[2] in (2.0, 3.0)]
        if len(bonds) != len(all_non_single):
            raise UnsupportedStructure(
                "a bond order other than single, double, or triple is not "
                "supported (see P-31.1.1.1)"
            )
        n_carbons = tuple(neighbor.GetIdx() for neighbor in neighbors)
        stereo = specified_stereocenters(mol)
        amine_name = _name_acyclic_secondary_tertiary_amine(mol, nitrogen.GetIdx(), n_carbons, bonds, stereo)
        return amine_name[:-1] + "ium"

    neutral_rw = Chem.RWMol(mol)
    neutral_nitrogen = neutral_rw.GetAtomWithIdx(nitrogen.GetIdx())
    neutral_nitrogen.SetFormalCharge(0)
    neutral_nitrogen.SetNoImplicit(True)
    neutral_nitrogen.SetNumExplicitHs(3 - degree)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    amine_name = name_amine(neutral_mol)
    return amine_name[:-1] + "ium"
