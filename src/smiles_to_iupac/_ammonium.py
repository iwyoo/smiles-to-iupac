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
- P-91.3/P-92 stereocenters: the
  1-3-substituted paths delegate to `name_amine` on the neutralized
  molecule, so they inherited full R/S support automatically once
  `_amine.py` gained it. The quaternary (4-substituted) path calls
  `_name_acyclic_secondary_tertiary_amine` directly and now passes
  `specified_stereo_elements(mol)` through the same way, e.g.
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

import re

from rdkit import Chem

from ._amine import CATIONIC_AMINES, _name_acyclic_secondary_tertiary_amine, name_amine
from ._common import UnsupportedStructure, non_single_bonds, specified_stereo_elements
from ._hetero_prefixes import is_functional_carbon
from ._onium_prefixes import onium_name


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


def has_polyammonium_shape(mol) -> bool:
    """Two or more +1 ammonium nitrogens (each bonded to carbon only, singly, with hydrogens making up four bonds) and
    no other charged atom."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) < 2 or any(a.GetAtomicNum() != 7 or a.GetFormalCharge() != 1 or a.GetIsotope() for a in charged):
        return False
    return all(
        a.GetDegree() + a.GetTotalNumHs() == 4
        and not a.IsInRing()
        and not a.GetIsAromatic()
        and all(
            n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in a.GetNeighbors()
        )
        for a in charged
    )


_MULTIPLIED_AMINIUM = {"di": "bis", "tri": "tris", "tetra": "tetrakis", "penta": "pentakis", "hexa": "hexakis"}


def name_polyammonium(mol) -> str:
    """Several ammonium nitrogens on one parent are named as a multiplied 'aminium' suffix, 'bis(aminium)', 'tris(...)'
    (P-73.1.2.1): each hydrogen-bearing cation is neutralised for the amine namer, which numbers the N-substituents."""
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetAtomicNum() == 6 and is_functional_carbon(mol, a.GetIdx()) for a in mol.GetAtoms()):
        raise UnsupportedStructure("a carbonyl-type group beside the cations is cited as a prefix by the chain engine")
    neutral = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() and atom.GetTotalNumHs():
            copy = neutral.GetAtomWithIdx(atom.GetIdx())
            copy.SetFormalCharge(0)
            copy.SetNoImplicit(True)
            copy.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    token = CATIONIC_AMINES.set(True)
    try:
        name = name_amine(base)
    finally:
        CATIONIC_AMINES.reset(token)
    match = re.search(r"(di|tri|tetra|penta|hexa)(amine|aniline)$", name)
    if match is None:
        raise UnsupportedStructure("the polycation is not named as a multiple amine parent")
    return f"{name[:match.start()]}{_MULTIPLIED_AMINIUM[match.group(1)]}(aminium)"


def _heteroatom_acyl(mol, carbon) -> bool:
    """A carbonyl carbon whose other neighbour is a heteroatom (carboxy, alkoxycarbonyl, carbamoyl, halocarbonyl): an acid
    derivative on the cationic nitrogen, cited as a prefix of 'azanium' because no amidium suffix exists for it."""
    if carbon.GetAtomicNum() != 6 or carbon.IsInRing():
        return False
    bonds = {n.GetIdx(): mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in carbon.GetNeighbors()}
    oxo = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 8 and bonds[n.GetIdx()] == 2.0 and n.GetDegree() == 1]
    hetero = [
        n
        for n in carbon.GetNeighbors()
        if n.GetAtomicNum() in (7, 8, 9, 16, 17, 35, 53) and bonds[n.GetIdx()] == 1.0 and not n.GetFormalCharge()
    ]
    return len(oxo) == 1 and carbon.GetDegree() == 3 and len(hetero) == 1


def name_ammonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() > 0
    ]
    if len(charged_nitrogens) != 1:
        raise UnsupportedStructure("several cationic nitrogens are named as a polyammonium, not here")
    (nitrogen,) = charged_nitrogens
    if nitrogen.GetFormalCharge() != 1 or nitrogen.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "ammonium nitrogen is supported (P-73.1.1.2)"
        )
    if any(
        atom.GetIdx() != nitrogen.GetIdx() and (atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0)
        for atom in mol.GetAtoms()
    ):
        # A second charged atom elsewhere (e.g. a carbanion) makes this a
        # dipolar/ylide species (P-74.2), not a plain ammonium salt cation
        # -- `has_ammonium_shape` only inspects the nitrogen's own local
        # bonding, so it matches both shapes identically; without this
        # check the ylide's second charge center was silently dropped.
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")

    if any(_heteroatom_acyl(mol, n) for n in nitrogen.GetNeighbors()) and not nitrogen.IsInRing():
        return onium_name(mol, nitrogen, "azanium")

    if any(atom.GetAtomicNum() == 6 and is_functional_carbon(mol, atom.GetIdx()) for atom in mol.GetAtoms()):
        raise UnsupportedStructure("a carbonyl-type group beside the cation is cited as a prefix by the chain engine")

    if nitrogen.IsInRing():
        raise UnsupportedStructure("a ring nitrogen cation is named on its ring parent hydride, not as an acyclic ammonium")

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

    if degree == 4 and mol.GetRingInfo().NumRings():
        raise UnsupportedStructure("a quaternary ammonium on or beside a ring is named by the ring-aware amine engine")
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
        stereo = specified_stereo_elements(mol)
        amine_name = _name_acyclic_secondary_tertiary_amine(mol, nitrogen.GetIdx(), n_carbons, bonds, stereo)
        return amine_name[:-1] + "ium"

    neutral_rw = Chem.RWMol(mol)
    neutral_nitrogen = neutral_rw.GetAtomWithIdx(nitrogen.GetIdx())
    neutral_nitrogen.SetFormalCharge(0)
    neutral_nitrogen.SetNoImplicit(True)
    neutral_nitrogen.SetNumExplicitHs(3 - degree)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    if sum(a.GetAtomicNum() == 7 for a in mol.GetAtoms()) > 1:
        raise UnsupportedStructure("a neutral amino group beside the ammonium group is a prefix of the aminium name")
    if len(neighbors) > 1 and not mol.GetRingInfo().NumRings() and _forced_prefixes_active():
        bonds = [b for b in non_single_bonds(neutral_mol) if b[2] in (2.0, 3.0)]
        amine_name = _name_acyclic_secondary_tertiary_amine(
            neutral_mol, nitrogen.GetIdx(), tuple(n.GetIdx() for n in neighbors), bonds, specified_stereo_elements(mol)
        )
        return amine_name[:-1] + "ium"
    try:
        amine_name = name_amine(neutral_mol)
        if not _cites_forced_prefixes(amine_name):
            raise UnsupportedStructure("the amine namer dropped a prefix supplied for a further onium group")
    except UnsupportedStructure:
        from ._polyfunctional import name_polyfunctional

        amine_name = name_polyfunctional(neutral_mol)
        if not amine_name.endswith("amine") or not _cites_forced_prefixes(amine_name):
            raise
    return amine_name[:-1] + "ium"


def _forced_prefixes_active():
    from ._substituents import FORCED_BRANCH_NAMES

    return FORCED_BRANCH_NAMES.get() is not None


def _cites_forced_prefixes(name):
    from ._substituents import FORCED_BRANCH_NAMES

    forced = FORCED_BRANCH_NAMES.get()
    return forced is None or all(text in name for text, _ in forced[1].values())
