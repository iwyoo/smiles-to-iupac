"""Naming of an acyclic carboxylic acid chain bearing a skeletal-carbon
isotope (12C/13C/14C) and/or an isotopically labeled carboxyl oxygen
(17O/18O), combined with the existing '-oic acid' suffix, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-82.5.1: numbering is unchanged from the isotopically-unmodified
  compound, so the winning chain numbering is exactly whatever
  `_carboxylic_acid.py`'s own acyclic path would pick for the
  isotopically-unmodified structure, reused directly via
  `_best_acyclic_carboxylic_acid_candidate` on a neutralized (isotope-
  stripped) copy of the molecule (mirrors `_isotope_alcohol.py`/
  `_isotope_ketone.py`'s identical approach).
- This case has a direct, exact-match confirmed worked example (not
  analogy from a different suffix): `tmp/bluebook/P8.txt` line 327,
  '(1-14C)pentan(3H)oic acid (PIN)' -- a carbon isotope on the carboxyl
  carbon itself (C1) placed at the front of the whole name, combined with
  a non-carbon isotope (tritium, on the acid's own exchangeable proton)
  placed immediately before "oic acid", as two separate parenthetical
  groups in the same name, not merged into one. This resolves both the
  single-isotope placements and the both-together case (unlike
  `_isotope_ketone.py`'s carbonyl case, which had no matching worked
  example for both together and defers that combination entirely).

Scope, deliberately narrow (P82-WS2 M1 step 3, mirrors #789/#803's own
pilot shape): an acyclic carboxylic acid chain, with halogen substituents
allowed (mirrors `_carboxylic_acid.py`'s own scope), bearing a skeletal
carbon isotope and/or a carboxyl-oxygen isotope. Explicitly out of scope
(raise `UnsupportedStructure`): any ring, any chain unsaturation,
deuterium or any other heteroatom besides the carboxyl oxygens and
halogens, more than one carboxylic acid group, a specified stereocenter,
mixing different carbon-isotope nuclides, and a carbon isotope on a
branch off the principal chain.
"""

import re

from rdkit import Chem

from ._carboxylic_acid import _best_acyclic_carboxylic_acid_candidate, _validate_and_collect_carboxyls
from ._common import UnsupportedStructure, non_single_bonds, specified_stereocenters
from ._isotope import _CARBON_ISOTOPES

_OXYGEN_ISOTOPES = {17, 18}
_RETAINED_STEM = re.compile(r"(?:formic|acetic) acid$")


def has_isotope_carboxylic_acid_shape(mol) -> bool:
    """True if `mol` contains some isotope label and looks like it might
    be an acyclic carboxylic acid (a carbon bearing both a doubly-bonded
    and a singly-bonded, one-H oxygen) -- deliberately loose (see module
    docstring for the real scope); `name_isotope_carboxylic_acid` raises
    for every detail this doesn't cover. Used by `core.py` to route here
    before `_isotope.py`'s own plain chain/methane path."""
    if not any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
        return False
    if mol.GetRingInfo().NumRings() > 0:
        return False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        has_carbonyl = any(
            o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for o in oxygens
        )
        has_hydroxyl = any(
            o.GetDegree() == 1
            and o.GetTotalNumHs() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            for o in oxygens
        )
        if has_carbonyl and has_hydroxyl:
            return True
    return False


def name_isotope_carboxylic_acid(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("a ring alongside an isotope-labeled carboxylic acid is not supported yet")
    other_unsaturation = [
        (a, b, order)
        for a, b, order in non_single_bonds(mol)
        if mol.GetAtomWithIdx(a).GetAtomicNum() != 8 and mol.GetAtomWithIdx(b).GetAtomicNum() != 8
    ]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation combined with an isotope-labeled carboxylic acid is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an isotope-labeled carboxylic acid is not supported yet"
        )

    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_isotope_positions = [c.GetIdx() for c in carbons if c.GetIsotope() != 0]
    carbon_isotopes = {mol.GetAtomWithIdx(c).GetIsotope() for c in carbon_isotope_positions}
    for isotope in carbon_isotopes:
        if isotope not in _CARBON_ISOTOPES:
            raise UnsupportedStructure(
                "only 12C/13C/14C skeletal carbon isotopes are supported (P-82.2.1)"
            )
    if len(carbon_isotopes) > 1:
        raise UnsupportedStructure(
            "mixing different skeletal carbon isotope nuclides (e.g. 13C "
            "and 14C together) is not supported yet"
        )
    has_carbon_isotope = bool(carbon_isotope_positions)

    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetIsotope(0)
    neutral_mol = rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    carboxyl_carbons, carboxyl_oxygens, extra_hydroxyls = _validate_and_collect_carboxyls(neutral_mol)
    if extra_hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside an isotope-labeled carboxylic acid is not supported yet"
        )
    if len(carboxyl_carbons) != 1:
        raise UnsupportedStructure(
            "more than one carboxylic acid group alongside an isotope label is not supported yet"
        )

    labeled_carboxyl_oxygens = [o for o in carboxyl_oxygens if mol.GetAtomWithIdx(o).GetIsotope() != 0]
    oxygen_isotopes = {mol.GetAtomWithIdx(o).GetIsotope() for o in labeled_carboxyl_oxygens}
    for isotope in oxygen_isotopes:
        if isotope not in _OXYGEN_ISOTOPES:
            raise UnsupportedStructure("only 17O/18O oxygen isotopes are supported (P-82.2.1)")
    if len(labeled_carboxyl_oxygens) > 1:
        raise UnsupportedStructure(
            "more than one isotopically labeled carboxyl oxygen at once is not supported yet"
        )
    has_oxygen_isotope = bool(labeled_carboxyl_oxygens)

    stray_isotopes = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetIsotope() != 0
        and atom.GetIdx() not in carbon_isotope_positions
        and atom.GetIdx() not in labeled_carboxyl_oxygens
    ]
    if stray_isotopes:
        raise UnsupportedStructure(
            "an isotope label outside the carboxyl oxygens and skeletal "
            "carbons is not supported yet"
        )
    if not has_carbon_isotope and not has_oxygen_isotope:
        raise UnsupportedStructure("no isotopically labeled atom found")

    best_name, best_position_of = _best_acyclic_carboxylic_acid_candidate(
        neutral_mol, carboxyl_carbons, carboxyl_oxygens, set(), ()
    )

    name = best_name
    descriptors = []
    if has_carbon_isotope:
        if any(c not in best_position_of for c in carbon_isotope_positions):
            raise UnsupportedStructure(
                "a carbon isotope on a substituent branch rather than the "
                "principal chain is not supported yet"
            )
        carbon_isotope_locants = sorted(best_position_of[c] for c in carbon_isotope_positions)
        nuclide = next(iter(carbon_isotopes))
        symbol = f"{nuclide}C" + (str(len(carbon_isotope_positions)) if len(carbon_isotope_positions) > 1 else "")
        locants_str = ",".join(str(loc) for loc in carbon_isotope_locants)
        descriptors.append(f"{locants_str}-{symbol}")
    retained = _RETAINED_STEM.search(name)
    if retained:
        if has_oxygen_isotope:
            (oxygen_isotope,) = oxygen_isotopes
            descriptors.append(f"{oxygen_isotope}O")
        return name[: retained.start()] + f"({','.join(descriptors)})" + name[retained.start() :]
    if descriptors:
        name = f"({descriptors[0]}){name}"
    if has_oxygen_isotope:
        (oxygen_isotope,) = oxygen_isotopes
        name = name[:-8] + f"({oxygen_isotope}O)oic acid"

    return name
