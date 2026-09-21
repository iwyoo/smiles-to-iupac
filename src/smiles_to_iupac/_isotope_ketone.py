"""Naming of an acyclic ketone chain bearing an isotopically labeled
carbonyl oxygen (17O/18O) XOR a skeletal-carbon isotope (12C/13C/14C),
combined with the existing '-one' suffix, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-82.5.1: numbering is unchanged from the isotopically-unmodified
  compound (isotopic substitution is considered last among P-14.4's
  numbering-priority criteria), so the winning chain numbering is exactly
  whatever `_ketone.py`'s own acyclic path would pick for the
  isotopically-unmodified structure, reused directly via
  `_best_acyclic_ketone_candidate` on a neutralized (isotope-stripped)
  copy of the molecule rather than re-derived here (mirrors
  `_isotope_alcohol.py`'s identical approach).
- Placement depends on which atom carries the isotope, confirmed by
  `tmp/bluebook/P8.txt`'s own worked examples:
  - A skeletal-carbon isotope (not the suffix-defining atom) goes at the
    *front of the whole name*, per P-82.2.5's own general rule ("In a
    name consisting of one word, the isotopic descriptor is placed
    before the name ... preferred to ... placing the descriptor before
    the implied name of the characteristic group"): '1-phenyl(1,2-13C2)-
    ethan-1-one (PIN)' (line 155).
  - An isotope on the suffix-defining atom itself (the carbonyl oxygen
    here) goes immediately before the suffix word, after the suffix's
    own locant, by direct analogy with the identical alcohol case:
    '1-(aminomethyl)cyclopentan-1-(18O)ol (PIN)' (line 187) -- the
    oxygen IS the suffix-defining atom there, unlike the ethanol/
    ethanone carbon-isotope examples above.
  - No confirmed worked example exists for *both* together in one name
    with an internal locant already present (`methan(2H,18O)ol (PIN)`
    combines both, but only for a locant-free mononuclear case where
    front-of-name and before-suffix are indistinguishable) -- this
    module raises `UnsupportedStructure` for that combination rather
    than guessing, matching `_isotope.py`'s own established practice.

Scope, deliberately narrow (P82-WS2 M1 step 2, mirrors #789/
`_isotope_alcohol.py`'s own pilot shape): an acyclic ketone chain, with
halogen substituents allowed (mirrors `_ketone.py`'s own scope), bearing
a carbonyl-oxygen isotope XOR a skeletal-carbon isotope (not both).
Explicitly out of scope (raise `UnsupportedStructure`): any ring, any
chain unsaturation, deuterium or any other heteroatom besides the
carbonyl oxygen and halogens, more than one ketone group, a specified
stereocenter, mixing different carbon-isotope nuclides, a carbon isotope
on a branch off the principal chain, and both isotope kinds present
together.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, non_single_bonds, specified_stereocenters
from ._isotope import _CARBON_ISOTOPES
from ._ketone import _best_acyclic_ketone_candidate, _validate_and_collect_ketones

_OXYGEN_ISOTOPES = {17, 18}


def has_isotope_ketone_shape(mol) -> bool:
    """True if `mol` contains some isotope label and looks like it might be
    an acyclic ketone -- deliberately loose (see module docstring for the
    real scope); `name_isotope_ketone` raises for every detail this
    doesn't cover. Used by `core.py` to route here before `_isotope.py`'s
    own plain chain/methane path, which rejects any oxygen outright.

    Excludes a -COOH carbon (a carbonyl carbon that also carries a
    single-bonded hydroxyl) so this doesn't shadow
    `_isotope_carboxylic_acid.py`'s own routing, mirroring
    `_isotope_alcohol.py`'s identical carboxylic-acid exclusion."""
    if not any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
        return False
    if mol.GetRingInfo().NumRings() > 0:
        return False
    return any(
        atom.GetAtomicNum() == 6
        and any(
            n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
        and not any(
            n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and n.GetTotalNumHs() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in atom.GetNeighbors()
        )
        for atom in mol.GetAtoms()
    )


def name_isotope_ketone(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("a ring alongside an isotope-labeled ketone is not supported yet")
    other_unsaturation = [
        (a, b, order)
        for a, b, order in non_single_bonds(mol)
        if mol.GetAtomWithIdx(a).GetAtomicNum() != 8 and mol.GetAtomWithIdx(b).GetAtomicNum() != 8
    ]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation combined with an isotope-labeled ketone is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an isotope-labeled ketone is not supported yet"
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

    ketones, hydroxyls = _validate_and_collect_ketones(neutral_mol)
    if hydroxyls:
        raise UnsupportedStructure(
            "a coexisting hydroxyl alongside an isotope-labeled ketone is not supported yet"
        )
    if len(ketones) != 1:
        raise UnsupportedStructure(
            "more than one ketone group alongside an isotope label is not supported yet"
        )
    (carbonyl_oxygen_idx,) = ketones

    carbonyl_isotope = mol.GetAtomWithIdx(carbonyl_oxygen_idx).GetIsotope()
    if carbonyl_isotope != 0 and carbonyl_isotope not in _OXYGEN_ISOTOPES:
        raise UnsupportedStructure("only 17O/18O oxygen isotopes are supported (P-82.2.1)")
    has_oxygen_isotope = carbonyl_isotope != 0

    stray_isotopes = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetIsotope() != 0
        and atom.GetIdx() not in carbon_isotope_positions
        and atom.GetIdx() != carbonyl_oxygen_idx
    ]
    if stray_isotopes:
        raise UnsupportedStructure(
            "an isotope label outside the carbonyl oxygen and skeletal "
            "carbons is not supported yet"
        )
    if not has_carbon_isotope and not has_oxygen_isotope:
        raise UnsupportedStructure("no isotopically labeled atom found")
    if has_carbon_isotope and has_oxygen_isotope:
        raise UnsupportedStructure(
            "a skeletal-carbon isotope combined with a carbonyl-oxygen "
            "isotope in the same name has no confirmed worked example yet "
            "(see module docstring)"
        )

    best_name, best_position_of = _best_acyclic_ketone_candidate(neutral_mol, ketones, hydroxyls, ())

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
        return f"({locants_str}-{symbol}){best_name}"

    return best_name[:-3] + f"({carbonyl_isotope}O)one"
