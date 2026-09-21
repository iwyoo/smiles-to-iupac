"""Naming of an acyclic alcohol chain bearing an isotopically labeled
hydroxyl oxygen (17O/18O) and/or a skeletal-carbon isotope (12C/13C/14C),
combined with the existing '-ol' suffix, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-82.5.1 (`tmp/bluebook/P8.txt` ~449-453): "Numbering of an
  isotopically substituted compound is not changed from that of an
  isotopically unmodified compound ... the presence of nuclides is
  considered last" among P-14.4's numbering-priority criteria -- so the
  winning chain numbering is exactly whatever `_alcohol.py`'s own acyclic
  path would already pick for the isotopically-unmodified structure,
  reused directly via `_best_acyclic_alcohol_candidate` on a neutralized
  (isotope-stripped) copy of the molecule rather than re-derived here.
- The isotope descriptor is inserted directly before the '-ol' suffix it
  modifies, confirmed by two worked examples already documented in
  `_isotope.py`'s own module docstring: 'methan(2H,18O)ol (PIN)'
  (mononuclear parent, no locant at all) and 'ethan(2H)ol (PIN)' (a
  2-carbon chain, whose own -OH locant is already omitted per
  `_alcohol.py`'s established P-44.4.1.8 rule regardless of isotopes).
  Since the descriptor always appears directly in front of the literal,
  unmodified 'ol' tail of `_alcohol.py`'s own plain name (elision only
  ever touches the *stem*'s trailing letter, never the bare suffix word
  itself), this module splices the descriptor in immediately before that
  tail rather than re-deriving the stem/elision/prefix logic itself.
- A skeletal carbon isotope's own locant is cited exactly like an
  ordinary substituent's (never omitted, even when it lands on the same
  carbon whose own suffix locant *is* omitted -- confirmed by
  `_alcohol.py`'s own established precedent 'CC(Cl)S(=O)(=O)O' ->
  '1-chloroethanesulfonic acid', where the chlorine's locant is cited
  despite sitting on the very carbon whose sulfonic-acid suffix locant is
  omitted). This resolves `_isotope.py`'s own "2-carbon chain, single
  isotope, no halogen" ambiguity for the bare-alkane case: an alcohol's
  own -OH already breaks the 2-carbon symmetry (P-44.4.1.8 fixes which
  carbon is C1), so the carbon-isotope's locant is never actually
  ambiguous here, unlike the plain, suffix-less 2-carbon chain case
  `_isotope.py` itself must defer on.
- P-82.3.1's alphabetical nuclide-symbol order ('C' before 'O') puts the
  carbon-isotope part first when both are present in one descriptor
  group, mirroring `_isotope.py`'s own carbon-before-hydrogen ordering.

Scope, deliberately narrow (P82-WS2 M1 step 1 pilot, chain length 1-2
only -- matching the two confirmed worked examples above exactly): a
1- or 2-carbon acyclic alcohol chain, with halogen substituents allowed
(mirrors `_alcohol.py`'s own scope), bearing a skeletal carbon isotope
and/or a hydroxyl-oxygen isotope. A longer chain, where the -OH's own
suffix locant would itself be cited (P-44.4.1.8), has no confirmed worked
example for where the isotope descriptor then splices in relative to
that locant, so it is deferred to a later step rather than guessed at.
Explicitly out of scope (raise `UnsupportedStructure`): any ring, any
chain unsaturation, any ether, deuterium or any other heteroatom besides
the hydroxyl oxygen and halogens, more than one hydroxyl, a specified
stereocenter, mixing different carbon-isotope nuclides, a carbon isotope
on a branch off the 1-2-carbon chain, and a chain length of 3 or more
carrying any isotope label at all (see above).
"""

from rdkit import Chem

from ._alcohol import _best_acyclic_alcohol_candidate, _validate_and_collect_hydroxyls
from ._common import UnsupportedStructure, non_single_bonds, specified_stereocenters
from ._isotope import _CARBON_ISOTOPES

_OXYGEN_ISOTOPES = {17, 18}


def has_isotope_alcohol_shape(mol) -> bool:
    """True if `mol` contains some isotope label and looks like it might
    be an acyclic alcohol -- deliberately loose (see module docstring for
    the real scope); `name_isotope_alcohol` raises for every detail this
    doesn't cover. Used by `core.py` to route here before `_isotope.py`'s
    own plain chain/methane path, which rejects any oxygen outright.

    Excludes a carboxylic acid's own -OH (a carbon bearing this oxygen
    that's also doubly-bonded to a second oxygen is a -COOH, not a plain
    alcohol) so this doesn't shadow `_isotope_carboxylic_acid.py`'s own
    routing -- confirmed necessary (2026-09-21): a plain -OH-shape check
    alone also matches a carboxylic acid's hydroxyl once any isotope is
    present anywhere in the molecule, and `core.py` checks this module
    first."""
    if not any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
        return False
    if mol.GetRingInfo().NumRings() > 0:
        return False
    return any(
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and not any(
            o.GetAtomicNum() == 8 and o.GetIdx() != atom.GetIdx() and mol.GetBondBetweenAtoms(n.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
            for o in n.GetNeighbors()
        )
        for atom in mol.GetAtoms()
    )


def name_isotope_alcohol(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("a ring alongside an isotope-labeled alcohol is not supported yet")
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "chain unsaturation combined with an isotope-labeled alcohol is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an isotope-labeled alcohol is not supported yet"
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

    hydroxyls, ethers = _validate_and_collect_hydroxyls(neutral_mol)
    if ethers:
        raise UnsupportedStructure(
            "an alkoxy ether alongside an isotope-labeled alcohol is not supported yet"
        )
    if len(hydroxyls) != 1:
        raise UnsupportedStructure(
            "more than one hydroxyl alongside an isotope label is not supported yet"
        )
    (oh_oxygen_idx,) = hydroxyls

    oxygen_isotope = mol.GetAtomWithIdx(oh_oxygen_idx).GetIsotope()
    if oxygen_isotope != 0 and oxygen_isotope not in _OXYGEN_ISOTOPES:
        raise UnsupportedStructure("only 17O/18O oxygen isotopes are supported (P-82.2.1)")
    has_oxygen_isotope = oxygen_isotope != 0

    stray_isotopes = [
        atom
        for atom in mol.GetAtoms()
        if atom.GetIsotope() != 0 and atom.GetIdx() not in carbon_isotope_positions and atom.GetIdx() != oh_oxygen_idx
    ]
    if stray_isotopes:
        raise UnsupportedStructure(
            "an isotope label outside the hydroxyl oxygen and skeletal "
            "carbons is not supported yet"
        )
    if not has_carbon_isotope and not has_oxygen_isotope:
        raise UnsupportedStructure("no isotopically labeled atom found")

    best_name, best_position_of = _best_acyclic_alcohol_candidate(neutral_mol, hydroxyls, ())
    chain_length = len(best_position_of)
    if chain_length > 2:
        raise UnsupportedStructure(
            "an isotope-labeled alcohol chain longer than 2 carbons is not "
            "covered by a confirmed worked example yet (where the isotope "
            "descriptor splices in relative to the -OH's own cited "
            "locant, P-44.4.1.8, is unconfirmed)"
        )

    descriptor_parts = []
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
        descriptor_parts.append(f"{locants_str}-{symbol}")
    if has_oxygen_isotope:
        descriptor_parts.append(f"{oxygen_isotope}O")

    descriptor = ",".join(descriptor_parts)
    return best_name[:-2] + f"({descriptor})ol"
