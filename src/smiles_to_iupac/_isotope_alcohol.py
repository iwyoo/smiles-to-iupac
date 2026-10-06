"""Naming of an acyclic alcohol chain bearing an isotopically labeled
hydroxyl oxygen (17O/18O) and/or a skeletal-carbon isotope (12C/13C/14C),
combined with the existing '-ol' suffix, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-82.5.1 (the Blue Book): "Numbering of an
  isotopically substituted compound is not changed from that of an
  isotopically unmodified compound ... the presence of nuclides is
  considered last" among P-14.4's numbering-priority criteria -- so the
  winning chain numbering is exactly whatever `_alcohol.py`'s own acyclic
  path would already pick for the isotopically-unmodified structure,
  reused directly via `_best_acyclic_alcohol_candidate` on a neutralized
  (isotope-stripped) copy of the molecule rather than re-derived here.
- Placement depends on which atom carries the isotope, corrected
  2026-09-21 (#804) after cross-checking more of the Blue Book's
  own worked examples than the original pass did:
  - An isotope on the hydroxyl oxygen itself (the suffix-defining atom)
    is inserted directly before the '-ol' suffix: 'methan(18O)ol'-style,
    confirmed by '1-(aminomethyl)cyclopentan-1-(18O)ol (PIN)' (line 187)
    -- the oxygen IS that molecule's own suffix-defining atom too.
  - A skeletal-carbon isotope (not the suffix-defining atom) instead goes
    at the front of the whole name, per P-82.2.5's own general rule ("In
    a name consisting of one word, the isotopic descriptor is placed
    before the name ... preferred to ... placing the descriptor before
    the implied name of the characteristic group") and confirmed
    directly on this exact ethanol shape: '(2-13C)ethan-1-ol (PIN)'
    (line ~178) -- **not** 'ethan(2-13C)ol' as an earlier pass of this
    module produced. A 2-carbon chain's own -OH locant, normally omitted
    (P-44.4.1.8, plain 'ethanol' never says 'ethan-1-ol'), is restored
    and explicitly cited once this front-of-name descriptor is present,
    exactly as the confirmed PIN shows. A 1-carbon (methanol) parent
    never cites a locant at all, mirroring P-82.2.1's own '(14C)methane'
    (not '(1-14C)methane').
  - No confirmed worked example was found for a skeletal-carbon isotope
    and a hydroxyl-oxygen isotope combined in the same name once an
    internal locant is actually possible (i.e. a 2-carbon chain) --
    'methan(2H,18O)ol (PIN)' merges both into one descriptor, but only
    for the locant-free 1-carbon case, which doesn't resolve how (or
    whether) the two would split for a 2-carbon chain. This module raises
    `UnsupportedStructure` for that combination rather than guessing,
    matching this project's established practice (see e.g.
    `_isotope_ketone.py`'s identical deferral for its own combined case).

Scope, deliberately narrow (P82-WS2 M1 step 1 pilot, chain length 1-2
only -- matching the confirmed worked examples above exactly): a 1- or
2-carbon acyclic alcohol chain, with halogen substituents allowed
(mirrors `_alcohol.py`'s own scope), bearing a skeletal carbon isotope
XOR a hydroxyl-oxygen isotope (not both, see above). A longer chain, where
the -OH's own suffix locant would already be cited regardless of any
isotope, has no confirmed worked example for exactly how it interacts
with a skeletal-carbon isotope's own front-of-name placement, so it is
deferred to a later step rather than guessed at. Explicitly out of scope
(raise `UnsupportedStructure`): any ring, any chain unsaturation, any
ether, deuterium or any other heteroatom besides the hydroxyl oxygen and
halogens, more than one hydroxyl, a specified stereocenter, mixing
different carbon-isotope nuclides, a carbon isotope on a branch off the
1-2-carbon chain, a chain length of 3 or more carrying any isotope label
at all, and both isotope kinds present together (see above).
"""

from rdkit import Chem

from ._alcohol import _best_acyclic_alcohol_candidate, _validate_and_collect_hydroxyls
from ._common import UnsupportedStructure, adjacency, non_single_bonds, specified_stereocenters
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

    if has_carbon_isotope and has_oxygen_isotope:
        raise UnsupportedStructure(
            "a skeletal-carbon isotope combined with a hydroxyl-oxygen "
            "isotope in the same name has no confirmed worked example yet "
            "(see module docstring)"
        )

    if has_carbon_isotope:
        if any(c not in best_position_of for c in carbon_isotope_positions):
            raise UnsupportedStructure(
                "a carbon isotope on a substituent branch rather than the "
                "principal chain is not supported yet"
            )
        nuclide = next(iter(carbon_isotopes))
        symbol = f"{nuclide}C" + (str(len(carbon_isotope_positions)) if len(carbon_isotope_positions) > 1 else "")
        if chain_length == 1:
            # A mononuclear parent never cites a locant (P-82.2.1's own
            # '(14C)methane', not '(1-14C)methane').
            return f"({symbol}){best_name}"
        carbon_isotope_locants = sorted(best_position_of[c] for c in carbon_isotope_positions)
        locants_str = ",".join(str(loc) for loc in carbon_isotope_locants)
        (oh_carbon,) = adjacency(neutral_mol)[oh_oxygen_idx]
        oh_locant = best_position_of[oh_carbon]
        stem = best_name[:-2]
        return f"({locants_str}-{symbol}){stem}-{oh_locant}-ol"

    return best_name[:-2] + f"({oxygen_isotope}O)ol"
