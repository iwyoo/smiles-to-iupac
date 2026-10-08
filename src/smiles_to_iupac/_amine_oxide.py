"""Naming of secondary/tertiary amine N-oxides (P-62.5, Chapter P-6,
https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):

Method (1) of P-62.5: a molecule with exactly one amine/imine oxide is
named by functional class nomenclature, appending ' N-oxide' to the name
of the underlying amine, e.g. (CH3)3N+-O- -> 'N,N-dimethylmethanamine
N-oxide (PIN)' (trimethylamine N-oxide, confirmed against the primary
source's own worked example).

This module accepts a primary, secondary or tertiary amine N-oxide, N-sulfide, N-selenide or N-telluride: a single
nitrogen, formal charge +1, bonded to exactly one terminal chalcogenide atom
(formal charge -1, single bond) and 1-3 carbon substituents shaped like
whatever `_amine.py`'s existing secondary/tertiary amine logic already
supports (unbranched, saturated, acyclic alkyl chains). The underlying
amine name is produced by literally calling `name_amine` on a version of
the molecule with the oxide oxygen removed and the nitrogen's charge reset
to neutral -- not a separate, parallel implementation.

Explicitly out of scope (the molecule simply isn't matched by
`has_amine_oxide_shape`, so it falls through to whatever other module or
rejection applies -- or `name_amine` itself rejects the reduced molecule):
- An imine oxide (P-62.5's other named class) -- `_imine.py`'s territory,
  not attempted here.
- More than one amine/imine oxide, or a coexisting separate amino group
  elsewhere in the molecule -- P-62.5's Method (1) restricts functional
  class nomenclature to a single such group; anything else needs the
  'amino'-prefix treatment the source describes, out of scope here.
- Anything `_amine.py` itself would reject for the reduced (neutral)
  molecule -- a ring, a branched/unsaturated N-substituent, a halogen
  substituent coexisting with the secondary/tertiary nitrogen, etc.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure


_CHALCOGEN_CLASS = {8: "oxide", 16: "sulfide", 34: "selenide", 52: "telluride"}


def has_amine_oxide_shape(mol) -> bool:
    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1
    ]
    if len(charged_nitrogens) != 1:
        return False
    nitrogen = charged_nitrogens[0]
    if nitrogen.GetIsotope() != 0 or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() not in (3, 4):
        return False

    oxide_neighbors = [
        n
        for n in nitrogen.GetNeighbors()
        if n.GetAtomicNum() in _CHALCOGEN_CLASS and n.GetFormalCharge() == -1 and n.GetDegree() == 1
    ]
    if len(oxide_neighbors) != 1:
        return False
    (oxide_oxygen,) = oxide_neighbors
    if oxide_oxygen.GetIsotope() != 0:
        return False
    if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), oxide_oxygen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False

    other_neighbors = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != oxide_oxygen.GetIdx()]
    if len(other_neighbors) not in (1, 2, 3):
        return False
    if any(n.GetAtomicNum() != 6 for n in other_neighbors):
        return False
    return all(
        mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in other_neighbors
    )


def name_amine_oxide(mol) -> str:
    from ._amine import name_amine

    nitrogen = next(
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1
    )
    oxide_oxygen = next(
        n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in _CHALCOGEN_CLASS and n.GetFormalCharge() == -1
    )
    class_word = _CHALCOGEN_CLASS[oxide_oxygen.GetAtomicNum()]

    reduced = Chem.RWMol(mol)
    reduced.GetAtomWithIdx(nitrogen.GetIdx()).SetFormalCharge(0)
    reduced.RemoveAtom(oxide_oxygen.GetIdx())
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)

    amine_nitrogens = [a for a in reduced_mol.GetAtoms() if _is_amine_nitrogen(a)]
    other_nitrogens = sum(a.GetAtomicNum() == 7 for a in reduced_mol.GetAtoms()) - len(amine_nitrogens)
    if len(amine_nitrogens) == 1 and other_nitrogens:
        return f"{_name_with_nitrile_prefix(reduced_mol)} N-{class_word}"
    if len(amine_nitrogens) > 1:
        return f"{_name_with_oxidized_parent(reduced_mol, nitrogen.GetIdx(), oxide_oxygen.GetIdx())} N-{class_word}"
    from ._common import heteroatom_stereo_prefix, specified_stereocenters

    centres = specified_stereocenters(mol) or []
    stereo_prefix = heteroatom_stereo_prefix(mol, nitrogen.GetIdx()) or "" if any(i == nitrogen.GetIdx() for i, _ in centres) else ""
    try:
        base_name = name_amine(reduced_mol)
    except UnsupportedStructure:
        if nitrogen.IsInRing():
            raise
        base_name = _name_with_nitrile_prefix(reduced_mol)
    return f"{stereo_prefix}{base_name} N-{class_word}"


def _is_amine_nitrogen(atom):
    """A nitrogen with single bonds to carbon or hydrogen only, none of the carbons an acyl or nitrile carbon."""
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge():
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if bond.GetBondTypeAsDouble() != 1.0 or other.GetAtomicNum() != 6:
            return False
        if any(b.GetBondTypeAsDouble() > 1.0 and b.GetOtherAtom(other).GetAtomicNum() != 6 for b in other.GetBonds()):
            return False
    return True


def _name_with_nitrile_prefix(reduced_mol):
    """P-62.5, P-67.1.6: the amine oxide is the parent, so a nitrile is cited as the prefix 'cyano' ('cyano-N,N-
    dimethylmethanamine N-oxide') and a carboxy group beside it as 'carboxy'."""
    from ._polyfunctional import FORCED_PRINCIPAL, name_polyfunctional

    token = FORCED_PRINCIPAL.set("amine")
    try:
        name = name_polyfunctional(reduced_mol)
    finally:
        FORCED_PRINCIPAL.reset(token)
    # P-14.3.4.4: the oxidized nitrogen carries no hydrogen, so the only position left for the prefix is carbon 1
    return re.sub(r"^1-(?=[a-z(\[{])", "", name) if name.endswith("methanamine") else name


def _name_with_oxidized_parent(reduced_mol, nitrogen_idx, oxide_idx):
    """P-62.5: the oxidized nitrogen is the amine suffix nitrogen of the parent, so every other amino group is cited
    as a prefix ('5-(dimethylamino)-N,N-dimethylpentan-1-amine N-oxide')."""
    from .core import _name_mol
    from ._common import UnsupportedStructure
    from ._hetero_chain import contract_hetero_groups_candidates

    position = nitrogen_idx - (1 if oxide_idx < nitrogen_idx else 0)
    names = []
    for contracted in contract_hetero_groups_candidates(reduced_mol, {position}):
        try:
            names.append(_name_mol(contracted))
        except UnsupportedStructure:
            continue
    if not names:
        raise UnsupportedStructure("no amine parent carries the oxidized nitrogen of this polyamine N-oxide")
    return min(names)
