"""Naming of simple ammonium cations (the '-aminium'/'azanium' suffix family, P-73.1.1.2).

A cation formed by adding a hydron to a parent hydride changes the terminal 'e' of its name to 'ium': 'methanamine' ->
'methanaminium'; NH4+ is 'azanium'. The 'aminium' derivation, not 'tetramethylazanium', is the PIN up to a quaternary
nitrogen (Table 7.3 lists 'N,N,N-trimethylmethanaminium'). A secondary or tertiary cation is named through its
neutralized amine; a quaternary one and a cation with supplied onium prefixes go to the chain engine directly. Anything
else is refused: a ring nitrogen, a carbonyl-type group beside the cation, several charged atoms, an isotope.
"""

import re

from rdkit import Chem

from ._common import CATIONIC_AMINES, UnsupportedStructure
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
        neighbor.GetAtomicNum() in (6, 9, 17, 35, 53)
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
_MULTIPLIER_OF = dict(enumerate(_MULTIPLIED_AMINIUM, start=2))


def name_polyammonium(mol) -> str:
    """Several ammonium nitrogens on one parent are named as a multiplied 'aminium' suffix, 'bis(aminium)', 'tris(...)'
    (P-73.1.2.1): each hydrogen-bearing cation is neutralised for the amine namer, which numbers the N-substituents."""
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    functional = any(a.GetAtomicNum() == 6 and is_functional_carbon(mol, a.GetIdx()) for a in mol.GetAtoms())
    neutral = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() and atom.GetTotalNumHs():
            copy = neutral.GetAtomWithIdx(atom.GetIdx())
            copy.SetFormalCharge(0)
            copy.SetNoImplicit(True)
            copy.SetNumExplicitHs(atom.GetTotalNumHs() - 1)
            copy.SetBoolProp("_cationic_amine", True)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    from ._multiplicative import name_if_multiplicative
    from ._polyfunctional import AMINIUM, cationic_multiplied_units, name_polyfunctional

    token = CATIONIC_AMINES.set(True)
    aminium_token = AMINIUM.set(True)
    try:
        try:
            name = None if functional else name_polyfunctional(base)
        except UnsupportedStructure:
            name = None
        match = re.search(r"(di|tri|tetra|penta|hexa)amine$", name or "")
        cations = sum(1 for a in mol.GetAtoms() if a.GetFormalCharge())
        if match is not None and match.group(1) == _MULTIPLIER_OF.get(cations):
            return f"{name[:match.start()]}{_MULTIPLIED_AMINIUM[match.group(1)]}(aminium)"
        multiplicative = name_if_multiplicative(base)
    finally:
        AMINIUM.reset(aminium_token)
        CATIONIC_AMINES.reset(token)
    named = cationic_multiplied_units(multiplicative, sum(1 for a in mol.GetAtoms() if a.GetFormalCharge())) if multiplicative else None
    if named is None and functional:
        raise UnsupportedStructure("a carbonyl-type group beside the cations is cited as a prefix by the chain engine")
    if named is None:
        return name_polyfunctional(mol)
    return named


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
    from ._polyfunctional import name_polyfunctional

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

    if (
        any(_heteroatom_acyl(mol, n) or n.GetAtomicNum() != 6 for n in nitrogen.GetNeighbors()) and not nitrogen.IsInRing()
    ):
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
        return name_polyfunctional(mol)

    neutral_rw = Chem.RWMol(mol)
    neutral_nitrogen = neutral_rw.GetAtomWithIdx(nitrogen.GetIdx())
    neutral_nitrogen.SetFormalCharge(0)
    neutral_nitrogen.SetNoImplicit(True)
    neutral_nitrogen.SetNumExplicitHs(3 - degree)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    if sum(a.GetAtomicNum() == 7 for a in mol.GetAtoms()) > 1:
        raise UnsupportedStructure("a neutral amino group beside the ammonium group is a prefix of the aminium name")
    amine_name = name_polyfunctional(neutral_mol)
    if not amine_name.endswith(("amine", "aniline")) or not _cites_forced_prefixes(amine_name):
        raise UnsupportedStructure("the amine name of this cation does not cite every prefix supplied for a further onium group")
    return amine_name[:-1] + "ium"


def _cites_forced_prefixes(name):
    from ._substituents import FORCED_BRANCH_NAMES

    forced = FORCED_BRANCH_NAMES.get()
    return forced is None or all(text in name for text, _ in forced[1].values())
