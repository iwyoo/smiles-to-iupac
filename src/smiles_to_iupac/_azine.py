"""Naming of symmetric azines (R2C=N-N=CR2, both C=N carbon substituents
structurally identical) on acyclic saturated carbon chains, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-68.3.1.2.3 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  an azine, like `_hydrazone.py`'s hydrazone, is not given its own suffix
  here -- the PIN treats one C=N side as an ordinary `_imine.py`-named
  imine parent and the other side as an "N-(...ylideneamino)" substituent
  prefix on that imine's nitrogen, reusing both modules' existing carbon-
  side chain search unchanged (no new algorithm). Confirmed via PubChem
  PUG REST (GET fails with a 400 for a SMILES containing '/', so these
  were queried via POST): CID 136329 (`C=NN=C`) ->
  "N-(methylideneamino)methanimine", CID 68964 (`CC=NN=CC`) ->
  "N-(ethylideneamino)ethanimine", CID 123355 (`CCC=NN=CCC`) ->
  "N-(propylideneamino)propan-1-imine" (the imine side keeps `_imine.py`'s
  own "cite the locant even at C1 for a 3+-carbon chain" rule), CID 79085
  (`CC(C)=NN=C(C)C`) -> "N-(propan-2-ylideneamino)propan-2-imine" (a
  ketone-shaped attachment on both sides).
- Which side becomes the imine parent and which becomes the
  "N-(...ylideneamino)" prefix is unverified for an *asymmetric* azine (two
  different C=N substituents) -- this module is deliberately restricted to
  a *symmetric* azine (both sides structurally identical), where the
  choice is arbitrary and doesn't affect the resulting name: this module
  names both sides with `_hydrazone.py`'s own `_name_hydrazone_carbon`
  ("-ylidene" form) and requires the two results to match textually before
  proceeding: a real, independent confirmation of "same structure" (not
  just "same atom count"), since `_name_hydrazone_carbon` already accounts
  for chain length, branching, and locant placement. If they match, one
  side (arbitrarily, the first C=N bond RDKit reports) is named as the
  imine parent via `_imine.py`'s own `_name_acyclic_imine`, and the other
  side's already-computed "-ylidene" name gets "amino" appended and wrapped
  as an "N-(...)" prefix, mirroring `_imine.py`'s own N-substituent
  formatting (P-62.3.1.1).

Explicitly out of scope (raise `UnsupportedStructure`):
- An asymmetric azine (the two "-ylidene" names differ) -- unverified,
  separate follow-up.
- Everything `_hydrazone.py` itself excludes on each C=N side: any ring,
  aromatic atom, unsaturation besides the azine's own two C=N bonds, any
  heteroatom other than the azine's two nitrogens and halogen substituents,
  and charged or isotopically modified atoms.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, specified_stereocenters
from ._hydrazone import _name_hydrazone_carbon
from ._imine import _name_acyclic_imine

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _find_cn_bonds(mol):
    return [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
    ]


def _find_azine_nitrogens(mol, cn_bonds):
    """Given exactly two C=N double bonds, return (c1, n1, c2, n2) if their
    two nitrogens are joined by a single bond and each nitrogen has no
    other neighbor (the R2C=N-N=CR2 azine shape); else None."""
    if len(cn_bonds) != 2:
        return None
    carbons, nitrogens = [], []
    for bond in cn_bonds:
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        carbons.append(a if a.GetAtomicNum() == 6 else b)
        nitrogens.append(a if a.GetAtomicNum() == 7 else b)
    n1, n2 = nitrogens
    if n1.GetIdx() == n2.GetIdx():
        return None
    nn_bond = mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx())
    if nn_bond is None or nn_bond.GetBondTypeAsDouble() != 1.0:
        return None
    if n1.GetDegree() != 2 or n2.GetDegree() != 2:
        return None
    return carbons[0].GetIdx(), n1.GetIdx(), carbons[1].GetIdx(), n2.GetIdx()


def has_azine_shape(mol) -> bool:
    return _find_azine_nitrogens(mol, _find_cn_bonds(mol)) is not None


def _validate_and_find_azine(mol):
    found = _find_azine_nitrogens(mol, _find_cn_bonds(mol))
    if found is None:
        raise UnsupportedStructure(
            "no azine (R2C=N-N=CR2) skeleton found; this module only "
            "handles azines"
        )
    c1, n1, c2, n2 = found

    for atom in (mol.GetAtomWithIdx(n1), mol.GetAtomWithIdx(n2)):
        if atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("charged atoms are not supported yet")
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the azine's own two nitrogens and "
                "halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic carbon skeleton is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    cn_bond_ids = {bond.GetIdx() for bond in _find_cn_bonds(mol)}
    other_non_single = [
        bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() != 1.0 and bond.GetIdx() not in cn_bond_ids
    ]
    if other_non_single:
        raise UnsupportedStructure(
            "a bond order other than single (besides the azine's own two "
            "C=N bonds) is not supported yet"
        )

    for carbon_idx in (c1, c2):
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) not in (0, 1, 2):
            raise UnsupportedStructure(
                "an azine carbon must have zero, one, or two carbon neighbors"
            )

    return c1, n1, c2, n2


def name_azine(mol) -> str:
    c1, n1, c2, n2 = _validate_and_find_azine(mol)
    if specified_stereocenters(mol) is not None:
        # Both of this molecule's C=N bonds are always flagged by RDKit's
        # `Chem.FindPotentialStereo` as unspecified potential Bond_Double
        # stereo elements, regardless of substituents -- same conclusion
        # as `_imine.py`/`_hydrazone.py` (this module reuses their
        # internal chain-assembly helpers directly, so it carries the
        # identical restriction independently rather than inheriting
        # their own entry-point checks). Any specified chain stereocenter
        # always coexists with at least one of them, so
        # `specified_stereocenters` correctly rejects rather than
        # silently dropping it.
        raise UnsupportedStructure(
            "a specified stereocenter alongside this azine's own "
            "always-unspecified C=N bonds is not supported yet (see "
            "P-92/P-93)"
        )

    ylidene1 = _name_hydrazone_carbon(mol, c1, n1)
    ylidene2 = _name_hydrazone_carbon(mol, c2, n2)
    if ylidene1 != ylidene2:
        raise UnsupportedStructure(
            "an asymmetric azine (two different C=N substituents) is not "
            "supported yet"
        )

    imine_name = _name_acyclic_imine(mol, c1, {n1})
    separator = "-" if imine_name[0].isdigit() else ""
    return f"N-({ylidene1}amino){separator}{imine_name}"
