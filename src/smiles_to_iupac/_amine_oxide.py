"""Naming of secondary/tertiary amine N-oxides (P-62.5, Chapter P-6,
https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):

Method (1) of P-62.5: a molecule with exactly one amine/imine oxide is
named by functional class nomenclature, appending ' N-oxide' to the name
of the underlying amine, e.g. (CH3)3N+-O- -> 'N,N-dimethylmethanamine
N-oxide (PIN)' (trimethylamine N-oxide, confirmed against the primary
source's own worked example).

This module only accepts a secondary or tertiary amine N-oxide: a single
nitrogen, formal charge +1, bonded to exactly one terminal oxide oxygen
(formal charge -1, single bond) and 2-3 carbon substituents shaped like
whatever `_amine.py`'s existing secondary/tertiary amine logic already
supports (unbranched, saturated, acyclic alkyl chains). The underlying
amine name is produced by literally calling `name_amine` on a version of
the molecule with the oxide oxygen removed and the nitrogen's charge reset
to neutral -- not a separate, parallel implementation.

Explicitly out of scope (the molecule simply isn't matched by
`has_amine_oxide_shape`, so it falls through to whatever other module or
rejection applies -- or `name_amine` itself rejects the reduced molecule):
- A primary amine N-oxide (nitrogen bonded to exactly one carbon) --
  PubChem itself doesn't name this with the same pattern (it computes an
  ammonium-salt-style name instead, P-73 territory), so this project has
  no verified example to build from yet.
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

from rdkit import Chem


def has_amine_oxide_shape(mol) -> bool:
    charged_nitrogens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1
    ]
    if len(charged_nitrogens) != 1:
        return False
    nitrogen = charged_nitrogens[0]
    if nitrogen.GetIsotope() != 0 or nitrogen.GetDegree() not in (3, 4):
        return False

    oxide_neighbors = [
        n
        for n in nitrogen.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 and n.GetDegree() == 1
    ]
    if len(oxide_neighbors) != 1:
        return False
    (oxide_oxygen,) = oxide_neighbors
    if oxide_oxygen.GetIsotope() != 0:
        return False
    if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), oxide_oxygen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False

    other_neighbors = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != oxide_oxygen.GetIdx()]
    if len(other_neighbors) not in (2, 3):
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
    oxide_oxygen = next(n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1)

    reduced = Chem.RWMol(mol)
    reduced.GetAtomWithIdx(nitrogen.GetIdx()).SetFormalCharge(0)
    reduced.RemoveAtom(oxide_oxygen.GetIdx())
    reduced_mol = reduced.GetMol()
    Chem.SanitizeMol(reduced_mol)

    base_name = name_amine(reduced_mol)
    return f"{base_name} N-oxide"
