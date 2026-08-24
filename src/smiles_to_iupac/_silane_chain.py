"""Naming of unbranched, unsubstituted acyclic silicon hydride chains
(silane, disilane, trisilane, ...), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-21.2.1/P-21.2.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  a chain of n silicon atoms (all single Si-Si bonds, each terminal Si
  bearing three H, each internal Si bearing two H -- i.e. SinH(2n+2), the
  silicon analogue of an unbranched alkane) is named with the P-14.2.1
  numerical-term multiplying prefix ('di', 'tri', 'tetra', ...) plus
  'silane'; n=1 (SiH4) is simply 'silane', with no 'mono-' prefix (the
  same pattern `_numerals.alkane_name` uses for 'methane'). Cross-checked
  against the well-known compounds disilane (Si2H6) and trisilane
  (Si3H8).
- This module is deliberately scoped to a *straight* chain only -- no
  branching (P-21.2.1's general branched-chain substitutive naming for
  silanes is structurally analogous to `_acyclic.py`'s alkane handling
  but is a separate, not-yet-implemented piece of work), no substituents,
  and no other heteroatom skeletal replacement (an -O-/-N-/-S- atom
  interrupting a chain, i.e. genuine P-21.2's 'a'-prefix skeletal
  replacement nomenclature for acyclic chains, is a materially different
  and still-open question -- see the task file this module was scoped
  from for why: it isn't obviously distinguishable in general from what
  `_ether.py`/a future `sulfide` module already cover via substituent-
  prefix naming instead).

Explicitly out of scope (raise `UnsupportedStructure`): any branching, any
substituent (halogen included), any atom other than silicon, any bond
order other than single, and any ring.
"""

from ._common import UnsupportedStructure, adjacency
from ._numerals import numerical_term

_SILICON = 14


def has_silane_chain_shape(mol) -> bool:
    """True iff every atom in `mol` is silicon (the shape this module
    accepts -- everything else is left to `core.py`'s other branches)."""
    atoms = mol.GetAtoms()
    return len(atoms) > 0 and all(atom.GetAtomicNum() == _SILICON for atom in atoms)


def name_silane_chain(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "a non-single Si-Si bond is not supported yet (see P-21.2.1)"
            )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("cyclic silicon skeletons are not supported yet (see P-22)")

    graph = adjacency(mol)
    n = len(graph)
    if n == 1:
        return "silane"
    degrees = sorted(len(neighbors) for neighbors in graph.values())
    if degrees != [1, 1] + [2] * (n - 2):
        raise UnsupportedStructure(
            "a branched silicon chain is not supported yet (see P-21.2.1's "
            "general branched-chain substitutive naming, not yet implemented "
            "for silanes)"
        )
    # Unlike '-ane' (P-21.2.1's alkane suffix, which starts with a vowel
    # and so triggers elision, e.g. 'hexa' + 'ane' -> 'hexane'), 'silane'
    # starts with a consonant, so the numerical term's terminal 'a' is
    # kept unchanged: 'tetra' + 'silane' -> 'tetrasilane', not
    # 'tetrsilane'.
    return numerical_term(n) + "silane"
