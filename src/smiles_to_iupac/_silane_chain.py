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

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    longest_branched_chain,
    substituent_locant_set_and_citation,
)
from ._numerals import numerical_term
from ._substituents import format_substituent_prefixes, name_branch

_SILICON = 14


def has_silane_chain_shape(mol) -> bool:
    """True iff every atom in `mol` is silicon (the shape this module
    accepts -- everything else is left to `core.py`'s other branches)."""
    atoms = mol.GetAtoms()
    return len(atoms) > 0 and all(atom.GetAtomicNum() == _SILICON for atom in atoms)


def name_silane_chain(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
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
        return _name_branched_silane(mol, graph)
    # Unlike '-ane' (P-21.2.1's alkane suffix, which starts with a vowel
    # and so triggers elision, e.g. 'hexa' + 'ane' -> 'hexane'), 'silane'
    # starts with a consonant, so the numerical term's terminal 'a' is
    # kept unchanged: 'tetra' + 'silane' -> 'tetrasilane', not
    # 'tetrsilane'.
    return numerical_term(n) + "silane"


def _name_branched_silane(mol, graph) -> str:
    """P-44.3 chain selection (longest, then most prefixes, then lowest locants) on the Si skeleton;
    every branch is a silyl-type substituent (P-29.4.1)."""
    best = None
    for leaf in (a for a, neighbors in graph.items() if len(neighbors) == 1):
        chain, branches = longest_branched_chain(graph, leaf, frozenset())
        substituents = {
            position: [name_branch(graph, root, chain[position - 1], {}, frozenset(), mol=mol) for root in roots]
            for position, roots in branches.items()
        }
        grouped = group_substituents(substituents)
        locant_set, total, citation = substituent_locant_set_and_citation(grouped)
        name = format_substituent_prefixes(grouped) + numerical_term(len(chain)) + "silane"
        key = (-len(chain), -total, locant_set, citation, name)
        if best is None or key < best:
            best = key
    return best[-1]
