"""Naming of telluroic acids (-C(=O)-TeH or -C(=Te)-OH), the tellurium
analogue of a thioic/selenoic acid, restricted to a single such group on an
acyclic unbranched saturated hydrocarbon chain, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-65.1.5 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Replacement of oxygen atom(s) of a carboxylic acid group by another
  chalcogen is indicated by the affixes 'thio', 'seleno', and 'telluro'" --
  i.e. exactly `_selenoic_acid.py`'s substitutive rule and letter-labeling
  mechanism (see that module's docstring), with 'telluro' standing in for
  'seleno' and the label 'Te' standing in for 'Se'.
- Confirmed via PubChem for the Te-acid-direction tautomer (the acidic
  hydrogen on the chalcogen, C=O retained), formic acid's analogue:
  '[TeH]C=O' resolves to CID 173268937, auto-generated name 'methanetelluroic
  Te-acid'.
- The longer-chain Te-acid case ('CC(=O)[TeH]') has no PubChem entry at all
  (CID 0), and the O-acid-direction tautomer ('CC(=[Te])O') has a
  PubChem-registered CID (85820967) but no computed name -- both are a
  reviewed, not directly PubChem-verified, mechanical extension of the same
  labeling rule already confirmed twice over (`_thioic_acid.py`,
  `_selenoic_acid.py`) for this identical chalcogen-replacement pattern.
- Formic acid's analogue (chain length 1, R = H) uses the 'methane' stem,
  mirroring `_thioic_acid.py`/`_selenoic_acid.py`'s identical treatment.

This module completes P-65.1.5's chalcogen set (S/Se/Te) for a single
monocarboxylic-acid-analogue group on an acyclic unbranched saturated
chain.

Explicitly out of scope (raise `UnsupportedStructure`): same list as
`_selenoic_acid.py` -- a branched or unsaturated R group, more than one
telluroic acid group, a ring anywhere in the molecule, any other
heteroatom, halogen substituent, charge, or isotopic label.
"""

from ._common import UnsupportedStructure, non_single_bonds
from ._numerals import alkane_name


def _telluroic_acid_carbons(mol):
    """Carbons shaped like a telluroic-acid group: bonded to exactly one
    double-bonded (terminal) O or Te, and exactly one single-bonded,
    one-H, terminal O or Te of the other element. Returns a list of
    (carbon, label, double_atom, single_atom) tuples, where label is 'O'
    or 'Te' -- whichever chalcogen carries the acidic hydrogen."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, 52)]
        if len(chalcogens) != 2:
            continue
        double_bonded = []
        single_bonded = []
        for chalcogen in chalcogens:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), chalcogen.GetIdx())
            if chalcogen.GetDegree() != 1:
                double_bonded = single_bonded = None
                break
            if bond.GetBondTypeAsDouble() == 2.0:
                double_bonded.append(chalcogen)
            elif bond.GetBondTypeAsDouble() == 1.0 and chalcogen.GetTotalNumHs() == 1:
                single_bonded.append(chalcogen)
            else:
                double_bonded = single_bonded = None
                break
        if not double_bonded or len(double_bonded) != 1 or len(single_bonded) != 1:
            continue
        (double_atom,) = double_bonded
        (single_atom,) = single_bonded
        if {double_atom.GetAtomicNum(), single_atom.GetAtomicNum()} != {8, 52}:
            continue
        label = "O" if single_atom.GetAtomicNum() == 8 else "Te"
        matches.append((atom, label, double_atom, single_atom))
    return matches


def has_telluroic_acid_shape(mol) -> bool:
    return bool(_telluroic_acid_carbons(mol))


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point. Mirrors
    `_selenoic_acid.py`'s identical helper."""
    length = 0
    previous = exclude_idx
    current = root_idx
    while True:
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        length += 1
        if not neighbors:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]


def name_telluroic_acid(mol) -> str:
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 52):
            raise UnsupportedStructure(
                "heteroatoms other than the telluroic acid's own chalcogens "
                "are not supported yet (P-65.1.5 is restricted to a plain "
                "acyclic telluroic acid here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    matches = _telluroic_acid_carbons(mol)
    if len(matches) != 1:
        raise UnsupportedStructure("more than one telluroic acid group is out of scope for this module")
    (acid_carbon, label, double_atom, single_atom) = matches[0]
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    acid_atom_idxs = {acid_carbon.GetIdx(), double_atom.GetIdx(), single_atom.GetIdx()}
    if any(a not in acid_atom_idxs and b not in acid_atom_idxs for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-65.1.5's scope "
            "here is limited to a single saturated chain)"
        )

    chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(chain_neighbors) > 1:
        raise UnsupportedStructure("a telluroic-acid carbon with more than one carbon neighbor is not valid")

    if chain_neighbors:
        (chain_root,) = chain_neighbors
        chain_length = _unbranched_chain_length(mol, chain_root.GetIdx(), acid_carbon.GetIdx())
        if chain_length is None:
            raise UnsupportedStructure(
                "a branched R group is out of scope for this module (see "
                "module docstring)"
            )
        chain_length += 1
    else:
        chain_length = 1

    return f"{alkane_name(chain_length)}telluroic {label}-acid"
