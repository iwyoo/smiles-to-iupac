"""Naming of thioic acids (-C(=O)-SH or -C(=S)-OH), the chalcogen analogue
of a carboxylic acid with one oxygen replaced by sulfur, restricted to a
single such group on an acyclic unbranched saturated hydrocarbon chain, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-65.1.5 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  functional replacement of one -COOH oxygen by sulfur is named with the
  'thioic acid' suffix. Because -CO-SH and -CS-OH are distinct tautomers,
  the element symbol of whichever atom carries the acidic hydrogen is
  prefixed (italicized in print) to 'acid': 'thioic S-acid' for -CO-SH (the
  sulfur carries the -H), 'thioic O-acid' for -CS-OH (the oxygen carries
  the -H).
- Although P-65.1.5.1 says these letter locants are "normally omitted...
  because the exact position... is not known or important", every worked
  example for a fully-determined structure (as any RDKit-parsed SMILES
  necessarily is) DOES cite the letter -- confirmed directly:
  'hexanethioic O-acid (PIN)' for CH3(CH2)4-CS-OH, 'ethanethioic O-acid
  (PIN)' for CH3-CS-OH, and 'methanethioic S-acid (PIN)' for HCO-SH
  (P-65.1.5.2). This module therefore always cites the letter.
- PubChem's own structure-lookup service was checked and found unusable as
  a cross-check here: querying it with 'CC(=O)S' and 'CC(=S)O' (the two
  distinct tautomers above) both resolve to the *same* CID (10484) and the
  same single auto-generated name -- PubChem's structure normalization
  collapses the very distinction this module exists to preserve. Verification
  instead relies directly on the three primary-source worked examples above.
- Formic acid's analogue (chain length 1, R = H) is named with the 'methane'
  stem, not a 'carbo-' prefix -- 'methanethioic S-acid (PIN)', mirroring
  `_carboxylic_acid.py`'s identical treatment of formic acid itself
  ('methanoic acid').

Explicitly out of scope (raise `UnsupportedStructure`):
- A branched or unsaturated R group (unlike `_carboxylic_acid.py`, this
  first cut mirrors `_sulfoxide.py`'s minimal-scope style).
- More than one thioic acid group (the dicarboxylic-acid analogue,
  'bis(thioic acid)'/'dithioic acid', uses a different construction per
  P-65.1.5.1 and is out of scope).
- A ring anywhere in the molecule (P-65.1.1.2's 'carbothioic acid'
  construction, out of scope here), except for one narrow case: a thioic
  acid's chain hanging off a single plain, unsubstituted benzene ring with
  no other ring substituent (`_name_phenyl_chain_thioic_acid`, e.g.
  '2-phenylethanethioic S-acid'), mirroring `_aldehyde.py`'s identical
  benzene-ring-substituent path. A thioic acid directly on the ring stays
  out of scope for this chain-parent module.
- The selenium/tellurium analogues (selenoic/telluroic acid) -- the
  selenium analogue is handled separately by `_selenoic_acid.py`; tellurium
  is a separate follow-up module.
- Any other heteroatom, halogen substituent, charge, or isotopic label.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ordered_chain,
    ring_chain_attachment,
)
from ._numerals import alkane_name


def _thioic_acid_carbons(mol):
    """Carbons shaped like a thioic-acid group: bonded to exactly one
    double-bonded (terminal) O or S, and exactly one single-bonded,
    one-H, terminal O or S of the other element. Returns a list of
    (carbon, label) pairs, where label is 'O' or 'S' -- whichever chalcogen
    carries the acidic hydrogen."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, 16)]
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
        if {double_atom.GetAtomicNum(), single_atom.GetAtomicNum()} != {8, 16}:
            continue
        label = "O" if single_atom.GetAtomicNum() == 8 else "S"
        matches.append((atom, label, double_atom, single_atom))
    return matches


def has_thioic_acid_shape(mol) -> bool:
    return bool(_thioic_acid_carbons(mol))


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point. Mirrors
    `_sulfoxide.py`'s identical helper."""
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


def _validate_and_collect_thioic_acid(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_thioic_acid`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_thioic_acid`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged. Returns
    (acid_carbon, label, acid_atom_idxs)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 16):
            raise UnsupportedStructure(
                "heteroatoms other than the thioic acid's own chalcogens are "
                "not supported yet (P-65.1.5 is restricted to a plain "
                "acyclic thioic acid here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure(
                "aromatic rings are out of scope for this module (see the "
                "separate aromatic-ring module)"
            )
    matches = _thioic_acid_carbons(mol)
    if len(matches) != 1:
        raise UnsupportedStructure("more than one thioic acid group is out of scope for this module")
    (acid_carbon, label, double_atom, single_atom) = matches[0]
    acid_atom_idxs = {acid_carbon.GetIdx(), double_atom.GetIdx(), single_atom.GetIdx()}
    if any(
        a not in acid_atom_idxs
        and b not in acid_atom_idxs
        and not (a in aromatic_ring_atoms and b in aromatic_ring_atoms)
        for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (P-65.1.5's scope "
            "here is limited to a single saturated chain)"
        )
    return acid_carbon, label, acid_atom_idxs


def _name_phenyl_chain_thioic_acid(mol, ring_atoms):
    """Name a thioic acid whose -C(=O)SH/-C(=S)OH lies entirely on a
    single unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. '2-phenylethanethioic S-acid'. The
    ring is cited as a 'phenyl' substituent prefix on the chain, which is
    the parent hydride, mirroring `_aldehyde.py`'s
    `_name_phenyl_chain_aldehyde`. This module never supports any
    substituent besides the ring itself (module docstring), so no locant
    tie-break is needed: the acid carbon is always C1 and the ring is
    always at the chain's far terminus (see
    tasks/phenyl-substituent-on-thioic-acid-chain.md's scope note)."""
    acid_carbon, label, acid_atom_idxs = _validate_and_collect_thioic_acid(mol, aromatic_ring_atoms=ring_atoms)
    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain thioic acid is not supported yet"
        )
    ring_atom, chain_root = attachment
    hetero_idxs = acid_atom_idxs - {acid_carbon.GetIdx()}
    chain = ordered_chain(graph, chain_root, ring_atom, hetero_idxs)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "thioic acid is not supported yet"
        )
    if chain[-1] != acid_carbon.GetIdx():
        raise UnsupportedStructure(
            "the thioic acid carbon must be the chain's far terminus from "
            "the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a thioic acid directly attached to the benzene ring uses a "
            "separate naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    return f"{chain_length}-phenyl{alkane_name(chain_length)}thioic {label}-acid"


def name_thioic_acid(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_thioic_acid(mol, ring_atoms)
    acid_carbon, label, acid_atom_idxs = _validate_and_collect_thioic_acid(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(chain_neighbors) > 1:
        raise UnsupportedStructure("a thioic-acid carbon with more than one carbon neighbor is not valid")

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

    return f"{alkane_name(chain_length)}thioic {label}-acid"
