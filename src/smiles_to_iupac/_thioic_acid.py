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
- A branched or unsaturated R group on the plain acyclic path (unlike
  `_carboxylic_acid.py`, this first cut mirrors `_sulfoxide.py`'s
  minimal-scope style) -- the benzene-ring-substituent path below is the
  one exception, where a branched chain is supported.
- More than one thioic acid group (the dicarboxylic-acid analogue,
  'bis(thioic acid)'/'dithioic acid', uses a different construction per
  P-65.1.5.1 and is out of scope).
- A ring anywhere in the molecule (P-65.1.1.2's 'carbothioic acid'
  construction, out of scope here), except for one narrow case: a thioic
  acid's branched chain hanging off a single benzene ring bearing only
  halogen/plain-alkyl substituents besides the chain itself
  (`_name_phenyl_chain_thioic_acid`, e.g. '2-phenylethanethioic S-acid'),
  mirroring `_carboxylic_acid.py`'s identical benzene-ring-substituent
  path (P-44.3.2, `longest_branched_chain`). A thioic acid directly on
  the ring stays out of scope for this chain-parent module.
- The selenium/tellurium analogues (selenoic/telluroic acid) -- the
  selenium analogue is handled separately by `_selenoic_acid.py`; tellurium
  is a separate follow-up module.
- Any other heteroatom, charge, or isotopic label. A halogen substituent
  (P-35.2.1) is allowed only on the benzene ring in the ring-substituent
  path above, not elsewhere.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    non_single_bonds,
    ring_chain_attachment_with_halogens,
    ring_chain_attachments_with_halogens,
    separate_aromatic_monocycles,
    unbranched_chain_length,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch, plain_alkyl_ring_substituents

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}


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
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the thioic acid's own chalcogens "
                "and halogen substituents (P-35.2.1) are not supported yet "
                "(P-65.1.5 is restricted to a plain acyclic thioic acid "
                "here)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure("a halogen atom must be a monovalent substituent (P-35.2.1)")
        if atomic_num == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
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
    """Name a thioic acid whose -C(=O)SH/-C(=S)OH lies on a (possibly
    branched) chain hanging off one atom of a benzene ring that otherwise
    bears only halogen/plain-alkyl substituents -- e.g.
    '2-phenylethanethioic S-acid'. The ring is cited as a 'phenyl' (or
    e.g. '4-chlorophenyl') substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`
    (P-44.3.2, `longest_branched_chain` absorbs a branch into the parent
    chain whenever that makes it longer). The acid carbon is always C1 --
    unlike a sulfonic/carboxylic acid group elsewhere in a chain, a
    thioic acid group is by definition a chain terminus, so no locant
    tie-break between chain directions is needed."""
    acid_carbon, label, acid_atom_idxs = _validate_and_collect_thioic_acid(mol, aromatic_ring_atoms=ring_atoms)
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_chain_attachments_with_halogens(graph, rings, set(), halogens)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain thioic acid is not "
            "supported yet"
        )
    acid_carbon_idx = acid_carbon.GetIdx()
    hetero_idxs = acid_atom_idxs - {acid_carbon_idx}
    chain, branches = longest_branched_chain(graph, acid_carbon_idx, ring_atoms, hetero_idxs, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a thioic acid directly attached to the benzene ring uses a "
            "separate naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{alkane_name(chain_length)}thioic {label}-acid"


def name_thioic_acid(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        return _name_phenyl_chain_thioic_acid(mol, set().union(*aromatic_rings))
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
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
        chain_length = unbranched_chain_length(mol, chain_root.GetIdx(), acid_carbon.GetIdx())
        if chain_length is None:
            raise UnsupportedStructure(
                "a branched R group is out of scope for this module (see "
                "module docstring)"
            )
        chain_length += 1
    else:
        chain_length = 1

    return f"{alkane_name(chain_length)}thioic {label}-acid"
