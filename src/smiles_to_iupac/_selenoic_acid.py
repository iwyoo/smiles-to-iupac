"""Naming of selenoic acids (-C(=O)-SeH or -C(=Se)-OH), the selenium
analogue of a thioic acid, restricted to a single such group on an acyclic
unbranched saturated hydrocarbon chain, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-65.1.5 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Replacement of oxygen atom(s) of a carboxylic acid group by another
  chalcogen is indicated by the affixes 'thio', 'seleno', and 'telluro'" --
  i.e. exactly `_thioic_acid.py`'s substitutive rule and letter-labeling
  mechanism (see that module's docstring), with 'seleno' standing in for
  'thio' and the label 'Se' standing in for 'S'.
- Confirmed via PubChem for the S-acid-direction tautomer (the acidic
  hydrogen on the chalcogen, C=O retained): 'CC(=O)[SeH]' resolves to CID
  57348413, auto-generated name 'ethaneselenoic Se-acid'.
- The O-acid-direction tautomer (acidic hydrogen on oxygen, C=Se) has no
  PubChem-computed name (its CID, 12189546, carries only a bare
  connectivity record) -- this direction is a reviewed, not directly
  PubChem-verified, mechanical extension of the same labeling rule
  `_thioic_acid.py` already established from direct Blue Book quotations
  (confirmed generally applicable to every chalcogen by P-65.1.5.1's own
  wording, and used identically elsewhere in the same section, e.g.
  'carbamoselenoic Se-acid (PIN)').
- Formic acid's analogue (chain length 1, R = H) uses the 'methane' stem,
  mirroring `_thioic_acid.py`'s identical treatment.

Explicitly out of scope (raise `UnsupportedStructure`): same list as
`_thioic_acid.py` -- a branched or unsaturated R group on the plain
acyclic path, more than one selenoic acid group, a ring anywhere in the
molecule, any other heteroatom, charge, or isotopic label, except for one
narrow case: a selenoic acid's branched chain hanging off a single
benzene ring bearing only halogen/plain-alkyl substituents besides the
chain itself (`_name_phenyl_chain_selenoic_acid`, e.g.
'2-phenylethaneselenoic Se-acid'), mirroring `_thioic_acid.py`'s
identical benzene-ring-substituent path. A selenoic acid directly on the
ring stays out of scope for this chain-parent module. The tellurium
analogue (telluroic acid) is handled separately by `_telluroic_acid.py`.
A halogen substituent (P-35.2.1) is allowed only on the benzene ring in
the ring-substituent path, not elsewhere.
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
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch, plain_alkyl_ring_substituents

_ALLOWED_ATOMIC_NUMS = {6, 8, 34, *HALOGEN_PREFIXES}


def _selenoic_acid_carbons(mol):
    """Carbons shaped like a selenoic-acid group: bonded to exactly one
    double-bonded (terminal) O or Se, and exactly one single-bonded,
    one-H, terminal O or Se of the other element. Returns a list of
    (carbon, label, double_atom, single_atom) tuples, where label is 'O'
    or 'Se' -- whichever chalcogen carries the acidic hydrogen."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, 34)]
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
        if {double_atom.GetAtomicNum(), single_atom.GetAtomicNum()} != {8, 34}:
            continue
        label = "O" if single_atom.GetAtomicNum() == 8 else "Se"
        matches.append((atom, label, double_atom, single_atom))
    return matches


def has_selenoic_acid_shape(mol) -> bool:
    return bool(_selenoic_acid_carbons(mol))


def _unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None if
    the chain branches, rings, or leaves carbon at any point. Mirrors
    `_thioic_acid.py`'s identical helper."""
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


def _validate_and_collect_selenoic_acid(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_selenoic_acid`'s benzene-
    ring-substituent path (see `_name_phenyl_chain_selenoic_acid`) can
    reuse this same validation for the rest of the molecule. Empty by
    default, so every other caller's behavior is unchanged. Returns
    (acid_carbon, label, acid_atom_idxs)."""
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the selenoic acid's own chalcogens "
                "and halogen substituents (P-35.2.1) are not supported yet "
                "(P-65.1.5 is restricted to a plain acyclic selenoic acid "
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
    matches = _selenoic_acid_carbons(mol)
    if len(matches) != 1:
        raise UnsupportedStructure("more than one selenoic acid group is out of scope for this module")
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


def _name_phenyl_chain_selenoic_acid(mol, ring_atoms):
    """Name a selenoic acid whose -C(=O)SeH/-C(=Se)OH lies on a (possibly
    branched) chain hanging off one atom of a benzene ring that otherwise
    bears only halogen/plain-alkyl substituents -- e.g.
    '2-phenylethaneselenoic Se-acid'. The ring is cited as a 'phenyl' (or
    e.g. '4-chlorophenyl') substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_thioic_acid.py`'s `_name_phenyl_chain_thioic_acid`
    (P-44.3.2, `longest_branched_chain`). The acid carbon is always C1,
    a chain terminus by definition, so no locant tie-break between chain
    directions is needed."""
    acid_carbon, label, acid_atom_idxs = _validate_and_collect_selenoic_acid(mol, aromatic_ring_atoms=ring_atoms)
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain selenoic acid is not "
            "supported yet"
        )
    acid_carbon_idx = acid_carbon.GetIdx()
    hetero_idxs = acid_atom_idxs - {acid_carbon_idx}
    chain, branches = longest_branched_chain(graph, acid_carbon_idx, ring_atoms, hetero_idxs, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a selenoic acid directly attached to the benzene ring uses a "
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
    return f"{prefix}{alkane_name(chain_length)}selenoic {label}-acid"


def name_selenoic_acid(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_selenoic_acid(mol, ring_atoms)
    acid_carbon, label, acid_atom_idxs = _validate_and_collect_selenoic_acid(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("rings are not supported by this module yet")

    chain_neighbors = [n for n in acid_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(chain_neighbors) > 1:
        raise UnsupportedStructure("a selenoic-acid carbon with more than one carbon neighbor is not valid")

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

    return f"{alkane_name(chain_length)}selenoic {label}-acid"
