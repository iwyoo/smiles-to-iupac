"""Naming of a diester on a plain unbranched diol chain, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-65.6.3.3.3.1: when both acyl groups are **identical**, the PIN cites
  the diol backbone as a multivalent organyl group ('...diyl') followed
  by the multiplied anion name -- the chapter's own worked example is
  'ethane-1,2-diyl diacetate (PIN)' (`tmp/bluebook/P6.txt` ~line 7118).
  This project's systematic-name convention (see `_ester.py`'s own
  'methyl methanoate'/'methyl ethanoate', not 'methyl formate'/'methyl
  acetate') carries over here too, so this module produces 'ethane-1,2-
  diyl diethanoate'-style names instead of PubChem's retained-name
  'diacetate'.
- P-65.6.3.3.2 (method 2): a **differing**-acid diester instead cites one
  ester as the principal characteristic group (suffix parent, "...yl
  <acid>oate") and the other acyl group as an '<acid>yloxy' substituent
  prefix on the alcohol chain, e.g. 'CC(=O)OCCOC(=O)CC' (mixed acetate/
  propanoate diester) -> '2-acetyloxyethyl propanoate' in PubChem's own
  (retained-name) style, '2-methanoyloxyethyl ethanoate'-style in this
  project's systematic convention -- this is only sanctioned as general
  nomenclature, not the PIN, but there's no PIN alternative implemented
  for the differing-acid case yet (P-65.6.3.3.3.2's method-1 form is a
  separate, not-yet-built mechanism), so it remains this module's output
  there.
- Which acid wins the suffix-parent slot in the differing-acid case:
  confirmed via PubChem PUG REST across three differing-length pairs
  (methanoic/ethanoic, ethanoic/propanoic, ethanoic/butanoic acid) that
  the longer acyl chain is always the suffix parent and the shorter one
  is demoted to the acyloxy prefix -- e.g. 'CCCC(=O)OCCOC(=O)C' (butanoic
  + ethanoic) -> '2-<acid>yloxyethyl butanoate', never the reverse.
- Scope of this first pass, deliberately narrow: exactly two ester
  groups, both acyl chains plain (unbranched, saturated, no
  halogens/stereocenters), and the alcohol backbone connecting the two
  ester oxygens a single plain unbranched saturated carbon chain with no
  other substituent anywhere on it (a real alkyl branch on the backbone,
  e.g. a methyl group alongside the second ester attachment, needs the
  full `name_branch`-style longest-chain/lowest-locant machinery and is
  deferred). A branched backbone, a branched or unsaturated/halogenated
  acyl chain, three or more ester groups, and any ring are all out of
  scope and raise `UnsupportedStructure`.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds, ordered_chain
from ._numerals import alkane_name, alkyl_name, multiplying_prefix

_ALLOWED_ATOMIC_NUMS = {6, 8}


def _find_ester_carbons(mol):
    """Every carbon shaped like an ester acyl carbon (a carbonyl oxygen plus
    a second, carbon-bonded ester oxygen), as (acyl_carbon, carbonyl_oxygen,
    ester_oxygen, alcohol_carbon) tuples."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if len(carbonyls) == 1 and len(ester_oxygens) == 1:
            alcohol_carbon = next(n for n in ester_oxygens[0].GetNeighbors() if n.GetIdx() != atom.GetIdx())
            matches.append((atom, carbonyls[0], ester_oxygens[0], alcohol_carbon))
    return matches


def has_diester_shape(mol) -> bool:
    return len(_find_ester_carbons(mol)) == 2


def _acyl_chain_length(mol, acyl_carbon, ester_oxygen_idx, carbonyl_oxygen_idx):
    graph = adjacency(mol)
    chain = ordered_chain(graph, acyl_carbon.GetIdx(), ester_oxygen_idx, {carbonyl_oxygen_idx})
    if chain is None or any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in chain):
        raise UnsupportedStructure(
            "a branched, unsaturated, or halogen-bearing acyl chain in a "
            "diester is not supported yet (see module docstring)"
        )
    return len(chain)


def _acid_stem(chain_length: int) -> str:
    return alkane_name(chain_length)[:-1]


def name_diester_acyloxy(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the two esters' own oxygens are not "
                "supported yet for a diester (see module docstring)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic ring is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")

    matches = _find_ester_carbons(mol)
    if len(matches) != 2:
        raise UnsupportedStructure("this module only handles exactly two ester groups")
    (acyl_1, carbonyl_1, ester_o_1, alcohol_1), (acyl_2, carbonyl_2, ester_o_2, alcohol_2) = matches

    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 4:
        raise UnsupportedStructure(
            "an oxygen outside the two esters' own carbonyl/ester pairs is "
            "out of scope for this module"
        )
    carbonyl_bonds = {
        frozenset((acyl_1.GetIdx(), carbonyl_1.GetIdx())),
        frozenset((acyl_2.GetIdx(), carbonyl_2.GetIdx())),
    }
    if any(frozenset((a, b)) not in carbonyl_bonds for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure("unsaturation alongside a diester's ester groups is not supported yet")

    graph = adjacency(mol)
    backbone = ordered_chain(graph, alcohol_1.GetIdx(), ester_o_1.GetIdx(), {ester_o_2.GetIdx()})
    if (
        backbone is None
        or backbone[-1] != alcohol_2.GetIdx()
        or any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in backbone)
    ):
        raise UnsupportedStructure(
            "the two esters must share a single plain, unbranched, "
            "saturated carbon chain with no other substituent (see module "
            "docstring)"
        )

    length_1 = _acyl_chain_length(mol, acyl_1, ester_o_1.GetIdx(), carbonyl_1.GetIdx())
    length_2 = _acyl_chain_length(mol, acyl_2, ester_o_2.GetIdx(), carbonyl_2.GetIdx())

    if length_1 == length_2:
        diyl_name = f"{alkane_name(len(backbone))}-1,{len(backbone)}-diyl"
        acid_name = f"{multiplying_prefix(2)}{_acid_stem(length_1)}oate"
        return f"{diyl_name} {acid_name}"

    if length_2 > length_1:
        suffix_length, acyloxy_length = length_2, length_1
        backbone = list(reversed(backbone))
    else:
        suffix_length, acyloxy_length = length_1, length_2

    acyloxy_locant = len(backbone)
    alcohol_chain_name = alkyl_name(len(backbone))
    acyloxy_prefix = f"{acyloxy_locant}-{_acid_stem(acyloxy_length)}oyloxy"
    suffix_acid_name = _acid_stem(suffix_length) + "oate"

    return f"{acyloxy_prefix}{alcohol_chain_name} {suffix_acid_name}"
