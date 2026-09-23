"""Naming of a polyester on a plain unbranched polyol chain, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-65.6.3.3.3.1: when all the acyl groups are **identical** (two or
  more), the PIN cites the polyol backbone as a multivalent organyl
  group ('...diyl'/'...triyl'/'...tetrayl'/...) followed by the
  multiplied anion name -- the chapter's own worked examples are
  'ethane-1,2-diyl diacetate (PIN)' and 'propane-1,2,3-triyl triacetate
  (PIN)' (`tmp/bluebook/P6.txt` ~line 7118/7256). This project's
  systematic-name convention (see `_ester.py`'s own 'methyl methanoate'/
  'methyl ethanoate', not 'methyl formate'/'methyl acetate') carries over
  here too, so this module produces 'ethane-1,2-diyl diethanoate'/
  'propane-1,2,3-triyl triethanoate'-style names instead of PubChem's
  retained-name 'diacetate'/'triacetate'.
- P-65.6.3.3.3.2 (method 1): when the acyl groups **differ** (two or
  more), the PIN instead cites each distinct acid as its own anion name,
  in alphanumerical order, each with its own locant set and multiplying
  prefix (`di`/`tri`/...) if it covers more than one position -- e.g.
  'propane-1,2,3-triyl 1,2-diacetate 3-propanoate (PIN)'. For exactly two
  esters, no acid locant is cited at all: the backbone's own two-fold
  symmetry (either numbering direction describes the identical molecule)
  makes one redundant, confirmed directly from the Blue Book's own
  'methylene acetate formate (PIN)' and '1,4-phenylene acetate
  dichloroacetate (PIN)' examples, neither of which cites one. Three or
  more positions lose that symmetry (an interior position is
  constitutionally distinct from the ends), so every group's locant(s)
  are always cited there. Method 2 (citing one ester as the suffix parent
  and the rest as '<acid>yloxy' prefixes) is sanctioned general
  nomenclature only, not implemented here -- this module always returns
  the method-1 PIN for the differing-acid case.
- Scope, deliberately narrow: two or more ester groups, all acyl chains
  plain (unbranched, saturated, no halogens/stereocenters), and every
  ester's alcohol-side carbon lying on a single plain unbranched
  saturated carbon chain with no other substituent anywhere on it (a
  real alkyl branch on the backbone, e.g. a methyl group alongside an
  ester attachment, needs the full `name_branch`-style longest-chain/
  lowest-locant machinery and is deferred). A backbone with a branch
  point/quaternary carbon (e.g. pentaerythritol's tetrahedral center) is
  its own separate, not-yet-scoped shape -- not attempted here even for
  identical acids, since it isn't a single chain at all. A branched
  backbone, a branched or unsaturated/halogenated acyl chain, and any
  ring are all out of scope and raise `UnsupportedStructure`.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds, ordered_chain
from ._numerals import alkane_name, multiplying_prefix

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
    return len(_find_ester_carbons(mol)) >= 2


def _find_full_chain(mol, graph, matches):
    """The single plain, unbranched, saturated carbon chain that carries
    every one of `matches`' alcohol-side carbons, tried from each match in
    turn as the starting end -- None if no such chain exists (a genuine
    branch point, a heteroatom on the path, or the esters simply aren't
    all on one chain). Excluding every *other* match's own ester oxygen
    from the walk (not just the two-ester case's single "other" oxygen)
    is what lets an interior attachment point (e.g. glycerol's middle
    carbon) pass through without `ordered_chain` mistaking its own ester
    oxygen for a second, branching neighbor."""
    ester_oxygens = {m[2].GetIdx() for m in matches}
    alcohol_atoms = {m[3].GetIdx() for m in matches}
    for _acyl, _carbonyl, start_oxygen, start_alcohol in matches:
        excluded = ester_oxygens - {start_oxygen.GetIdx()}
        chain = ordered_chain(graph, start_alcohol.GetIdx(), start_oxygen.GetIdx(), excluded)
        if chain is None:
            continue
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in chain):
            continue
        if chain[-1] not in alcohol_atoms:
            # The chain must terminate exactly at an ester's own alcohol
            # carbon, not run past it into a further plain substituent --
            # that plain tail is a real alkyl branch off the ester
            # backbone (out of scope, see module docstring), not part of
            # a longer multivalent-organyl parent chain.
            continue
        if {idx for idx in chain if idx in alcohol_atoms} != alcohol_atoms:
            continue
        return chain
    return None


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
                "heteroatoms other than the esters' own oxygens are not "
                "supported yet for a polyester (see module docstring)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic ring is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")

    matches = _find_ester_carbons(mol)
    if len(matches) < 2:
        raise UnsupportedStructure("this module only handles two or more ester groups")

    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    if total_oxygens != 2 * len(matches):
        raise UnsupportedStructure(
            "an oxygen outside the esters' own carbonyl/ester pairs is "
            "out of scope for this module"
        )
    carbonyl_bonds = {frozenset((acyl.GetIdx(), carbonyl.GetIdx())) for acyl, carbonyl, _, _ in matches}
    if any(frozenset((a, b)) not in carbonyl_bonds for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure("unsaturation alongside a polyester's ester groups is not supported yet")

    graph = adjacency(mol)
    backbone = _find_full_chain(mol, graph, matches)
    if backbone is None:
        raise UnsupportedStructure(
            "the esters must share a single plain, unbranched, saturated "
            "carbon chain with no other substituent (see module docstring)"
        )

    lengths = [_acyl_chain_length(mol, acyl, ester_o.GetIdx(), carbonyl.GetIdx()) for acyl, carbonyl, ester_o, _ in matches]

    alcohol_atoms = {alcohol.GetIdx() for _, _, _, alcohol in matches}
    best_locants = None
    best_position_of = None
    for candidate in (backbone, list(reversed(backbone))):
        position_of = {idx: i + 1 for i, idx in enumerate(candidate)}
        locants = sorted(position_of[idx] for idx in alcohol_atoms)
        if best_locants is None or locants < best_locants:
            best_locants, best_position_of = locants, position_of

    yl_prefix = multiplying_prefix(len(matches))
    locant_str = ",".join(str(loc) for loc in best_locants)
    group_name = f"{alkane_name(len(backbone))}-{locant_str}-{yl_prefix}yl"

    if len(set(lengths)) == 1:
        acid_name = f"{yl_prefix}{_acid_stem(lengths[0])}oate"
        return f"{group_name} {acid_name}"

    # P-65.6.3.3.3.2 method 1: differing acyl groups are cited as separate
    # anion names, each with its own locant set and multiplying prefix (if
    # repeated), in alphanumerical order of the acid name -- the PIN, in
    # place of the acyloxy-prefix method 2 this module used to return here
    # (general nomenclature only, still `_diester_acyloxy.py`'s territory
    # for reference in the module docstring, not this function's output).
    groups = {}
    for length, (_, _, _, alcohol) in zip(lengths, matches):
        groups.setdefault(length, []).append(best_position_of[alcohol.GetIdx()])

    # For exactly two esters, the backbone's own two-fold symmetry (either
    # numbering direction describes the identical molecule) means no
    # locant is needed to know which acid sits where -- confirmed against
    # the Blue Book's own 'methylene acetate formate (PIN)' and
    # '1,4-phenylene acetate dichloroacetate (PIN)' examples, neither of
    # which cites an acid locant. Three or more positions lose that
    # symmetry (an interior position is constitutionally distinct from
    # the ends), so every group's locant(s) must be cited there, per the
    # 'propane-1,2,3-triyl 1,2-diacetate 3-propanoate (PIN)' example.
    cite_locants = len(matches) != 2

    acid_parts = []
    for length in sorted(groups, key=lambda length: _acid_stem(length)):
        locants = sorted(groups[length])
        count = len(locants)
        prefix = multiplying_prefix(count) if count > 1 else ""
        acid_name = f"{prefix}{_acid_stem(length)}oate"
        if cite_locants:
            loc_str = ",".join(str(loc) for loc in locants)
            acid_name = f"{loc_str}-{acid_name}"
        acid_parts.append(acid_name)

    return f"{group_name} " + " ".join(acid_parts)
