"""Naming of carbamate esters (R-O-C(=O)-NHR'), per the IUPAC 2013
Recommendations ("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'carbamic
  acid' (H2N-COOH) is a retained name that is itself the preferred IUPAC
  name -- not a systematic '...oic acid' derivative -- so its esters are
  named the same way any other retained-acid ester is: 'R carbamate', two
  words, the alkyl group R cited first exactly as `_ester.py`'s alcohol
  part is, e.g. 'ethyl carbamate' for CH3CH2-O-CO-NH2. Unlike
  `_ester.py`'s acyl part (an open-ended acid chain, 'oate'), the acid
  side here is fixed: it is always the single word 'carbamate', with no
  chain-length logic of its own.
- This module's scope mirrors `_ester.py`'s own first pass for its alcohol
  part: R is restricted to a plain, unbranched, unsubstituted, saturated
  alkyl group attached at its own chain terminus (e.g. 'methyl', 'ethyl',
  'propyl'); a branched, substituted, unsaturated, or ring-bearing R is
  deferred.
- The amide nitrogen is unsubstituted -NH2, mono-substituted, or
  N,N-disubstituted with one or two plain, unsubstituted, saturated,
  acyclic alkyl groups (branched or unbranched), each cited as its own
  'N-'-prefixed substituent directly ahead of 'carbamate' (P-16.3.3/
  P-66.1, mirroring how an amide's N-substituents are cited), in
  alphanumerical order (P-14.5.2, ignoring italicized prefixes like
  'tert-' -- via `alpha_sort_key`), with a 'di' multiplying prefix (and a
  single shared 'N,N-' locant pair) when both substituents are identical:
  'methyl N-methylcarbamate' for CH3-NH-CO-O-CH3 (PubChem CID 81151),
  'methyl N,N-dimethylcarbamate' for (CH3)2N-CO-O-CH3, 'methyl
  N-ethyl-N-methylcarbamate' for CH3(C2H5)N-CO-O-CH3, 'methyl
  N-propan-2-ylcarbamate' (CID 568334), 'methyl N-tert-butylcarbamate'
  (CID 575684), 'methyl N-tert-butyl-N-ethylcarbamate' (CID 87089877,
  alphabetized as 'b' before 'e', ignoring 'tert-') -- structures
  confirmed via PubChem, which returns these exact names. A compound
  identical-pair 'di' name is parenthesized to avoid ambiguity ('methyl
  N,N-di(propan-2-yl)carbamate', CID 568417) while a non-compound
  (retained-name) one is not ('methyl N,N-ditert-butylcarbamate', CID
  12567954) -- a single N-substituent is never parenthesized either way.
  A cyclic/unsaturated N-substituent and ring-attached amide nitrogens
  are still deferred.
- P-29.3.2.1: both R's and R''s names are built with `name_branch` (P-29
  PIN style, fixed project-wide by PR #237 -- previously this module
  avoided `name_branch` for exactly this reason, but that blocker no
  longer applies), e.g. 'propan-2-yl' for R = isopropyl, matching the
  verified PIN 'propan-2-yl carbamate' (PubChem CID 15628) and
  'tert-butyl carbamate' (CID 77922). R itself is never parenthesized
  regardless of `name_branch`'s `is_compound` flag -- the 'R carbamate'
  two-word pattern (mirroring `_ester.py`'s alcohol part) has no
  nested-prefix ambiguity to guard against; each individual N-substituent
  is likewise never parenthesized -- only the multiplied 'N,N-di(...)'
  form needs it (see above).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule.
- An N-substituent that is unsaturated or ring-bearing (a branched but
  otherwise plain saturated acyclic N-substituent is supported, see
  above).
- An unsaturated or cyclic R (a branched but otherwise plain saturated
  acyclic R is supported, see above).
- More than one carbamate group, or any other heteroatom/oxygen not part
  of this single carbamate group (an ether, alcohol, or second carbonyl
  elsewhere) -- including free carbamic acid itself (H2N-COOH, R = H),
  which is a distinct parent-hydride-less special case deferred entirely.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, non_single_bonds
from ._substituents import alpha_sort_key, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8}


def _carbamate_cores(mol):
    """List of (carbamate_carbon, carbonyl_o, ester_o, alkyl_carbon,
    amide_n, n_alkyl_carbons) for every R-O-C(=O)-N(R'')R' pattern: a carbon
    carrying exactly one doubly-bonded (terminal) oxygen, one singly-bonded
    oxygen itself bonded to a second carbon, and one singly-bonded
    nitrogen that is unsubstituted (degree 1, two H), mono-substituted
    (degree 2, one H), or N,N-disubstituted (degree 3, no H) with carbon --
    `n_alkyl_carbons` is a tuple of those carbons' indices (0, 1, or 2 of
    them)."""
    cores = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
        if len(oxygens) != 2 or len(nitrogens) != 1:
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
        if len(carbonyls) != 1 or len(ester_oxygens) != 1:
            continue
        (amide_n,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), amide_n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        n_substituents = [n for n in amide_n.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        n_alkyl_carbons = tuple(n.GetIdx() for n in n_substituents if n.GetAtomicNum() == 6)
        if len(n_alkyl_carbons) != len(n_substituents):
            continue
        expected_hs = 2 - len(n_alkyl_carbons)
        if amide_n.GetDegree() != 1 + len(n_alkyl_carbons) or amide_n.GetTotalNumHs() != expected_hs:
            continue
        ester_oxygen = ester_oxygens[0]
        alkyl_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        cores.append((
            atom.GetIdx(),
            carbonyls[0].GetIdx(),
            ester_oxygen.GetIdx(),
            alkyl_carbon.GetIdx(),
            amide_n.GetIdx(),
            n_alkyl_carbons,
        ))
    return cores


def has_carbamate_shape(mol) -> bool:
    return bool(_carbamate_cores(mol))


def name_carbamate(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-attached carbamate is out of scope for this "
            "acyclic-only module"
        )
    cores = _carbamate_cores(mol)
    if len(cores) != 1:
        raise UnsupportedStructure(
            "exactly one carbamate group is required; zero or multiple "
            "carbamate groups are not supported yet"
        )
    carbamate_c, carbonyl_o, ester_o, alkyl_c, amide_n, n_alkyl_cs = cores[0]

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the carbamate's own oxygens/"
                "nitrogen are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num == 7:
            if atom.GetIdx() != amide_n:
                raise UnsupportedStructure(
                    "a nitrogen other than the carbamate's own -NH2 needs "
                    "Table 3.3 seniority handling not yet implemented here"
                )
        elif atom.GetIdx() not in (carbonyl_o, ester_o):
            raise UnsupportedStructure(
                "an oxygen other than the carbamate's own two oxygens (a "
                "coexisting ether/alcohol/second carbonyl) is out of scope "
                "for this module"
            )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon alkyl "
            "group to name (free carbamic acid, R = H, is out of scope)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    other_non_single = [b for b in non_single_bonds(mol) if carbamate_c not in (b[0], b[1])]
    if other_non_single:
        raise UnsupportedStructure(
            "unsaturation in the R group is not supported yet"
        )

    full_graph = adjacency(mol)
    r_name, _ = name_branch(full_graph, alkyl_c, ester_o, {})
    if not n_alkyl_cs:
        return f"{r_name} carbamate"

    n_entries = [name_branch(full_graph, n_alkyl_c, amide_n, {}) for n_alkyl_c in n_alkyl_cs]

    if len(n_entries) == 2 and n_entries[0][0] == n_entries[1][0]:
        name, is_compound = n_entries[0]
        di_name = f"({name})" if is_compound else name
        n_prefix = f"N,N-di{di_name}"
    else:
        n_prefix = "-".join(f"N-{name}" for name, _ in sorted(n_entries, key=lambda e: alpha_sort_key(e[0])))
    return f"{r_name} {n_prefix}carbamate"
