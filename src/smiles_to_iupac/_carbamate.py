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
- The amide nitrogen is either unsubstituted -NH2, or mono-substituted
  with a single plain, unbranched, unsubstituted, saturated alkyl group
  R' cited as an 'N-'-prefixed substituent directly ahead of 'carbamate'
  (P-16.3.3/P-66.1, mirroring how an amide's N-substituent is cited):
  'methyl N-methylcarbamate' for CH3-NH-CO-O-CH3 (PubChem CID 81151).
  N,N-disubstitution, a branched/cyclic/unsaturated N-substituent, and
  ring-attached amide nitrogens are all deferred.
- P-29.3.2.1: both R's and R''s names are plain alkyl substituent-group
  names, built with `alkyl_name` directly since both are restricted to an
  unbranched chain here -- deliberately not using `name_branch` for R,
  since that function still produces pre-PIN branched-substituent names
  (e.g. '1-methylethyl' instead of 'propan-2-yl'; see the P-29 blocker
  tracked in the roadmap) that would mismatch a verified PIN like
  'propan-2-yl carbamate' or 'tert-butyl carbamate'.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule.
- N,N-disubstitution, or an N-substituent that is branched, unsaturated,
  or ring-bearing.
- A branched, substituted, unsaturated, or cyclic R -- only a plain
  unbranched saturated alkyl R is supported in this first pass.
- More than one carbamate group, or any other heteroatom/oxygen not part
  of this single carbamate group (an ether, alcohol, or second carbonyl
  elsewhere) -- including free carbamic acid itself (H2N-COOH, R = H),
  which is a distinct parent-hydride-less special case deferred entirely.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    carbon_adjacency,
    linear_branch,
    non_single_bonds,
)
from ._numerals import alkyl_name

_ALLOWED_ATOMIC_NUMS = {6, 7, 8}


def _carbamate_cores(mol):
    """List of (carbamate_carbon, carbonyl_o, ester_o, alkyl_carbon,
    amide_n, n_alkyl_carbon) for every R-O-C(=O)-NHR' pattern: a carbon
    carrying exactly one doubly-bonded (terminal) oxygen, one singly-bonded
    oxygen itself bonded to a second carbon, and one singly-bonded
    nitrogen that is either unsubstituted (degree 1, two H) or
    mono-substituted with a single carbon (degree 2, one H) -- `n_alkyl_carbon`
    is that carbon's index, or None for the unsubstituted -NH2 case."""
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
        n_alkyl_carbon = None
        if amide_n.GetDegree() == 1 and amide_n.GetTotalNumHs() == 2:
            pass
        elif amide_n.GetDegree() == 2 and amide_n.GetTotalNumHs() == 1:
            (n_alkyl_carbon,) = (
                n for n in amide_n.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom.GetIdx()
            )
        else:
            continue
        ester_oxygen = ester_oxygens[0]
        alkyl_carbon = next(n for n in ester_oxygen.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        cores.append((
            atom.GetIdx(),
            carbonyls[0].GetIdx(),
            ester_oxygen.GetIdx(),
            alkyl_carbon.GetIdx(),
            amide_n.GetIdx(),
            n_alkyl_carbon.GetIdx() if n_alkyl_carbon is not None else None,
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
    carbamate_c, carbonyl_o, ester_o, alkyl_c, amide_n, n_alkyl_c = cores[0]

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

    carbon_graph = carbon_adjacency(mol)
    length = linear_branch(carbon_graph, alkyl_c, None)
    if length is None:
        raise UnsupportedStructure(
            "a branched R group is not supported yet"
        )
    if n_alkyl_c is None:
        return f"{alkyl_name(length)} carbamate"

    n_length = linear_branch(carbon_graph, n_alkyl_c, None)
    if n_length is None:
        raise UnsupportedStructure(
            "a branched N-substituent is not supported yet"
        )
    return f"{alkyl_name(length)} N-{alkyl_name(n_length)}carbamate"
