"""Naming of simple phosphine oxides ("phosphanones", R-P(=O)< bearing
1-3 unbranched, saturated alkyl and/or plain phenyl substituents on the
phosphorus), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-68.3.2.3.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  a P=O group on a phosphane parent is named substitutively with the
  '-one' suffix -- the same suffix ketones use, just on 'phosphane'
  instead of an alkane, eliding the final 'e' before the vowel-starting
  suffix (P-16.3.3): 'phosphane' + 'one' -> 'phosphanone'. This outranks
  functional class nomenclature ('phosphane oxide'), confirmed by the
  Blue Book's own direct worked examples: 'C6H5-P=O' ->
  'phenylphosphanone (PIN)' [explicitly "not oxo(phenyl)phosphane"], and
  '(C6H5)3P=O' -> 'triphenyl-λ5-phosphanone (PIN)' [explicitly "not
  oxotriphenyl-λ5-phosphane"; 'triphenylphosphane oxide' is cited
  alongside as an acceptable but non-preferred functional-class name].
  PubChem's own auto-generated names for this shape (e.g.
  'phosphorosomethane' for `CP=O`) use a completely different, non-PIN
  retained-name pattern and are not usable as a cross-check here -- the
  same kind of PubChem-autoname-vs-PIN mismatch already documented
  elsewhere in this project (see e.g. `_isocyanate.py`), but this one is
  settled by a *direct* Blue Book worked-example quotation rather than an
  analogy to a different group.
- The λ-convention (P-14.1): normal trivalent phosphorus (the P=O double
  bond plus at most one more single-bonded substituent, total bond order
  <= 3) needs no λ label -- confirmed by the phenyl (single-substituent)
  example above. Three single-bonded substituents plus the P=O double
  bond gives a total bond order of 5, exceeding normal valence, so the
  nonstandard valence is flagged with a 'λ5' label inserted (hyphens on
  both sides) directly before 'phosphanone' -- confirmed by the
  triphenyl example above. Two single-bonded substituents plus P=O (total
  bond order 4; RDKit fills the remaining valence with one implicit H to
  reach 5, per `Chem.MolFromSmiles`'s own default valence model for P) is
  treated the same way by direct analogy -- reviewed, not independently
  confirmed by a worked example of its own, since the mechanism (cite
  whatever real substituents are present, silently default everything
  else to implicit H, same as `_phosphane.py` already does for 0-3
  substituents on a plain trivalent phosphane) doesn't depend on which
  exact count triggers it.
- Substituent prefix assembly (identical/different substituents,
  P-16.5.1.3.1 parenthesization) reuses `_phosphane.py`'s own
  `format_mononuclear_prefixes` unchanged -- the alkyl-substituent
  collection logic (linear, unbranched, saturated alkyl chains) is copied
  from that module too.
- A plain, unsubstituted benzene ring bonded directly to the phosphorus
  is cited as a 'phenyl' substituent -- the Blue Book's own two worked
  examples above are exactly this shape (a lone phenyl, and three
  identical phenyls), so it's supported outright rather than deferred.
  Detection mirrors `_carbamate.py`'s `plain_phenyl_substituent_atoms`
  (`is_plain_benzene_ring` + `ring_chain_attachment`), adapted for
  phosphorus's own direct-neighbor substituent roots instead of a
  chain-root carbon one hop away from the parent heteroatom.

Explicitly out of scope (raise `UnsupportedStructure`):
- A halogen substituent on phosphorus (unverified for this shape;
  `_phosphane.py`'s own halogen support is not reused here), or an
  aromatic substituent other than a plain, unsubstituted phenyl (a
  substituted or heteroaromatic ring, e.g.).
- Zero substituents (bare 'phosphanone', an exotic/unconfirmed shape).
- More than one phosphorus atom, more than one P=O group, any other
  heteroatom, a branched or unsaturated (non-aromatic) substituent, a
  ring other than a plain phenyl substituent, or charged/isotopically
  modified atoms.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    linear_branch,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
from ._numerals import alkyl_name
from ._substituents import format_mononuclear_prefixes


def _find_phosphanone_phosphorus(mol):
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return None
    (phosphorus,) = phosphorus_atoms
    if phosphorus.IsInRing():
        return None
    oxo_bonds = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0 and phosphorus.GetIdx() in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    ]
    if len(oxo_bonds) != 1:
        return None
    (bond,) = oxo_bonds
    oxo = bond.GetOtherAtom(phosphorus)
    if oxo.GetAtomicNum() != 8 or oxo.GetDegree() != 1 or oxo.GetFormalCharge() != 0:
        return None
    return phosphorus, oxo


def has_phosphanone_shape(mol) -> bool:
    return _find_phosphanone_phosphorus(mol) is not None


def _validate_and_collect_substituents(mol):
    found = _find_phosphanone_phosphorus(mol)
    if found is None:
        raise UnsupportedStructure(
            "no plain phosphine oxide (R-P(=O)<) shape found; this module "
            "only handles phosphanones"
        )
    phosphorus, oxo = found
    graph = adjacency(mol)
    roots = {r for r in graph[phosphorus.GetIdx()] if r != oxo.GetIdx()}
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots)

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (15, 6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the phosphanone's own phosphorus "
                "and oxo oxygen are not supported yet"
            )
        if atom.GetAtomicNum() == 8 and atom.GetIdx() != oxo.GetIdx():
            raise UnsupportedStructure(
                "an oxygen atom not shaped like the phosphanone's own P=O "
                "is out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and atom.GetIdx() not in phenyl_atoms:
            raise UnsupportedStructure(
                "an aromatic substituent other than a plain, unsubstituted "
                "phenyl group is out of scope for this module"
            )
    all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if all_ring_atoms - phenyl_atoms:
        raise UnsupportedStructure(
            "a ring other than a plain phenyl substituent directly on "
            "phosphorus is out of scope for this module"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in phenyl_atoms and b[1] not in phenyl_atoms
    ]
    if len(non_ring_unsaturation) != 1:
        raise UnsupportedStructure(
            "unsaturation other than the phosphanone's own P=O (and a "
            "plain phenyl substituent's own aromaticity) is out of scope "
            "for this module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    substituent_names = []
    for root in roots:
        if root in phenyl_atoms:
            substituent_names.append(("phenyl", False))
            continue
        length = linear_branch(graph, root, phosphorus.GetIdx())
        if length is None:
            raise UnsupportedStructure("a branched substituent is out of scope for this module")
        substituent_names.append((alkyl_name(length), False))
    if not substituent_names:
        raise UnsupportedStructure(
            "a phosphanone with no substituents on phosphorus (bare "
            "'phosphanone') is not supported yet"
        )
    return substituent_names


def name_phosphanone(mol) -> str:
    substituent_names = _validate_and_collect_substituents(mol)
    prefix = format_mononuclear_prefixes(substituent_names)
    if len(substituent_names) >= 2:
        return f"{prefix}-λ5-phosphanone"
    return prefix + "phosphanone"
