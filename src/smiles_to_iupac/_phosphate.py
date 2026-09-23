"""Naming of trialkyl phosphate esters (P(=O)(OR)3, R groups identical or
mixed), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group(s) as separate word(s), in alphanumeric
  order if more than one, followed by the acid's anion name. Confirmed
  worked example: `P(O-CH3)3` -> phosphite... no wait, the P=O analogue
  `P(O)(O-CH3)3` -> 'trimethyl phosphate (PIN)' (`tmp/bluebook/P6a.txt`
  line ~4053-4054). For three identical R groups this is a single
  multiplying-prefixed word ('trimethyl'), not three repeated words -- the
  ordinary P-14.2.1/P-14.2.2 multiplying-prefix convention already used
  throughout this project (`_numerals.multiplying_prefix`), not a
  substituent-prefix-on-a-parent-hydride construction.
- Mixed-alkyl esters (P-14.5.2): each distinct R name is its own word
  (multiplying-prefixed if it appears more than once), words cited in
  alphanumeric order via `alpha_sort_key` -- e.g. `ethyl methyl phenyl
  phosphate` (PubChem CID 12494385/12494386, all three distinct) and
  `diethyl methyl phosphate` (PubChem CID 120420, two identical + one
  distinct). Both PubChem-computed names are directly usable here (no
  von-Baeyer ambiguity, unlike a ring parent).

Scope: a single phosphorus atom shaped like a phosphate ester -- one P=O
double bond, three P-O-R single bonds, each R named via `name_branch` (a
plain alkyl chain, a branched chain, or a plain benzene ring, and their
halogenated variants, exactly like `_phosphonic_acid.py`'s R), identical
or mixed freely. Explicitly out of scope (raise `UnsupportedStructure`):
partial ("hydrogen") esters (a remaining P-OH), any chalcogen-replacement
analogue, and any other heteroatom.

The word-assembly step (group identical R names, multiplying-prefix +
enclosure per group, alphanumeric word order) is shared with the
phosphite ester case (no P=O, P(III) instead of P(V)) via
`format_ester_words`, exported for `_phosphite.py` to reuse.
"""

from rdkit import Chem

from collections import Counter

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._numerals import multiplying_prefix
from ._substituents import name_branch

_PHOSPHORUS = 15


def _phosphate_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a fully-esterified phosphate: bonded
    to exactly one double-bonded (terminal) oxygen and exactly three
    single-bonded ester oxygens (each degree 2 -- P plus one R carbon;
    this excludes a hydroxyl oxygen, degree 1, i.e. a partial ester)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _PHOSPHORUS or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygens) != 4:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 1 or len(ester_os) != 3:
            continue
        (double_o,) = double_os
        if double_o.GetDegree() != 1:
            continue
        if any(o.GetDegree() != 2 for o in ester_os):
            continue
        matches.append(atom)
    return matches


def has_phosphate_shape(mol) -> bool:
    return bool(_phosphate_phosphorus_atoms(mol))


def name_phosphate(mol) -> str:
    phosphorus_atoms = _phosphate_phosphorus_atoms(mol)
    if len(phosphorus_atoms) != 1:
        raise UnsupportedStructure("more than one phosphate group is not supported yet")
    (phosphorus,) = phosphorus_atoms

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _PHOSPHORUS and atom.GetIdx() != phosphorus.GetIdx():
            raise UnsupportedStructure("more than one phosphorus atom is not supported yet")
        if atomic_num not in (1, 6, 8, _PHOSPHORUS, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the phosphate's own phosphorus/"
                "oxygens and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetIdx() not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the phosphate's own "
                "P(=O)(OR)3 group is out of scope for this module"
            )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    double_o = next(
        n
        for n in phosphorus.GetNeighbors()
        if n.GetAtomicNum() == 8
        and mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    )
    ester_oxygens = [idx for idx in group_oxygens if idx != double_o.GetIdx()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != phosphorus.GetIdx()]
        roots.append((root, oxygen_idx))

    names = [name_branch(graph, root, coming_from, halogens, aromatic_atoms, mol=mol)[0] for root, coming_from in roots]
    return format_ester_words(names) + " phosphate"


def format_ester_words(names) -> str:
    """P-67.1.3.2's own worked examples confirm this ester word-citation
    style is NOT the same as the substituent-prefix-on-a-parent-hydride
    convention (`format_mononuclear_prefixes`'s P-16.5.1.3.1 rule, which
    would parenthesize even a lone compound name): PubChem's own
    PIN-matching name for a branched-but-internally-locanted R gives
    'tripropan-2-yl phosphate' (no parens, plain 'tri'), while a R whose
    own name starts with a locant digit ('2-chloroethyl') gives
    'tris(2-chloroethyl) phosphate' -- enclosing marks are needed there
    purely to keep the leading digit from reading as part of the
    multiplying term itself, confirmed against real PubChem structures for
    both shapes. A mixed-alkyl ester (P-14.5.2) groups identical R names
    under one multiplied word each, the words then cited in alphanumeric
    order (e.g. 'diethyl methyl phosphate', PubChem CID 120420) -- three
    identical names collapse to the same single-word case. Shared with
    `_phosphite.py` (same P-67.1.3.2 citation style, a different anion
    word appended by the caller)."""
    counts = Counter(names)
    words = []
    for name in sorted(counts, key=alpha_sort_key):
        count = counts[name]
        needs_enclosure = name[0].isdigit()
        group = f"({name})" if needs_enclosure else name
        prefix = multiplying_prefix(count, compound=needs_enclosure) if count > 1 else ""
        words.append(f"{prefix}{group}")
    return " ".join(words)
