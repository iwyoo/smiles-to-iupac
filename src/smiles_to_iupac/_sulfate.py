"""Naming of dialkyl sulfate esters (S(=O)(=O)(OR)2, R groups identical or
mixed), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group(s) as separate word(s), in alphanumeric
  order if more than one, followed by the acid's anion name -- the exact
  same rule and citation style `_phosphate.py`/`_phosphite.py` already
  established, just with sulfur's own two double-bonded oxygens instead
  of phosphorus's one, and the anion word "sulfate".
- Word-assembly logic (grouping identical R names, multiplying-prefix +
  enclosure per group, alphanumeric word order) is shared verbatim with
  `_phosphate.py` via its exported `format_ester_words`.
- Partial esters (P-67.1.3.2): sulfate is dibasic, so a partial ester
  always has exactly one R group and one remaining S-OH, cited by
  inserting the word "hydrogen" (never "dihydrogen") between the
  R-group word and "sulfate" -- mirrors `_phosphate.py`'s own
  partial-ester citation, confirmed worked example "CH3-O-SO2-OH ->
  methyl hydrogen sulfate (PIN)" (the Blue Book).
- Salts of partial esters (P-67.1.3.2): a single deprotonated S-O^-
  balanced by one +1 monoatomic cation (`_salt.py`'s
  `_MONOATOMIC_CATION_NAMES`) cites the cation's name before the R-group
  word, e.g. "sodium methyl sulfate" (PubChem CID 2735086) -- mirrors
  `_phosphate.py`'s own salt-of-partial-ester citation.

Scope: a single sulfur atom shaped like a sulfate ester -- two S=O double
bonds, and two more S-O positions, each either an ester (S-O-R), a plain
S-OH partial-ester remainder, or a single deprotonated S-O^- salt
position, each R named via `name_branch` (a plain alkyl chain, a
branched chain, or a plain benzene ring, and their halogenated variants,
exactly like `_phosphate.py`'s R), identical or mixed freely. Explicitly
out of scope (raise `UnsupportedStructure`): a 2+/3+ or ammonium cation,
sulfite (one fewer double-bonded oxygen -- a separate module), any
chalcogen-replacement analogue, and any other heteroatom.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import format_ester_words
from ._salt import _MONOATOMIC_CATION_NAMES
from ._substituents import name_branch

_SULFUR = 16


def _sulfate_sulfur_atoms(mol):
    """Sulfur atoms shaped like a fully- or partially-esterified sulfate:
    bonded to exactly two double-bonded (terminal) oxygens and exactly
    two more single-bonded oxygens, each either an ester oxygen (degree
    2 -- S plus one R carbon), a plain uncharged hydroxyl (degree 1, a
    neutral partial ester's remainder), or -- P-67.1.3.2, a salt of a
    partial ester -- a single degree-1 oxygen with formal charge -1 (its
    counter-cation is validated later, by `name_sulfate`, since that
    needs the whole multi-fragment molecule) -- at least one of the two
    must be an ester oxygen (both being hydroxyl/charged would be
    sulfuric acid itself or its bare anion, not an ester, out of scope
    here)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SULFUR or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygens) != 4:
            continue
        double_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        single_os = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 2 or len(single_os) != 2:
            continue
        if any(o.GetDegree() != 1 or o.GetFormalCharge() != 0 for o in double_os):
            continue
        charged_os = [o for o in single_os if o.GetFormalCharge() == -1]
        neutral_os = [o for o in single_os if o.GetFormalCharge() == 0]
        if len(charged_os) + len(neutral_os) != 2 or len(charged_os) > 1:
            continue
        if any(o.GetDegree() != 1 for o in charged_os):
            continue
        ester_os = [o for o in neutral_os if o.GetDegree() == 2]
        hydroxyl_os = [o for o in neutral_os if o.GetDegree() == 1]
        if len(ester_os) + len(hydroxyl_os) + len(charged_os) != 2 or not ester_os:
            continue
        matches.append(atom)
    return matches


def has_sulfate_shape(mol) -> bool:
    return bool(_sulfate_sulfur_atoms(mol))


def name_sulfate(mol) -> str:
    sulfur_atoms = _sulfate_sulfur_atoms(mol)
    if len(sulfur_atoms) != 1:
        raise UnsupportedStructure("more than one sulfate group is not supported yet")
    (sulfur,) = sulfur_atoms

    group_oxygens = {n.GetIdx() for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 8}
    charged_oxygens = {idx for idx in group_oxygens if mol.GetAtomWithIdx(idx).GetFormalCharge() == -1}

    frags = Chem.GetMolFrags(mol)
    anion_frag = next(frag for frag in frags if sulfur.GetIdx() in frag)
    other_frags = [frag for frag in frags if frag is not anion_frag]

    cation_name = None
    if charged_oxygens:
        # P-67.1.3.2: a salt of a partial ester cites the cation's name
        # before the R-group word -- scope limited to a single +1
        # monoatomic cation balancing the one deprotonated position above.
        if len(other_frags) != 1 or len(other_frags[0]) != 1:
            raise UnsupportedStructure("a salt with other than one monoatomic counter-ion is not supported yet")
        (cation_idx,) = other_frags[0]
        cation_atom = mol.GetAtomWithIdx(cation_idx)
        cation_name = _MONOATOMIC_CATION_NAMES.get((cation_atom.GetSymbol(), cation_atom.GetFormalCharge()))
        if cation_name is None or cation_atom.GetFormalCharge() != 1:
            raise UnsupportedStructure("only a single +1 monoatomic cation is supported for a partial-ester salt yet")
    elif other_frags:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    for idx in anion_frag:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetIsotope() != 0 or (atom.GetFormalCharge() != 0 and idx not in charged_oxygens):
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        atomic_num = atom.GetAtomicNum()
        if atomic_num == _SULFUR and idx != sulfur.GetIdx():
            raise UnsupportedStructure("more than one sulfur atom is not supported yet")
        if atomic_num not in (1, 6, 8, _SULFUR, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the sulfate's own sulfur/"
                "oxygens and a halogen substituent are not supported yet"
            )
        if atomic_num == 8 and idx not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the sulfate's own "
                "S(=O)(=O)(OR)2 group is out of scope for this module"
            )

    double_os = {
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 8
        and mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    }
    single_oxygens = [idx for idx in group_oxygens if idx not in double_os]
    ester_oxygens = [idx for idx in single_oxygens if mol.GetAtomWithIdx(idx).GetDegree() == 2]
    has_hydroxyl = len(single_oxygens) - len(ester_oxygens) - len(charged_oxygens) > 0

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != sulfur.GetIdx()]
        roots.append((root, oxygen_idx))

    names = [name_branch(graph, root, coming_from, halogens, aromatic_atoms, mol=mol)[0] for root, coming_from in roots]
    ester_words = format_ester_words(names)
    # P-67.1.3.2: a partial ester of a dibasic acid inserts the word
    # "hydrogen" (never "dihydrogen" -- sulfate has only two acidic
    # S-OH positions, so a partial ester always has exactly one R group
    # and one remaining hydroxyl) between the R-group word and "sulfate"
    # -- mirrors `_phosphate.py`'s own partial-ester citation, confirmed
    # worked example "CH3-O-SO2-OH -> methyl hydrogen sulfate (PIN)"
    # (the Blue Book).
    anion_name = f"{ester_words} hydrogen sulfate" if has_hydroxyl else ester_words + " sulfate"
    return f"{cation_name} {anion_name}" if cation_name else anion_name
