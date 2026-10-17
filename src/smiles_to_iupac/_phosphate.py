"""Naming of trialkyl phosphate esters (P(=O)(OR)3, R groups identical or
mixed), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group(s) as separate word(s), in alphanumeric
  order if more than one, followed by the acid's anion name. Confirmed
  worked example: `P(O-CH3)3` -> phosphite... no wait, the P=O analogue
  `P(O)(O-CH3)3` -> 'trimethyl phosphate (PIN)' (the Blue Book). For three identical R groups this is a single
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
- Partial esters (P-67.1.3.2, the Blue Book): 1 or 2
  of the three acidic P-OH positions may be left unesterified, cited by
  inserting the word "hydrogen" (2 R groups, 1 remaining OH) or
  "dihydrogen" (1 R group, 2 remaining OH) between the R-group word(s)
  and "phosphate" -- confirmed worked example `P(O)(O-CH3)(OH)2` ->
  "methyl dihydrogen phosphate (PIN)", and by direct analogy (same rule,
  one more R group) the 2-R "hydrogen" case, both confirmed real via
  PubChem (CIDs 13130/74190 for the 1-R case, 13134/654 for the 2-R
  case).

Scope: a single phosphorus atom shaped like a phosphate ester -- one P=O
double bond, and three more P-O positions, each either an ester (P-O-R),
a plain P-OH partial-ester remainder, or (P-67.1.3.2, salts of partial
acid esters) a single deprotonated P-O^- balanced by one +1 monoatomic
cation (`_salt.py`'s `_MONOATOMIC_CATION_NAMES`), each R named via
`name_branch` with the groups that rank below the ester as prefixes (P-41),
identical or mixed freely. The cation's name is cited before the
R-group word(s), e.g. "sodium methyl hydrogen phosphate" (PubChem CID
23690691). Explicitly out of scope (raise `UnsupportedStructure`): a 2+/3+
or ammonium cation, more than one remaining deprotonated/acidic position,
any chalcogen-replacement analogue, and an acid, anhydride or other ester in R.

The word-assembly step (group identical R names, multiplying-prefix +
enclosure per group, alphanumeric word order) is shared with the
phosphite ester case (no P=O, P(III) instead of P(V)) via
`format_ester_words`, exported for `_phosphite.py` to reuse.
"""

import contextlib
import re

from rdkit import Chem

from collections import Counter

from ._common import assign_cip_labels, match_cip_labels
from ._multiplicative_text import enclose
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents, is_nitro_nitrogen
from ._numerals import multiplying_prefix
from ._salt import _MONOATOMIC_CATION_NAMES
from ._substituents import BRANCH_STEREO, name_branch, prefix_multiplier

_PHOSPHORUS = 15


def _phosphate_phosphorus_atoms(mol):
    """Phosphorus atoms shaped like a fully- or partially-esterified
    phosphate: bonded to exactly one double-bonded (terminal) oxygen and
    exactly three more single-bonded oxygens, each either an ester
    oxygen (degree 2 -- P plus one R carbon), a plain uncharged hydroxyl
    (degree 1, a neutral partial ester's remainder), or -- P-67.1.3.2, a
    salt of a partial ester -- a single degree-1 oxygen with formal
    charge -1 (its counter-cation is validated later, by `name_phosphate`,
    since that needs the whole multi-fragment molecule, not just this
    P-shaped group) -- at least one of the three must be an ester oxygen
    (all three being hydroxyl/charged would be phosphoric acid itself or
    its bare anion, not an ester, out of scope here). A phosphorus in a ring is a ring heteroatom of a heterocycle
    (P-65.6.3.5 pseudoketone, P-22.2), not the centre of an acyclic phosphate ester."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _PHOSPHORUS or atom.GetDegree() != 4 or atom.GetFormalCharge() != 0:
            continue
        if atom.IsInRing():
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
        if len(double_os) != 1 or len(single_os) != 3:
            continue
        (double_o,) = double_os
        if double_o.GetDegree() != 1 or double_o.GetFormalCharge() != 0:
            continue
        charged_os = [o for o in single_os if o.GetFormalCharge() == -1]
        neutral_os = [o for o in single_os if o.GetFormalCharge() == 0]
        if len(charged_os) + len(neutral_os) != 3 or len(charged_os) > 1:
            continue
        if any(o.GetDegree() != 1 for o in charged_os):
            continue
        ester_os = [o for o in neutral_os if o.GetDegree() == 2]
        hydroxyl_os = [o for o in neutral_os if o.GetDegree() == 1]
        if len(ester_os) + len(hydroxyl_os) + len(charged_os) != 3 or not ester_os:
            continue
        matches.append(atom)
    return matches


def _is_hydroxy_on_carbon(oxygen) -> bool:
    return oxygen.GetDegree() == 1 and oxygen.GetNeighbors()[0].GetAtomicNum() == 6 and oxygen.GetTotalNumHs() == 1


_ACETATE_ENDING = re.compile(r"(di|tri|tetra|penta|hexa|hepta|octa)acetate$")


def name_polyester_through_acetates(mol, centers, word, prefix):
    """P-65.6.3.2.2: 'propane-1,2,3-triyl trinitrate'. The organyl group joining identical nitrate or nitrite groups is
    named as in the acetate of the same polyol, whose acetyloxy groups stand in for the nitrogen ester groups; `word`
    replaces 'acetate' and `prefix` the 'acetyloxy' groups cited on the organyl group."""
    from .core import smiles_to_iupac

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    proxy = Chem.RWMol(mol)
    for center in centers:
        terminals = sorted(
            (n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors() if n.GetDegree() == 1),
            key=lambda idx: mol.GetBondBetweenAtoms(center, idx).GetBondTypeAsDouble() != 2.0,
        )
        proxy.GetAtomWithIdx(center).SetAtomicNum(6)
        proxy.GetAtomWithIdx(center).SetFormalCharge(0)
        proxy.GetAtomWithIdx(terminals[0]).SetFormalCharge(0)
        proxy.GetBondBetweenAtoms(center, terminals[0]).SetBondType(Chem.BondType.DOUBLE)
        if len(terminals) == 2:
            proxy.GetAtomWithIdx(terminals[1]).SetAtomicNum(6)
            proxy.GetAtomWithIdx(terminals[1]).SetFormalCharge(0)
            proxy.GetAtomWithIdx(terminals[1]).SetNoImplicit(False)
            proxy.GetBondBetweenAtoms(center, terminals[1]).SetBondType(Chem.BondType.SINGLE)
        else:
            proxy.AddBond(center, proxy.AddAtom(Chem.Atom(6)), Chem.BondType.SINGLE)
    proxy = proxy.GetMol()
    Chem.SanitizeMol(proxy)
    name = smiles_to_iupac(Chem.MolToSmiles(proxy))
    matched = Chem.Mol(proxy)
    match_cip_labels(mol, matched)
    if Chem.MolToSmiles(matched) != Chem.MolToSmiles(proxy):
        # a changed CIP order rewrites R/S descriptors, but a configurational name (a sugar) follows the geometry
        cited = smiles_to_iupac(Chem.MolToSmiles(matched))
        if re.search(r"\b\d+[a-z]?[RSEZ]\b", cited):
            name = cited
    ending = _ACETATE_ENDING.search(name)
    if ending is None:
        raise UnsupportedStructure("the acetate of this polyol is not a multiplied ester name")
    name = (name[: ending.start()] + ending.group(1) + word).replace("acetyloxy", prefix)
    if "acetate" in name or "acetyloxy" in name:
        raise UnsupportedStructure("this polyol keeps an acetate group in its name")
    return name


_CHALCOGENS = (8, 16, 34, 52)


def _outranks_ester(mol, atom) -> bool:
    """An acid, an anhydride or an ester (P-41): a boron atom with a hydroxy group, or an atom doubly bonded to a chalcogen
    that carries a hydroxy-type group or a chalcogen bonded to a second atom outside the ring of the atom."""
    if atom.GetAtomicNum() == 5:
        return any(n.GetAtomicNum() in _CHALCOGENS and n.GetDegree() == 1 and n.GetTotalNumHs() == 1 for n in atom.GetNeighbors())
    if not any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in _CHALCOGENS for b in atom.GetBonds()):
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() not in _CHALCOGENS or bond.GetBondTypeAsDouble() != 1.0 or bond.IsInRing():
            continue
        if other.GetDegree() == 2 or other.GetTotalNumHs() or other.GetFormalCharge():
            return True
    return False


def check_ester_organyl(mol, group_atoms, anionic=False, charged=frozenset(), atoms=None):
    """Raise unless every atom of `atoms` (all by default) outside the ester group `group_atoms` belongs to an organyl group
    whose substituents rank below the ester (P-41): an acid, another ester or an anhydride outranks it, except beside an
    anionic ester. The atoms of `charged` may carry a charge."""
    nitro = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in group_atoms and is_nitro_nitrogen(mol, atom.GetIdx()):
            nitro |= {atom.GetIdx(), *(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)}
    for atom in mol.GetAtoms() if atoms is None else (mol.GetAtomWithIdx(i) for i in atoms):
        if atom.GetIdx() in group_atoms or atom.GetIdx() in nitro:
            continue
        if atom.GetIsotope() or (atom.GetFormalCharge() and atom.GetIdx() not in charged):
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if not anionic and _outranks_ester(mol, atom):
            raise UnsupportedStructure("an acid or another ester outranks the ester of the noncarbon acid")


def has_phosphate_shape(mol) -> bool:
    return bool(_phosphate_phosphorus_atoms(mol))


def name_phosphate(mol) -> str:
    phosphorus_atoms = _phosphate_phosphorus_atoms(mol)
    if len(phosphorus_atoms) != 1:
        raise UnsupportedStructure("more than one phosphate group is not supported yet")
    (phosphorus,) = phosphorus_atoms

    group_oxygens = {n.GetIdx() for n in phosphorus.GetNeighbors() if n.GetAtomicNum() == 8}
    charged_oxygens = {idx for idx in group_oxygens if mol.GetAtomWithIdx(idx).GetFormalCharge() == -1}

    frags = Chem.GetMolFrags(mol)
    anion_frag = next(frag for frag in frags if phosphorus.GetIdx() in frag)
    other_frags = [frag for frag in frags if frag is not anion_frag]

    cation_name = None
    positive = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() > 0]
    zwitterion = bool(charged_oxygens) and not other_frags and sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0
    if zwitterion and not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in positive):
        raise UnsupportedStructure("a zwitterion with other than ammonium centres is not supported yet")
    anion_only = bool(charged_oxygens) and not other_frags and not positive
    if charged_oxygens and not zwitterion and not anion_only:
        # P-67.1.3.2: a salt of a partial ester cites the cation's name
        # before the R-group word(s) -- scope limited to a single +1
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

    check_ester_organyl(
        mol, group_oxygens | {phosphorus.GetIdx()}, bool(charged_oxygens), set(positive) if zwitterion else frozenset(), anion_frag
    )

    double_o = next(
        n
        for n in phosphorus.GetNeighbors()
        if n.GetAtomicNum() == 8
        and mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    )
    single_oxygens = [idx for idx in group_oxygens if idx != double_o.GetIdx()]
    ester_oxygens = [idx for idx in single_oxygens if mol.GetAtomWithIdx(idx).GetDegree() == 2]
    hydroxyl_count = len(single_oxygens) - len(ester_oxygens) - len(charged_oxygens)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    roots = []
    for oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[oxygen_idx] if n != phosphorus.GetIdx()]
        roots.append((root, oxygen_idx))

    entries = cited_ester_groups(mol, graph, roots, halogens, aromatic_atoms)
    ester_words = format_ester_entries(entries)
    # P-67.1.3.2: a partial ester of a polybasic acid inserts the word
    # "hydrogen" (with a multiplying prefix if more than one remaining
    # acidic P-OH) between the R-group word(s) and the anion name --
    # confirmed worked examples "methyl dihydrogen phosphate (PIN)" and
    # "dimethyl hydrogen phosphate (PIN)" (the Blue Book).
    if hydroxyl_count:
        hydrogen_word = multiplying_prefix(hydroxyl_count) + "hydrogen" if hydroxyl_count > 1 else "hydrogen"
        anion_name = f"{ester_words} {hydrogen_word} phosphate"
    else:
        anion_name = ester_words + " phosphate"
    return f"{cation_name} {anion_name}" if cation_name else anion_name


@contextlib.contextmanager
def citing_stereo(mol):
    """Within the block `name_branch` cites the CIP descriptors of the stereo elements of `mol` (P-93.5); leaving the block
    raises when one of them was not cited by any name."""
    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    context = {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.HasProp("_CIPCode")},
        "bonds": {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode") for b in probe.GetBonds() if b.HasProp("_CIPCode")},
        "used": set(),
    }
    token = BRANCH_STEREO.set(context)
    try:
        yield
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of the ester groups is not cited by any supported name")


def cited_ester_groups(mol, graph, roots, halogens, aromatic_atoms):
    """[(name, is compound)] of the organyl groups of an ester, `roots` being (root, ester oxygen) pairs, each with the CIP
    descriptors of its stereo elements cited in front of the group name (P-93.5)."""
    if any(mol.GetAtomWithIdx(root).GetAtomicNum() in (7, 8) for root, _ in roots):
        raise UnsupportedStructure("an acid group on nitrogen or oxygen is a hydroxylamine or peroxy acid derivative, not an ester")
    with citing_stereo(mol):
        return [name_branch(graph, root, via, halogens, aromatic_atoms, mol=mol) for root, via in roots]


def format_ester_entries(entries) -> str:
    """`format_ester_words` for (name, is compound) pairs as `name_branch` returns them: a substituted group is cited
    with 'bis', 'tris' (P-16.3.5)."""
    return format_ester_words([name for name, _ in entries], frozenset(name for name, compound in entries if compound))


def format_ester_words(names, compound_names=frozenset()) -> str:
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
        if count > 1 and name in compound_names:
            multiplier, enclosed = prefix_multiplier(count, name, True)
            words.append(multiplier + (enclose(name) if enclosed else name))
            continue
        needs_enclosure = (name[0].isdigit() or name[0] == "(" or name in compound_names) and count > 1
        group = enclose(name) if needs_enclosure else name
        prefix = multiplying_prefix(count, compound=needs_enclosure) if count > 1 else ""
        words.append(f"{prefix}{group}")
    return " ".join(words)
