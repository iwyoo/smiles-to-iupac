"""Linear peptides named by the acyl groups of their N-terminal residues and the C-terminal amino acid (P-103.3.2, P-103.3.4).

The molecule is cut at every peptide bond; each piece is named as the free amino acid and the names of all but the
C-terminal one change their ending to 'yl' (P-103.2.5). An isopeptide bond to the side-chain carboxyl of aspartic or
glutamic acid makes a beta-aspartyl or gamma-glutamyl residue (P-103.2.5), a substituent of a residue amino group stays
in front of the residue name, and a C-terminal amide or ester ends the name as 'alaninamide' or 'methyl ... alaninate'
(P-103.2.6, P-103.2.7).
"""

import re
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure

_AMINO = "[NX3;!$(N=*)]-[CX4]-[CX3]=O"
_BONDS = {
    "alpha": Chem.MolFromSmarts(f"[CX3](=O)(-[CX4]-[NX3])-{_AMINO}"),
    "beta": Chem.MolFromSmarts(f"[CX3](=O)(-[CH2]-[CX4](-[NX3])-[CX3]=O)-{_AMINO}"),
    "gamma": Chem.MolFromSmarts(f"[CX3](=O)(-[CH2]-[CH2]-[CX4](-[NX3])-[CX3]=O)-{_AMINO}"),
}
_SPECIAL_ACYL = {
    "cysteine": "cysteinyl",
    "asparagine": "asparaginyl",
    "glutamine": "glutaminyl",
    "tryptophan": "tryptophyl",
    "aspartic acid": "aspartyl",
    "glutamic acid": "glutamyl",
}
_SIDE_CHAIN = {("aspartic acid", "beta"): "β-aspartyl", ("glutamic acid", "gamma"): "γ-glutamyl"}
_RESIDUES = {
    "alanine", "arginine", "histidine", "isoleucine", "leucine", "lysine", "methionine", "phenylalanine",
    "proline", "serine", "threonine", "valine", "tyrosine", "glycine", *_SPECIAL_ACYL,
}
_BASE = "|".join(sorted(_RESIDUES, key=len, reverse=True))
_DESCRIPTOR = re.compile(r"^(?P<pre>.*?)(?:(?<![A-Za-z])(?P<descriptor>[LD])-)?(?P<allo>allo)?$")


def _acyl(base):
    return _SPECIAL_ACYL.get(base) or base[: -len("ine")] + "yl"


def _last_forms(base):
    """The names a C-terminal residue can take: the amino acid and its ester anion (P-103.2.6)."""
    return {base: base, base[:-1] + "ate": base[:-1] + "ate"} if base.endswith("ine") else {base: base}


def acyl_word(residue, last=False, xi=True, side_chain=None):
    """The acyl (or, for the C-terminal residue, the amino acid) word of a residue name such as 'L-alanine' or
    'N-methyl-L-alanine'; None if not a common amino acid."""
    words = _words(residue, last)
    if words is None:
        return None
    pre, descriptor, allo, base, ending = words
    if descriptor is None and base != "glycine" and xi:
        descriptor = "ξ"
    if last:
        word = ending
    elif side_chain and (base, side_chain) in _SIDE_CHAIN:
        word = _SIDE_CHAIN[(base, side_chain)]
    else:
        word = _acyl(base)
    return pre + (f"{descriptor}-" if descriptor else "") + (allo or "") + word


def _words(residue, last):
    for base in sorted(_RESIDUES, key=len, reverse=True):
        endings = _last_forms(base).values() if last else (base,)
        for ending in endings:
            if not residue.endswith(ending):
                continue
            head = _DESCRIPTOR.match(residue[: -len(ending)])
            if head is not None and (not head.group("pre") or head.group("pre").endswith("-")):
                return head.group("pre"), head.group("descriptor"), head.group("allo"), base, ending
    return None


def _cut(mol):
    cuts = {}
    for kind, query in _BONDS.items():
        for match in mol.GetSubstructMatches(query):
            carbon = match[0]
            nitrogen = next(i for i in match[1:] if mol.GetAtomWithIdx(i).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(carbon, i) is not None)
            cuts.setdefault((carbon, nitrogen), kind)
    editable = Chem.RWMol(mol)
    for carbon, nitrogen in cuts:
        editable.RemoveBond(carbon, nitrogen)
        editable.AddBond(carbon, editable.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    pieces = editable.GetMol()
    Chem.SanitizeMol(pieces)
    return cuts, pieces


def _sequence(cuts, fragments):
    """Fragment indices from the N-terminal residue to the C-terminal one, None unless they form one chain."""
    owner = {a: i for i, atoms in enumerate(fragments) for a in atoms}
    following = {owner[c]: owner[n] for c, n in cuts}
    if len(following) != len(cuts) or len(set(following.values())) != len(cuts):
        return None
    first = [i for i in range(len(fragments)) if i not in following.values()]
    if len(first) != 1:
        return None
    order = [first[0]]
    while order[-1] in following:
        order.append(following[order[-1]])
    return order if len(order) == len(fragments) else None


@lru_cache(maxsize=256)
def _name(smiles):
    from .core import smiles_to_iupac

    mol = Chem.MolFromSmiles(smiles)
    try:
        cuts, pieces = _cut(mol)
    except Chem.MolSanitizeException:
        return None
    fragments = Chem.GetMolFrags(pieces)
    if not cuts or len(fragments) != len(cuts) + 1:
        return None
    order = _sequence(cuts, fragments)
    if order is None:
        return None
    owner = {a: i for i, atoms in enumerate(fragments) for a in atoms}
    kind_of = {owner[c]: kind for (c, _), kind in cuts.items()}
    molecules = Chem.GetMolFrags(pieces, asMols=True)
    words, ester = [], None
    for position, index in enumerate(order):
        last = position == len(order) - 1
        try:
            residue = smiles_to_iupac(Chem.MolToSmiles(molecules[index]))
        except (UnsupportedStructure, ValueError):
            return None
        if last:
            alkyl = re.match(r"^(?P<alkyl>[^ ]+) (?P<rest>.+ate)$", residue)
            if alkyl is not None:
                ester, residue = alkyl.group("alkyl"), alkyl.group("rest")
        word = acyl_word(residue, last, side_chain=kind_of.get(index))
        if word is None:
            return None
        words.append(word)
    text = words[0] + "".join(("-" if w[1:2] == "-" else "") + w for w in words[1:])
    return f"{ester} {text}" if ester else text


def has_peptide_shape(mol) -> bool:
    return any(mol.HasSubstructMatch(query) for query in _BONDS.values()) and _name(Chem.MolToSmiles(mol)) is not None


def name_peptide(mol) -> str:
    return _name(Chem.MolToSmiles(mol))
