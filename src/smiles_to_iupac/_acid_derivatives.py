"""Functional class names of esters, pseudoesters, anhydrides and acyl halides (P-65.5-P-65.7): the molecule is cut at
the junction bonds, each acid piece is named as the free acid by the substitutive engine and turned into an anion,
acid-word or acyl name, and organyl pieces are cited as substituent groups (which also expresses polyester prefixes).
"""

import itertools
import re
from dataclasses import dataclass

from rdkit import Chem

from ._common import assign_cip_labels, has_acyloxy_amine_group, match_cip_labels
from ._acid_lexicon import CHALCOGENS
from ._amino_acid import SYSTEMATIC_ACID_PROBE
from ._amino_acid_retained import has_retained_amino_acid_shape
from ._carbonic_family import pseudohalide_at
from ._chalcogenourea import has_group_senior_to_urea, is_urea_substituent_root
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, citation_order_key, halogen_substituents, is_halide_ligand
from ._locant_omission import CITE_SUBSTITUENT_LOCANTS
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._isotope_labels import restore_labels
from ._hetero_prefixes import EXTENDED_PREFIXES
from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS, name_branch, prefix_multiplier

_SYMBOL = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_CENTERS = {6: "C", 16: "S", 34: "Se", 52: "Te"}
_INORGANIC_CENTERS = {5, 15, 33, 51}
_CHALCOGEN_CENTERS = {16, 34, 52}
_HEAVY_PNICTOGEN_STEMS = {33: "ars", 51: "stib"}
_FREE_INORGANIC_ACIDS = {"OB(O)O": "boric acid", "OP(O)O": "phosphorous acid",
    "O=P(O)(O)O": "phosphoric acid",
    "O=[PH](O)O": "phosphonic acid",
    "O=[PH2]O": "phosphinic acid",
    "OPO": "phosphonous acid",
    "PO": "phosphinous acid",
}
_HALIDES = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_PSEUDOESTER_ELEMENTS = {5, 7, 13, 14, 15, 31, 32, 33, 49, 50, 51, 81, 82, 83}
_RING_ONLY_PSEUDOESTER = {5, 15, 33, 51, 83, 7}
_ANHYDRIDE_WORDS = {
    ("O",): "anhydride",
    ("S",): "thioanhydride",
    ("Se",): "selenoanhydride",
    ("Te",): "telluroanhydride",
    ("O", "O"): "peroxyanhydride",
    ("O", "S"): "thioperoxyanhydride",
    ("S", "O"): "thioperoxyanhydride",
    ("S", "S"): "dithioperoxyanhydride",
}


@dataclass
class Link:
    kind: str
    center: int
    chain: tuple
    far: int


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _oxo_slots(mol, center):
    atom = mol.GetAtomWithIdx(center)
    slots = []
    for n in atom.GetNeighbors():
        if _bond(mol, center, n.GetIdx()) != 2.0 or n.GetFormalCharge():
            continue
        if n.GetAtomicNum() in _SYMBOL and n.GetDegree() == 1:
            slots.append(_SYMBOL[n.GetAtomicNum()])
        elif n.GetAtomicNum() == 7 and n.GetTotalNumHs() == 1 and n.GetDegree() == 1:
            slots.append("NH")
        elif n.GetAtomicNum() == 7 and not n.GetTotalNumHs():
            far = [m for m in n.GetNeighbors() if m.GetIdx() != center]
            if len(far) == 1 and far[0].GetAtomicNum() == 7 and far[0].GetDegree() == 1 and far[0].GetTotalNumHs() == 2:
                slots.append("NNH2")
            elif far and all(_bond(mol, n.GetIdx(), m.GetIdx()) == 1.0 for m in far):
                slots.append("NH")
    return slots


def _role(mol, center, atom):
    z = atom.GetAtomicNum()
    if z in _HALIDES and atom.GetTotalDegree() == 1 and is_halide_ligand(atom):
        return "hal"
    if pseudohalide_at(mol, atom.GetIdx(), center) is not None:
        return "pseudo"
    if z in _SYMBOL:
        return "Y"
    return {6: "C", 7: "N"}.get(z, "other")


def _inorganic_center(mol, idx):
    """'inorganic' kind for a boron or phosphorus acid centre (borates, phosphates, phosphonates, phosphinates):
    only P=O, hydroxy/oxy and organyl neighbours, at least one oxy position."""
    atom = mol.GetAtomWithIdx(idx)
    oxo = [n for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 2.0]
    chalcogen = atom.GetAtomicNum() in _CHALCOGEN_CENTERS
    if len(oxo) > (2 if chalcogen else 1) or any(n.GetAtomicNum() != 8 or n.GetDegree() != 1 for n in oxo):
        return None
    rest = [n for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 1.0]
    direct = [n for n in rest if n.GetAtomicNum() == atom.GetAtomicNum()]
    rest = [n for n in rest if n not in direct]
    if len(oxo) + len(rest) + len(direct) != atom.GetDegree() or any(n.GetAtomicNum() not in (6, 8, 16, 34, 52) for n in rest):
        return None
    positions = [(n.GetIdx(), "Y") for n in rest if n.GetAtomicNum() in (8, 16, 34, 52)]
    if not positions or atom.GetAtomicNum() == 5 and oxo:
        return None
    return "inorganic", ["O"] * len(oxo), positions


def _acylated_oxoacid_center(mol, atom):
    """A sulfur, selenium or tellurium oxoacid centre (no organyl neighbour) with an acyloxy, nitrooxy or sulfonyloxy neighbour."""
    if atom.GetFormalCharge() or atom.GetIsotope() or atom.IsInRing():
        return False
    neighbors = atom.GetNeighbors()
    if any(n.GetAtomicNum() == 6 for n in neighbors):
        return False
    return any(
        n.GetAtomicNum() == 8
        and any(
            c.GetIdx() != atom.GetIdx()
            and (
                (c.GetAtomicNum() == 6 and _oxo_slots(mol, c.GetIdx()))
                or _nitrogen_oxoacid_centre(mol, c, n)
                or (c.GetAtomicNum() in _CHALCOGEN_CENTERS and _oxo_slots(mol, c.GetIdx()) and any(m.GetAtomicNum() == 6 for m in c.GetNeighbors()))
            )
            for c in n.GetNeighbors()
        )
        for n in neighbors
    )


def _halogen_oxoacid_center(mol, idx):
    """'inorganic' kind for the halogen of a hypohalous, halous, halic or perhalic acid joined through one oxygen
    (P-67.1.3.3): X-O-R with up to three terminal oxygens on X."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() not in _HALIDES or atom.GetIsotope() or atom.IsInRing():
        return None
    bridges, oxo = [], []
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != 8:
            return None
        double = _bond(mol, idx, n.GetIdx()) == 2.0 and not n.GetFormalCharge()
        anionic = _bond(mol, idx, n.GetIdx()) == 1.0 and n.GetFormalCharge() == -1 and n.GetDegree() == 1
        if n.GetDegree() == 1 and (double or anionic):
            oxo.append(n)
        elif n.GetDegree() == 2 and _bond(mol, idx, n.GetIdx()) == 1.0 and not n.GetFormalCharge():
            bridges.append(n)
        else:
            return None
    if len(bridges) != 1 or len(oxo) > 3 or atom.GetFormalCharge() != sum(-o.GetFormalCharge() for o in oxo):
        return None
    if atom.GetTotalValence() != 1 + 2 * len(oxo) + atom.GetFormalCharge():
        return None
    return "inorganic", ["O"] * len(oxo), [(bridges[0].GetIdx(), "Y")]


def _ring_nitrogen_parent(atom):
    """A neutral ring nitrogen that is the parent-hydride position carrying a sulfonyl-type group (P-65.3.1)."""
    return atom.GetAtomicNum() == 7 and atom.IsInRing() and not atom.GetFormalCharge()


def _center(mol, idx):
    """(kind, oxo slots, [(neighbor, role)] of the Y positions) for an acyl-type ('acyl'), carbonic
    ('carbonic') or cyanic ('cyanic') centre, else None."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() in _HALIDES:
        return _halogen_oxoacid_center(mol, idx)
    if atom.GetAtomicNum() in _INORGANIC_CENTERS and not (atom.GetFormalCharge() or atom.GetIsotope() or atom.IsInRing()):
        return _inorganic_center(mol, idx)
    if atom.GetAtomicNum() in _CHALCOGEN_CENTERS and _acylated_oxoacid_center(mol, atom):
        return _inorganic_center(mol, idx)
    symbol = _CENTERS.get(atom.GetAtomicNum())
    if symbol is None or atom.GetFormalCharge() or atom.GetIsotope() or atom.IsInRing():
        return None
    if symbol == "C":
        triple = [n for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 3.0 and n.GetAtomicNum() == 7 and n.GetDegree() == 1]
        if triple:
            rest = [n for n in atom.GetNeighbors() if n.GetIdx() != triple[0].GetIdx()]
            if len(rest) == 1 and _role(mol, idx, rest[0]) in ("Y", "hal", "pseudo"):
                return "cyanic", [], [(rest[0].GetIdx(), _role(mol, idx, rest[0]))]
            return None
    oxo = _oxo_slots(mol, idx)
    if symbol == "C" and len(oxo) != 1 or symbol != "C" and len(oxo) not in (1, 2):
        return None
    rest = [n for n in atom.GetNeighbors() if _bond(mol, idx, n.GetIdx()) == 1.0]
    if atom.GetDegree() != len(oxo) + len(rest):
        return None
    roles = [(n.GetIdx(), _role(mol, idx, n)) for n in rest]
    organics = [r for r in roles if r[1] == "C" or symbol != "C" and r[1] == "N" and _ring_nitrogen_parent(mol.GetAtomWithIdx(r[0]))]
    positions = [r for r in roles if r not in organics]
    if symbol != "C":
        kind = "acyl" if len(organics) == 1 and len(positions) == 1 else None
    elif len(organics) == 1 and len(positions) == 1:
        kind = "acyl"
    elif not organics and atom.GetTotalNumHs() == 1 and len(positions) == 1:
        kind = "acyl"
    elif not organics and not atom.GetTotalNumHs() and len(positions) == 2:
        kind = "carbonic"
    else:
        kind = None
    if kind is None or any(role == "other" for _, role in positions):
        return None
    if kind == "acyl" and positions[0][1] == "N":
        return None
    return kind, oxo, positions


def _inorganic_anhydride(near, far):
    """Whether a link between two noncarbon acid centres is an anhydride (P-67.1.3.3): one is an acid of one position
    with an oxo group, as in phosphinic acid. Between two -ous centres the compound is an alternating heteroatom chain
    and between two polybasic acids a polynuclear acid (P-68.1.1.2.2, P-67.1.2)."""
    return (len(near[2]) == 1 or len(far[2]) == 1) and bool(near[1] or far[1])


def _nitrogen_oxoacid_centre(mol, atom, bridge):
    """The nitrogen of a nitrate (N(+)(=O)O(-)) or nitrite (N=O) group joined through the oxygen `bridge`: the acid
    component of a mixed anhydride of nitric or nitrous acid (P-65.7.2)."""
    if atom.GetAtomicNum() != 7 or atom.GetIsotope() or atom.IsInRing() or atom.GetDegree() != 3 - (atom.GetFormalCharge() == 0):
        return False
    oxygens = [n for n in atom.GetNeighbors() if n.GetIdx() != bridge.GetIdx()]
    if atom.GetFormalCharge() == 0:
        return len(oxygens) == 1 and oxygens[0].GetAtomicNum() == 8 and _bond(mol, atom.GetIdx(), oxygens[0].GetIdx()) == 2.0
    kinds = sorted((_bond(mol, atom.GetIdx(), o.GetIdx()), o.GetFormalCharge()) for o in oxygens if o.GetAtomicNum() == 8 and o.GetDegree() == 1)
    return atom.GetFormalCharge() == 1 and kinds == [(1.0, -1), (2.0, 0)]


def find_links(mol):
    """Every junction at an acyl, carbonic or cyanic centre: (links, free acid present)."""
    links, acid = [], []
    for atom in mol.GetAtoms():
        found = _center(mol, atom.GetIdx())
        if found is None:
            continue
        idx = atom.GetIdx()
        for position, role in found[2]:
            if role in ("hal", "pseudo"):
                links.append(Link("halide" if role == "hal" else "pseudohalide", idx, (), position))
                continue
            if role != "Y":
                continue
            y = mol.GetAtomWithIdx(position)
            if y.GetFormalCharge():
                if y.GetFormalCharge() == -1 and y.GetDegree() == 1:
                    acid.append(idx)
                continue
            chain = [position]
            beyond = [n for n in y.GetNeighbors() if n.GetIdx() != idx]
            if (
                len(beyond) == 1
                and beyond[0].GetAtomicNum() in _SYMBOL
                and not beyond[0].GetFormalCharge()
                and _center(mol, beyond[0].GetIdx()) is None
            ):
                chain.append(beyond[0].GetIdx())
                beyond = [n for n in beyond[0].GetNeighbors() if n.GetIdx() != y.GetIdx()]
            last = mol.GetAtomWithIdx(chain[-1])
            if not beyond:
                if last.GetTotalNumHs() == 1:
                    acid.append(idx)
                continue
            if len(beyond) != 1 or last.GetTotalNumHs():
                continue
            far = beyond[0]
            if mol.GetBondBetweenAtoms(chain[-1], far.GetIdx()).IsInRing():
                raise UnsupportedStructure("a cyclic ester or anhydride is named as a heterocycle")
            if found[0] in ("acyl", "inorganic") and _nitrogen_oxoacid_centre(mol, far, mol.GetAtomWithIdx(chain[-1])):
                links.append(Link("anhydride", idx, tuple(chain), far.GetIdx()))
            elif _center(mol, far.GetIdx()) is not None:
                if found[0] == "carbonic" and _center(mol, far.GetIdx())[0] == "carbonic" and len(chain) == 1:
                    continue
                if found[0] == "inorganic" and _center(mol, far.GetIdx())[0] == "inorganic" and len(chain) == 1:
                    if not _inorganic_anhydride(found, _center(mol, far.GetIdx())):
                        continue
                links.append(Link("anhydride", idx, tuple(chain), far.GetIdx()))
            elif far.GetAtomicNum() == 6 or (
                far.GetAtomicNum() in _PSEUDOESTER_ELEMENTS
                and (far.GetAtomicNum() not in _RING_ONLY_PSEUDOESTER or far.IsInRing())
            ):
                links.append(Link("ester", idx, tuple(chain), far.GetIdx()))
    return links, acid


def _restore_nuclides(editable, mol, keep, skip=()):
    """The nuclide labels of the kept atoms go back on an acid cut from the unlabelled `mol` (P-82.2.4): the acid is
    named, and its descriptor cited, as a word of its own."""
    context = ISOTOPE_LABELS.get()
    if not context or context.get("mol") is not mol:
        return
    kept = set(keep) - set(skip)
    restore_labels(editable, mol, context["labels"], skip=set(range(mol.GetNumAtoms())) - kept)
    context["consumed"].update(kept & set(context["labels"]))


def _capped_acid(mol, keep, caps, skip=()):
    """The free acid made of the atoms `keep` of `mol` with a hydroxy group on every atom of `caps`, and the index
    of each kept atom in it."""
    editable = Chem.RWMol(mol)
    for idx in caps:
        atom = editable.GetAtomWithIdx(idx)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
    _restore_nuclides(editable, mol, keep, skip)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    return acid, {idx: i for i, idx in enumerate(sorted(keep))}


def _acid_name(mol, keep, caps, skip=()):
    """Name of the free acid made of the atoms `keep` of `mol` with a hydroxy group on every atom of `caps`."""
    acid, index = _capped_acid(mol, keep, caps, skip)
    match_cip_labels(mol, acid, {i: idx for idx, i in index.items()})
    for atom in acid.GetAtoms():
        if atom.GetAtomicNum() in (15, 16):
            # the descriptor of a stereogenic acid centre is cited in front of the anion, not inside the acid name
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return _free_acid_name(acid)


def _free_acid_name(acid):
    from ._phosphanyl_group import PREFIX_PROP
    from .core import _name_mol, smiles_to_iupac

    if any(atom.HasProp(PREFIX_PROP) for atom in acid.GetAtoms()):
        return _name_mol(acid)
    smiles = Chem.MolToSmiles(acid)
    if smiles in _FREE_INORGANIC_ACIDS:
        return _FREE_INORGANIC_ACIDS[smiles]
    heavy = [
        a
        for a in acid.GetAtoms()
        if a.GetAtomicNum() in _HEAVY_PNICTOGEN_STEMS
        and any(
            n.GetAtomicNum() in (8, 16, 34, 52) and n.GetDegree() == 1 and acid.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in a.GetNeighbors()
        )
    ]
    if heavy:
        return _heavy_pnictogen_acid_name(acid, heavy)
    if has_retained_amino_acid_shape(acid):
        # the retained name of a Table 10.5 acid stands for the free acid only (P-62.2.3)
        token = SYSTEMATIC_ACID_PROBE.set(True)
        try:
            return smiles_to_iupac(smiles)
        finally:
            SYSTEMATIC_ACID_PROBE.reset(token)
    return _phosphorous_acid_name(acid) or smiles_to_iupac(smiles)


def _heavy_pnictogen_acid_name(acid, heavy):
    """'methylarsonic acid', 'stiboric acid' (P-67.1.1.1): the acid is named as its phosphorus analogue and the
    stem 'phosph' becomes 'ars' or 'stib'."""
    if len(heavy) != 1 or any(a.GetAtomicNum() == 15 for a in acid.GetAtoms()):
        raise UnsupportedStructure("several group 15 acid centres in one acid piece")
    stem = _HEAVY_PNICTOGEN_STEMS[heavy[0].GetAtomicNum()]
    analogue = Chem.RWMol(acid)
    analogue.GetAtomWithIdx(heavy[0].GetIdx()).SetAtomicNum(15)
    phosphorus = analogue.GetMol()
    Chem.SanitizeMol(phosphorus)
    name = _free_acid_name(phosphorus)
    renamed = re.sub(r"phosph(?=(?:on|in|or)(?:ic|ous) acid$)", stem, name)
    if renamed == name:
        raise UnsupportedStructure("the acid is not a simple oxoacid of a heavier group 15 element")
    return renamed


def _phosphorous_acid_name(acid):
    """'dimethylphosphinous acid' / 'phenylphosphonous acid' (P-67.1.1.1): the P(III) acid is named from its
    P=O analogue, whose retained 'phosphinic'/'phosphonic' name becomes '-ous'; None for any other piece."""
    from .core import smiles_to_iupac

    def hydroxy_of(center):
        return [
            n
            for n in center.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and not n.GetFormalCharge()
            and acid.GetBondBetweenAtoms(center.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]

    phosphorus = [a for a in acid.GetAtoms() if a.GetAtomicNum() == 15 and not a.GetFormalCharge() and hydroxy_of(a)]
    if len(phosphorus) != 1 or phosphorus[0].GetDegree() + phosphorus[0].GetTotalNumHs() != 3:
        return None
    center = phosphorus[0]
    hydroxy = hydroxy_of(center)
    organics = [n for n in center.GetNeighbors() if n.GetAtomicNum() in (6, 15)]
    if len(hydroxy) + len(organics) == center.GetDegree() and (len(hydroxy), len(organics)) in ((2, 1), (1, 2)):
        oxidised = Chem.RWMol(acid)
        oxo = oxidised.AddAtom(Chem.Atom(8))
        oxidised.AddBond(center.GetIdx(), oxo, Chem.BondType.DOUBLE)
        name = smiles_to_iupac(Chem.MolToSmiles(oxidised.GetMol()))
        for oxo_word, ous_word in (("phosphonic", "phosphonous"), ("phosphinic", "phosphinous")):
            if name.endswith(oxo_word + " acid"):
                return name[: -len(oxo_word + " acid")] + ous_word + " acid"
    return None


def name_phosphorous_acid(mol):
    """P-67.3.1: a P(III) acid centre carrying phosphanyl groups outranks the phosphane chain ('the acid is senior to
    the heterol')."""
    if len(Chem.GetMolFrags(mol)) > 1 or sum(a.GetAtomicNum() == 15 for a in mol.GetAtoms()) < 2:
        raise UnsupportedStructure("not a phosphorous acid with phosphorus substituents")
    name = _phosphorous_acid_name(mol)
    if name is None:
        raise UnsupportedStructure("not a phosphorous acid with phosphorus substituents")
    return name


_ACID_ENDING = re.compile(r"(?: (?:[A-Za-z]+\d*(?:,[A-Za-z]+\d*)*)-acid| acid)(\)?)$")


def anion_name(acid_name, letters_in_anion=True):
    """P-65.6.1: the '-ic acid' ending becomes '-ate'. The element letters of a peroxy acid stay in the anion term
    unless an ester cites them in front of its organyl group (P-65.6.3.3.7.1), which also moves 'carbo' into the
    parentheses: 'benzene(carbothioperoxoate)'."""
    from ._acid_lexicon import _LETTERED_TAIL

    lettered = _LETTERED_TAIL.search(acid_name)
    if lettered:
        stem = acid_name[: lettered.start()]
        if letters_in_anion:
            return stem + f"({lettered.group(2)}-{lettered.group(1)}ate)"
        if stem.endswith("carbo"):
            return stem[: -len("carbo")] + f"(carbo{lettered.group(1)}ate)"
        return stem + f"({lettered.group(1)}ate)"
    match = _ACID_ENDING.search(acid_name)
    if match is None:
        raise UnsupportedStructure("the acid part has no 'acid' ending to turn into an anion name")
    stem = acid_name[: match.start()]
    close = match.group(1)
    if stem.endswith("phosphoric"):
        return stem[: -len("phosphoric")] + "phosphate" + close
    if stem.endswith("ic"):
        return stem[:-2] + "ate" + close
    if stem.endswith("ic)"):
        return stem[:-3] + "ate)" + close
    if stem.endswith("phosphorous"):
        return stem[: -len("phosphorous")] + "phosphite" + close
    if stem.endswith("ous"):
        return stem[:-3] + "ite" + close
    raise UnsupportedStructure("the acid part has an unexpected acid name ending")


def acid_word(acid_name):
    """The acid name without its class word, for anhydrides (P-65.7.1)."""
    match = _ACID_ENDING.search(acid_name)
    if match is None:
        raise UnsupportedStructure("the acid part has no 'acid' ending")
    return acid_name[: match.start()] + match.group(1) if match.group(1) else acid_name[: match.start()]


def _position_lengths(mol, center):
    """The numbers of chalcogen atoms (1 for -OH, 2 for -OOH) of the single-bonded positions of an acid centre."""
    lengths = set()
    for position, role in _center(mol, center)[2]:
        if role != "Y":
            continue
        y = mol.GetAtomWithIdx(position)
        onward = [n for n in y.GetNeighbors() if n.GetIdx() != center and n.GetAtomicNum() in _SYMBOL and _center(mol, n.GetIdx()) is None]
        lengths.add(2 if onward else 1)
    return lengths


def _free_groups(mol, center):
    """The atoms of each free -OH (or -SH, ...) and -OOH (or -SSH, ...) position of an acid centre, nearest atom first."""
    groups = []
    for n in mol.GetAtomWithIdx(center).GetNeighbors():
        if n.GetAtomicNum() not in _SYMBOL or n.GetFormalCharge() or _bond(mol, center, n.GetIdx()) != 1.0:
            continue
        if n.GetDegree() == 1 and n.GetTotalNumHs() == 1:
            groups.append([n.GetIdx()])
            continue
        onward = [m for m in n.GetNeighbors() if m.GetIdx() != center]
        if (
            not n.GetTotalNumHs()
            and len(onward) == 1
            and onward[0].GetAtomicNum() in _SYMBOL
            and onward[0].GetDegree() == 1
            and onward[0].GetTotalNumHs() == 1
            and not onward[0].GetFormalCharge()
        ):
            groups.append([n.GetIdx(), onward[0].GetIdx()])
    return groups


def _letters(mol, link):
    """Element letters of the atom chain bearing the organyl group when several chalcogens occur at the centre (P-65.6.3.3.7)."""
    found = _center(mol, link.center)
    elements = {x for x in _oxo_slots(mol, link.center) if x in CHALCOGENS}
    for position, role in found[2]:
        if role == "Y":
            y = mol.GetAtomWithIdx(position)
            elements.add(_SYMBOL[y.GetAtomicNum()])
            for n in y.GetNeighbors():
                if n.GetIdx() != link.center and n.GetAtomicNum() in _SYMBOL and position not in link.chain:
                    elements.add(_SYMBOL[n.GetAtomicNum()])
    for a in link.chain:
        elements.add(_SYMBOL[mol.GetAtomWithIdx(a).GetAtomicNum()])
    if len(elements) < 2 and len(_position_lengths(mol, link.center)) < 2:
        return ""
    return "".join(_SYMBOL[mol.GetAtomWithIdx(a).GetAtomicNum()] for a in link.chain)


def _fragments(mol, cut_bonds):
    editable = Chem.RWMol(mol)
    for a, b in cut_bonds:
        editable.RemoveBond(a, b)
    frags = Chem.GetMolFrags(editable)
    owner = {}
    for i, atoms in enumerate(frags):
        for a in atoms:
            owner[a] = i
    return frags, owner


def _stereo_context(mol, removed):
    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    return {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in removed and a.HasProp("_CIPCode")},
        "bonds": {
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
            for b in probe.GetBonds()
            if b.GetBeginAtomIdx() in removed and b.GetEndAtomIdx() in removed and b.HasProp("_CIPCode")
        },
        "used": set(),
    }


def _isocyanide_carbon(atom):
    return (
        atom.GetAtomicNum() == 6
        and atom.GetFormalCharge() == -1
        and atom.GetDegree() == 1
        and atom.GetNeighbors()[0].GetAtomicNum() == 7
        and atom.GetNeighbors()[0].GetFormalCharge() == 1
    )


def in_acid_derivative_scope(mol):
    """Whether the acid derivative engine names `mol` as it stands: a structure with a nuclide, a radical or an ion
    centre away from nitrogen and oxygen is named with that modification, by the namers of those classes."""
    return not any(
        atom.GetIsotope()
        or atom.GetNumRadicalElectrons()
        or atom.GetFormalCharge() and atom.GetAtomicNum() not in (7, 8) and not _isocyanide_carbon(atom)
        for atom in mol.GetAtoms()
    )


def _multiplied(entries):
    """Citation of the organyl groups in alphanumerical order with multipliers and element letters;
    entries are (name, compound, count, letters) tuples."""
    parts = []
    for name, compound, count, letters in sorted(entries, key=lambda e: citation_order_key(e[0])):
        if count > 1:
            multiplier, enclosed = prefix_multiplier(count, name, compound and not _plain(name))
            text = multiplier + (enclose(name) if enclosed or compound else name)
        else:
            text = enclose(name) if compound and any(letters) else name
        parts.append((f"{','.join(sorted(letters, key=lambda t: (len(t.split('-')[0]), t)))}-" if any(letters) else "") + text)
    return " ".join(parts)


def _plain(name):
    from ._substituents import is_plain_stem_prefix

    return is_plain_stem_prefix(name)


def _piece_key(mol, atoms):
    """Canonical text of the piece `atoms`, each cut bond ending in a dummy atom, so that equivalent pieces compare
    equal whatever the atom order."""
    inside = set(atoms)
    editable = Chem.RWMol(mol)
    cut = {n.GetIdx() for i in inside for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in inside}
    for idx in cut:
        editable.GetAtomWithIdx(idx).SetAtomicNum(0)
    for first, second in [(a, b) for a in cut for b in cut if a < b and mol.GetBondBetweenAtoms(a, b) is not None]:
        editable.RemoveBond(first, second)
    for idx in sorted(set(range(mol.GetNumAtoms())) - inside - cut, reverse=True):
        editable.RemoveAtom(idx)
    piece = editable.GetMol()
    piece.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(piece)
    return Chem.MolToSmiles(piece, isomericSmiles=False)


_DIYL = re.compile(r"^(?P<prefix>.*?)benzene-(?P<locants>[\d,]+)-(?P<count>di|tri)ol$")


def _diyl_name(mol, atoms, attachments):
    """The multivalent organyl name of the piece `atoms` joined through the `attachments` atoms, from its polyol name."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    hydroxy = [
        i for i in atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 8 and mol.GetAtomWithIdx(i).GetDegree() == 1 and mol.GetAtomWithIdx(i).GetTotalNumHs() == 1
    ]
    if hydroxy and any(mol.GetAtomWithIdx(i).GetAtomicNum() == 53 for i in atoms):
        return None
    for idx in hydroxy:
        placeholder = editable.GetAtomWithIdx(idx)
        placeholder.SetAtomicNum(53)
        placeholder.SetNumExplicitHs(0)
        placeholder.SetNoImplicit(True)
    for idx in attachments:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(idx, oxygen, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    polyol = editable.GetMol()
    Chem.SanitizeMol(polyol)
    name = smiles_to_iupac(Chem.MolToSmiles(polyol)).replace("iodo", "hydroxy")
    match = _DIYL.match(name)
    if match:
        return f"{match.group('prefix')}{match.group('locants')}-phenylene" if match.group("count") == "di" else None
    if name.endswith("diol"):
        return name[: -len("ol")] + "yl"
    if name.endswith("triol"):
        return name[: -len("ol")] + "yl"
    return _diyl_from_acetates(mol, atoms, attachments)


def _diyl_from_acetates(mol, atoms, attachments):
    """The organyl group as it is cited in the acetate of the same polyol (P-65.6.3.2.2), for a group whose polyol name
    carries the hydroxy groups on a compound multiplying group ('oxydi(ethane-2,1-diyl)')."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for idx in attachments:
        oxygen = editable.AddAtom(Chem.Atom(8))
        carbon = editable.AddAtom(Chem.Atom(6))
        editable.AddBond(idx, oxygen, Chem.BondType.SINGLE)
        editable.AddBond(oxygen, carbon, Chem.BondType.SINGLE)
        editable.AddBond(carbon, editable.AddAtom(Chem.Atom(8)), Chem.BondType.DOUBLE)
        editable.AddBond(carbon, editable.AddAtom(Chem.Atom(6)), Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    acetate = editable.GetMol()
    Chem.SanitizeMol(acetate)
    name = smiles_to_iupac(Chem.MolToSmiles(acetate))
    tail = _ACETATE_TAIL.search(name)
    if tail is None or tail.group(1) != multiplying_prefix(len(attachments)):
        return None
    return name[: tail.start()].rstrip()


def _pendant_text(mol, graph, links):
    halogens = halogen_substituents(mol)
    counts = {}
    for l in links:
        name = name_branch(graph, l.far, l.chain[-1], halogens, mol=mol)[0]
        counts[name] = counts.get(name, 0) + 1
    return _multiplied([(n, not n.isalpha(), k, []) for n, k in counts.items()])


def _multiplicative_ester(mol, graph, esters, frags, owner, acid_pieces, r_pieces, suffixed):
    """'dimethyl ethane-1,2-diyl dibutanedioate' (P-65.6.3.3.4.1): one organyl piece joining identical acid pieces
    whose remaining ester groups carry organyl groups; None when the molecule is not of that shape."""

    hubs = [i for i, ls in r_pieces.items() if len(ls) >= 2 and i not in acid_pieces]
    if len(hubs) == 2 and len(acid_pieces) == 3:
        return _bridged_multiplicative_ester(mol, graph, frags, owner, acid_pieces, r_pieces, hubs)
    if len(hubs) != 1:
        return None
    hub = hubs[0]
    hub_links = r_pieces[hub]
    centers = [owner[l.center] for l in hub_links]
    if len(set(centers)) != len(centers) or len(acid_pieces) != len(centers):
        return None
    if any(l not in suffixed[owner[l.center]] for l in hub_links):
        return None
    acid_atoms = {}
    for c in centers:
        component = _fragments(mol, [(l.chain[-1], l.far) for l in suffixed[c]])
        acid_atoms[c] = component[0][component[1][suffixed[c][0].center]]
    if len({_piece_key(mol, acid_atoms[c]) for c in centers}) != 1:
        return None
    pendants = [l for c in centers for l in suffixed[c] if owner[l.far] != hub]
    if any(len(r_pieces[owner[l.far]]) != 1 for l in pendants):
        return None
    modified = any(_chalcogen_modified(mol, l) for l in hub_links)
    if not pendants and (
        (any(mol.GetAtomWithIdx(l.center).GetAtomicNum() == 6 for l in hub_links) and not modified)
        or any(_free_acid_hydroxyl(mol, a) for c in centers for a in acid_atoms[c])
    ):
        return None
    hydroxyls = {sum(_free_acid_hydroxyl(mol, a) for a in acid_atoms[c]) for c in centers}
    if len(hydroxyls) != 1:
        return None
    first = suffixed[centers[0]]
    acid = _acid_name(mol, acid_atoms[centers[0]], [l.chain[-1] for l in first])
    anion = anion_name(acid, letters_in_anion=False)
    if hydroxyls != {0}:
        anion = f"{multiplying_prefix(next(iter(hydroxyls)))}hydrogen {anion}".replace("monohydrogen", "hydrogen")
    diyl = _diyl_name(mol, frags[hub], [l.far for l in hub_links])
    if diyl is None:
        return None
    n_substituted = any(
        n.GetAtomicNum() == 7 and n.GetDegree() > 1 for n in mol.GetAtomWithIdx(first[0].center).GetNeighbors()
    )
    compound = not anion.isalpha() or n_substituted
    multiplier = multiplying_prefix(len(centers), compound=compound)
    anion_text = f"{multiplier}{enclose(anion)}" if compound else multiplier + anion
    if modified:
        symbols = ["".join(_SYMBOL[mol.GetAtomWithIdx(a).GetAtomicNum()] for a in l.chain) for l in hub_links]
        diyl = ",".join(f"{symbol}{_PRIMES[i]}" for i, symbol in enumerate(sorted(symbols))) + "-" + diyl
    return " ".join(part for part in (_pendant_text(mol, graph, pendants), diyl, anion_text) if part)


_PRIMES = ("", "\u2032", "\u2033", "\u2034", "\u2057")
_ACETATE_TAIL = re.compile(r"(di|tri|tetra|penta|hexa)acetate$")


def _cyanic_polyester(mol, esters, acid_pieces):
    """P-65.6.3.2.2, P-65.6.3.3.7.2.1: 'propane-1,3-diyl dithiocyanate'. The organyl group joining identical cyanic acid
    anions is named as in the acetate of the same polyol, whose cyano-chalcogen groups stand in for the acetyloxy ones; the
    anion word then replaces 'acetate'."""
    from .core import smiles_to_iupac

    if any((_center(mol, l.center) or ("",))[0] != "cyanic" or len(l.chain) != 1 for l in esters):
        return None
    if len({mol.GetAtomWithIdx(l.chain[0]).GetAtomicNum() for l in esters}) != 1:
        return None
    proxy = Chem.RWMol(mol)
    for l in esters:
        proxy.GetAtomWithIdx(l.chain[0]).SetAtomicNum(8)
        nitrogen = next(n.GetIdx() for n in mol.GetAtomWithIdx(l.center).GetNeighbors() if n.GetAtomicNum() == 7)
        proxy.GetAtomWithIdx(nitrogen).SetAtomicNum(8)
        proxy.GetBondBetweenAtoms(l.center, nitrogen).SetBondType(Chem.BondType.DOUBLE)
        proxy.AddBond(l.center, proxy.AddAtom(Chem.Atom(6)), Chem.BondType.SINGLE)
    proxy = proxy.GetMol()
    Chem.SanitizeMol(proxy)
    match_cip_labels(mol, proxy)
    name = smiles_to_iupac(Chem.MolToSmiles(proxy))
    tail = _ACETATE_TAIL.search(name)
    if tail is None:
        return None
    link = esters[0]
    pieces, owner = _fragments(mol, [(link.chain[-1], link.far)])
    acid = _acid_name(mol, pieces[owner[link.center]], [link.chain[-1]])
    return name[: tail.start()] + tail.group(1) + anion_name(acid, letters_in_anion=False)

def _ester_pieces(mol, esters):
    frags, owner = _fragments(mol, [(l.chain[-1], l.far) for l in esters])
    acid_pieces, r_pieces = {}, {}
    for l in esters:
        acid_pieces.setdefault(owner[l.center], []).append(l)
        r_pieces.setdefault(owner[l.far], []).append(l)
    return frags, owner, acid_pieces, r_pieces


def _suffixed_links(mol, frags, acid_pieces):
    return {c: _suffix_links(mol, ls, frags[c], lambda l: l.chain[-1]) if len(ls) > 1 else ls for c, ls in acid_pieces.items()}


_ORGANIC_ELEMENTS_OF_LINKERS = frozenset({1, 6, 7, 8, 9, 15, 16, 17, 34, 35, 52, 53})
_METALLOID_ELEMENTS_OF_LINKERS = frozenset({5, 14, 32, 33, 50, 51, 82, 83})


def _linker_elements_allowed(mol, esters):
    """Boron and the metalloids of groups 14 and 15 may join the organyl groups of the esters; bonded to an ester oxygen
    they make an organometallic ester that is named otherwise."""
    allowed = _ORGANIC_ELEMENTS_OF_LINKERS
    if all(mol.GetAtomWithIdx(l.far).GetAtomicNum() == 6 for l in esters):
        allowed = allowed | _METALLOID_ELEMENTS_OF_LINKERS
    return all(a.GetAtomicNum() in allowed for a in mol.GetAtoms())


def _chalcogen_modified(mol, link):
    """Whether the acid group of `link` has an oxygen replaced by sulfur, selenium or tellurium (P-65.1.3)."""
    return any(x != "O" for x in _oxo_slots(mol, link.center)) or any(
        mol.GetAtomWithIdx(a).GetAtomicNum() != 8 for a in link.chain
    )


def _chalcogen_modified_carboxylic(mol, link):
    """A thio-, seleno- or tellurocarboxylic acid group, which ranks below the carboxylic acid itself (P-41)."""
    return (
        mol.GetAtomWithIdx(link.center).GetAtomicNum() == 6
        and (_center(mol, link.center) or ("",))[0] == "acyl"
        and _chalcogen_modified(mol, link)
    )


def _acid_piece_rank(mol, frags, acid_pieces, r_pieces, piece):
    carboxylic = any(
        mol.GetAtomWithIdx(l.center).GetAtomicNum() == 6
        and (_center(mol, l.center) or ("",))[0] == "acyl"
        and not _chalcogen_modified_carboxylic(mol, l)
        for l in acid_pieces[piece]
    )
    return carboxylic, len(acid_pieces[piece]), _acid_rank(mol, frags[piece], acid_pieces[piece], r_pieces.get(piece, []))


def _polyester_with_acyloxy_groups(mol, links):
    """P-65.6.3.3.4.2: '2-(acetyloxy)ethyl methyl ethane-1,2-diyl dibutanedioate': an acyl group that is the only ester
    of its acid and ranks below the others is an acyloxy prefix of an organyl group, not an anion of the name."""
    if any(l.kind != "ester" for l in links):
        return None
    frags, owner, acid_pieces, r_pieces = _ester_pieces(mol, links)
    if len(acid_pieces) < 3:
        return None
    graph = adjacency(mol)
    esters = links
    ranks = {i: _acid_piece_rank(mol, frags, acid_pieces, r_pieces, i) for i in acid_pieces}
    top = max(ranks.values())
    junior = {i for i, links in acid_pieces.items() if len(links) == 1 and ranks[i] < top}
    if not junior:
        return None
    kept = [l for l in esters if owner[l.center] not in junior]
    frags, owner, acid_pieces, r_pieces = _ester_pieces(mol, kept)
    if len(acid_pieces) != 2:  # three or more anion residues are named by skeletal replacement (P-65.6.3.3.4.3)
        return None
    return _multiplicative_ester(mol, graph, kept, frags, owner, acid_pieces, r_pieces, _suffixed_links(mol, frags, acid_pieces))


def _free_acid_hydroxyl(mol, idx):
    """A hydroxy group on an acid centre, cited as 'hydrogen' in the anion word."""
    atom = mol.GetAtomWithIdx(idx)
    return (
        atom.GetAtomicNum() in _SYMBOL
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and _center(mol, atom.GetNeighbors()[0].GetIdx()) is not None
    )


def _bridged_multiplicative_ester(mol, graph, frags, owner, acid_pieces, r_pieces, hubs):
    """'dimethyl butanedioylbis(oxy-2,1-phenylene)dibutanedioate' (P-65.6.3.3.2.2.1): two identical terminal acids joined
    through identical organyl pieces and a central diacyl group; None for any other arrangement."""
    from ._chain_multiplicative import _make_context
    from ._multiplicative_linker import DecompositionRejected, name_component

    inner = [c for c, ls in acid_pieces.items() if len(ls) == 2 and all(owner[l.far] in hubs for l in ls)]
    if len(inner) != 1:
        return None
    inner = inner[0]
    terminals = [c for c in acid_pieces if c != inner]
    if any(len(acid_pieces[t]) < len(acid_pieces[inner]) for t in terminals):
        return None
    if len({_piece_key(mol, frags[c]) for c in terminals}) != 1 or len({_piece_key(mol, frags[h]) for h in hubs}) != 1:
        return None
    pendants, arms = [], []
    for t in terminals:
        toward = [l for l in acid_pieces[t] if owner[l.far] in hubs]
        if len(toward) != 1:
            return None
        pendants += [l for l in acid_pieces[t] if l is not toward[0]]
        hub = owner[toward[0].far]
        center_link = next(l for l in r_pieces[hub] if owner[l.center] == inner)
        arms.append((hub, toward[0], center_link))
    if len({a[0] for a in arms}) != 2 or any(len(r_pieces[owner[l.far]]) != 1 for l in pendants):
        return None
    ctx = _make_context(mol, graph)
    arm_texts = []
    for hub, unit_link, center_link in arms:
        atoms = list(frags[hub])
        kind = "ring" if any(mol.GetAtomWithIdx(a).IsInRing() for a in atoms) else "carbon"
        attachments = [(unit_link.far, unit_link.chain[-1], 1), (center_link.far, center_link.chain[-1], 1)]
        try:
            part = name_component(mol, kind, atoms, attachments, ctx, (unit_link.far, center_link.far))
        except (DecompositionRejected, UnsupportedStructure):
            return None
        arm_texts.append("oxy" + (enclose(part.text) if part.has_prefix else ("-" if part.text[0].isdigit() else "") + part.text))
    if arm_texts[0] != arm_texts[1]:
        return None
    acyl = acyl_name(_acid_name(mol, frags[inner], [l.chain[-1] for l in acid_pieces[inner]]))
    anion = anion_name(_acid_name(mol, frags[terminals[0]], [l.chain[-1] for l in acid_pieces[terminals[0]]]), letters_in_anion=False)
    multiplier = multiplying_prefix(2, compound=not anion.isalpha())
    anion_text = multiplier + anion if anion.isalpha() else f"{multiplier}{enclose(anion)}"
    linker = f"{acyl}bis{enclose(arm_texts[0])}"
    return " ".join(part for part in (_pendant_text(mol, graph, pendants), linker + anion_text) if part)


def _suffix_links(mol, links, keep, cap_of):
    """The `links` whose acid centre the polyfunctional engine cites as a suffix group when every link of `links`
    is turned into a free acid (P-44.1.1): when one parent cannot carry them all, the other links stay as ester or
    halide groups and are cited as prefixes of that acid. All of `links` when the parent cannot be determined."""
    from ._polyfunctional import _group_of, _principal_class, _ring_occurrences, _try_select

    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetIntProp("_origin", atom.GetIdx())
    for l in links:
        atom = editable.GetAtomWithIdx(cap_of(l))
        atom.SetAtomicNum(8)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    token = ISOTOPE_LABELS.set(None)
    try:
        Chem.SanitizeMol(acid)
        selected = _try_select(acid)
    except Exception:
        return list(links)
    finally:
        ISOTOPE_LABELS.reset(token)
    if selected is None:
        return list(links)
    position_of = selected[2][4]
    origin = {a.GetIdx(): a.GetIntProp("_origin") for a in acid.GetAtoms()}
    parent = {origin[i] for i in position_of if i in origin}
    occurrences = _ring_occurrences(acid)
    class_at = {}
    for a in acid.GetAtoms():
        found = _group_of(acid, a.GetIdx())
        if found is not None:
            class_at[origin[a.GetIdx()]] = found[0]
    for cls, _, owned in occurrences:
        for i in owned:
            class_at.setdefault(origin[i], cls)
    principal = _principal_class(set(class_at.values()))
    kept = []
    for l in links:
        neighbours = {n.GetIdx() for n in mol.GetAtomWithIdx(l.center).GetNeighbors() if n.GetIdx() in keep}
        in_parent = l.center in parent or neighbours & parent
        cls = class_at.get(l.center) or next((class_at[n] for n in neighbours if n in class_at), None)
        if in_parent and (cls is None or principal is None or cls == principal):
            kept.append(l)
    return kept or list(links)


def _acid_rank(mol, atoms, links, entering):
    """(ring parent, parent size, substituents cited as prefixes) of the acid component made of `atoms`, in the order of
    P-44.1.2.2, P-44.3 and P-45.2.1; an ester oxygen entering the component counts as an acyloxy prefix."""
    inside = set(atoms)
    centers = {l.center for l in links}
    seeds = set(centers)
    for c in centers:
        seeds |= {n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() in inside and n.GetAtomicNum() == 6}
    ringed = {i for i in seeds if mol.GetAtomWithIdx(i).IsInRing()}
    parent = set(ringed or seeds)
    stack = list(parent)
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            j = n.GetIdx()
            if j in inside and j not in parent and (n.IsInRing() if ringed else n.GetAtomicNum() == 6 and not n.IsInRing()):
                parent.add(j)
                stack.append(j)
    substituents = len(entering)
    for idx in parent:
        for n in mol.GetAtomWithIdx(idx).GetNeighbors():
            j = n.GetIdx()
            if j in inside and j not in parent and not (idx in centers and n.GetDegree() == 1):
                substituents += 1
    return bool(ringed), len(parent), substituents


_ACID_CLASS_RANK = {"acyl": 0, "carbonic": 0, "cyanic": 0, "inorganic": 1}


_CENTRE_ORDER = (15, 33, 51, 14, 5, 16, 34, 52)


def _acid_kind_rank(mol, idx):
    """Position of the noncarbon acid on the centre `idx` in the order of P-42.3 and P-42.4: acids with a substitutable
    hydrogen or a carbon group on the centre (phosphonic) before those without (phosphoric), then the central atom (P > As >
    Sb > Si > B > S > Se > Te), then the acid with an oxo group (phosphonous after phosphonic), then the one with more acid
    positions (boronic before borinic); () for any other centre."""
    found = _center(mol, idx)
    if found is None or found[0] != "inorganic":
        return ()
    atom = mol.GetAtomWithIdx(idx)
    substituted = atom.GetTotalNumHs() or any(n.GetAtomicNum() not in _SYMBOL for n in atom.GetNeighbors())
    z = atom.GetAtomicNum()
    return (
        0 if substituted else 1,
        _CENTRE_ORDER.index(z) if z in _CENTRE_ORDER else len(_CENTRE_ORDER),
        0 if found[1] else 1,
        -len(found[2]),
    )


def _acid_kind_key(mol, idx, acidic):
    """Kinds of noncarbon acid with a free acid group (`acidic`, a set of ranks) come before the kinds that only occur as
    esters (P-41 class 7 before class 9), each in the order of `_acid_kind_rank`."""
    rank = _acid_kind_rank(mol, idx)
    return (0 if rank in acidic else 1, rank) if rank else ()


def _acidic_kinds(mol):
    return {
        _acid_kind_rank(mol, a.GetIdx())
        for a in mol.GetAtoms()
        if _acid_kind_rank(mol, a.GetIdx()) and _free_groups(mol, a.GetIdx())
    }


def _senior_acid_esters(mol, links):
    """`links` without the esters of the junior acids: an ester takes the rank of its acid (P-41 classes 7a, 7b before 7c;
    P-42.3 within the noncarbon acids), and the esters of the junior acid are prefixes of the organyl group."""
    acidic = _acidic_kinds(mol)
    rank = {
        id(l): (_ACID_CLASS_RANK.get(_center(mol, l.center)[0], 3), _acid_kind_key(mol, l.center, acidic))
        for l in links
        if l.kind == "ester"
    }
    top = min(rank.values(), default=None)
    return [l for l in links if l.kind != "ester" or rank[id(l)] == top]


def _senior_kind_has_no_ester(mol, links):
    """A kind of noncarbon acid that ranks above every one carrying an ester is a free acid, named as the parent elsewhere."""
    acidic = _acidic_kinds(mol)
    chosen = min((_acid_kind_key(mol, l.center, acidic) for l in links if l.kind == "ester"), default=())
    present = [_acid_kind_key(mol, a.GetIdx(), acidic) for a in mol.GetAtoms()]
    return bool(chosen) and any(key and key < chosen for key in present)


def _ester_oxygen_nuclides(mol, links):
    """{ester oxygen: label entry} of the labelled oxygens that join an acid centre to its organyl groups; they are cited
    with the organyl word and in the anion, not as oxygens of the acid (P-82.6.4)."""
    context = ISOTOPE_LABELS.get()
    if not context or context.get("mol") is not mol:
        return {}
    found = {}
    for link in links:
        if len(link.chain) > 2 or any(mol.GetAtomWithIdx(a).GetAtomicNum() != 8 for a in link.chain):
            continue
        for oxygen in link.chain:
            entry = context["labels"].get(oxygen)
            if entry and entry["skeleton"] and not entry["H"]:
                found[oxygen] = entry
    context["consumed"].update(found)
    return found


def _ester_letters(mol, link, ester_nuclides):
    """The locant letters of the organyl group of `link`, a modified oxygen taking its nuclide symbol in place of 'O' (P-82.6.4)."""
    letters = _letters(mol, link)
    if len(link.chain) == 1:
        return ester_nuclides[link.chain[0]]["skeleton"] if link.chain[0] in ester_nuclides else letters
    if letters:
        return "".join(ester_nuclides[a]["skeleton"] if a in ester_nuclides else "O" for a in link.chain)
    return ester_nuclides[link.chain[-1]]["skeleton"] if link.chain[-1] in ester_nuclides else ""


def _cite_ester_oxygens(anion, ester_nuclides):
    from ._ester_by_parts import _insert_before_ending

    counts = {}
    for entry in ester_nuclides.values():
        counts[entry["skeleton"]] = counts.get(entry["skeleton"], 0) + 1
    return _insert_before_ending(anion, "(" + ",".join(f"{n}{c}" for n, c in sorted(counts.items())) + ")")


def name_ester(mol, links):
    if has_acyloxy_amine_group(mol):
        raise UnsupportedStructure("an acyclic nitrogen gives no traditional ester: the group is an acyloxy prefix (P-65.6.3.4.1)")
    links = _senior_acid_esters(mol, links)
    if _senior_kind_has_no_ester(mol, links):
        raise UnsupportedStructure("a free acid of a senior kind outranks the ester (P-41, P-42.3)")
    kinds = {_center(mol, l.center)[0] == "inorganic" for l in links if l.kind == "ester"}
    if kinds == {True} and any(a.GetFormalCharge() for a in mol.GetAtoms()):
        raise UnsupportedStructure("the seniority of this phosphorus or boron ester among the other groups is not decided here")
    esters = [l for l in links if l.kind == "ester"]
    graph = adjacency(mol)
    frags, owner, acid_pieces, r_pieces = _ester_pieces(mol, esters)
    if len(acid_pieces) > 1:
        cyanic = _cyanic_polyester(mol, esters, acid_pieces)
        if cyanic is not None:
            return cyanic
        suffixed = _suffixed_links(mol, frags, acid_pieces)
        multiplicative = _multiplicative_ester(mol, graph, esters, frags, owner, acid_pieces, r_pieces, suffixed)
        if multiplicative is not None:
            return multiplicative
        if len({owner[l.far] for ls in suffixed.values() for l in ls}) == 1:
            from ._diester_acyloxy import name_diester_acyloxy

            if not any(_chalcogen_modified_carboxylic(mol, l) for ls in suffixed.values() for l in ls):
                try:
                    return name_diester_acyloxy(mol)
                except UnsupportedStructure:
                    if all(_ACID_CLASS_RANK.get(_center(mol, l.center)[0], 3) == 0 for l in esters) and _linker_elements_allowed(mol, esters):
                        from ._chain_multiplicative import ester_multiplicative_name

                        multiplicative = ester_multiplicative_name(mol, [l.chain[-1] for l in esters])
                        if multiplicative is not None:
                            return multiplicative
                    if len({_piece_key(mol, frags[i]) for i in acid_pieces}) == 1:
                        raise
    principal = max(acid_pieces, key=lambda i: _acid_piece_rank(mol, frags, acid_pieces, r_pieces, i))
    chosen = acid_pieces[principal]
    component = _fragments(mol, [(l.chain[-1], l.far) for l in chosen])
    acid_atoms = component[0][component[1][chosen[0].center]]
    if len(chosen) > 1:
        suffixed = _suffix_links(mol, chosen, acid_atoms, lambda l: l.chain[-1])
        if len(suffixed) < len(chosen):
            chosen = suffixed
            component = _fragments(mol, [(l.chain[-1], l.far) for l in chosen])
            acid_atoms = component[0][component[1][chosen[0].center]]
    removed = set()
    for l in chosen:
        far_side = _far_side(graph, l.far, l.chain[-1])
        if far_side & set(acid_atoms):
            raise UnsupportedStructure("a ring-closing ester is not named part-wise")
        removed |= far_side
    ester_nuclides = _ester_oxygen_nuclides(mol, chosen)
    acid = _acid_name(mol, acid_atoms, [l.chain[-1] for l in chosen], ester_nuclides)
    if _other_inorganic_group_is_prefix(mol, acid_atoms, chosen, acid):
        raise UnsupportedStructure("an acid group of the same kind on another part of the skeleton is a prefix of the acid name")
    anion = anion_name(acid, letters_in_anion=False)
    if ester_nuclides:
        anion = _cite_ester_oxygens(anion, ester_nuclides)
    context = _stereo_context(mol, removed)
    halogens = halogen_substituents(mol)
    named = {}
    per_link = []
    token = BRANCH_STEREO.set(context)
    prefixes = EXTENDED_PREFIXES.set(True)
    try:
        for l in chosen:
            name, compound = name_branch(graph, l.far, l.chain[-1], halogens, mol=mol)
            entry = named.setdefault(name, [0, compound, []])
            entry[0] += 1
            entry[2].append(_ester_letters(mol, l, ester_nuclides))
            per_link.append((l, name))
    finally:
        EXTENDED_PREFIXES.reset(prefixes)
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of the organyl part is not cited by any supported name")
    acid_class = _center(mol, chosen[0].center)[0]
    chosen_rank = _acid_kind_rank(mol, chosen[0].center)
    free_slots = [
        (a, group[0])
        for a in acid_atoms
        if (_center(mol, a) or (None,))[0] == acid_class and _acid_kind_rank(mol, a) == chosen_rank
        for group in _free_groups(mol, a)
    ]
    free_centers = [center for center, _ in free_slots]
    free = len(free_centers)
    hydrogen_positions = None
    polycarbonic = _polycarbonic_locants(mol, acid_atoms, per_link, free_slots) if len(chosen) + free > 1 else None
    if polycarbonic is not None:
        cited, hydrogen_positions = polycarbonic
    else:
        cited = _ester_locants(mol, acid_atoms, per_link, free_centers) if (len(named) > 1 or free) and len(chosen) + free > 1 else None
        if cited is not None:
            token = CITE_SUBSTITUENT_LOCANTS.set(True)
            try:
                anion = anion_name(_acid_name(mol, acid_atoms, [l.chain[-1] for l in chosen], ester_nuclides), letters_in_anion=False)
                if ester_nuclides:
                    anion = _cite_ester_oxygens(anion, ester_nuclides)
            finally:
                CITE_SUBSTITUENT_LOCANTS.reset(token)
        else:
            cited = _diphosphate_locants(mol, acid_atoms, per_link, free_centers)
    if not cited and len({l.center for l, _ in per_link}) > 1 and any(any(letters) for _, _, letters in named.values()):
        raise UnsupportedStructure("the element letters of several acid centres need primes and centre locants")
    entries = [
        (name, compound, count, cited[name] if cited else letters) for name, (count, compound, letters) in named.items()
    ]
    hydrogen = "" if free == 0 else multiplying_prefix(free) + "hydrogen " if free > 1 else "hydrogen "
    if hydrogen_positions:
        hydrogen = f"{','.join(hydrogen_positions)}-{hydrogen}"
    return _multiplied(entries) + " " + hydrogen + anion


_ACID_GROUP_PREFIX = re.compile(
    r"(?<![a-z])(?:borono|phosphono|arsono|stibono)"
    r"(?!us\b|thio|seleno|telluro|peroxo|chlorid|bromid|fluorid|iodid|amid|azid|cyan|isocyan|isothio|nitrid|imid|hydrazid)"
)


def _other_inorganic_group_is_prefix(mol, acid_atoms, chosen, acid_name):
    """Another acid centre of the kind of the chosen ones, free or with chosen esters, that the free acid name cites as a
    borono, phosphono, arsono or stibono prefix: the acid suffix of the ester name then holds fewer positions than the esters
    to be multiplied, whose other groups would be lost (P-65.6.3.2.3 names them as prefixes of the organyl group)."""
    kind = _acid_kind_rank(mol, chosen[0].center)
    if not kind or _center(mol, chosen[0].center)[0] != "inorganic" or not _ACID_GROUP_PREFIX.search(acid_name):
        return False
    named = {l.center for l in chosen}
    return sum(_acid_kind_rank(mol, a) == kind and (a in named or bool(_free_groups(mol, a))) for a in acid_atoms) > 1


def _ester_locants(mol, acid_atoms, per_link, free_centers):
    """{organyl name: acid locants} when the acid positions are not equivalent (P-65.6.3.2.2: '16-ethyl 18-methyl
    yohimban-16,18-dicarboxylate'), else None."""
    ends = {a for l, _ in per_link for a in l.chain}
    centers = sorted({l.center for l, _ in per_link} | set(free_centers))
    free_groups = {a for c in free_centers for group in _free_groups(mol, c) for a in group}
    core = [a for a in acid_atoms if a not in ends and a not in free_groups]
    locants, symmetric = _locants_of(mol, core, centers, with_symmetry=True)
    if symmetric:
        return None
    cited = {}
    for l, name in per_link:
        cited.setdefault(name, []).append(str(locants[l.center]))
    return {name: sorted(locs, key=lambda t: (len(t), t)) for name, locs in cited.items()}


def _polycarbonic_locants(mol, acid_atoms, per_link, free):
    """P-65.6.3.3.7.2.1: ({organyl name: positions}, positions of the hydrogen atoms) for an ester of a di- to
    tetracarbonic acid, each position a locant of the acid chain followed by the element symbol when the carbonic unit
    carries more than one chalcogen ('1-S-methyl'); (None, None) when every free position is equivalent; None when the
    acid is not a polycarbonic acid."""
    from ._polycarbonic import end_numberings

    acid, index = _capped_acid(mol, acid_atoms, [l.chain[-1] for l, _ in per_link])
    try:
        numberings = end_numberings(acid)
    except UnsupportedStructure:
        return None
    slots = [(l.center, l.chain[-1], l, name) for l, name in per_link]
    slots += [(center, oxygen, Link("ester", center, (oxygen,), -1), "") for center, oxygen in free]
    ranks = Chem.CanonicalRankAtoms(acid, breakTies=False)
    if len({ranks[index[oxygen]] for _, oxygen, _, _ in slots}) < 2:
        return None, None
    slots.sort(key=lambda slot: (not slot[3], alpha_sort_key(slot[3])))
    best = min(tuple(locants[index[center]] for center, _, _, _ in slots) for _, locants in numberings)
    cited, hydrogen = {}, []
    for locant, (_, _, link, name) in zip(best, slots):
        symbol = _letters(mol, link)
        position = f"{locant}-{symbol}" if symbol else str(locant)
        (cited.setdefault(name, []) if name else hydrogen).append(position)
    return cited, hydrogen


def _diphosphate_locants(mol, acid_atoms, per_link, free_centers):
    """{organyl name: locants} for the esters of a diphosphate P1-O2-P3 when its organyl groups and hydrogen atoms can sit
    on the two phosphorus atoms in more than one way, which the bare name cannot tell apart (P-65.6.3.3.2.1: '1,1-diethyl
    3-methyl butane-1,1,3-tricarboxylate'); None when the arrangement is determined."""
    phosphorus = [a for a in acid_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 15]
    if len(phosphorus) != 2:
        return None
    around = [{n.GetIdx() for n in mol.GetAtomWithIdx(p).GetNeighbors()} for p in phosphorus]
    if len(around[0] & around[1]) != 1 or any(mol.GetAtomWithIdx(b).GetAtomicNum() != 8 for b in around[0] & around[1]):
        return None
    slots = {p: [] for p in phosphorus}
    for link, name in per_link:
        slots[link.center].append(name)
    for center in free_centers:
        slots[center].append("")
    if any(len(entries) != 2 for entries in slots.values()):
        return None
    items = [name for entries in slots.values() for name in entries]
    arrangements = {
        tuple(sorted([tuple(sorted(order[:2])), tuple(sorted(order[2:]))])) for order in itertools.permutations(items)
    }
    if len(arrangements) == 1:
        return None
    best = None
    for first, second in (phosphorus, phosphorus[::-1]):
        locant = {first: "1", second: "3"}
        cited = {}
        for link, name in per_link:
            cited.setdefault(name, []).append(locant[link.center])
        cited = {name: sorted(found) for name, found in cited.items()}
        key = (
            sorted(x for found in cited.values() for x in found),
            [x for name in sorted(cited, key=alpha_sort_key) for x in cited[name]],
        )
        if best is None or key < best[0]:
            best = (key, cited)
    return best[1]


def _far_side(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _component_acid(mol, atoms, centers):
    """Name of the free acid made of `atoms` with a hydroxy group added on every centre in `centers`
    (the acyl centres that were joined by an anhydride bridge)."""
    editable = Chem.RWMol(mol)
    for center in centers:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(center, oxygen, Chem.BondType.SINGLE)
    _restore_nuclides(editable, mol, atoms)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    match_cip_labels(mol, acid, dict(enumerate(sorted(atoms))))
    for atom in acid.GetAtoms():
        if atom.GetAtomicNum() in (15, 16):
            # the descriptor of a stereogenic acid centre is cited in front of the anion, not inside the acid name
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    # P-14.3.4.1 lists no anhydride, so the substituent locants of the acid piece stay cited
    token = CITE_SUBSTITUENT_LOCANTS.set(True)
    systematic = SYSTEMATIC_ACID_PROBE.set(True)
    try:
        return _free_acid_name(acid)
    finally:
        SYSTEMATIC_ACID_PROBE.reset(systematic)
        CITE_SUBSTITUENT_LOCANTS.reset(token)


def _bridged_components(mol, links):
    """({bridge: (centres, chain)}, fragments, owner) for the anhydride bridges of the molecule."""
    bridges = {}
    for l in links:
        if l.kind == "anhydride":
            bridges.setdefault(frozenset(l.chain), [l.chain, l.center, l.far])
    if not bridges:
        raise UnsupportedStructure("no anhydride bridge")
    cuts = []
    for chain, first, second in bridges.values():
        cuts += [(first, chain[0]), (second, chain[-1])]
    frags, owner = _fragments(mol, cuts)
    return bridges, frags, owner


class _UnlocatedGroups(UnsupportedStructure):
    def __init__(self, message, centers):
        super().__init__(message)
        self.centers = centers


def _locants_of(mol, atoms, centers, with_symmetry=False):
    """{centre: locant} of the acid groups `centers` in the substitutive name of the acid made of `atoms`;
    with `with_symmetry` also whether all the centres are equivalent in that acid (then no locants are needed)."""
    from ._polyfunctional import _check_scope, _try_select

    editable = Chem.RWMol(mol)
    for center in centers:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(center, oxygen, Chem.BondType.SINGLE)
    for atom in editable.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    if with_symmetry:
        plain = Chem.Mol(acid)
        mapped = {a.GetAtomMapNum() - 1: a.GetIdx() for a in acid.GetAtoms()}
        for atom in plain.GetAtoms():
            atom.SetAtomMapNum(0)
        ranks = Chem.CanonicalRankAtoms(plain, breakTies=False)
        if len({ranks[mapped[c]] for c in centers if c in mapped}) == 1:
            return {}, True
    selected = _try_select(acid, stereo=_check_scope(acid))
    if selected is None:
        raise _UnlocatedGroups("the locants of the anhydride groups are not available", list(centers))
    parts = selected[2]
    position = {acid.GetAtomWithIdx(i).GetAtomMapNum() - 1: loc for i, loc in parts[4].items()}
    found, unlocated = {}, []
    for center in centers:
        if center in position:
            found[center] = position[center]
            continue
        attach = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors() if n.GetIdx() in position]
        if len(attach) != 1:
            unlocated.append(center)
            continue
        found[center] = position[attach[0]]
    if unlocated:
        raise _UnlocatedGroups("the locants of the anhydride groups are not available", unlocated)
    return (found, False) if with_symmetry else found


def _located_acid_word(mol, atoms, centers):
    """The acid word of the component made of `atoms` with every locant cited, for a core whose acid groups are
    told apart by locants (P-65.7.6.2)."""
    from ._polyfunctional import FORCE_LOCANTS, _try_select

    editable = Chem.RWMol(mol)
    for center in centers:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(center, oxygen, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    token = FORCE_LOCANTS.set(True)
    try:
        selected = _try_select(acid)
    finally:
        FORCE_LOCANTS.reset(token)
    if selected is None:
        raise UnsupportedStructure("the acid groups of this core are told apart by locants that the parent hydride does not give")
    return acid_word(selected[1])


def _multiplied_word(word, count):
    if count == 1:
        return word
    return multiplying_prefix(count, compound=not word.isalpha()) + (word if word.isalpha() else enclose(word))


def _component_word(name):
    """The acid word of an anhydride component; a component that is itself an ester is cited whole in parentheses
    (P-67.1.3.3), as in '(methyl dihydrogen phosphate)'."""
    if _ACID_ENDING.search(name):
        return acid_word(name)
    if " " not in name:
        raise UnsupportedStructure("an anhydride component is an acid or an ester")
    return f"({name})"


_CHALCOGEN_SENIORITY = ("O", "S", "Se", "Te")


def _senior_bridge_links(mol, links):
    """P-65.7.6.4.2: when single-chalcogen linkages of different elements occur, the most senior (O > S > Se > Te)
    names the anhydride and the others stay in the acid pieces as substituent groups."""
    anhydrides = [l for l in links if l.kind == "anhydride"]
    kinds = {id(l): _SYMBOL[mol.GetAtomWithIdx(l.chain[0]).GetAtomicNum()] for l in anhydrides if len(l.chain) == 1}
    if len(kinds) != len(anhydrides) or len(set(kinds.values())) < 2:
        return links
    senior = min(kinds.values(), key=_CHALCOGEN_SENIORITY.index)
    return [l for l in links if l.kind != "anhydride" or kinds[id(l)] == senior]


def _acid_group_count(mol, atoms, centers):
    """Acid groups of an anhydride component: its hydroxy groups on oxo-bearing atoms plus the bridged ones (P-65.7.2)."""
    hydroxy = sum(
        1
        for atom in atoms
        for n in mol.GetAtomWithIdx(atom).GetNeighbors()
        if n.GetAtomicNum() == 8
        and n.GetDegree() == 1
        and n.GetTotalNumHs() == 1
        and n.GetIdx() in atoms
        and (atom in centers or any(_bond(mol, atom, m.GetIdx()) == 2.0 for m in mol.GetAtomWithIdx(atom).GetNeighbors()))
    )
    return hydroxy + len(centers)


def _mononuclear_noncarbon_center(mol, centers):
    return len(set(centers)) == 1 and mol.GetAtomWithIdx(centers[0]).GetAtomicNum() != 6


def _alphabetical_acid_words(counts, class_word):
    """P-67.1.3.3: the acids of an anhydride of a mononuclear noncarbon oxoacid are cited alphabetically; one kind of
    acyl group is cited once, mixed kinds keep their multipliers so that the name stays unambiguous."""
    cited = [_multiplied_word(w, counts[w]) for w in sorted(counts, key=citation_order_key)]
    return f"{' '.join(cited)} {class_word}"


def _bridge_nuclide(mol, bridges, owner, words):
    """The descriptor of the labelled oxygens of an -O- or -O-O- bridge, cited in front of the class word that denotes
    them (P-82.2.1); the oxygens of a peroxy bridge are numbered from the acid component cited first (P-82.5.2). The
    oxygens of the acid components carry theirs in the acid words."""
    from ._isotope_labels import descriptor

    context = ISOTOPE_LABELS.get()
    if not context or context.get("mol") is not mol or len(bridges) != 1:
        return ""
    (chain, first, _), = bridges.values()
    labelled = {a: context["labels"][a] for a in chain if a in context["labels"]}
    if not labelled or any(e["H"] or not e["skeleton"] for e in labelled.values()):
        return ""
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 8 for a in chain):
        return ""
    cited_first = min(words, key=lambda i: citation_order_key(words[i]))
    forward = chain if owner[first] == cited_first else chain[::-1]
    orders = (forward, forward[::-1]) if len(set(words.values())) == 1 else (forward,)
    positions = min(({a: i + 1 for i, a in enumerate(o)} for o in orders), key=lambda p: sorted(p[a] for a in labelled))
    context["consumed"].update(labelled)
    return descriptor(labelled, positions, len(chain) == 1)


def name_anhydride(mol, links):
    links = _senior_bridge_links(mol, links)
    bridges, frags, owner = _bridged_components(mol, links)
    kinds = {tuple(_SYMBOL[mol.GetAtomWithIdx(a).GetAtomicNum()] for a in chain) for chain, _, _ in bridges.values()}
    if len(kinds) != 1 or next(iter(kinds)) not in _ANHYDRIDE_WORDS:
        raise UnsupportedStructure("the anhydride bridges are not all of one kind named by P-65.7")
    word = _ANHYDRIDE_WORDS[next(iter(kinds))]
    bridge_atoms = set().union(*(set(chain) for chain, _, _ in bridges.values()))
    components = [i for i, atoms in enumerate(frags) if not set(atoms) <= bridge_atoms]
    adjacency_of = {i: [] for i in components}
    centers = {i: [] for i in components}
    for chain, first, second in bridges.values():
        a, b = owner[first], owner[second]
        if a == b:
            raise UnsupportedStructure("a cyclic anhydride is named as a heterocycle")
        adjacency_of[a].append(b)
        adjacency_of[b].append(a)
        centers[a].append(first)
        centers[b].append(second)
    if len(components) != len(bridges) + 1:
        raise UnsupportedStructure("the acid components are not joined in a tree")
    count = len(bridges)
    class_word = word if count == 1 else (
        multiplying_prefix(count) + "anhydride" if word == "anhydride" else f"{'bis' if count == 2 else 'tris'}({word})"
    )
    degrees = {i: len(adjacency_of[i]) for i in components}
    word_of = {
        i: _component_word(_component_acid(mol, [a for a in frags[i] if a not in bridge_atoms], centers[i]))
        for i in components
    }
    if count == 1:
        words = list(word_of.values())
        text = words[0] if words[0] == words[1] else " ".join(sorted(words, key=citation_order_key))
        if word == "anhydride" and any(
            _acid_group_count(mol, [a for a in frags[i] if a not in bridge_atoms], centers[i]) > 1 for i in components
        ):
            class_word = "monoanhydride"
        return f"{text} {_bridge_nuclide(mol, bridges, owner, word_of)}{class_word}"
    hubs = [i for i in components if degrees[i] > 2]
    if hubs:
        hub = min(
            hubs,
            key=lambda i: (-degrees[i], -sum(1 for a in frags[i] if a not in bridge_atoms), citation_order_key(word_of[i])),
        )
        if any(degrees[i] != 1 for i in components if i != hub):
            hub_keys = {key for key, (_, first, second) in bridges.items() if hub in (owner[first], owner[second])}
            return name_anhydride(
                mol, [l for l in links if l.kind != "anhydride" or frozenset(l.chain) in hub_keys]
            )
        hub_keep = [a for a in frags[hub] if a not in bridge_atoms]
        carbon_hub = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in hub_keep) and all(
            mol.GetAtomWithIdx(c).GetAtomicNum() == 6 for c in centers[hub]
        )
        try:
            locants = _locants_of(mol, hub_keep, centers[hub]) if carbon_hub else {}
        except _UnlocatedGroups as unlocated:
            # P-65.7.6.3: an anhydride group of a substituent is cited by its acyloxy prefix, not by the class name
            kept = [
                l for l in links
                if l.kind != "anhydride" or not {l.center, l.far} & set(unlocated.centers)
            ]
            if len(kept) == len(links):
                raise
            return name_anhydride(mol, kept)
        groups = {}
        for chain, first, second in bridges.values():
            leaf_center, hub_center = (first, second) if owner[first] != hub else (second, first)
            leaf = owner[leaf_center]
            groups.setdefault(word_of[leaf], []).append(str(locants.get(hub_center, "")))
        if not carbon_hub and len(groups) > 1 and len(set(centers[hub])) > 1:
            raise UnsupportedStructure("different acyl groups on the acid centres of a hub need centre locants")
        if _mononuclear_noncarbon_center(mol, centers[hub]):
            counts = {w: len(g) if len(groups) > 1 else 1 for w, g in groups.items()}
            return _alphabetical_acid_words({**counts, word_of[hub]: 1}, class_word)
        cited = []
        for leaf_word in sorted(groups, key=citation_order_key):
            group = sorted(groups[leaf_word], key=lambda t: (len(t), t))
            body = _multiplied_word(leaf_word, len(group))
            if carbon_hub and len(group) == 1 and not leaf_word.isalpha():
                body = enclose(leaf_word)
            cited.append(f"{','.join(group)}-{body}" if carbon_hub else body)
        return f"{' '.join(cited)} {word_of[hub]} {class_word}"
    ends = [i for i in components if degrees[i] == 1]
    order = [min(ends, key=lambda i: citation_order_key(word_of[i]))]
    while len(order) < len(components):
        order.append(next(n for n in adjacency_of[order[-1]] if n not in order))
    if len(order) == 3:
        w = [word_of[i] for i in order]
        if _mononuclear_noncarbon_center(mol, centers[order[1]]):
            return _alphabetical_acid_words(dict.fromkeys(w, 1), class_word)
        if w[0] == w[2]:
            return f"{_multiplied_word(w[0], 2)} {w[1]} {class_word}"
        return f"{' '.join(w)} {class_word}"
    if len(order) == 4:
        senior = {
            i: (
                min(mol.GetAtomWithIdx(c).GetAtomicNum() for c in centers[i]),
                -sum(1 for a in frags[i] if a not in bridge_atoms),
            )
            for i in order[1:3]
        }
        if senior[order[2]] < senior[order[1]]:
            order = order[::-1]
    tail = order[2:]
    tail_atoms = set()
    for i in tail:
        tail_atoms |= set(frags[i])

    def bridge_between(first_component, second_component):
        return next(
            (chain, first, second)
            for chain, first, second in bridges.values()
            if {owner[first], owner[second]} == {first_component, second_component}
        )

    def center_in(bridge, component):
        return bridge[1] if owner[bridge[1]] == component else bridge[2]

    tail_center = center_in(bridge_between(order[1], order[2]), order[2])
    for chain, first, second in bridges.values():
        if owner[first] in tail and owner[second] in tail:
            tail_atoms |= set(chain)
    merged = acid_word(_component_acid(mol, tail_atoms, [tail_center]))
    core = word_of[order[1]]
    items = [[word_of[order[0]], ""], [core, ""], [merged, ""]]
    core_atoms = [a for a in frags[order[1]] if a not in bridge_atoms]
    if core_atoms and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in centers[order[1]]):
        locants, symmetric = _locants_of(mol, core_atoms, centers[order[1]], with_symmetry=True)
        if not symmetric:
            items[1][0] = _located_acid_word(mol, core_atoms, centers[order[1]])
            items[0][1] = str(locants[center_in(bridge_between(order[0], order[1]), order[1])])
            items[2][1] = str(locants[center_in(bridge_between(order[1], order[2]), order[1])])
    cited = [
        f"{locant}-{enclose(word) if re.search(r'[-(]', word) else word}" if locant else word
        for word, locant in sorted(items, key=lambda item: citation_order_key(item[0]))
    ]
    return f"{' '.join(cited)} dianhydride"


_RETAINED_ACYL = {
    "oxamic": "oxamoyl",
    "carbamic": "carbamoyl",
    "carbonic": "carbonyl",
    "carbamimidic": "carbamimidoyl",
    "carbonimidic": "carbonimidoyl",
}


def acyl_name(acid_name):
    """P-65.1.7: an acid name turned into its acyl group name."""
    match = _ACID_ENDING.search(acid_name)
    if match is None:
        raise UnsupportedStructure("the acid has no 'acid' ending to turn into an acyl name")
    stem = acid_name[: match.start()]
    close = match.group(1)
    for word, acyl in _RETAINED_ACYL.items():
        if stem.endswith(word):
            return stem[: -len(word)] + acyl + close
    if stem.endswith("carboxylic"):
        return stem[: -len("carboxylic")] + "carbonyl" + close
    if stem.endswith("oic"):
        return stem[:-3] + "oyl" + close
    for ending, replacement in (("imidic", "imidoyl"), ("hydrazonic", "hydrazonoyl")):
        if stem.endswith(ending):
            return stem[: -len(ending)] + replacement + close
    if stem.endswith("idic"):
        return stem[:-2] + "oyl" + close
    if stem.endswith("ic"):
        return stem[:-2] + "yl" + close
    raise UnsupportedStructure("the acid name has an unexpected ending")


_PSEUDOHALIDE_WORDS = {
    "N3": "azide",
    "CN": "cyanide",
    "NC": "isocyanide",
    "NCO": "isocyanate",
    "NCS": "isothiocyanate",
    "NCSe": "isoselenocyanate",
    "NCTe": "isotellurocyanate",
}


def _class_words(mol, links):
    """The halide words of an acyl or sulfonyl halide, a modified halogen cited in front of the word that holds it and
    before the unmodified ones (P-82.2.4)."""
    counts = {}
    for l in links:
        far = mol.GetAtomWithIdx(l.far)
        if l.kind == "halide":
            word = _HALIDES[far.GetAtomicNum()]
        else:
            word = _PSEUDOHALIDE_WORDS[pseudohalide_at(mol, l.far, l.center)]
        nuclide = f"{far.GetIsotope()}{far.GetSymbol()}" if far.GetIsotope() else ""
        counts[(word, nuclide)] = counts.get((word, nuclide), 0) + 1
    words = []
    for word, nuclide in sorted(counts, key=lambda key: (alpha_sort_key(key[0]), not key[1])):
        count = counts[(word, nuclide)]
        text = multiplying_prefix(count) + word if count > 1 else word
        words.append(f"({nuclide}{count if count > 1 else ''}){text}" if nuclide else text)
    return " ".join(words)


def _halocarbonic_amide(mol, graph):
    """'carbonochloridic amide' for X-C(=O)-NR2 (P-66.1.1.1.2.2): the amide of the halogenated carbonic acid."""
    from ._substituents import format_substituent_prefixes

    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6 or carbon.GetDegree() != 3 or carbon.IsInRing() or carbon.GetFormalCharge():
            continue
        halide = [n for n in carbon.GetNeighbors() if is_halide_ligand(n)]
        oxo = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and _bond(mol, carbon.GetIdx(), n.GetIdx()) == 2.0]
        amino = [
            n
            for n in carbon.GetNeighbors()
            if n.GetAtomicNum() == 7
            and not n.GetFormalCharge()
            and not n.IsInRing()
            and not any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() != 6 for b in n.GetBonds())
        ]
        if len(halide) != 1 or len(oxo) != 1 or len(amino) != 1:
            continue
        nitrogen = amino[0]
        substituents = [n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
        halogens = halogen_substituents(mol)
        grouped = {}
        arms = set().union(*(_far_side(graph, root, nitrogen.GetIdx()) for root in substituents)) if substituents else set()
        context = _stereo_context(mol, arms)
        token = BRANCH_STEREO.set(context)
        try:
            for root in substituents:
                name, compound = name_branch(graph, root, nitrogen.GetIdx(), halogens, mol=mol, unsaturated=True)
                grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append("N")
        finally:
            BRANCH_STEREO.reset(token)
        if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
            ("bond", b) not in context["used"] for b in context["bonds"]
        ):
            continue
        covered = {carbon.GetIdx(), halide[0].GetIdx(), oxo[0].GetIdx(), nitrogen.GetIdx()}
        if has_group_senior_to_urea(mol, covered) or any(
            a.GetIdx() not in covered and not is_urea_substituent_root(mol, a) for a in mol.GetAtoms()
        ):
            continue
        acid = f"carbono{_HALIDES[halide[0].GetAtomicNum()][:-3]}idic"
        return (format_substituent_prefixes(grouped) if grouped else "") + f"{acid} amide"
    return None


def name_acyl_halide(mol, links):
    amide = _halocarbonic_amide(mol, adjacency(mol))
    if amide is not None:
        return amide
    halides = [l for l in links if l.kind == "halide"] or [l for l in links if l.kind == "pseudohalide"]
    if not halides:
        raise UnsupportedStructure("no acyl halide group")
    graph = adjacency(mol)
    centers = {l.center for l in halides}
    if len(centers) == 1 and _center(mol, next(iter(centers)))[0] == "cyanic":
        return "carbononitridic " + _class_words(mol, halides)
    if len(halides) > 1:
        beyond = set().union(*(_far_side(graph, l.far, l.center) - {l.far} for l in halides if l.kind == "pseudohalide"))
        suffixed = _suffix_links(mol, halides, set(range(mol.GetNumAtoms())) - beyond, lambda l: l.far)
        if len(suffixed) < len(halides):
            halides = suffixed
    editable = Chem.RWMol(mol)
    drop = set()
    for l in halides:
        if l.kind == "pseudohalide":
            drop |= _far_side(graph, l.far, l.center)
        else:
            atom = editable.GetAtomWithIdx(l.far)
            atom.SetAtomicNum(8)
            atom.SetNumExplicitHs(1)
            atom.SetNoImplicit(True)
    for l in halides:
        if l.kind == "pseudohalide":
            oxygen = editable.AddAtom(Chem.Atom(8))
            editable.AddBond(l.center, oxygen, Chem.BondType.SINGLE)
    for idx in sorted(drop, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    survivors = [i for i in range(mol.GetNumAtoms()) if i not in drop]
    match_cip_labels(mol, acid, dict(enumerate(survivors)))
    from .core import smiles_to_iupac

    acid_text = smiles_to_iupac(Chem.MolToSmiles(acid))
    if not _ACID_ENDING.search(acid_text):
        token = SYSTEMATIC_ACID_PROBE.set(True)
        try:
            acid_text = smiles_to_iupac(Chem.MolToSmiles(acid))
        finally:
            SYSTEMATIC_ACID_PROBE.reset(token)
    if re.fullmatch(r"[a-z0-9,\-]*(?:di|tri|tetra)carbonic acid", acid_text):
        return f"{acid_text[: -len(' acid')]} {_class_words(mol, halides)}"
    return f"{acyl_name(acid_text)} {_class_words(mol, halides)}"


_LEADING_HYDROGEN_NUCLIDES = re.compile(r"^\((?:[\d,]+-)?(\d+H)(\d*)\)")


def has_isotope_acyl_halide_shape(mol):
    """A single acyl or sulfonyl halide group whose isotopic labels all sit on carbon or hydrogen atoms, the centre of the
    group, its doubly bonded oxygens or its halogens."""
    if not any(a.GetIsotope() for a in mol.GetAtoms()) or len(Chem.GetMolFrags(mol)) != 1:
        return False
    links, acid = find_links(_without_nuclides(mol))
    if not links or acid or not all(l.kind in ("halide", "pseudohalide") for l in links):
        return False
    carbonyl_oxygens = {n.GetIdx() for l in links for n in mol.GetAtomWithIdx(l.center).GetNeighbors() if n.GetAtomicNum() == 8}
    centres = {l.center for l in links}
    halogens = {l.far for l in links}
    return all(
        a.GetAtomicNum() in (1, 6) or a.GetIdx() in carbonyl_oxygens or a.GetIdx() in centres or a.GetIdx() in halogens
        or (a.GetAtomicNum() in _HALIDES and a.GetDegree() == 1)
        for a in mol.GetAtoms()
        if a.GetIsotope()
    )


def _without_nuclides(mol):
    plain = Chem.RWMol(mol)
    for atom in plain.GetAtoms():
        atom.SetIsotope(0)
    return plain.GetMol()


def name_isotope_acyl_halide(mol):
    """P-82.2.1, P-82.6.1.1: the nuclide descriptor of an acyl halide is that of its acid, whose hydroxy group stands for the
    halogen; deuterium on the only hydrogen-bearing atom needs no locant, because the acid's O-H is not there."""
    from .core import nuclides_cited, smiles_to_iupac

    plain = _without_nuclides(mol)
    amide = _halocarbonic_amide(plain, adjacency(plain))
    if amide is not None:
        labelled = _labelled_halocarbonic_amide(mol, amide)
        if labelled is not None:
            return labelled
    links, _ = find_links(plain)
    editable = Chem.RWMol(mol)
    graph = adjacency(mol)
    dropped = set()
    for l in links:
        atom = editable.GetAtomWithIdx(l.far)
        atom.SetAtomicNum(8)
        atom.SetIsotope(0)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(0)
        dropped |= _far_side(graph, l.far, l.center) - {l.far}
    for idx in sorted(dropped, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    kept = [i for i in range(mol.GetNumAtoms()) if i not in dropped]
    match_cip_labels(mol, acid, dict(enumerate(kept)))
    name = f"{acyl_name(smiles_to_iupac(Chem.MolToSmiles(acid)))} {_class_words(mol, links)}"
    hydrogen_bearing = {
        (a.GetIdx() if a.GetAtomicNum() != 1 else a.GetNeighbors()[0].GetIdx())
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 1 or a.GetTotalNumHs()
    }
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8) == 1:
        name = re.sub(r"^\(O-(\d+O)\)", r"(\1)", name)
    found = _LEADING_HYDROGEN_NUCLIDES.match(name)
    only_hydrogen = all(a.GetAtomicNum() == 1 for a in mol.GetAtoms() if a.GetIsotope())
    if found and only_hydrogen and len(hydrogen_bearing) == 1:
        name = f"({found.group(1)}{found.group(2)})" + name[found.end():]
    return name if nuclides_cited(mol, name) else None


def _labelled_halocarbonic_amide(mol, amide):
    """'carbono(36Cl)chloridic amide', 'N-methyl(13C)carbonochloridic amide': the nuclide of the halogen stands in front of
    the infix that names it and that of the carbon in front of the acid name (P-82.2.1); None when another atom is modified."""
    from .core import nuclides_cited

    halogens = [a for a in mol.GetAtoms() if is_halide_ligand(a) and a.GetNeighbors()[0].GetAtomicNum() == 6]
    modified = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsotope()}
    for halogen in halogens:
        carbon = halogen.GetNeighbors()[0]
        if modified == {i for i in (halogen.GetIdx(), carbon.GetIdx()) if mol.GetAtomWithIdx(i).GetIsotope()} and halogen.GetIsotope():
            break
    else:
        return None
    infix = re.search(r"carbono(\w+?id)ic amide$", amide)
    carbon_nuclide = f"({carbon.GetIsotope()}C)" if carbon.GetIsotope() else ""
    name = amide[: infix.start()] + carbon_nuclide + f"carbono({halogen.GetIsotope()}{halogen.GetSymbol()}){infix.group(1)}ic amide"
    return name if nuclides_cited(mol, name) else None


_OXOACID_STEMS = {5: ("boron", "borin"), 15: ("phosphon", "phosphin"), 33: ("arson", "arsin"), 51: ("stibon", "stibin")}


def _acyloxy_oxoacid(mol, acid):
    """'(acetyloxy)phosphonic acid', '(propanoyloxy)boronic acid' (P-67.1.3.3): an acidic anhydride of a phosphorus,
    boron, arsenic or antimony oxoacid is named with the oxoacid as parent and the acyl-oxy group as prefix."""
    from ._common import group_substituents
    from ._substituents import format_substituent_prefixes

    if len(set(acid)) != 1:
        raise UnsupportedStructure("an anhydride of several inorganic acid groups is not named by a single parent")
    center = mol.GetAtomWithIdx(acid[0])
    oxo = [n for n in center.GetNeighbors() if _bond(mol, center.GetIdx(), n.GetIdx()) == 2.0]
    hydroxy = [n for n in center.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1]
    taken = {n.GetIdx() for n in oxo + hydroxy}
    others = [n for n in center.GetNeighbors() if n.GetIdx() not in taken]
    boron = center.GetAtomicNum() == 5
    if len(hydroxy) not in (1, 2) or len(others) != 3 - len(hydroxy) or len(oxo) > (0 if boron else 1):
        raise UnsupportedStructure("this oxoacid anhydride is not named by a prefix on the acid parent")
    for n in others:
        if n.GetAtomicNum() == 8 and not any(m.GetIdx() != center.GetIdx() and _center(mol, m.GetIdx()) for m in n.GetNeighbors()):
            raise UnsupportedStructure("an ester group on an acidic anhydride of an oxoacid is not named by a single parent")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    entries = {1: [name_branch(graph, n.GetIdx(), center.GetIdx(), halogens, mol=mol) for n in others]}
    ending = "ic acid" if boron or oxo else "ous acid"
    parent = _OXOACID_STEMS[center.GetAtomicNum()][2 - len(hydroxy)] + ending
    return format_substituent_prefixes(group_substituents(entries), omit_locants=True) + parent


_ORGANIC_ELEMENTS = {1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 31, 32, 33, 34, 35, 49, 50, 51, 52, 53, 81, 82, 83}


def _free_anion(atom):
    """An anionic oxygen, sulfur or nitrogen that is not part of a zwitterionic group such as nitro or an N-oxide."""
    return (
        atom.GetFormalCharge() == -1
        and atom.GetAtomicNum() in (7, 8, 16)
        and not any(n.GetFormalCharge() > 0 for n in atom.GetNeighbors())
    )


def _free_cation(atom):
    """A cationic centre that is not part of a zwitterionic group such as nitro or an N-oxide."""
    return atom.GetFormalCharge() > 0 and not any(n.GetFormalCharge() < 0 for n in atom.GetNeighbors())


def _nitrogen_oxide_centre(atom):
    """The positive nitrogen of an amine or imine oxide, a zwitterionic group; a nitro group has a doubly bonded oxygen."""
    return (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == 1
        and any(n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 and n.GetDegree() == 1 for n in atom.GetNeighbors())
        and not any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8 for b in atom.GetBonds()
        )
    )


def name_acid_derivative(mol):
    """Name of an ester, anhydride or acyl halide with no free acid group, else raises."""
    if any(atom.GetAtomicNum() not in _ORGANIC_ELEMENTS for atom in mol.GetAtoms()):
        raise UnsupportedStructure("elements outside the organic acid derivative engine")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments are named as salts or adducts")
    links, acid = find_links(mol)
    if not links:
        raise UnsupportedStructure("no acid derivative group")
    if not acid and any(_free_anion(a) for a in mol.GetAtoms()):
        raise UnsupportedStructure("an anionic group outranks the ester, which is then a prefix (P-41)")
    if not acid and any(_free_cation(a) for a in mol.GetAtoms()):
        raise UnsupportedStructure("a cationic centre outranks the ester, which is then a prefix (P-41)")
    if not acid and any(_nitrogen_oxide_centre(a) for a in mol.GetAtoms()):
        raise UnsupportedStructure("an amine or imine oxide outranks the ester, which is then a prefix (P-41, P-62.5)")
    from ._heteroacyclic import name_heteroacyclic_ester

    if not acid:
        acyloxy = _polyester_with_acyloxy_groups(mol, links)
        if acyloxy is not None:
            return acyloxy
    skeletal = name_heteroacyclic_ester(mol)
    if skeletal is not None:
        return skeletal
    from ._appendix3_skeletons import name_appendix3_skeleton

    if name_appendix3_skeleton(mol) is not None:
        raise UnsupportedStructure("an acid derivative on an Appendix 3 retained parent is named on that parent")
    if acid:
        esters = [l for l in links if l.kind == "ester"]
        if all(_center(mol, c)[0] in ("carbonic", "inorganic") for c in acid) and esters and all(
            _center(mol, l.center)[0] in ("carbonic", "inorganic") for l in esters
        ) and all(l.kind in ("ester", "halide") for l in links):
            return name_ester(mol, links)
        if (
            not any(a.GetFormalCharge() for a in mol.GetAtoms())
            and all(_center(mol, c)[0] == "inorganic" and c in {l.center for l in esters} for c in acid)
            and any(_ACID_CLASS_RANK.get(_center(mol, l.center)[0], 3) == 0 for l in esters)
        ):
            try:
                return name_ester(mol, links)
            except UnsupportedStructure:
                pass
        if all(_center(mol, c)[0] == "carbonic" for c in acid) and all(l.kind == "anhydride" for l in links):
            return name_anhydride(mol, links)
        if all(_center(mol, c)[0] == "inorganic" for c in acid) and any(l.kind == "anhydride" for l in links):
            if any("cyanic" in (_center(mol, l.center)[0], (_center(mol, l.far) or ("inorganic",))[0]) for l in links if l.kind == "anhydride"):
                # P-67.1.2.4.1.3: an -OCN group on the acid is cited by functional replacement, the acid outranking the anhydride
                return _acyloxy_oxoacid(mol, acid)
            try:
                return name_anhydride(mol, links)
            except UnsupportedStructure:
                return _acyloxy_oxoacid(mol, acid)
        from ._multiplicative import name_if_multiplicative
        from ._polyfunctional import name_polyfunctional

        failure = None
        try:
            name = name_if_multiplicative(mol) or name_polyfunctional(mol)
        except UnsupportedStructure as error:
            name, failure = None, error
        if name is not None:
            return name
        # an ester of a senior acid (P-41 class 7a, 7b) is the parent beside the hydrogen esters of a noncarbon acid, but
        # a free acid group, an anion or a cation outranks every ester
        hydrogen_esters = {l.center for l in esters}
        if (
            esters
            and all(c in hydrogen_esters and _center(mol, c)[0] in ("carbonic", "inorganic") for c in acid)
            and all(l.kind in ("ester", "halide") for l in links)
            and not any(a.GetFormalCharge() for a in mol.GetAtoms())
        ):
            return name_ester(mol, links)
        if failure is not None:
            raise failure
        raise UnsupportedStructure("an acid group outranks the esters and no supported parent carries them")
    kinds = {l.kind for l in links}
    if "anhydride" in kinds:
        return name_anhydride(mol, links)
    if "ester" in kinds:
        return _stereogenic_phosphorus_ester(mol, name_ester(mol, links))
    return name_acyl_halide(mol, links)


def name_sulfur_acid_ester(mol):
    """The ester name of a sulfonic or sulfinic acid ester, whose labelled oxygens are cited as those of the acid or of the
    ester linkage (P-82.6.4)."""
    links, _ = find_links(mol)
    esters = [l for l in links if l.kind == "ester"]
    if not esters or any(mol.GetAtomWithIdx(l.center).GetAtomicNum() != 16 for l in esters):
        raise UnsupportedStructure("no sulfonic or sulfinic acid ester")
    return name_acid_derivative(mol)


def name_labelled_anhydride(mol):
    """The anhydride name with the nuclides of its acid components and bridge cited, for a labelled structure whose
    other acid derivatives keep their own part-wise names."""
    links, _ = find_links(mol)
    if not any(l.kind == "anhydride" for l in links):
        raise UnsupportedStructure("no anhydride bridge")
    return name_acid_derivative(mol)


def _stereogenic_phosphorus_ester(mol, name):
    """P-93.2.4, P-93.3.4.1: the descriptor of a stereogenic phosphorus or sulfur centre precedes the anion, bracketed for
    phosphorus: 'methyl (S)-[methyl(phenyl)phosphinate]', 'ethyl (R)-4-nitrobenzene-1-sulfinate'."""
    from ._common import specified_stereo_elements, with_unlocated_descriptor

    centres = [
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (15, 16) and a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ]
    if len(centres) != 1 or centres[0] not in {l.center for l in find_links(mol)[0]}:
        return name
    code = next((c for kind, idx, c in specified_stereo_elements(mol) or () if kind == "atom" and idx == centres[0]), None)
    if code is None:
        return name
    cations, _, anion = name.rpartition(" ")
    centre = mol.GetAtomWithIdx(centres[0])
    if centre.GetAtomicNum() == 16 and centre.GetDegree() == 3:
        return f"{cations} {with_unlocated_descriptor(anion, code)}"
    return f"{cations} ({code})-{enclose(anion)}"
