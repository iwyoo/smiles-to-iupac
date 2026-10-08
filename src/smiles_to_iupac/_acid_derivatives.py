"""Functional class names of esters, pseudoesters, anhydrides and acyl halides (P-65.5-P-65.7): the molecule is cut at
the junction bonds, each acid piece is named as the free acid by the substitutive engine and turned into an anion,
acid-word or acyl name, and organyl pieces are cited as substituent groups (which also expresses polyester prefixes).
"""

import itertools
import re
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._acid_lexicon import CHALCOGENS
from ._carbonic_family import pseudohalide_at
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._substituents import BRANCH_STEREO, name_branch

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
    if z in _HALIDES and atom.GetDegree() == 1:
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
    if len(oxo) + len(rest) + len(direct) != atom.GetDegree() or any(n.GetAtomicNum() not in (6, 8) for n in rest):
        return None
    positions = [(n.GetIdx(), "Y") for n in rest if n.GetAtomicNum() == 8]
    if not positions or atom.GetAtomicNum() == 5 and oxo:
        return None
    return "inorganic", ["O"] * len(oxo), positions


def _acylated_oxoacid_center(mol, atom):
    """A sulfur, selenium or tellurium oxoacid centre (no organyl neighbour) with an acyloxy neighbour."""
    if atom.GetFormalCharge() or atom.GetIsotope() or atom.IsInRing():
        return False
    neighbors = atom.GetNeighbors()
    if any(n.GetAtomicNum() == 6 for n in neighbors):
        return False
    return any(
        n.GetAtomicNum() == 8
        and any(c.GetIdx() != atom.GetIdx() and c.GetAtomicNum() == 6 and _oxo_slots(mol, c.GetIdx()) for c in n.GetNeighbors())
        for n in neighbors
    )


def _center(mol, idx):
    """(kind, oxo slots, [(neighbor, role)] of the Y positions) for an acyl-type ('acyl'), carbonic
    ('carbonic') or cyanic ('cyanic') centre, else None."""
    atom = mol.GetAtomWithIdx(idx)
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
    organics = [r for r in roles if r[1] == "C"]
    positions = [r for r in roles if r[1] != "C"]
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
            if _center(mol, far.GetIdx()) is not None:
                if found[0] == "carbonic" and _center(mol, far.GetIdx())[0] == "carbonic" and len(chain) == 1:
                    continue
                if found[0] == "inorganic" and _center(mol, far.GetIdx())[0] == "inorganic":
                    continue
                links.append(Link("anhydride", idx, tuple(chain), far.GetIdx()))
            elif far.GetAtomicNum() == 6 or (
                far.GetAtomicNum() in _PSEUDOESTER_ELEMENTS
                and (far.GetAtomicNum() not in _RING_ONLY_PSEUDOESTER or far.IsInRing())
            ):
                links.append(Link("ester", idx, tuple(chain), far.GetIdx()))
    return links, acid


def _acid_name(mol, keep, caps):
    """Name of the free acid made of the atoms `keep` of `mol` with a hydroxy group on every atom of `caps`."""
    editable = Chem.RWMol(mol)
    for idx in caps:
        atom = editable.GetAtomWithIdx(idx)
        atom.SetNumExplicitHs(1)
        atom.SetNoImplicit(True)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    for atom in acid.GetAtoms():
        if atom.GetAtomicNum() in (15, 16):
            # the descriptor of a stereogenic acid centre is cited in front of the anion, not inside the acid name
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return _free_acid_name(acid)


def _free_acid_name(acid):
    from .core import smiles_to_iupac

    smiles = Chem.MolToSmiles(acid)
    if smiles in _FREE_INORGANIC_ACIDS:
        return _FREE_INORGANIC_ACIDS[smiles]
    heavy = [a for a in acid.GetAtoms() if a.GetAtomicNum() in _HEAVY_PNICTOGEN_STEMS]
    if heavy:
        return _heavy_pnictogen_acid_name(acid, heavy)
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

    phosphorus = [
        a
        for a in acid.GetAtoms()
        if a.GetAtomicNum() == 15 and any(n.GetAtomicNum() == 8 and n.GetDegree() == 1 for n in a.GetNeighbors())
    ]
    if len(phosphorus) != 1 or phosphorus[0].GetDegree() + phosphorus[0].GetTotalNumHs() != 3:
        return None
    center = phosphorus[0]
    hydroxy = [n for n in center.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
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


_ACID_ENDING = re.compile(r"(?: (?:[A-Za-z]+(?:,[A-Za-z]+)*)-acid| acid)(\)?)$")


def anion_name(acid_name):
    """P-65.6.1: the '-ic acid' ending becomes '-ate' (italic acid letters dropped)."""
    from ._acid_lexicon import _LETTERED_TAIL

    lettered = _LETTERED_TAIL.search(acid_name)
    if lettered:
        return acid_name[: lettered.start()] + f"({lettered.group(2)}-{lettered.group(1)}ate)"
    match = _ACID_ENDING.search(acid_name)
    if match is None:
        raise UnsupportedStructure("the acid part has no 'acid' ending to turn into an anion name")
    stem = acid_name[: match.start()]
    close = match.group(1)
    if stem.endswith("phosphoric"):
        return stem[: -len("phosphoric")] + "phosphate" + close
    if stem.endswith("ic"):
        return stem[:-2] + "ate" + close
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
    if len(elements) < 2:
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
    rdCIPLabeler.AssignCIPLabels(probe)
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


def _reject_unsupported(mol):
    for atom in mol.GetAtoms():
        if atom.GetIsotope() or atom.GetNumRadicalElectrons() or atom.GetFormalCharge() and atom.GetAtomicNum() not in (7, 8) and not _isocyanide_carbon(atom):
            raise UnsupportedStructure("isotopes, radicals and charges are not supported by the acid derivative engine")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments are named as salts or adducts")


def _multiplied(entries):
    """Citation of the organyl groups in alphanumerical order with multipliers and element letters;
    entries are (name, compound, count, letters) tuples."""
    parts = []
    for name, compound, count, letters in sorted(entries, key=lambda e: alpha_sort_key(e[0])):
        if count > 1:
            multiplier = multiplying_prefix(count, compound=compound and not _plain(name))
            text = f"{multiplier}{enclose(name)}" if compound else f"{multiplier}{name}"
        else:
            text = enclose(name) if compound and any(letters) else name
        parts.append((f"{','.join(sorted(letters, key=lambda t: (len(t), t)))}-" if any(letters) else "") + text)
    return " ".join(parts)


def _plain(name):
    from ._substituents import is_plain_stem_prefix

    return is_plain_stem_prefix(name)


def _piece_key(mol, atoms):
    return Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(atoms), isomericSmiles=False)


_DIYL = re.compile(r"^(?P<prefix>.*?)benzene-(?P<locants>[\d,]+)-(?P<count>di|tri)ol$")


def _diyl_name(mol, atoms, attachments):
    """The multivalent organyl name of the piece `atoms` joined through the `attachments` atoms, from its polyol name."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for idx in attachments:
        oxygen = editable.AddAtom(Chem.Atom(8))
        editable.AddBond(idx, oxygen, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    polyol = editable.GetMol()
    Chem.SanitizeMol(polyol)
    name = smiles_to_iupac(Chem.MolToSmiles(polyol))
    match = _DIYL.match(name)
    if match:
        return f"{match.group('prefix')}{match.group('locants')}-phenylene" if match.group("count") == "di" else None
    if name.endswith("diol"):
        return name[: -len("ol")] + "yl"
    if name.endswith("triol"):
        return name[: -len("ol")] + "yl"
    return None


def _pendant_text(mol, graph, links):
    halogens = halogen_substituents(mol)
    counts = {}
    for l in links:
        name = name_branch(graph, l.far, l.chain[-1], halogens, mol=mol)[0]
        counts[name] = counts.get(name, 0) + 1
    return _multiplied([(n, not n.isalpha(), k, []) for n, k in counts.items()])


def _multiplicative_ester(mol, graph, esters, frags, owner, acid_pieces, r_pieces):
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
    if len({_piece_key(mol, frags[c]) for c in centers}) != 1:
        return None
    pendants = [l for ls in acid_pieces.values() for l in ls if owner[l.far] != hub]
    if not pendants or any(len(r_pieces[owner[l.far]]) != 1 for l in pendants):
        return None
    first = acid_pieces[centers[0]]
    acid = _acid_name(mol, frags[centers[0]], [l.chain[-1] for l in first])
    anion = anion_name(acid)
    diyl = _diyl_name(mol, frags[hub], [l.far for l in hub_links])
    if diyl is None:
        return None
    multiplier = multiplying_prefix(len(centers), compound=not anion.isalpha())
    anion_text = multiplier + anion if anion.isalpha() else f"{multiplier}{enclose(anion)}"
    return " ".join(part for part in (_pendant_text(mol, graph, pendants), diyl, anion_text) if part)


def _bridged_multiplicative_ester(mol, graph, frags, owner, acid_pieces, r_pieces, hubs):
    """'dimethyl butanedioylbis[oxy(2,1-phenylene)] dibutanedioate' (P-13.6.2): two identical terminal acids joined
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
        arm_texts.append("oxy" + enclose(part.text))
    if arm_texts[0] != arm_texts[1]:
        return None
    acyl = acyl_name(_acid_name(mol, frags[inner], [l.chain[-1] for l in acid_pieces[inner]]))
    anion = anion_name(_acid_name(mol, frags[terminals[0]], [l.chain[-1] for l in acid_pieces[terminals[0]]]))
    multiplier = multiplying_prefix(2, compound=not anion.isalpha())
    anion_text = multiplier + anion if anion.isalpha() else f"{multiplier}{enclose(anion)}"
    linker = f"{acyl}bis{enclose(arm_texts[0])}"
    return " ".join(part for part in (_pendant_text(mol, graph, pendants), linker, anion_text) if part)


def _suffix_links(mol, links, keep, cap_of):
    """The `links` whose acid centre the polyfunctional engine cites as a suffix group when every link of `links`
    is turned into a free acid (P-44.1.1): when one parent cannot carry them all, the other links stay as ester or
    halide groups and are cited as prefixes of that acid. All of `links` when the parent cannot be determined."""
    from ._polyfunctional import _group_of, _principal_class, _ring_occurrences, _select

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
    try:
        Chem.SanitizeMol(acid)
        position_of = _select(acid)[2][4]
    except Exception:
        return list(links)
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


def name_ester(mol, links):
    _reject_unsupported(mol)
    kinds = {_center(mol, l.center)[0] == "inorganic" for l in links}
    if len(kinds) > 1 or kinds == {True} and any(a.GetFormalCharge() for a in mol.GetAtoms()):
        raise UnsupportedStructure("the seniority of this phosphorus or boron ester among the other groups is not decided here")
    esters = [l for l in links if l.kind == "ester"]
    graph = adjacency(mol)
    cuts = [(l.chain[-1], l.far) for l in esters]
    frags, owner = _fragments(mol, cuts)
    acid_pieces = {}
    for l in esters:
        acid_pieces.setdefault(owner[l.center], []).append(l)
    r_pieces = {}
    for l in esters:
        r_pieces.setdefault(owner[l.far], []).append(l)
    if len(acid_pieces) > 1:
        multiplicative = _multiplicative_ester(mol, graph, esters, frags, owner, acid_pieces, r_pieces)
        if multiplicative is not None:
            return multiplicative
    if len(acid_pieces) > 1 and len(r_pieces) == 1:
        raise UnsupportedStructure("an organyl component with several acid components uses a multiplied name")
    def _carboxylic(i):
        return any(mol.GetAtomWithIdx(l.center).GetAtomicNum() == 6 and (_center(mol, l.center) or ("",))[0] == "acyl" for l in acid_pieces[i])

    principal = max(acid_pieces, key=lambda i: (_carboxylic(i), len(acid_pieces[i]), len(frags[i])))
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
    acid = _acid_name(mol, acid_atoms, [l.chain[-1] for l in chosen])
    anion = anion_name(acid)
    context = _stereo_context(mol, removed)
    halogens = halogen_substituents(mol)
    named = {}
    per_link = []
    token = BRANCH_STEREO.set(context)
    try:
        for l in chosen:
            name, compound = name_branch(graph, l.far, l.chain[-1], halogens, mol=mol)
            entry = named.setdefault(name, [0, compound, []])
            entry[0] += 1
            entry[2].append(_letters(mol, l))
            per_link.append((l, name))
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of the organyl part is not cited by any supported name")
    free_centers = [
        a
        for a in acid_atoms
        for n in mol.GetAtomWithIdx(a).GetNeighbors()
        if n.GetAtomicNum() in _SYMBOL and n.GetDegree() == 1 and n.GetTotalNumHs() == 1 and _center(mol, a) is not None
    ]
    free = len(free_centers)
    cited = _ester_locants(mol, acid_atoms, per_link, free_centers) if (len(named) > 1 or free) and len(chosen) + free > 1 else None
    if cited is None:
        cited = _diphosphate_locants(mol, acid_atoms, per_link, free_centers)
    entries = [
        (name, compound, count, cited[name] if cited else letters) for name, (count, compound, letters) in named.items()
    ]
    hydrogen = "" if free == 0 else multiplying_prefix(free) + "hydrogen " if free > 1 else "hydrogen "
    return _multiplied(entries) + " " + hydrogen + anion


def _ester_locants(mol, acid_atoms, per_link, free_centers):
    """{organyl name: acid locants} when the acid positions are not equivalent (P-65.6.3.2.2: '16-ethyl 18-methyl
    yohimban-16,18-dicarboxylate'), else None."""
    ends = {a for l, _ in per_link for a in l.chain}
    core = [
        a
        for a in acid_atoms
        if a not in ends
        and not (
            mol.GetAtomWithIdx(a).GetAtomicNum() in _SYMBOL
            and mol.GetAtomWithIdx(a).GetDegree() == 1
            and mol.GetAtomWithIdx(a).GetTotalNumHs() == 1
        )
    ]
    centers = sorted({l.center for l, _ in per_link} | set(free_centers))
    locants, symmetric = _locants_of(mol, core, centers, with_symmetry=True)
    if symmetric:
        return None
    cited = {}
    for l, name in per_link:
        cited.setdefault(name, []).append(str(locants[l.center]))
    return {name: sorted(locs, key=lambda t: (len(t), t)) for name, locs in cited.items()}


def _diphosphate_locants(mol, acid_atoms, per_link, free_centers):
    """{organyl name: locants} for the esters of a diphosphate P1-O2-P3 when its organyl groups and hydrogen atoms can sit
    on the two phosphorus atoms in more than one way, which the bare name cannot tell apart (P-65.6.3.3.2.1: '1,1-diethyl
    3-methyl butane-1,1,3-tricarboxylate'); None when the arrangement is determined."""
    phosphorus = [a for a in acid_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 15]
    if len(phosphorus) != 2:
        return None
    around = [{n.GetIdx() for n in mol.GetAtomWithIdx(p).GetNeighbors()} for p in phosphorus]
    if len(around[0] & around[1]) != 1:
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
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    for atom in acid.GetAtoms():
        if atom.GetAtomicNum() in (15, 16):
            # the descriptor of a stereogenic acid centre is cited in front of the anion, not inside the acid name
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return _free_acid_name(acid)


def _bridged_components(mol, links):
    """({bridge: (centres, chain)}, fragments, owner) for the anhydride bridges of the molecule."""
    bridges = {}
    for l in links:
        if l.kind == "anhydride":
            entry = bridges.setdefault(frozenset(l.chain), [l.chain, l.center, l.far])
    if not bridges:
        raise UnsupportedStructure("no anhydride bridge")
    cuts = []
    for chain, first, second in bridges.values():
        cuts += [(first, chain[0]), (second, chain[-1])]
    frags, owner = _fragments(mol, cuts)
    return bridges, frags, owner


def _locants_of(mol, atoms, centers, with_symmetry=False):
    """{centre: locant} of the acid groups `centers` in the substitutive name of the acid made of `atoms`;
    with `with_symmetry` also whether all the centres are equivalent in that acid (then no locants are needed)."""
    from ._polyfunctional import _group_of, _principal_class, _ring_occurrences, _select

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
    _, _, parts = _select(acid)
    position = {acid.GetAtomWithIdx(i).GetAtomMapNum() - 1: loc for i, loc in parts[4].items()}
    found = {}
    for center in centers:
        if center in position:
            found[center] = position[center]
            continue
        attach = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors() if n.GetIdx() in position]
        if len(attach) != 1:
            raise UnsupportedStructure("the locants of the anhydride groups are not available")
        found[center] = position[attach[0]]
    return (found, False) if with_symmetry else found


def _multiplied_word(word, count):
    if count == 1:
        return word
    return multiplying_prefix(count, compound=not word.isalpha()) + (word if word.isalpha() else enclose(word))


def _component_word(name):
    """The acid word of an anhydride component; a component that is itself an ester is cited whole in parentheses
    (P-67.1.3.3), as in '(methyl dihydrogen phosphate)'."""
    return acid_word(name) if _ACID_ENDING.search(name) else f"({name})"


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


def name_anhydride(mol, links):
    _reject_unsupported(mol)
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
        text = words[0] if words[0] == words[1] else " ".join(sorted(words, key=alpha_sort_key))
        return f"{text} {class_word}"
    hubs = [i for i in components if degrees[i] > 2]
    if hubs:
        if len(hubs) != 1:
            raise UnsupportedStructure("a branched polyanhydride beyond one polybasic centre is not handled")
        hub = hubs[0]
        if any(degrees[i] != 1 for i in components if i != hub):
            hub_keys = {key for key, (_, first, second) in bridges.items() if hub in (owner[first], owner[second])}
            return name_anhydride(
                mol, [l for l in links if l.kind != "anhydride" or frozenset(l.chain) in hub_keys]
            )
        hub_keep = [a for a in frags[hub] if a not in bridge_atoms]
        carbon_hub = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in hub_keep) and all(
            mol.GetAtomWithIdx(c).GetAtomicNum() == 6 for c in centers[hub]
        )
        locants = _locants_of(mol, hub_keep, centers[hub]) if carbon_hub else {}
        groups = {}
        for chain, first, second in bridges.values():
            leaf_center, hub_center = (first, second) if owner[first] != hub else (second, first)
            leaf = owner[leaf_center]
            groups.setdefault(word_of[leaf], []).append(str(locants.get(hub_center, "")))
        if not carbon_hub and len(groups) > 1 and len(set(centers[hub])) > 1:
            raise UnsupportedStructure("different acyl groups on the acid centres of a hub need centre locants")
        cited = []
        for leaf_word in sorted(groups, key=alpha_sort_key):
            group = sorted(groups[leaf_word], key=lambda t: (len(t), t))
            body = _multiplied_word(leaf_word, len(group))
            cited.append(f"{','.join(group)}-{body}" if carbon_hub else body)
        return f"{' '.join(cited)} {word_of[hub]} {class_word}"
    ends = [i for i in components if degrees[i] == 1]
    order = [min(ends, key=lambda i: alpha_sort_key(word_of[i]))]
    while len(order) < len(components):
        order.append(next(n for n in adjacency_of[order[-1]] if n not in order))
    if len(order) == 3:
        w = [word_of[i] for i in order]
        if w[0] == w[2]:
            return f"{_multiplied_word(w[0], 2)} {w[1]} {class_word}"
        return f"{' '.join(w)} {class_word}"
    tail = order[2:]
    tail_atoms = set()
    for i in tail:
        tail_atoms |= set(frags[i])
    join_bridge = next(
        (chain, first, second)
        for chain, first, second in bridges.values()
        if {owner[first], owner[second]} == {order[1], order[2]}
    )
    tail_center = join_bridge[1] if owner[join_bridge[1]] == order[2] else join_bridge[2]
    for chain, first, second in bridges.values():
        if owner[first] in tail and owner[second] in tail:
            tail_atoms |= set(chain)
    merged = acid_word(_component_acid(mol, tail_atoms, [tail_center]))
    return f"{word_of[order[0]]} {word_of[order[1]]} {merged} dianhydride"


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
    counts = {}
    for l in links:
        if l.kind == "halide":
            word = _HALIDES[mol.GetAtomWithIdx(l.far).GetAtomicNum()]
        else:
            word = _PSEUDOHALIDE_WORDS[pseudohalide_at(mol, l.far, l.center)]
        counts[word] = counts.get(word, 0) + 1
    words = []
    for word in sorted(counts, key=alpha_sort_key):
        words.append(multiplying_prefix(counts[word]) + word if counts[word] > 1 else word)
    return " ".join(words)


def _halocarbonic_amide(mol, graph):
    """'carbonochloridic amide' for X-C(=O)-NR2 (P-66.1.1.1.2.2): the amide of the halogenated carbonic acid."""
    from ._substituents import format_substituent_prefixes

    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6 or carbon.GetDegree() != 3 or carbon.IsInRing() or carbon.GetFormalCharge():
            continue
        halide = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() in _HALIDES and n.GetDegree() == 1]
        oxo = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and _bond(mol, carbon.GetIdx(), n.GetIdx()) == 2.0]
        amino = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7 and not n.GetFormalCharge() and not n.IsInRing() and _bond(mol, carbon.GetIdx(), n.GetIdx()) == 1.0]
        if len(halide) != 1 or len(oxo) != 1 or len(amino) != 1:
            continue
        nitrogen = amino[0]
        substituents = [n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
        halogens = halogen_substituents(mol)
        grouped = {}
        for root in substituents:
            name, compound = name_branch(graph, root, nitrogen.GetIdx(), halogens, mol=mol, unsaturated=True)
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append("N")
        covered = {carbon.GetIdx(), halide[0].GetIdx(), oxo[0].GetIdx(), nitrogen.GetIdx()}
        if any(a.GetIdx() not in covered and a.GetAtomicNum() not in (1, 6, 9, 17, 35, 53) for a in mol.GetAtoms()):
            continue
        acid = f"carbono{_HALIDES[halide[0].GetAtomicNum()][:-3]}idic"
        return (format_substituent_prefixes(grouped) if grouped else "") + f"{acid} amide"
    return None


def name_acyl_halide(mol, links):
    _reject_unsupported(mol)
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
    from .core import smiles_to_iupac

    acid_text = smiles_to_iupac(Chem.MolToSmiles(acid))
    if re.fullmatch(r"[a-z0-9,\-]*(?:di|tri|tetra)carbonic acid", acid_text):
        return f"{acid_text[: -len(' acid')]} {_class_words(mol, halides)}"
    return f"{acyl_name(acid_text)} {_class_words(mol, halides)}"


_OXOACID_STEMS = {5: ("boron", "borin"), 15: ("phosphon", "phosphin"), 33: ("arson", "arsin"), 51: ("stibon", "stibin")}


def _acyloxy_oxoacid(mol, acid):
    """'(acetyloxy)phosphonic acid', '(propanoyloxy)boronic acid' (P-67.1.3.3): an acidic anhydride of a phosphorus,
    boron, arsenic or antimony oxoacid is named with the oxoacid as parent and the acyl-oxy group as prefix."""
    from ._common import group_substituents
    from ._substituents import format_substituent_prefixes

    _reject_unsupported(mol)
    if len(set(acid)) != 1:
        raise UnsupportedStructure("an anhydride of several inorganic acid groups is not named by a single parent")
    center = mol.GetAtomWithIdx(acid[0])
    oxo = [n for n in center.GetNeighbors() if _bond(mol, center.GetIdx(), n.GetIdx()) == 2.0]
    hydroxy = [n for n in center.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetTotalNumHs() == 1]
    taken = {n.GetIdx() for n in oxo + hydroxy}
    others = [n for n in center.GetNeighbors() if n.GetIdx() not in taken]
    boron = center.GetAtomicNum() == 5
    if len(hydroxy) not in (1, 2) or len(others) != 3 - len(hydroxy) or len(oxo) != (0 if boron else 1):
        raise UnsupportedStructure("this oxoacid anhydride is not named by a prefix on the acid parent")
    for n in others:
        if n.GetAtomicNum() == 8 and not any(m.GetIdx() != center.GetIdx() and _center(mol, m.GetIdx()) for m in n.GetNeighbors()):
            raise UnsupportedStructure("an ester group on an acidic anhydride of an oxoacid is not named by a single parent")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    entries = {1: [name_branch(graph, n.GetIdx(), center.GetIdx(), halogens, mol=mol) for n in others]}
    parent = _OXOACID_STEMS[center.GetAtomicNum()][2 - len(hydroxy)] + "ic acid"
    return format_substituent_prefixes(group_substituents(entries), omit_locants=True) + parent


_ORGANIC_ELEMENTS = {1, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 31, 32, 33, 34, 35, 49, 50, 51, 52, 53, 81, 82, 83}


def name_acid_derivative(mol):
    """Name of an ester, anhydride or acyl halide with no free acid group, else raises."""
    if any(atom.GetAtomicNum() not in _ORGANIC_ELEMENTS for atom in mol.GetAtoms()):
        raise UnsupportedStructure("elements outside the organic acid derivative engine")
    links, acid = find_links(mol)
    if not links:
        raise UnsupportedStructure("no acid derivative group")
    from ._heteroacyclic import name_heteroacyclic_ester

    skeletal = name_heteroacyclic_ester(mol)
    if skeletal is not None:
        return skeletal
    from ._appendix3_skeletons import name_appendix3_skeleton

    if name_appendix3_skeleton(mol) is not None:
        raise UnsupportedStructure("an acid derivative on an Appendix 3 retained parent is named on that parent")
    if acid:
        if all(_center(mol, c)[0] in ("carbonic", "inorganic") for c in acid) and all(
            l.kind == "ester" and _center(mol, l.center)[0] in ("carbonic", "inorganic") for l in links
        ):
            return name_ester(mol, links)
        if all(_center(mol, c)[0] == "inorganic" for c in acid) and any(l.kind == "anhydride" for l in links):
            return _acyloxy_oxoacid(mol, acid)
        from ._multiplicative import name_if_multiplicative
        from ._polyfunctional import name_polyfunctional

        return name_if_multiplicative(mol) or name_polyfunctional(mol)
    kinds = {l.kind for l in links}
    if "anhydride" in kinds:
        return name_anhydride(mol, links)
    if "ester" in kinds:
        return _stereogenic_phosphorus_ester(mol, name_ester(mol, links))
    return name_acyl_halide(mol, links)


def _stereogenic_phosphorus_ester(mol, name):
    """P-93.2.4, P-93.3.4.1: the descriptor of a stereogenic phosphorus or sulfur centre precedes the anion, bracketed for
    phosphorus: 'methyl (S)-[methyl(phenyl)phosphinate]', 'ethyl (R)-4-nitrobenzene-1-sulfinate'."""
    from ._common import heteroatom_stereo_prefix

    centres = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (15, 16)]
    if len(centres) != 1:
        return name
    descriptor = heteroatom_stereo_prefix(mol, centres[0])
    if descriptor is None:
        return name
    cations, _, anion = name.rpartition(" ")
    if mol.GetAtomWithIdx(centres[0]).GetAtomicNum() == 16:
        return f"{cations} {descriptor}{anion}"
    return f"{cations} {descriptor}[{anion}]"
