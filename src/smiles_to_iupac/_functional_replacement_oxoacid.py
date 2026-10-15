"""Mononuclear oxoacids of P, As, Sb, B and N modified by functional replacement (P-67.1.2.1, P-67.1.2.3, P-67.1.4.1.1.4):
the acid stem follows the number of chalcogen ligands on the centre (one 'phosphin-', two 'phosphon-', three 'phosphor-'),
'ic' for the =X form of the pentavalent centre and 'ous' for the trivalent one, carbon groups replace its hydrogens as
prefixes, the replacing S, Se, Te and NH2 are infixes ('phosphinothious acid', 'boronamidic acid'), and the
hydrogen-bearing elements are cited when the ligands differ ('phosphorothioic O,O,O-acid')."""

from rdkit import Chem

from ._hetero_prefixes import pentavalent_nitrogen_oxo
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, group_substituents, halogen_substituents, is_halide_ligand
from ._numerals import numerical_term
from ._substituents import cited_stereo_around, format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_STEMS = {
    15: ("phosphin", "phosphon", "phosphor"),
    33: ("arsin", "arson", "arsor"),
    51: ("stibin", "stibon", "stibor"),
    5: ("borin", "boron", "bor"),
    7: ("azin", "azon", "azor"),
}
_PENTAVALENT = {7, 15, 33, 51}
_LOCANT = {5: "B", 15: "P", 33: "As", 51: "Sb"}
_VOWELS = "aeiou"
_CHALCOGENS = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_ORDER = {8: 0, 16: 1, 34: 2, 52: 3}
_HALIDE_INFIX = {9: "fluorid", 17: "chlorid", 35: "bromid", 53: "iodid"}
# pseudohalide groups bonded to the acid centre by their first atom, with the infix that replaces an OH (P-67.1.2.4.1.3)
_PSEUDOHALIDES = (
    ("[O;D2]-C#N", "cyanatid"),
    ("[N;D2]=C=O", "isocyanatid"),
    ("[S;D2]-C#N", "thiocyanatid"),
    ("[N;D2]=C=S", "isothiocyanatid"),
    ("[Se;D2]-C#N", "selenocyanatid"),
    ("[C;D2]#N", "cyanid"),
    ("[N;D2]#[C;D1]", "isocyanid"),
    ("[N;D2]=[N+]=[N-]", "azid"),
)
_PSEUDOHALIDE_PATTERNS = [(Chem.MolFromSmarts(smarts), infix) for smarts, infix in _PSEUDOHALIDES]


def _replacer(mol, center, neighbor):
    """(infix, atoms) of a halogen or pseudohalogen group bonded to the acid centre, else None."""
    if neighbor.GetAtomicNum() in _HALIDE_INFIX and is_halide_ligand(neighbor):
        return _HALIDE_INFIX[neighbor.GetAtomicNum()], {neighbor.GetIdx()}
    for pattern, infix in _PSEUDOHALIDE_PATTERNS:
        for match in mol.GetSubstructMatches(pattern):
            if match[0] == neighbor.GetIdx():
                return infix, set(match)
    return None


def _amide_nitrogen(mol, atom, center):
    """An acyclic NH2, NHR or NR2 group of carbon substituents on the acid centre: the amid infix (P-67.1.2.4.1)."""
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != center.GetIdx()]
    return (
        not atom.IsInRing()
        and not atom.GetIsAromatic()
        and not atom.GetFormalCharge()
        and not atom.GetIsotope()
        and all(n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0 for n in others)
        and len(others) + atom.GetTotalNumHs() == 2
    )


def _ring_nitrogen_group(atom):
    """A ring nitrogen bonded to the acid centre is a substituent group (piperidin-1-yl), not an amid infix, and counts as
    an atom of a parent hydride for the acid stem (P-67.1.2.4.2)."""
    return atom.GetAtomicNum() == 7 and atom.IsInRing() and not atom.GetFormalCharge() and not atom.GetIsotope()


def _ring_system(mol, root):
    """The atoms of the ring system that contains `root`."""
    system = set()
    pending = [set(ring) for ring in mol.GetRingInfo().AtomRings() if root in ring]
    while pending:
        ring = pending.pop()
        if ring <= system:
            continue
        system |= ring
        pending += [set(other) for other in mol.GetRingInfo().AtomRings() if system & set(other) and not set(other) <= system]
    return system


def _ligands(mol, center):
    """([=X atoms], [XH atoms], [carbon neighbours], [halide groups], [amide nitrogens]) of the acid centre, or None if it
    has any other neighbour."""
    oxo, hydroxy, carbon, replacers, amides = [], [], [], [], []
    for neighbor in center.GetNeighbors():
        z = neighbor.GetAtomicNum()
        order = mol.GetBondBetweenAtoms(center.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble()
        if center.GetAtomicNum() == 7 and pentavalent_nitrogen_oxo(mol, center.GetIdx(), neighbor.GetIdx()):
            oxo.append(neighbor)
        elif z in _CHALCOGENS and neighbor.GetDegree() == 1 and not neighbor.GetFormalCharge() and not neighbor.GetIsotope():
            if order == 2.0 and neighbor.GetTotalNumHs() == 0:
                oxo.append(neighbor)
            elif order == 1.0 and neighbor.GetTotalNumHs() == 1:
                hydroxy.append(neighbor)
            else:
                return None
        elif z == 6 and order == 1.0:
            carbon.append(neighbor)
        elif order == 1.0 and _replacer(mol, center, neighbor) is not None:
            replacers.append(neighbor)
        elif z == 7 and order == 1.0 and _amide_nitrogen(mol, neighbor, center):
            amides.append(neighbor)
        elif order == 1.0 and _ring_nitrogen_group(neighbor):
            carbon.append(neighbor)
        else:
            return None
    return oxo, hydroxy, carbon, replacers, amides


def _acid_parts(mol):
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _STEMS and a.GetAtomicNum() != 7]
    centers = centers or [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    if len(centers) != 1 or len(Chem.GetMolFrags(mol)) != 1:
        return None
    center = centers[0]
    if (center.GetFormalCharge() and center.GetAtomicNum() != 7) or center.GetIsotope() or center.IsInRing():
        return None
    ligands = _ligands(mol, center)
    if ligands is None:
        return None
    oxo, hydroxy, carbon, replacers, amides = ligands
    z = center.GetAtomicNum()
    pentavalent = z in _PENTAVALENT and len(oxo) == 1
    if len(oxo) > 1 or (len(oxo) == 1 and not pentavalent) or not 1 <= len(hydroxy) <= 3:
        return None
    if len(hydroxy) + len(replacers) + len(amides) + len(carbon) + center.GetTotalNumHs() != 3:
        return None
    if z == 7 and (amides or len(hydroxy) + len(replacers) < (1 if pentavalent else 2)):
        return None
    chalcogens = [*oxo, *hydroxy]
    if all(a.GetAtomicNum() == 8 for a in chalcogens) and pentavalent and not replacers and not amides:
        return None
    ligand_atoms = {x.GetIdx() for x in chalcogens} | {center.GetIdx()} | {n.GetIdx() for n in amides}
    for neighbor in replacers:
        ligand_atoms |= _replacer(mol, center, neighbor)[1]
    ring_hetero = set().union(*(_ring_system(mol, n.GetIdx()) for n in carbon if n.GetAtomicNum() == 7))
    if any(
        a.GetAtomicNum() not in (6, *HALOGEN_PREFIXES) and a.GetIdx() not in ring_hetero
        for a in mol.GetAtoms()
        if a.GetIdx() not in ligand_atoms
    ):
        return None
    return center, oxo, hydroxy, carbon, pentavalent, replacers, amides


def has_functional_replacement_oxoacid_shape(mol) -> bool:
    return _acid_parts(mol) is not None


def name_functional_replacement_oxoacid(mol) -> str:
    parts = _acid_parts(mol)
    if parts is None:
        raise UnsupportedStructure("this is not a mononuclear oxoacid modified by functional replacement")
    center, oxo, hydroxy, carbon, pentavalent, replacers, amides = parts
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    with cited_stereo_around(mol, center.GetIdx()):
        entries = [name_branch(graph, c.GetIdx(), center.GetIdx(), halogens, mol=mol) for c in carbon]
        arms = sorted(
            (
                [name_branch(graph, c.GetIdx(), n.GetIdx(), halogens, mol=mol) for c in n.GetNeighbors() if c.GetIdx() != center.GetIdx()]
                for n in amides
            ),
            key=lambda found: (-len(found), sorted(name for name, _ in found)),
        )

    chalcogens = [*oxo, *hydroxy]
    replaced = sorted((a.GetAtomicNum() for a in chalcogens if a.GetAtomicNum() != 8), key=lambda z: _INFIX[z])
    pieces = []
    for z in sorted(set(replaced), key=lambda z: _INFIX[z]):
        count = replaced.count(z)
        pieces.append(((numerical_term(count) if count > 1 else "") + _INFIX[z], _INFIX[z]))
    infixes = [_replacer(mol, center, n)[0] for n in replacers] + ["amid"] * len(amides)
    for text in sorted(set(infixes)):
        count = infixes.count(text)
        pieces.append(((numerical_term(count) if count > 1 else "") + text, text))
    pieces.sort(key=lambda piece: piece[1])
    z_center = center.GetAtomicNum()
    if z_center == 7 and amides:
        raise UnsupportedStructure("the locants of the amide nitrogen of a nitrogen acid are not defined by the rule")
    stem = _STEMS[z_center][len(hydroxy) + len(replacers) + len(amides) - 1]
    if z_center == 7 and pentavalent and len(hydroxy) + len(replacers) == 3:
        stem = "nitror"
    word = stem
    for text, _ in pieces:
        word += ("o" if text[0] not in _VOWELS and not word.endswith("o") else "") + text
    if z_center == 5 or pentavalent:
        word += "ic acid"
    else:
        word = (word[:-1] if word.endswith("o") and pieces else word) + "ous acid"
    if len({a.GetAtomicNum() for a in chalcogens}) > 1 and (len(chalcogens) >= 3 or pentavalent):
        # the replacing atoms that carry hydrogen, else the oxygens that do (P-67.1.2.4.1)
        cited = [a.GetAtomicNum() for a in hydroxy if a.GetAtomicNum() != 8] or [a.GetAtomicNum() for a in hydroxy]
        word = word[: -len("acid")] + ",".join(_CHALCOGENS[z] for z in sorted(cited, key=lambda z: _ORDER[z])) + "-acid"
    return _substituent_prefix(center, entries, arms) + word


def _substituent_prefix(center, entries, arms):
    """Carbon groups on the centre and on its amid nitrogens, with the letter locants B, P, ... and N, N' once an amide
    nitrogen is substituted (P-67.1.2.4.1)."""
    if not any(arms):
        return format_mononuclear_prefixes(entries)
    located = {}
    for level, found in enumerate(arms):
        located.setdefault("N" + "'" * level, []).extend(found)
    if entries:
        located[_LOCANT[center.GetAtomicNum()]] = entries
    return format_substituent_prefixes(group_substituents(located))
