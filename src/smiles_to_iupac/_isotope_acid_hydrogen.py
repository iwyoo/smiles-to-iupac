"""A nuclide on the O-H, S-H or N-H of an oxoacid or of the hydroxy, sulfanyl or amino group of a silane, germane, stannane or
alumane (P-82.2.1, P-82.6.2): the compound is named without it and the descriptor, with the letter locant of the modified
atom, is inserted before the parent, after its prefixes: 'methyl(O-2H)phosphonic acid'; the hydrogen of the centre itself is cited
with it; the hydrogen of a partial ester is the word 'hydrogen' that carries it, the modified one cited first (P-82.2.2.1):
'methyl (2H)hydrogen hydrogen phosphate'."""

import re
from functools import lru_cache

from rdkit import Chem

from ._common import UnsupportedStructure
from ._isotope_labels import restore_labels, split_isotopes

_ACID_PARENT = re.compile(
    r"(?:(?:phosph|ars|stib)[a-z]*(?:ic|ous)|sulfur\w*(?:ic|ous)|sulfuro\(peroxoic\)|sulfamic|bor(?:on|in)?ic|carbamic|carbonic|oxalic|cyanic"
    r"|nitr(?:ic|ous)|(?:(?:hypo|per)?(?:chlor|brom|iod))(?:ic|ous)|tellur(?:ic|ous)|selen(?:ic|ous)|silicic) acid$|cyanamide$"
    r"|(?:sil|germ|stann|plumb|alum)ane?(?:di|tri|tetra)?(?:thiol|amine|ol)$"
)
_CENTRE_WITHOUT_OXO = {5, 7, 13, 14, 15, 17, 32, 33, 34, 35, 50, 51, 52, 53}
_POLY_STEM = re.compile(r"(?<![a-z])(?:di|tri|tetra)$")
_ACID_HYDROGEN_WORD = re.compile(r"(?<![a-z])(di)?hydrogen(?= )")


def _is_oxoacid_centre(atom):
    """Whether a hydroxy, sulfanyl or amino group on `atom` is cited by an acid or hydride parent: a Group 13 to 17 atom
    other than carbon and sulfur, or a carbon or sulfur atom with a terminal =O (carbon also with a triple bond to nitrogen)."""
    z = atom.GetAtomicNum()
    if z in _CENTRE_WITHOUT_OXO:
        return True
    if z not in (6, 16):
        return False
    return any(
        (b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8 and b.GetOtherAtom(atom).GetDegree() == 1)
        or (z == 6 and b.GetBondTypeAsDouble() == 3.0 and b.GetOtherAtom(atom).GetAtomicNum() == 7)
        for b in atom.GetBonds()
    )


def _group_centre(atom):
    """The one centre that holds the hydroxy, sulfanyl or amino group at `atom`, through the oxygen of a peroxy group
    too; None when it has none or several."""
    centres = set()
    for n in atom.GetNeighbors():
        if _is_oxoacid_centre(n):
            centres.add(n.GetIdx())
        elif atom.GetAtomicNum() == 8 and n.GetAtomicNum() == 8 and n.GetDegree() == 2 and n.GetTotalNumHs() == 0:
            centres.update(m.GetIdx() for m in n.GetNeighbors() if m.GetIdx() != atom.GetIdx() and _is_oxoacid_centre(m))
    return next(iter(centres)) if len(centres) == 1 else None


def _is_carbamic_or_cyanamide_centre(carbon):
    """A carbon that bears the N-H of carbamic acid (a hydroxy group beside =O) or of cyanamide (C#N); an amide is no acid."""
    if carbon.GetAtomicNum() != 6:
        return False
    return any(
        (b.GetBondTypeAsDouble() == 3.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 7)
        or (b.GetBondTypeAsDouble() == 1.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 8 and b.GetOtherAtom(carbon).GetTotalNumHs() > 0)
        for b in carbon.GetBonds()
    )


def _labelled_group_atom(mol):
    """(mol without the nuclide of the acid hydrogen, its atom, the label entry, the atoms with the same label, the label
    entries of the hydrogens of the one centre) for a structure whose only acid hydrogens modified are one nuclide on
    O-H or S-H groups (or one N-H); the other nuclides stay where they are."""
    split = split_isotopes(mol)
    if split is None or len(Chem.GetMolFrags(mol)) != 1:
        return None
    clean, labels, _ = split
    acid = {}
    for index, entry in labels.items():
        atom = clean.GetAtomWithIdx(index)
        if entry["skeleton"] or len(entry["H"]) != 1 or atom.GetAtomicNum() not in (7, 8, 16):
            continue
        if sum(entry["H"].values()) != 1 and atom.GetAtomicNum() != 7:
            continue
        centre = _group_centre(atom)
        if centre is None:
            continue
        centre_atom = clean.GetAtomWithIdx(centre)
        if centre_atom.GetAtomicNum() == 6 and (atom.GetAtomicNum() != 8 and not _is_carbamic_or_cyanamide_centre(centre_atom)):
            continue
        acid[index] = entry
    if not acid or len({tuple(entry["H"].items()) for entry in acid.values()}) != 1:
        return None
    if len({clean.GetAtomWithIdx(i).GetAtomicNum() for i in acid}) != 1:
        return None
    first = min(acid)
    centres = {_group_centre(clean.GetAtomWithIdx(i)) for i in acid}
    centre_labels = {}
    if len(centres) == 1:
        (centre,) = centres
        entry = labels.get(centre)
        if entry and not entry["skeleton"] and clean.GetAtomWithIdx(centre).GetAtomicNum() != 6 and len(_centre_chain(clean, centre)) == 1:
            centre_labels[centre] = entry
    editable = Chem.RWMol(clean)
    restore_labels(editable, clean, labels, skip=set(acid) | set(centre_labels))
    return editable.GetMol(), clean.GetAtomWithIdx(first), acid[first], sorted(acid), centre_labels


@lru_cache(maxsize=256)
def _unlabelled_name(smiles):
    from .core import smiles_to_iupac

    try:
        return smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None


def _named_group(mol):
    """The `_labelled_group_atom` result with the name of the structure without the acid hydrogen label, or None when
    that name holds neither an acid parent nor the word 'hydrogen' that can carry the label."""
    found = _labelled_group_atom(mol)
    if found is None:
        return None
    name = _unlabelled_name(Chem.MolToSmiles(found[0]))
    if name is None or not (_ACID_PARENT.search(name) or _ACID_HYDROGEN_WORD.search(name)):
        return None
    return found, name


def has_isotope_acid_hydrogen_shape(mol) -> bool:
    return _named_group(mol) is not None


def name_isotope_acid_hydrogen(mol) -> str:
    named = _named_group(mol)
    if named is None:
        raise UnsupportedStructure("no labelled hydrogen on an oxoacid group")
    (mol_clean, atom, entry, modified, centre_labels), name = named
    (nuclide,) = entry["H"]
    match = _ACID_PARENT.search(name)
    if match is None:
        return _hydrogen_word_citation(name, nuclide, len(modified), atom)
    stem = _POLY_STEM.search(name[: match.start()])
    groups = {nuclide: [_letter_locants(mol_clean, modified), entry["H"][nuclide] * len(modified), atom.GetTotalNumHs() > 1]}
    for index, centre_entry in centre_labels.items():
        for centre_nuclide, count in centre_entry["H"].items():
            group = groups.setdefault(centre_nuclide, [[], 0, False])
            group[0].insert(0, "1")
            group[1] += count
            group[2] |= mol_clean.GetAtomWithIdx(index).GetTotalNumHs() > 1
    items = [
        f"{','.join(locants)}-{n}{count if count > 1 else '1' if repeatable else ''}"
        for n, (locants, count, repeatable) in sorted(groups.items(), key=lambda item: int(item[0][:-1]))
    ]
    start = stem.start() if stem else match.start()
    return f"{name[:start]}({','.join(items)}){name[start:]}"


def _letter_locants(mol, modified):
    """The letter locants of the modified atoms (P-82.6.2): the element symbol, primed for each further atom on one
    centre; on an acid of several centres joined through oxygen the symbol takes the number of its centre, the centres
    numbered 1, 3, 5 along the chain from the end that gives the modified atoms the lowest locants."""
    centre = {a: _group_centre(mol.GetAtomWithIdx(a)) for a in modified}
    chain = _centre_chain(mol, centre[modified[0]])
    options = [{c: 2 * i + 1 for i, c in enumerate(order)} for order in (chain, chain[::-1])] if len(chain) > 1 else [{}]
    number = min(options, key=lambda numbering: sorted(numbering.get(centre[a], 0) for a in modified))
    seen = {}
    located = []
    for a in sorted(modified, key=lambda a: (number.get(centre[a], 0), a)):
        primes = seen[centre[a]] = seen.get(centre[a], -1) + 1
        located.append(f"{mol.GetAtomWithIdx(a).GetSymbol()}{number.get(centre[a], '')}{"'" * primes}")
    return located


def _centre_chain(mol, first):
    """The acid centres of the same element joined through bridging oxygens to `first`, in order along the chain."""
    element = mol.GetAtomWithIdx(first).GetAtomicNum()

    def joined(c):
        return [
            other.GetIdx()
            for o in mol.GetAtomWithIdx(c).GetNeighbors()
            if o.GetAtomicNum() == 8 and o.GetDegree() == 2
            for other in o.GetNeighbors()
            if other.GetIdx() != c and other.GetAtomicNum() == element and _is_oxoacid_centre(other)
        ]

    component, queue = [first], [first]
    while queue:
        for other in joined(queue.pop()):
            if other not in component:
                component.append(other)
                queue.append(other)
    end = next((c for c in component if len(joined(c)) <= 1), first)
    order = [end]
    while len(order) < len(component):
        order.append(next(n for n in joined(order[-1]) if n not in order))
    return order


def _hydrogen_word_citation(name, nuclide, modified, atom):
    """The ester name with the nuclide in front of the word 'hydrogen' that holds the modified atom; a nuclide cited
    separately from the unmodified hydrogen, which follows it (P-82.2.2.1, P-82.2.4)."""
    word = _ACID_HYDROGEN_WORD.search(name)
    present = 2 if word.group(1) else 1
    if modified > present:
        raise UnsupportedStructure("more modified acid hydrogens than the acid hydrogens of the ester")
    if modified == present:
        citation = f"({nuclide}{modified if modified > 1 else ''})" + ("dihydrogen" if present == 2 else "hydrogen")
    else:
        citation = f"({nuclide})hydrogen hydrogen"
    return name[: word.start()] + citation + name[word.end():]
