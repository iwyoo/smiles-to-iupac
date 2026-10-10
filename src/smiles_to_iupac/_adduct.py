"""Adducts and solvates of neutral components (P-14.8, P-77.1.3(1)): the component names are joined by em dashes in
order of the seniority of classes (P-41, organic before inorganic) and the proportions follow as '(m/n)'. Components
of one class are cited alphanumerically."""

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key
from ._seniority import SUFFIX_CLASS_RANK

_NOBLE_GAS_NAMES = {2: "helium", 10: "neon", 18: "argon", 36: "krypton", 54: "xenon", 86: "radon"}
_MAIN_GROUP_NAMES = {
    1: "hydrogen", 5: "boron", 7: "nitrogen", 8: "oxygen", 9: "fluorine", 13: "aluminium", 14: "silicon",
    15: "phosphorus", 16: "sulfur", 17: "chlorine", 31: "gallium", 32: "germanium", 33: "arsenic", 34: "selenium",
    35: "bromine", 49: "indium", 50: "tin", 51: "antimony", 52: "tellurium", 53: "iodine", 81: "thallium",
    82: "lead", 83: "bismuth", 84: "polonium", 85: "astatine",
}

_CLASS_PATTERNS = (
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H1]"),
    ("ester", "[CX3](=O)[OX2][#6]"),
    ("amide", "[CX3](=O)[NX3]"),
    ("nitrile", "[CX2]#N"),
    ("aldehyde", "[CX3H1](=O)[#6]"),
    ("ketone", "[#6][CX3](=O)[#6]"),
    ("alcohol", "[#6][OX2H1]"),
    ("amine", "[NX3;!$(N[#6]=[O,S,N]);!$(N-[a])]"),
)
_JUNIOR_CARBON_PATTERNS = (
    (0.1, "[#6][OX2;!R][#6]"),
    (0.2, "[#6][SX2,SeX2,TeX2,SX3,SeX3,TeX3,SX4,SeX4,TeX4,OX2;!R][#6,O,S,Se,Te]"),
)
_DONOR_ATOMS = {7, 8, 15, 16, 33, 34, 51, 52}
_PARENT_ATOMS = (7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81, 8, 16, 34, 52)
_CARBON_RANK = max(SUFFIX_CLASS_RANK.values()) + 1
_INORGANIC_RANK = _CARBON_RANK + 1
_WATER_RANK = _INORGANIC_RANK + 1


def _is_bare_water(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 2


def _inorganic_name(frag):
    """Name of a bare atom, dihydrogen, hydron or hydride component (P-14.8.2, ref. 12 IR-5.5), else None."""
    from ._coordination import _METAL_NAMES

    atoms = list(frag.GetAtoms())
    if any(a.GetIsotope() for a in atoms):
        return None
    charge = sum(a.GetFormalCharge() for a in atoms)
    if all(a.GetAtomicNum() == 1 for a in atoms):
        hydrogens = len(atoms) + sum(a.GetTotalNumHs() for a in atoms)
        if charge == 1 and hydrogens == 1:
            return "hydron"
        if charge == -1 and hydrogens == 1:
            return "hydride"
        if charge == 0 and hydrogens == 2:
            return "dihydrogen"
        return None
    if len(atoms) != 1 or charge or atoms[0].GetTotalNumHs():
        return None
    z = atoms[0].GetAtomicNum()
    return _NOBLE_GAS_NAMES.get(z) or _METAL_NAMES.get(z) or _MAIN_GROUP_NAMES.get(z)


def _parent_atom_class(frag):
    """Index in P-44.1.2 of the most senior skeletal heteroatom (classes 21-39 of Table 4.1), or None for a carbon
    compound: heteroatoms of ethers, sulfides and nitro or nitroso groups belong to substituents."""
    found = []
    for atom in frag.GetAtoms():
        z = atom.GetAtomicNum()
        if z not in _PARENT_ATOMS:
            continue
        if not atom.IsInRing():
            if z in (8, 16, 34, 52) or (z == 7 and any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors())):
                continue
        found.append(_PARENT_ATOMS.index(z))
    return min(found, default=None)


def _parent_key(frag):
    """Seniority of the parent structure within one class: a ring system before a chain (P-44.1.2.2), ring systems by
    P-44.2.1."""
    from ._ring_system_seniority import general_key

    systems = []
    for ring in frag.GetRingInfo().AtomRings():
        joined = set(ring)
        rest = []
        for system in systems:
            if system & joined:
                joined |= system
            else:
                rest.append(system)
        systems = rest + [joined]
    if not systems:
        return (1,)
    return (0,) + min(general_key(frag, system) for system in systems)


def _class_rank(frag):
    if _is_bare_water(frag):
        return (_WATER_RANK,)
    if not any(a.GetAtomicNum() == 6 for a in frag.GetAtoms()):
        return (_INORGANIC_RANK,)
    ranks = [
        SUFFIX_CLASS_RANK[name]
        for name, smarts in _CLASS_PATTERNS
        if frag.HasSubstructMatch(Chem.MolFromSmarts(smarts))
    ]
    if ranks:
        return (min(ranks),) + _parent_key(frag)
    atom_class = _parent_atom_class(frag)
    if atom_class is not None:
        return (_CARBON_RANK - 1 + atom_class / 100,) + _parent_key(frag)
    junior = min(
        (offset for offset, smarts in _JUNIOR_CARBON_PATTERNS if frag.HasSubstructMatch(Chem.MolFromSmarts(smarts))),
        default=0,
    )
    return (_CARBON_RANK + junior,) + _parent_key(frag)


def _components(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) < 2:
        return None
    charged = [f for f in frags if sum(a.GetFormalCharge() for a in f.GetAtoms())]
    if charged and not (
        all(_inorganic_name(f) in ("hydron", "hydride") for f in charged)
        and any(a.GetAtomicNum() == 6 for f in frags for a in f.GetAtoms())
    ):
        return None
    counts = {}
    for frag in frags:
        key = Chem.MolToSmiles(frag)
        entry = counts.setdefault(key, [frag, 0])
        entry[1] += 1
    return list(counts.items())


def _donor_text(base, donor):
    """(donor symbol or N locant, component name citing its locants) for a donor atom that shares its element with other
    atoms of the component: the nitrogens of a urea or guanidine take the locants N, N', N'' (P-66.1.1.4.1.1,
    P-66.4.1.2.1.1)."""
    from ._guanidine import nitrogen_locants as guanidine_locants
    from ._urea import nitrogen_locants as urea_locants

    if donor.GetAtomicNum() != 7:
        return None
    for locants in (urea_locants, guanidine_locants):
        found = locants(base)
        if found is not None and donor.GetIdx() in found[0]:
            return found[0][donor.GetIdx()], found[1]
    return None


def _attachment(base, acid):
    """(donor text, acceptor symbol, component name or None) for a drawn donor-acceptor pair, cited when the donor
    component has several possible donor atoms (P-68.1.6.2); None otherwise. The donor text is the element symbol when
    that is unambiguous in the component, else the locant of the donor atom."""
    donors = [a for a in base.GetAtoms() if a.HasProp("_adduct_donor")]
    acceptors = [a for a in acid.GetAtoms() if a.HasProp("_adduct_acceptor")]
    if len(donors) != 1 or len(acceptors) != 1:
        return None
    donor, acceptor = donors[0], acceptors[0]
    candidates = [a.GetAtomicNum() for a in base.GetAtoms() if a.GetAtomicNum() in _DONOR_ATOMS]
    if len(candidates) < 2:
        return None
    if [a.GetAtomicNum() for a in acid.GetAtoms()].count(acceptor.GetAtomicNum()) != 1:
        return None
    if candidates.count(donor.GetAtomicNum()) == 1:
        return donor.GetSymbol(), acceptor.GetSymbol(), None
    located = _donor_text(base, donor)
    return None if located is None else (located[0], acceptor.GetSymbol(), located[1])


def has_bare_inorganic_component(mol) -> bool:
    from ._ocene import has_ocene_shape

    index_frags = Chem.GetMolFrags(mol)
    if len(index_frags) < 2 or not any(len(f) <= 2 for f in index_frags):
        return False
    try:
        frags = Chem.GetMolFrags(mol, asMols=True)
    except Chem.rdchem.MolSanitizeException:
        return False
    return any(_inorganic_name(f) for f in frags) and has_adduct_shape(mol) and not has_ocene_shape(mol)


def has_inorganic_component(mol) -> bool:
    return any(not any(a.GetAtomicNum() == 6 for a in f.GetAtoms()) for f in Chem.GetMolFrags(mol, asMols=True))


def has_adduct_shape(mol) -> bool:
    return _components(mol) is not None


def _is_lewis_acid(frag):
    return any(a.HasProp("_adduct_acceptor") for a in frag.GetAtoms())


def name_adduct(mol, namer) -> str:
    named = []
    for smiles, (frag, count) in _components(mol):
        name = "water" if _is_bare_water(frag) else _inorganic_name(frag) or namer(smiles)
        named.append((_class_rank(frag), alpha_sort_key(name), name, count, frag))
    if len(named) == 1:
        raise UnsupportedStructure("identical components form a repeated molecule, not an adduct")
    named.sort(key=lambda item: (item[0] == (_WATER_RANK,), _is_lewis_acid(item[4])) + item[:4])
    proportions = f"({'/'.join(str(item[3]) for item in named)})"
    if len(named) == 2 and all(item[3] == 1 for item in named):
        pair = _attachment(named[0][4], named[1][4])
        if pair is not None:
            return f"{pair[2] or named[0][2]}({pair[0]}\u2014{pair[1]}){named[1][2]} {proportions}"
    joined = "\u2014".join(item[2] for item in named)
    return f"{joined} {proportions}"
