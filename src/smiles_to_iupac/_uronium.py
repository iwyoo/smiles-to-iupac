"""Guanidinium, uronium and chalcogen uronium cations (P-73.1.2.1, P-73.1.2.2): a carbon bonded to three heteroatoms, two
nitrogens and a nitrogen, oxygen, sulfur, selenium or tellurium atom, one of which carries the positive charge.

The charge is delocalized over the core, so the substituents are located on N, N', N'' (guanidinium) or on N, N' and the
chalcogen symbol (uronium), the nitrogen with more substituents unprimed, then the one holding the alphanumerically first."""

from rdkit import Chem

from ._chalcogenourea import n_substituent_names
from ._common import UnsupportedStructure, adjacency, group_substituents
from ._substituents import alpha_sort_key, format_substituent_prefixes

_CHALCOGEN_WORD = {8: ("O", ""), 16: ("S", "thio"), 34: ("Se", "seleno"), 52: ("Te", "telluro")}


def _core(mol):
    """(carbon, nitrogens, chalcogen or None) of the cation core, or None."""
    if len(Chem.GetMolFrags(mol)) != 1 or sum(a.GetFormalCharge() for a in mol.GetAtoms()) != 1:
        return None
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    for carbon in mol.GetAtoms():
        if carbon.GetAtomicNum() != 6 or carbon.GetDegree() != 3 or carbon.GetFormalCharge() or carbon.IsInRing():
            continue
        neighbors = list(carbon.GetNeighbors())
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        others = [n for n in neighbors if n.GetAtomicNum() != 7]
        if len(nitrogens) == 3 and not others:
            chalcogen = None
        elif len(nitrogens) == 2 and len(others) == 1 and others[0].GetAtomicNum() in _CHALCOGEN_WORD:
            chalcogen = others[0]
        else:
            continue
        core = [*nitrogens, *([chalcogen] if chalcogen is not None else [])]
        if any(n.IsInRing() or n.GetIsAromatic() for n in core):
            continue
        orders = sorted(mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in core)
        charged = [n for n in core if n.GetFormalCharge()]
        if orders != [1.0, 1.0, 2.0] or len(charged) != 1 or charged[0].GetFormalCharge() != 1:
            continue
        if mol.GetBondBetweenAtoms(carbon.GetIdx(), charged[0].GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        return carbon, nitrogens, chalcogen
    return None


def has_uronium_shape(mol) -> bool:
    return _core(mol) is not None


def _order(names_by_atom):
    """Atoms ordered for the locants N, N', N'': more substituents first, then the alphanumerically first substituent."""
    key = lambda item: (-len(item[1]), min((alpha_sort_key(name) for name, _ in item[1]), default=""))
    return [atom for atom, _ in sorted(names_by_atom, key=key)]


def name_uronium(mol) -> str:
    found = _core(mol)
    if found is None:
        raise UnsupportedStructure("not a guanidinium or uronium cation")
    carbon, nitrogens, chalcogen = found
    members = [*nitrogens, *([chalcogen] if chalcogen is not None else [])]
    core_atoms = {carbon.GetIdx(), *(m.GetIdx() for m in members)}
    names = n_substituent_names(mol, adjacency(mol), core_atoms, tuple(m.GetIdx() for m in members), carbon.GetIdx(), junior_groups=True)
    by_atom = dict(zip((m.GetIdx() for m in members), names))
    amino = _order([(n.GetIdx(), by_atom[n.GetIdx()]) for n in nitrogens])
    locants = [("N", "N'", "N''")[i] for i in range(len(amino))]
    positions = {locant: by_atom[atom] for locant, atom in zip(locants, amino)}
    if chalcogen is None:
        parent = "guanidinium"
    else:
        symbol, word = _CHALCOGEN_WORD[chalcogen.GetAtomicNum()]
        positions[symbol] = by_atom[chalcogen.GetIdx()]
        parent = f"{word}uronium"
    grouped = group_substituents({locant: names for locant, names in positions.items() if names})
    return f"{format_substituent_prefixes(grouped)}{parent}"
