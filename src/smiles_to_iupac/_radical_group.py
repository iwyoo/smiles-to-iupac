"""A single radical centre on a skeleton that carries other groups (P-71.2.1.2, P-71.3.4, P-71.7): a radical
outranks every characteristic group (P-41), so the structure is the parent radical -- the substituent group
obtained by removing one hydrogen atom from the parent hydride -- with every other group cited as a prefix."""

import re

from rdkit import Chem

from ._common import YLO_MAP_NUMBER, UnsupportedStructure, adjacency, halogen_substituents, specified_stereo_elements
from ._hetero_prefixes import PEROXY_PREFIXES, is_functional_carbon
from ._numerals import alkyl_name
from ._substituents import name_branch

_MAX_ATOMS = 80
_CHALCOGENS = {8, 16, 34, 52}
_SENIORITY = (7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81, 8, 16, 34, 52, 6)


def _polyradical_centres(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if (
        len(radicals) > 1
        and all(a.GetNumRadicalElectrons() == 1 and a.GetAtomicNum() == 6 and not a.GetIsotope() for a in radicals)
        and not any(a.GetFormalCharge() for a in mol.GetAtoms())
    ):
        return radicals
    return None


def _multi_centres(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if (
        len(radicals) > 1
        and not any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms())
        and all(
            (a.GetAtomicNum() == 6 and a.GetNumRadicalElectrons() <= 3)
            or (a.GetAtomicNum() in (8, 16) and a.GetNumRadicalElectrons() == 1 and a.GetDegree() == 1 and a.GetTotalNumHs() == 0)
            or (a.GetAtomicNum() == 7 and a.GetNumRadicalElectrons() == 1 and a.GetDegree() <= 2 and not a.GetIsAromatic())
            for a in radicals
        )
    ):
        return radicals
    return None


def _centre(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(radicals) > 1 and (_polyradical_centres(mol) or _multi_centres(mol)):
        return radicals[0]
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() not in (1, 2, 3):
        return None
    centre = radicals[0]
    if (centre.GetIsotope() and centre.GetAtomicNum() != 6) or centre.GetFormalCharge():
        return None
    if centre.GetNumRadicalElectrons() > 1 and centre.GetAtomicNum() != 6:
        return None
    if centre.GetAtomicNum() == 6:
        return centre
    if centre.GetAtomicNum() == 7 and centre.GetDegree() <= 2 and not centre.GetIsAromatic() and centre.GetNumRadicalElectrons() == 1:
        return centre
    if (
        centre.GetAtomicNum() in _CHALCOGENS
        and centre.GetDegree() == 1
        and centre.GetTotalNumHs() == 0
        and centre.GetNumRadicalElectrons() == 1
        and centre.GetNeighbors()[0].GetAtomicNum() in (6, 8)
    ):
        return centre
    return None


def _hydride(mol, centre):
    """`mol` with the radical electrons of `centre` capped: one hydrogen for a 'yl' centre, a doubly or triply bonded
    methylene or methine for 'ylidene' and 'ylidyne'; (capped mol, capping atom index)."""
    order = centre.GetNumRadicalElectrons()
    capped = Chem.RWMol(mol)
    atom = capped.GetAtomWithIdx(centre.GetIdx())
    atom.SetNumRadicalElectrons(0)
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)
    cap = capped.AddAtom(Chem.Atom(1 if order == 1 else 6))
    capped.AddBond(centre.GetIdx(), cap, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[order])
    hydride = capped.GetMol()
    Chem.SanitizeMol(hydride)
    return hydride, cap


def has_radical_group_shape(mol) -> bool:
    centre = _centre(mol)
    return centre is not None and (mol.GetNumAtoms() > 1 or bool(centre.GetIsotope()))


def name_radical_group(mol) -> str:
    from ._isotope_labels import split_isotopes
    from ._substituents import ISOTOPE_LABELS

    labels = {}
    split = split_isotopes(mol)
    if split is not None:
        mol, labels, _ = split
    centre = _centre(mol)
    if centre is not None and _polyradical_centres(mol):
        try:
            return _name_polyradical(mol)
        except UnsupportedStructure:
            return _name_parent_radical(mol)
    if centre is not None and _multi_centres(mol):
        return _name_parent_radical(mol)
    if (
        centre is None
        or mol.GetNumAtoms() > _MAX_ATOMS
        or len(Chem.GetMolFrags(mol)) != 1
        or specified_stereo_elements(mol)
    ):
        raise UnsupportedStructure("this radical is not a single radical centre on a plain skeleton")
    if labels and centre.GetAtomicNum() != 6:
        raise UnsupportedStructure("isotopic modification of this radical is not supported yet")
    if centre.GetAtomicNum() == 7 and not centre.IsInRing():
        return _name_aminyl(mol, centre)
    if centre.GetAtomicNum() == 8 and centre.GetNeighbors()[0].GetAtomicNum() == 8:
        return _name_peroxyl(mol, centre)
    hydride, hydrogen = _hydride(mol, centre)
    graph = adjacency(hydride)
    aromatic = frozenset(a.GetIdx() for a in hydride.GetAtoms() if a.GetIsAromatic())
    context = {"labels": labels, "consumed": set(), "mol": hydride}
    token = ISOTOPE_LABELS.set(context if labels else None)
    peroxy = PEROXY_PREFIXES.set(True)
    try:
        name, _ = name_branch(
            graph, centre.GetIdx(), hydrogen, halogen_substituents(hydride), aromatic, mol=hydride, unsaturated=True
        )
    finally:
        PEROXY_PREFIXES.reset(peroxy)
        ISOTOPE_LABELS.reset(token)
    if set(labels) - context["consumed"]:
        raise UnsupportedStructure("an isotopically modified atom of the radical is not cited by any supported name")
    if centre.GetAtomicNum() == 8:
        if not name.endswith("oxy"):
            raise UnsupportedStructure("the oxygen radical has no 'oxy' prefix to turn into 'oxyl' (P-71.3.4)")
        return name + "l"
    return name


def _name_aminyl(mol, centre):
    """Aminyl and amidyl radicals: the amine or amide suffix with '-yl' added, every other group a prefix (P-71.3.2)."""
    from ._polyfunctional import FORCED_PRINCIPAL, _name_labelled

    hydride, _ = _hydride(mol, centre)
    carbons = [n for n in centre.GetNeighbors() if n.GetAtomicNum() == 6]
    if not carbons or any(n.GetAtomicNum() not in (1, 6) for n in centre.GetNeighbors()):
        raise UnsupportedStructure("this nitrogen radical is not an aminyl or amidyl radical")

    def acyl(carbon):
        return any(
            n.GetAtomicNum() == 8 and hydride.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in carbon.GetNeighbors()
        )

    amide = any(acyl(c) for c in carbons)
    hydride = Chem.RemoveHs(hydride)
    token = FORCED_PRINCIPAL.set("amide" if amide else "amine")
    try:
        name = _name_labelled(hydride, {})
    finally:
        FORCED_PRINCIPAL.reset(token)
    if amide and name.endswith("amide"):
        return name[:-1] + "yl"
    if not amide and name.endswith("amine"):
        return name[:-1] + "yl"
    if not amide and name.endswith("aniline"):
        return name[: -len("aniline")] + "benzenaminyl"
    raise UnsupportedStructure("the nitrogen radical is not named as an amine or amide parent")


def _name_peroxyl(mol, centre):
    hydride, hydrogen = _hydride(mol, centre)
    peroxide = centre.GetNeighbors()[0]
    carbons = [n for n in peroxide.GetNeighbors() if n.GetIdx() != centre.GetIdx()]
    if len(carbons) != 1 or carbons[0].GetAtomicNum() != 6:
        raise UnsupportedStructure("this oxygen radical is not a peroxyl radical")
    graph = adjacency(hydride)
    aromatic = frozenset(a.GetIdx() for a in hydride.GetAtoms() if a.GetIsAromatic())
    name, compound = name_branch(
        graph, carbons[0].GetIdx(), peroxide.GetIdx(), halogen_substituents(hydride), aromatic, mol=hydride, unsaturated=True
    )
    return (f"({name})" if compound else name) + "peroxyl"


def _name_polyradical(mol):
    """Several radical carbons (P-71.2.2.2, P-71.2.3): named as the polyanion of the same skeleton with 'ide' read as
    'yl', every other group a prefix."""
    from ._anion import name_anion

    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetIsotope() for a in mol.GetAtoms()) or specified_stereo_elements(mol):
        raise UnsupportedStructure("this polyradical is not a single plain skeleton")
    anionic = Chem.RWMol(mol)
    for atom in _polyradical_centres(mol):
        target = anionic.GetAtomWithIdx(atom.GetIdx())
        target.SetNumRadicalElectrons(0)
        target.SetFormalCharge(-1)
    anionic = anionic.GetMol()
    Chem.SanitizeMol(anionic)
    name = name_anion(anionic)
    for ending in ("diide", "triide", "tetraide"):
        if name.endswith(ending) and "-" in name[: -len(ending)]:
            return name[: -len("ide")] + "yl"
    raise UnsupportedStructure("the polyradical is not named as a polyanion parent")


def _carbenium_valence(atom):
    """A carbenium centre: three bonds in total, three single bonds or one double and one single bond (a vinyl or aryl
    cation, named through the carbanion of the same skeleton)."""
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    if atom.GetTotalNumHs() == 0 and sorted(orders) == [1.0, 2.0] and atom.GetBonds()[0].GetOtherAtom(atom).GetAtomicNum() == 6:
        return all(b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds())
    return atom.GetTotalNumHs() + atom.GetDegree() == 3 and all(order == 1.0 for order in orders)


def _group_cation_centre(mol):
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if (
        len(cations) != 1
        or cations[0].GetFormalCharge() != 1
        or cations[0].GetAtomicNum() != 6
        or cations[0].GetIsotope()
        or not _carbenium_valence(cations[0])
        or any(a.GetNumRadicalElectrons() or a.GetIsotope() for a in mol.GetAtoms())
        or len(Chem.GetMolFrags(mol)) != 1
        or mol.GetNumAtoms() < 2
    ):
        return None
    return cations[0]


def _is_unbranched_alkane_with_terminal_cation(mol, centre) -> bool:
    return (
        centre.GetDegree() == 1
        and all(a.GetAtomicNum() == 6 and not a.IsInRing() and a.GetDegree() <= 2 for a in mol.GetAtoms())
        and all(b.GetBondTypeAsDouble() == 1.0 for b in mol.GetBonds())
    )


def has_group_cation_shape(mol) -> bool:
    return _group_cation_centre(mol) is not None


def name_group_cation(mol) -> str:
    """A carbenium centre on a skeleton with other groups (P-73.2.2.1): the carbanion name of the same skeleton
    with 'ide' read as 'ylium', every other group a prefix."""
    from ._anion import name_anion

    centre = _group_cation_centre(mol)
    if centre is None or specified_stereo_elements(mol):
        raise UnsupportedStructure("this cation is not a single carbenium centre on a plain skeleton")
    if _is_unbranched_alkane_with_terminal_cation(mol, centre):
        # P-73.2.2.1.1: the specific method turns the 'ane' of the parent hydride into 'ylium'
        return alkyl_name(mol.GetNumAtoms()) + "ium"
    anionic = Chem.RWMol(mol)
    anionic.GetAtomWithIdx(centre.GetIdx()).SetFormalCharge(-1)
    anionic = anionic.GetMol()
    Chem.SanitizeMol(anionic)
    name = name_anion(anionic)
    if re.search(r"\d+H\)-ide$", name):
        if centre.GetDegree() == 2:
            raise UnsupportedStructure("an added-hydrogen anion name does not carry over to an aryl cation")
        # P-73.2.2.2: the added hydrogen of a mancude cation is cited as for the anion, 'anthracen-4a(2H)-ylium'
        return name[: -len("ide")] + "ylium"
    if name.endswith("anide"):
        return name[: -len("anide")] + "ylium"
    if name.endswith("benzenide"):
        return name[: -len("ide")] + "ylium"
    if name.endswith("ide") and name[:-3].endswith("-"):
        return name[:-3] + "ylium"
    raise UnsupportedStructure("the cation is not named as a carbanion-like parent")


def _with_ylo(mol, kept):
    """`mol` with every radical centre outside `kept` replaced by dummy halogens named 'ylo' (P-71.5)."""
    rw = Chem.RWMol(mol)
    for atom in mol.GetAtoms():
        electrons = atom.GetNumRadicalElectrons()
        if not electrons or atom.GetIdx() in kept:
            continue
        target = rw.GetAtomWithIdx(atom.GetIdx())
        target.SetNumRadicalElectrons(0)
        target.SetNoImplicit(True)
        target.SetNumExplicitHs(atom.GetTotalNumHs())
        for _ in range(electrons):
            dummy = Chem.Atom(53)
            dummy.SetProp("_named_prefix", "ylo")
            dummy.SetAtomMapNum(YLO_MAP_NUMBER)
            rw.AddBond(atom.GetIdx(), rw.AddAtom(dummy), Chem.BondType.SINGLE)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def _name_parent_radical(mol):
    """Several radical centres that do not fit one parent (P-71.7): the parent holds the most centres, then rings
    outrank chains; the other centres are cited with the nondetachable prefix 'ylo' (P-71.5)."""
    from itertools import combinations

    if (
        len(Chem.GetMolFrags(mol)) != 1
        or mol.GetNumAtoms() > _MAX_ATOMS
        or any(a.GetAtomicNum() == 53 for a in mol.GetAtoms())
        or specified_stereo_elements(mol)
    ):
        raise UnsupportedStructure("these radical centres are not on one plain skeleton")
    radicals = [a.GetIdx() for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    for size in range(len(radicals), 0, -1):
        found = []
        for kept in combinations(radicals, size):
            if any(mol.GetAtomWithIdx(i).GetNumRadicalElectrons() != 1 for i in kept):
                continue
            try:
                marked = _with_ylo(mol, set(kept))
                if size > 1 and not _polyradical_centres(marked):
                    continue
                name = name_radical_group(marked) if size == 1 else _name_polyradical(marked)
            except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
                continue
            if size == 1 and "ylo" not in name and "oxylcarbonyl" not in name and len(radicals) > 1:
                continue
            seniority = tuple(
                -sum(1 for i in kept if mol.GetAtomWithIdx(i).GetAtomicNum() == z) for z in _SENIORITY
            )
            acyl_centres = sum(
                1
                for i in kept
                for n in mol.GetAtomWithIdx(i).GetNeighbors()
                if n.GetAtomicNum() == 6 and is_functional_carbon(mol, n.GetIdx())
            )
            found.append((seniority, -acyl_centres, not all(mol.GetAtomWithIdx(i).IsInRing() for i in kept), name))
        if found:
            return min(found)[-1]
    raise UnsupportedStructure("no parent radical holds the radical centres")
