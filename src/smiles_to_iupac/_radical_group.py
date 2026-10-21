"""A single radical centre on a skeleton that carries other groups (P-71.2.1.2, P-71.3.4, P-71.7): a radical
outranks every characteristic group (P-41), so the structure is the parent radical -- the substituent group
obtained by removing one hydrogen atom from the parent hydride -- with every other group cited as a prefix."""

import re

from rdkit import Chem

from ._common import (
    YLO_MAP_NUMBER,
    UnsupportedStructure,
    adjacency,
    assign_cip_labels,
    halogen_substituents,
    match_cip_labels,
    phantom_capped,
    refuse_specified_stereo,
    specified_stereo_elements,
)
from ._hetero_prefixes import (
    DIPOLAR_GROUPS,
    EXTENDED_PREFIXES,
    MONONUCLEAR_HYDRIDES,
    RADICAL_CENTRE,
    RADICAL_PARENT,
    is_functional_carbon,
)
from ._numerals import alkyl_name
from ._radical_hydryl import hydryl_radical_name
from ._radical_poly import _is_acyl, polyradical_name
from ._substituents import name_branch

_MAX_ATOMS = 80
_IMINE = Chem.MolFromSmarts("[CX3]=[NX2]")
_CHALCOGENS = {8, 16, 34, 52}
_NOT_ATE_CENTRE = {1, 6, 7, 8, 9, 16, 17, 34, 35, 52, 53}
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


def _share_a_heteroatom_chain(mol, radicals):
    """Two radical centres of one element joined through atoms of that element alone lie in one parent hydride (P-71.7a)."""
    for atom in radicals:
        z = atom.GetAtomicNum()
        if z == 6:
            continue
        seen, stack = {atom.GetIdx()}, [atom]
        while stack:
            for n in stack.pop().GetNeighbors():
                if n.GetAtomicNum() == z and n.GetIdx() not in seen:
                    seen.add(n.GetIdx())
                    stack.append(n)
        if any(o.GetIdx() in seen and o.GetIdx() != atom.GetIdx() for o in radicals if o.GetAtomicNum() == z):
            return True
    return False


def _multi_centres(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if _share_a_heteroatom_chain(mol, radicals):
        return None
    if (
        len(radicals) > 1
        and not any((a.GetFormalCharge() and a.GetAtomicNum() != 6) or a.GetIsotope() for a in mol.GetAtoms())
        and all(
            (a.GetAtomicNum() == 6 and a.GetNumRadicalElectrons() <= 3)
            or (a.GetAtomicNum() in (8, 16) and a.GetNumRadicalElectrons() == 1 and a.GetDegree() == 1 and a.GetTotalNumHs() == 0)
            or (a.GetAtomicNum() == 7 and a.GetNumRadicalElectrons() == 1 and a.GetDegree() <= 2 and not a.GetIsAromatic())
            or (a.GetAtomicNum() in MONONUCLEAR_HYDRIDES and a.GetAtomicNum() not in _CHALCOGENS and a.GetNumRadicalElectrons() == 1)
            for a in radicals
        )
    ):
        return radicals
    return None


def _centre(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(radicals) > 1:
        ions_fit = not sum(a.GetFormalCharge() for a in mol.GetAtoms()) or all(
            a.GetAtomicNum() == 6 for a in mol.GetAtoms() if a.GetFormalCharge()
        )
        return radicals[0] if not any(a.GetIsotope() for a in mol.GetAtoms()) and ions_fit else None
    if len(radicals) != 1 or radicals[0].GetNumRadicalElectrons() not in (1, 2, 3):
        return None
    centre = radicals[0]
    if (centre.GetIsotope() and centre.GetAtomicNum() != 6) or centre.GetFormalCharge():
        return None
    if centre.GetNumRadicalElectrons() > 1 and centre.GetAtomicNum() not in (6, 7) and (
        centre.GetAtomicNum() not in MONONUCLEAR_HYDRIDES or centre.GetAtomicNum() in _CHALCOGENS
    ):
        return None
    if centre.GetAtomicNum() == 6:
        return centre
    if centre.GetAtomicNum() == 7 and centre.GetDegree() <= 2 and not centre.GetIsAromatic() and centre.GetNumRadicalElectrons() == 1:
        return centre
    if centre.GetAtomicNum() == 7 and centre.GetDegree() == 1 and not centre.GetIsAromatic() and centre.GetNumRadicalElectrons() == 2:
        return centre
    if (
        centre.GetAtomicNum() in MONONUCLEAR_HYDRIDES
        and centre.GetAtomicNum() not in _CHALCOGENS
        and centre.GetDegree() >= 1
        and not centre.GetIsAromatic()
    ):
        return centre
    if (
        centre.GetAtomicNum() in _CHALCOGENS
        and centre.GetDegree() == 1
        and centre.GetTotalNumHs() == 0
        and centre.GetNumRadicalElectrons() == 1
    ):
        return centre
    if centre.GetAtomicNum() in (16, 34, 52) and centre.GetNumRadicalElectrons() == 1 and _is_acyl(centre):
        return centre
    return None


def _sanitize(mol):
    """Sanitize, leaving the valence of an ate centre marked by `marked_neutral` unchecked (it is neutral there)."""
    if any(a.HasProp("_anion_word") for a in mol.GetAtoms()):
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        for atom in mol.GetAtoms():
            if atom.HasProp("_anion_word"):
                atom.SetNumRadicalElectrons(0)
    else:
        Chem.SanitizeMol(mol)


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
    _sanitize(hydride)
    return hydride, cap


def has_radical_shape(mol) -> bool:
    return any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())


def has_radical_group_shape(mol) -> bool:
    centre = _centre(mol)
    return centre is not None or (mol.GetNumAtoms() == 1 and sum(1 for a in mol.GetAtoms() if a.GetNumRadicalElectrons()) == 1)


def name_radical_group(mol) -> str:
    token = DIPOLAR_GROUPS.set(True)
    parent_token = RADICAL_PARENT.set(True)
    try:
        return _name_radical_group(mol)
    finally:
        RADICAL_PARENT.reset(parent_token)
        DIPOLAR_GROUPS.reset(token)


def _name_radical_group(mol) -> str:
    from ._isotope_labels import split_isotopes
    from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS

    labels = {}
    original = mol
    split = split_isotopes(mol)
    if split is not None:
        mol, labels, _ = split
    centre = _centre(mol)
    if centre is not None and sum(1 for a in mol.GetAtoms() if a.GetNumRadicalElectrons()) > 1:
        from ._radical_hetero import hetero_radical_name

        found = hetero_radical_name(mol) or polyradical_name(mol)
        if found is not None:
            return found
        if (
            len(Chem.GetMolFrags(mol)) != 1
            or mol.GetNumAtoms() > _MAX_ATOMS
            or any(a.GetAtomicNum() == 53 for a in mol.GetAtoms())
        ):
            raise UnsupportedStructure("these radical centres are not on one parent")
        if _polyradical_centres(mol):
            try:
                return _name_polyradical(mol)
            except UnsupportedStructure:
                return _name_parent_radical(mol)
        if _multi_centres(mol):
            return _name_parent_radical(mol)
        raise UnsupportedStructure("these radical centres are not on one parent")
    if centre is None or (labels and centre.GetAtomicNum() != 6):
        from ._radical_hetero import hetero_radical_name

        found = hetero_radical_name(original)
        if found is not None:
            return found
    if (
        centre is None
        or mol.GetNumAtoms() > _MAX_ATOMS
        or len(Chem.GetMolFrags(mol)) != 1
    ):
        raise UnsupportedStructure("this radical is not a single radical centre on a plain skeleton")
    if labels and centre.GetAtomicNum() != 6:
        raise UnsupportedStructure("isotopic modification of this radical is not supported yet")
    phantom = phantom_capped(mol)
    stereo = specified_stereo_elements(phantom)
    if stereo and labels:
        raise UnsupportedStructure("isotopic modification of a radical with stereo elements is not supported yet")
    if any(
        a.GetFormalCharge() < 0
        and (a.GetAtomicNum() in (6, 7) or (a.GetAtomicNum() not in _NOT_ATE_CENTRE and not a.IsInRing()))
        and not (a.GetAtomicNum() == 7 and (a.IsInRing() or a.GetFormalCharge() != -1))
        for a in mol.GetAtoms()
    ):
        from ._anion import marked_neutral

        mol = marked_neutral(mol)
        centre = mol.GetAtomWithIdx(centre.GetIdx())
    if centre.GetAtomicNum() == 7 and not centre.IsInRing() and _suffix_nitrogen_neighbours(centre):
        refuse_specified_stereo(mol)
        suffixed = _name_aminyl(mol, centre)
        if suffixed is not None:
            return suffixed
    if centre.GetAtomicNum() != 6 and centre.GetNumRadicalElectrons() > 1:
        return _name_multivalent_centre(mol, centre)
    if centre.GetAtomicNum() == 8 and centre.GetNeighbors()[0].GetAtomicNum() == 8:
        refuse_specified_stereo(mol)
        try:
            return _name_peroxyl(mol, centre)
        except UnsupportedStructure:
            from ._radical_hetero import hetero_radical_name

            found = hetero_radical_name(original)
            if found is None:
                raise
            return found
    hydride, hydrogen = _hydride(mol, centre)
    stereo_token = None
    if stereo:
        match_cip_labels(phantom, hydride)
        reference = Chem.Mol(phantom)
        assign_cip_labels(reference)
        cited = {
            "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in reference.GetAtoms() if a.HasProp("_CIPCode")},
            "bonds": {
                (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
                for b in reference.GetBonds()
                if b.HasProp("_CIPCode")
            },
            "used": set(),
        }
        stereo_token = BRANCH_STEREO.set(cited)
    graph = adjacency(hydride)
    aromatic = frozenset(a.GetIdx() for a in hydride.GetAtoms() if a.GetIsAromatic())
    context = {"labels": labels, "consumed": set(), "mol": hydride}
    token = ISOTOPE_LABELS.set(context if labels else None)
    extended = EXTENDED_PREFIXES.set(centre.GetAtomicNum() in (16, 34, 52) or mol.HasSubstructMatch(_IMINE))
    centre_token = RADICAL_CENTRE.set(centre.GetIdx())
    try:
        hydryl = None if labels else hydryl_radical_name(mol, centre)
        if hydryl is not None:
            return hydryl
        name, _ = name_branch(
            graph, centre.GetIdx(), hydrogen, halogen_substituents(hydride), aromatic, mol=hydride, unsaturated=True
        )
    finally:
        RADICAL_CENTRE.reset(centre_token)
        EXTENDED_PREFIXES.reset(extended)
        ISOTOPE_LABELS.reset(token)
        if stereo_token is not None:
            BRANCH_STEREO.reset(stereo_token)
    if stereo and (
        any(("atom", a) not in cited["used"] for a in cited["atoms"]) or any(("bond", b) not in cited["used"] for b in cited["bonds"])
    ):
        raise UnsupportedStructure("a stereo element of the radical is not cited by any supported name")
    if set(labels) - context["consumed"]:
        raise UnsupportedStructure("an isotopically modified atom of the radical is not cited by any supported name")
    if centre.GetAtomicNum() in (6, 16) and name in _ACYL_OF_HYDROXY_ACID:
        return _ACYL_OF_HYDROXY_ACID[name]
    if centre.GetAtomicNum() == 7 and name.endswith("amino"):
        return name[: -len("amino")] + "azanyl"
    if centre.GetAtomicNum() == 8:
        if not name.endswith("oxy"):
            raise UnsupportedStructure("the oxygen radical has no 'oxy' prefix to turn into 'oxyl' (P-71.3.4)")
        return _aminoxyl(name) + "l"
    return name


_ACYL_OF_HYDROXY_ACID = {"carboxy": "hydroxycarbonyl", "sulfo": "hydroxysulfonyl", "sulfino": "hydroxysulfinyl"}


def _name_multivalent_centre(mol, centre):
    """P-71.2.2.2: 'ylidene' or 'ylidyne' replaces the 'yl' of the monovalent radical whose extra radical electrons
    are filled with hydrogens, so the parent hydride and the locants are those of the 'yl' name."""
    capped = Chem.RWMol(mol)
    target = capped.GetAtomWithIdx(centre.GetIdx())
    extra = centre.GetNumRadicalElectrons() - 1
    target.SetNumRadicalElectrons(1)
    target.SetNoImplicit(True)
    target.SetNumExplicitHs(centre.GetTotalNumHs() + extra)
    monovalent = capped.GetMol()
    _sanitize(monovalent)
    monovalent.GetAtomWithIdx(centre.GetIdx()).SetNumRadicalElectrons(1)
    name = _name_radical_group(monovalent)
    if not name.endswith("yl"):
        raise UnsupportedStructure("the monovalent radical name does not end in 'yl'")
    return name[:-2] + {1: "ylidene", 2: "ylidyne"}[extra]


def _imidoyl_carbon(mol, carbon):
    return any(
        n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in carbon.GetNeighbors()
    )


def _aminoxyl(oxy_name):
    """P-71.3.4: 'aminooxy' contracts to 'aminoxy', also under substituents of the nitrogen."""
    match = re.fullmatch(r"[\[(](.*amino)[\])]oxy", oxy_name)
    if match:
        return match.group(1) + "xy"
    return oxy_name[: -len("aminooxy")] + "aminoxy" if oxy_name.endswith("aminooxy") else oxy_name


def _nitrogen_hydride(mol, centre):
    capped = Chem.RWMol(mol)
    target = capped.GetAtomWithIdx(centre.GetIdx())
    target.SetNumRadicalElectrons(0)
    target.SetNoImplicit(True)
    target.SetNumExplicitHs(target.GetNumExplicitHs() + centre.GetNumRadicalElectrons())
    hydride = capped.GetMol()
    _sanitize(hydride)
    return hydride


def _suffix_nitrogen_neighbours(centre):
    """Whether the radical nitrogen carries only hydrogen and carbon, or a sulfonyl-type chalcogen as in a sulfonamide,
    so that it is the nitrogen of an amine, imine or amide suffix; any other neighbour makes it the azanyl or
    hydrazinyl parent radical (P-71.3.2, P-71.3.3)."""
    neighbours = list(centre.GetNeighbors())
    if len(neighbours) == 1 and neighbours[0].GetAtomicNum() == 15:
        return True
    if not any(n.GetAtomicNum() == 6 or n.GetAtomicNum() in (16, 34, 52) for n in neighbours):
        return False
    return all(
        n.GetAtomicNum() in (1, 6)
        or (
            n.GetAtomicNum() in (16, 34, 52)
            and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() in (8, 7) for b in n.GetBonds())
        )
        for n in neighbours
    )


def _name_aminyl(mol, centre):
    """Aminyl and amidyl radicals: the amine or amide suffix with '-yl' (or '-ylidene' for a nitrene-like centre) added,
    every other group a prefix (P-71.3.2)."""
    from ._polyfunctional import _SULFONYL_AMIDE_ORDER, FORCED_PRINCIPAL, _name_labelled

    suffix = {1: "yl", 2: "ylidene"}[centre.GetNumRadicalElectrons()]
    hydride = _nitrogen_hydride(mol, centre)
    neighbours = list(centre.GetNeighbors())
    if len(neighbours) == 1 and neighbours[0].GetAtomicNum() == 15 and suffix == "yl":
        name = _plain_name(hydride)
        if name.endswith("imine"):
            return name[:-1] + "yl"
    carbons = [n for n in neighbours if n.GetAtomicNum() == 6]
    if suffix == "ylidene" and len(carbons) == 1 and _imidoyl_carbon(hydride, carbons[0]):
        graph = adjacency(hydride)
        token = EXTENDED_PREFIXES.set(True)
        try:
            group, compound = name_branch(
                graph, carbons[0].GetIdx(), centre.GetIdx(), halogen_substituents(hydride), frozenset(), mol=hydride
            )
        finally:
            EXTENDED_PREFIXES.reset(token)
        return (f"({group})" if compound else group) + "azanylidene"

    def acyl(carbon):
        return any(
            n.GetAtomicNum() == 8 and hydride.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in carbon.GetNeighbors()
        )

    sulfonyl = any(n.GetAtomicNum() in (16, 34, 52) for n in neighbours)
    amide = sulfonyl or any(acyl(c) for c in carbons)
    hydride = Chem.RemoveHs(hydride, sanitize=False)
    _sanitize(hydride)
    imine = any(b.GetBondTypeAsDouble() == 2.0 for b in centre.GetBonds())
    candidates = ("sulfonamide", *_SULFONYL_AMIDE_ORDER) if sulfonyl else ("amide" if amide else "imine" if imine else "amine",)
    for forced in candidates:
        token = FORCED_PRINCIPAL.set(forced)
        try:
            name = _name_labelled(hydride, {})
        except UnsupportedStructure:
            continue
        finally:
            FORCED_PRINCIPAL.reset(token)
        if name is None:
            continue
        if amide and name.endswith("amide"):
            return name[:-1] + suffix
        if not amide and name.endswith(("amine", "imine")):
            return name[:-1] + suffix
        if not amide and name.endswith("aniline"):
            return name[: -len("aniline")] + "benzenamin" + suffix
    return None


def _plain_name(mol):
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(mol))


def _name_peroxyl(mol, centre):
    hydride, hydrogen = _hydride(mol, centre)
    peroxide = centre.GetNeighbors()[0]
    carbons = [n for n in peroxide.GetNeighbors() if n.GetIdx() != centre.GetIdx()]
    if len(carbons) != 1 or (carbons[0].GetAtomicNum() != 6 and carbons[0].GetAtomicNum() not in MONONUCLEAR_HYDRIDES):
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

    anionic = Chem.RWMol(mol)
    for atom in _polyradical_centres(mol):
        target = anionic.GetAtomWithIdx(atom.GetIdx())
        target.SetNumRadicalElectrons(0)
        target.SetFormalCharge(-1)
    anionic = anionic.GetMol()
    Chem.SanitizeMol(anionic)
    match_cip_labels(phantom_capped(mol), anionic)
    name = name_anion(anionic)
    for ending in ("diide", "triide", "tetraide"):
        if name.endswith(ending) and "-" in name[: -len(ending)]:
            if re.search(r"-u?ide?-", name[: -len(ending)]):
                break
            return name[: -len("ide")] + "yl"
    raise UnsupportedStructure("the polyradical is not named as a polyanion parent")


def _carbenium_valence(atom):
    """A carbenium centre: three bonds in total, three single bonds or one double and one single bond (a vinyl or aryl
    cation, named through the carbanion of the same skeleton)."""
    orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
    if sorted(orders) in ([1.0, 2.0], [2.0]) and atom.GetTotalNumHs() + sum(orders) == 3 and (
        atom.GetBonds()[0].GetOtherAtom(atom).GetAtomicNum() == 6 or atom.IsInRing()
    ):
        return atom.IsInRing() or all(b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds())
    return atom.GetTotalNumHs() + atom.GetDegree() == 3 and all(order == 1.0 for order in orders)


def _pi_cation_kekule(mol):
    """`mol` with its aromatic bonds written as alternating ones when the carbenium centre is an aromatic atom with three
    sigma bonds (the empty p orbital of a tropylium-type cation), else `mol`."""
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(cations) != 1 or not cations[0].GetIsAromatic() or cations[0].GetTotalNumHs() + cations[0].GetDegree() != 3:
        return mol
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.rdchem.KekulizeException:
        return mol
    return kekule


def _group_cation_centre(mol):
    mol = _pi_cation_kekule(mol)
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
    from ._anion import HYDRO_PREFIXED_CENTRE, name_anion

    centre = _group_cation_centre(mol)
    if centre is None:
        raise UnsupportedStructure("this cation is not a single carbenium centre on a plain skeleton")
    mol = _pi_cation_kekule(mol)
    centre = mol.GetAtomWithIdx(centre.GetIdx())
    if _is_unbranched_alkane_with_terminal_cation(mol, centre):
        # P-73.2.2.1.1: the specific method turns the 'ane' of the parent hydride into 'ylium'
        return alkyl_name(mol.GetNumAtoms()) + "ium"
    anionic = Chem.RWMol(mol)
    anionic.GetAtomWithIdx(centre.GetIdx()).SetFormalCharge(-1)
    anionic = anionic.GetMol()
    Chem.SanitizeMol(anionic)
    match_cip_labels(mol, anionic)
    token = HYDRO_PREFIXED_CENTRE.set(True)
    try:
        name = name_anion(anionic)
    finally:
        HYDRO_PREFIXED_CENTRE.reset(token)
    if re.search(r"\d+H\)-ide$", name):
        if centre.GetDegree() == 2:
            raise UnsupportedStructure("an added-hydrogen anion name does not carry over to an aryl cation")
        # P-73.2.2.2: the added hydrogen of a mancude cation is cited as for the anion, 'anthracen-4a(2H)-ylium'
        return name[: -len("ide")] + "ylium"
    if name.endswith("anide"):
        return name[: -len("anide")] + "ylium"
    if name.endswith(("benzenide", "ethenide")):
        return name[: -len("ide")] + "ylium"
    if name.endswith("ide") and name[:-3].endswith("-"):
        return _delocalized_name(mol, centre, name[:-3] + "ylium")
    raise UnsupportedStructure("the cation is not named as a carbanion-like parent")


def _delocalized_name(mol, centre, name):
    """P-76: a cation whose centre lies in a fused ring system with no sp3 atom is delocalized over the whole system, so
    the mancude parent name takes 'ylium' without indicated hydrogen or locant (benzo[7]annulenylium)."""
    from ._multiplicative import _ring_systems

    match = re.fullmatch(r"\d+[a-z]?H-([^-]+)-\d+[a-z]?-ylium", name)
    system = next((s for s in _ring_systems(mol) if centre.GetIdx() in s), set())
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= system]
    sp2 = all(mol.GetAtomWithIdx(a).GetTotalNumHs() + mol.GetAtomWithIdx(a).GetDegree() <= 3 for a in system)
    return match.group(1) + "ylium" if match and len(rings) > 1 and sp2 else name


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
    match_cip_labels(mol, out)
    return out


def _name_parent_radical(mol):
    """Several radical centres that do not fit one parent (P-71.7): the parent holds the most centres, then rings
    outrank chains; the other centres are cited with the nondetachable prefix 'ylo' (P-71.5)."""
    from itertools import combinations

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


def has_labelled_group_cation_shape(mol) -> bool:
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    return (
        len(cations) == 1
        and cations[0].GetAtomicNum() == 6
        and cations[0].GetFormalCharge() == 1
        and not cations[0].GetIsAromatic()
        and _carbenium_valence(cations[0])
        and any(a.GetIsotope() for a in mol.GetAtoms())
        and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        and len(Chem.GetMolFrags(mol)) == 1
    )


def name_labelled_group_cation(mol) -> str:
    """A carbenium centre with nuclides (P-82.2.1): the name of the carbanion of the same labelled skeleton with 'ide' read
    as 'ylium' (P-73.2.2.1)."""
    from .core import smiles_to_iupac

    anionic = Chem.RWMol(mol)
    for atom in anionic.GetAtoms():
        if atom.GetFormalCharge():
            atom.SetFormalCharge(-1)
    anionic = anionic.GetMol()
    Chem.SanitizeMol(anionic)
    name = smiles_to_iupac(Chem.MolToSmiles(anionic))
    terminal = re.fullmatch(r"((?:\([^()]*\))?[a-z]+)an-1-ide", name)
    if terminal:
        return terminal.group(1) + "ylium"
    if name.endswith("anide"):
        return name[: -len("anide")] + "ylium"
    if re.search(r"-\d+-ide$", name):
        return name[: -len("ide")] + "ylium"
    raise UnsupportedStructure("the carbanion of this skeleton is not named by a parent hydride and the suffix 'ide'")
