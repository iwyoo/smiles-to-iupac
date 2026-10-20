"""Amides of the halogen oxoacids, R-NH-X, R-NH-XO and R2N-X (P-61.3.2.2, P-62.4, P-68.5.3): a halogen atom on nitrogen
makes an amide of an inorganic acid, which outranks the amine and halo classes of the 'N-halogenoamine' name, so the
PIN is '<substituents>hypochlorous amide' (bromous, chloric, perchloric amide with one to three oxygens on the halogen)
with every group junior to the amide (nitrile, aldehyde, ketone, hydroxy, sulfanyl, amino, ...) cited without locants;
a carboxylic or sulfonic acid, ester, acid halide or carboxamide elsewhere is senior and keeps the halogen a prefix."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency
from ._substituents import format_mononuclear_prefixes

_ACID_WORDS = {
    9: ("hypofluorous", "fluorous", "fluoric", "perfluoric"),
    17: ("hypochlorous", "chlorous", "chloric", "perchloric"),
    35: ("hypobromous", "bromous", "bromic", "perbromic"),
    53: ("hypoiodous", "iodous", "iodic", "periodic"),
}


def _nuclide(halogen):
    """'(18F)': the descriptor of a halogen of a given mass number, cited without a count (P-82.2.1, P-82.2.2)."""
    from ._isotope_labels import descriptor

    if not halogen.GetIsotope():
        return ""
    entry = {"skeleton": f"{halogen.GetIsotope()}{halogen.GetSymbol()}", "H": {}}
    return descriptor({halogen.GetIdx(): entry}, {halogen.GetIdx(): 1}, True)


def _oxo_oxygens(halogen, nitrogen):
    """The terminal oxygens (=O or -O(-)) of a halogen bonded to the nitrogen, None when it carries anything else."""
    found = []
    for bond in halogen.GetBonds():
        other = bond.GetOtherAtom(halogen)
        if other.GetIdx() == nitrogen.GetIdx():
            continue
        double = bond.GetBondTypeAsDouble() == 2.0 and not other.GetFormalCharge()
        anionic = bond.GetBondTypeAsDouble() == 1.0 and other.GetFormalCharge() == -1
        if other.GetAtomicNum() != 8 or other.GetDegree() != 1 or not (double or anionic):
            return None
        found.append(other)
    return found if len(found) <= 3 and sum(o.GetFormalCharge() for o in found) + halogen.GetFormalCharge() == 0 else None


def _halogen_amide_nitrogen(mol):
    found = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 7
        and sum(n.GetAtomicNum() in HALOGEN_PREFIXES and not n.GetAtomMapNum() for n in a.GetNeighbors()) in (1, 2, 3)
    ]
    return found[0] if len(found) == 1 else None


def has_halogen_amide_shape(mol):
    from ._chalcogenourea import has_group_senior_to_urea

    nitrogen = _halogen_amide_nitrogen(mol)
    if nitrogen is None:
        return False
    owned = {nitrogen.GetIdx()}
    for halogen in (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES and not n.GetAtomMapNum()):
        owned |= {halogen.GetIdx(), *(o.GetIdx() for o in _oxo_oxygens(halogen, nitrogen) or ())}
    return not has_group_senior_to_urea(mol, owned)


def name_halogen_amide(mol) -> str:
    from ._chalcogenourea import _HYDRIDE_GROUP_CENTRES, is_element_hydride_group, n_substituent_names

    nitrogen = _halogen_amide_nitrogen(mol)
    if nitrogen is None:
        raise UnsupportedStructure("no nitrogen bearing one to three halogen atoms")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    halogen, *second = sorted(
        (n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES and not n.GetAtomMapNum()),
        key=lambda n: -n.GetIsotope(),
    )
    oxygens = _oxo_oxygens(halogen, nitrogen)
    if oxygens is None:
        raise UnsupportedStructure("a halogen on nitrogen carries something other than terminal oxygens")
    owned = {halogen.GetIdx(), *(o.GetIdx() for o in oxygens)}
    prefixes = []
    for other in second:
        if other.GetAtomicNum() != halogen.GetAtomicNum() or oxygens or other.GetDegree() != 1:
            raise UnsupportedStructure("two different halogens on nitrogen need a decided choice of the acid")
        owned.add(other.GetIdx())
        prefixes.append((_nuclide(other) + HALOGEN_PREFIXES[other.GetAtomicNum()], False))
    if nitrogen.IsInRing() or any(
        a.GetAtomicNum() == 7 and not a.IsInRing() and a.GetIdx() in {n.GetIdx() for n in nitrogen.GetNeighbors()}
        for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("a ring nitrogen or a second nitrogen is not an amide of a halogen acid")
    if any(a.GetFormalCharge() or a.GetIsotope() for a in mol.GetAtoms() if a.GetIdx() not in owned) or any(
        o.GetIsotope() for o in oxygens
    ):
        raise UnsupportedStructure("a charged or labelled atom beside an amide of a halogen acid is not supported yet")
    if any(a.GetAtomicNum() in _HYDRIDE_GROUP_CENTRES and not is_element_hydride_group(a) for a in mol.GetAtoms()):
        raise UnsupportedStructure("an oxoacid or amide centre of another element outranks an amide of a halogen acid")
    if any(
        b.GetBondTypeAsDouble() != 1.0 for b in nitrogen.GetBonds() if b.GetOtherAtomIdx(nitrogen.GetIdx()) not in owned
    ):
        raise UnsupportedStructure("an amide of a halogen acid needs single bonds from the nitrogen to its groups")
    graph = adjacency(mol)
    (groups,) = n_substituent_names(mol, graph, owned | {nitrogen.GetIdx()}, (nitrogen.GetIdx(),), halogen.GetIdx(), junior_groups=True)
    entries = prefixes + groups
    return f"{format_mononuclear_prefixes(entries)}{_nuclide(halogen)}{_ACID_WORDS[halogen.GetAtomicNum()][len(oxygens)]} amide"
