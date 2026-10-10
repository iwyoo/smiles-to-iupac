"""Amides of the halogen oxoacids, R-NH-X, R-NH-XO and R2N-X (P-61.3.2.2, P-62.4, P-68.5.3): a halogen atom on nitrogen
makes an amide of an inorganic acid, which outranks the amine and halo classes of the 'N-halogenoamine' name, so the
PIN is '<substituents>hypochlorous amide' (bromous, chloric, perchloric amide with one to three oxygens on the halogen)
with the carbon groups on nitrogen cited without locants."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_ACID_WORDS = {
    9: ("hypofluorous", "fluorous", "fluoric", "perfluoric"),
    17: ("hypochlorous", "chlorous", "chloric", "perchloric"),
    35: ("hypobromous", "bromous", "bromic", "perbromic"),
    53: ("hypoiodous", "iodous", "iodic", "periodic"),
}


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
        and sum(n.GetAtomicNum() in HALOGEN_PREFIXES for n in a.GetNeighbors()) == 1
    ]
    return found[0] if len(found) == 1 else None


def has_halogen_amide_shape(mol):
    return _halogen_amide_nitrogen(mol) is not None


def name_halogen_amide(mol) -> str:
    nitrogen = _halogen_amide_nitrogen(mol)
    if nitrogen is None:
        raise UnsupportedStructure("no nitrogen bearing exactly one halogen atom")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    halogen = next(n for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES)
    oxygens = _oxo_oxygens(halogen, nitrogen)
    if oxygens is None:
        raise UnsupportedStructure("a halogen on nitrogen carries something other than terminal oxygens")
    owned = {halogen.GetIdx(), *(o.GetIdx() for o in oxygens)}
    if any(
        a.GetIdx() not in owned and (a.GetAtomicNum() not in (6, 7, *HALOGEN_PREFIXES) or a.GetFormalCharge() or a.GetIsotope())
        for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("only carbon groups and halogens may accompany an amide of a halogen acid")
    if nitrogen.IsInRing() or any(
        a.GetAtomicNum() == 7 and a.GetIdx() != nitrogen.GetIdx() and not a.IsInRing() for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("a ring nitrogen or a second nitrogen is not an amide of a halogen acid")
    carbons = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != halogen.GetIdx()]
    if not carbons or any(
        mol.GetBondBetweenAtoms(nitrogen.GetIdx(), c.GetIdx()).GetBondTypeAsDouble() != 1.0
        or not (c.GetAtomicNum() == 6 or (c.GetAtomicNum() == 7 and c.IsInRing()))
        for c in carbons
    ):
        raise UnsupportedStructure("an amide of a halogen acid needs single-bonded carbon or ring nitrogen groups on nitrogen")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, c.GetIdx(), nitrogen.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for c in carbons
    ]
    return f"{format_mononuclear_prefixes(entries)}{_ACID_WORDS[halogen.GetAtomicNum()][len(oxygens)]} amide"
