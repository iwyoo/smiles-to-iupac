"""Acyclic chain parents carrying any mix of supported groups (P-41, P-44.1.1,
P-44.3, P-45): the most senior class present becomes the suffix and every
other group, ring or branch is cited as a substituent prefix through
`name_branch` (P-29.3.3, P-29.4). Rings may only be substituents; a principal
group on a ring is not handled here.
"""

import contextlib
import contextvars
import itertools
from collections import Counter
from types import MappingProxyType
import re

import networkx as nx
from rdkit import Chem, rdBase
from rdkit.Chem import CanonicalRankAtoms

from ._common import KEEP_ASSEMBLY_ENDING
from ._common import assign_cip_labels
from ._common import (
    CITE_SKELETAL_LAMBDA,
    alphanumerical_name_key,
    citation_order_key,
    assembly_join,
    alpha_sort_key,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_prefix_placeholder,
    lowest_locant_set,
    multiplied_word,
    nonstandard_bonding,
    name_from_substituents,
    ring_cycle,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._amino_acyl_group import CATION_INSERT_PROP
from ._anion import ANION_PROP, anion_weight
from ._numerals import multiplying_prefix
from ._chalcogenourea import is_oxo_nitrogen
from ._functional_prefixes import is_nitro_nitrogen
from ._hetero_prefixes import (
    CATION_PARENT,
    EXTENDED_PREFIXES,
    IMINE_PARENT,
    LAMBDA_CENTRE_STEMS,
    MONONUCLEAR_HYDRIDES,
    _oxo_or_imino,
    _ring_outranks_centre,
    azoxy_atoms,
    chalcogen_link_is_parent,
    is_functional_carbon,
    is_halogen_oxo_part,
    azono_ligands,
)
from ._multiplicative import _bare_key
from ._multiplicative_text import JunctionBondLocant, PrimedLocant, enclose, unit_phrase
from ._multiplicative_ring import (
    NO_RETAINED_BENZENE,
    RingSpec,
    _RETAINED_BENZENE,
    _SUFFIX_WORDS,
    _citation_key,
    _join,
    _prefix_text,
    _suffix_text,
    monocycle_spec,
    multiple_locants,
    name_ring_component,
    numberings,
    spec_of,
)
from ._fused_numbering import HETERO_RANK as _HETERO_RANK
from ._ring_diyl_numbering import OMIT_SUFFIX_LOCANT, _exocyclic_oxo
from ._ring_system_seniority import ring_seniority_key
from ._acid_groups import acid_group_at
from ._acid_lexicon import SULFONYL_CENTRES, SULFONYL_GROUP_NAMES, carbo_suffix, chain_suffix, make_spec, rank_key, spec_from_key
from ._retained_acids import retained_chain_acid, single_site_prefixes, split_retained_prefix
from ._hydroperoxide import peroxol_suffix
from ._locant_omission import has_single_arrangement, omits_all_locants, omits_group_locants
from ._substituents import BRANCH_STEREO, _branch_stereo_rank, format_substituent_prefixes, name_branch

_AMIDRAZONE = ("hydrazonamide", "imidohydrazide", "hydrazonohydrazide")
_CHALCOGEN_AMIDE = {16: "thioamide", 34: "selenoamide", 52: "telluroamide"}
_CHALCOGEN_AMIDE_CLASSES = tuple(_CHALCOGEN_AMIDE.values())
_REPLACED_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}
_RING_ACID_CENTRES = frozenset({5, 13, 14, 15, 31, 32, 33, 49, 50, 51, 81, 82, 83})
_SULFONYL_AMIDE_ORDER = [n for (_, _, t), n in SULFONYL_GROUP_NAMES.items() if t == "amide" and n != "sulfonamide"]
_SULFONYL_HYDRAZIDE_ORDER = [n for (_, _, t), n in SULFONYL_GROUP_NAMES.items() if t == "hydrazide"]


def _sulfonyl_names(terminal, nitrogen):
    """Names of the groups with `terminal` ('amide' or 'hydrazide') whose =X positions hold no nitrogen group ('plain'),
    an imido group but no hydrazono group ('imide') or a hydrazono group ('hydrazone')."""
    kinds = {"NNH2": "hydrazone", "NH": "imide"}
    return {
        key: name
        for key, name in SULFONYL_GROUP_NAMES.items()
        if key[2] == terminal and next((kinds[x] for x in ("NNH2", "NH") if x in key[1]), "plain") == nitrogen
    }


_PEROXOL_ELEMENTS = (8, 16, 34, 52)
_PEROXOL_ANALOGUES = tuple(
    peroxol_suffix(a, t)
    for a, t in ((16, 8), (34, 8), (52, 8), (8, 16), (8, 34), (8, 52), (16, 16), (34, 16), (52, 16), (16, 34), (16, 52), (34, 34), (52, 34), (34, 52), (52, 52))
)
_PEROXOL_CLASSES = ("peroxol", *_PEROXOL_ANALOGUES)
_CHALCOGEN_HYDRAZIDE = _sulfonyl_names("hydrazide", "plain")
_CHALCOGEN_IMIDAMIDE = {**_sulfonyl_names("amide", "imide"), **_sulfonyl_names("hydrazide", "imide")}
_CHALCOGEN_HYDRAZIDINE = {**_sulfonyl_names("amide", "hydrazone"), **_sulfonyl_names("hydrazide", "hydrazone")}
_CHALCOGEN_SULFONAMIDE_CLASSES = tuple(name for name in _sulfonyl_names("amide", "plain").values() if name != "sulfonamide")
_SENIORITY = [
    "ide", "acid", "thioic", "peroxoic", "imidic", "sulfonic", "amide", *_CHALCOGEN_AMIDE_CLASSES, "amidine", *_AMIDRAZONE, "sulfonamide", *_SULFONYL_AMIDE_ORDER, "hydrazide", *_SULFONYL_HYDRAZIDE_ORDER, "nitrile", "aldehyde", "ketone", "thione", "selone", "tellone", "alcohol",
    "thiol", "selenol", "tellurol", "peroxol", *_PEROXOL_ANALOGUES, "amine", "imine",
]
_TERMINAL = {"acid", "thioic", "peroxoic", "imidic", "amide", *_CHALCOGEN_AMIDE_CLASSES, "amidine", *_AMIDRAZONE, "hydrazide", "nitrile", "aldehyde"}
_PLAIN_ACID = make_spec("C", ["O"], ["O"])
_PLAIN_SULFONIC = make_spec("S", ["O", "O"], ["O"])


def _is_variant(name):
    return bool(name) and name.startswith("acid:")


def _acid_rank(name):
    spec = {"acid": _PLAIN_ACID, "sulfonic": _PLAIN_SULFONIC}.get(name) or spec_from_key(name)
    return rank_key(spec)


AMINIUM = contextvars.ContextVar("aminium", default=False)
RING_CENTER = contextvars.ContextVar("ring_center", default=False)
CATION_RING_ATOM = contextvars.ContextVar("cation_ring_atom", default=None)
FORCE_LOCANTS = contextvars.ContextVar("force_locants", default=False)
LAST_POSITIONS = contextvars.ContextVar("last_positions", default=None)
FORCED_PRINCIPAL = contextvars.ContextVar("forced_principal", default=None)
SUBSTITUTED_AMINE_PREFIX = contextvars.ContextVar("substituted_amine_prefix", default=False)
PREFERRED_AMINE_ATOM = contextvars.ContextVar("preferred_amine_atom", default=None)


def _senior_anionic_centres(mol, groups, ring_groups, marked):
    """P-72.7(d): between a peroxolate and a chalcogenolate the anionic centre of the first element in O > S > Se > Te
    decides; with the same element the thiol class outranks the peroxol class (P-72.7(e))."""
    thiols = ("thiol", "selenol", "tellurol")

    def centre(owned):
        return min(mol.GetAtomWithIdx(a).GetAtomicNum() for a in owned & marked)

    peroxols = [centre(o) for cls, members in groups.items() if cls in _PEROXOL_CLASSES for o in members.values()]
    peroxols += [centre(owned) for cls, _, owned in ring_groups if cls in _PEROXOL_CLASSES]
    chalcogenols = [centre(o) for cls, members in groups.items() if cls in thiols for o in members.values()]
    chalcogenols += [centre(owned) for cls, _, owned in ring_groups if cls in thiols]
    if peroxols and chalcogenols and min(peroxols) < min(chalcogenols):
        groups = {cls: members for cls, members in groups.items() if cls not in thiols}
        ring_groups = [g for g in ring_groups if g[0] not in thiols]
    return groups, ring_groups


def _principal_class(classes):
    if FORCED_PRINCIPAL.get():
        return FORCED_PRINCIPAL.get() if FORCED_PRINCIPAL.get() in classes else None
    if RING_CENTER.get():
        return None
    if AMINIUM.get():
        parent = AMINIUM.get() if AMINIUM.get() in ("imine", "amide", "nitrile", "amidine") else "amine"
        return parent if parent in classes else None
    if "ide" in classes:
        return "ide"
    acids = [c for c in classes if c in ("acid", "sulfonic") or _is_variant(c)]
    if acids:
        return min(acids, key=_acid_rank)
    return next((name for name in _SENIORITY if name in classes), None)


def _anionic_group_atom(mol, atom):
    """An -O(-) or -S(-) atom that is the anionic end of an acid group (P-72.2.2.2.1)."""
    if atom.GetFormalCharge() != -1 or atom.GetAtomicNum() not in (8, 16, 34, 52):
        return False
    for n in atom.GetNeighbors():
        for start in (n, *n.GetNeighbors()):
            group = acid_group_at(mol, start.GetIdx())
            if group is not None and group.spec.anion and atom.GetIdx() in group.owned:
                return True
    return False


def _is_acid_family(name):
    return name in ("acid", "sulfonic") or _is_variant(name)


def _junior_end_group(mol, idx):
    """A carbon of an acid derivative (ester, amide, anhydride, thioester ...) at
    the end of a chain: it joins the chain and is cited as 'oxo' with the
    prefix of its other attachment (P-65.1.6.1, P-65.5.4, P-65.6.3.2.3)."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.IsInRing() or not is_functional_carbon(mol, idx):
        return False
    group = acid_group_at(mol, idx)
    if group is not None and group.spec.plain:
        return False
    if not any(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds()):
        return False
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6) == 1


def _is_terminal(name):
    return name in _TERMINAL or (_is_variant(name) and spec_from_key(name).center == "C")
_MAX_ATOMS = 80
IDE_EXTRA = contextvars.ContextVar("ide_extra", default=MappingProxyType({}))


def _double_oxygens(mol, carbon):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]


def _single_neighbors(mol, carbon, atomic_num):
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(carbon).GetNeighbors()
        if n.GetAtomicNum() == atomic_num and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
    ]


def _terminal_heteroatom(mol, idx, hydrogens):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetDegree() == 1 and atom.GetTotalNumHs() == hydrogens and not atom.GetFormalCharge()


_CHALCOGEN_KETONE_CLASS = {16: "thione", 34: "selone", 52: "tellone"}
_CARBOXYLIC_OR_SULFONIC = Chem.MolFromSmarts("[#6,#16](=[O,S,Se,Te])[O,S,Se,Te;!$(*C#N)]")
_SENIOR_TO_NITRILE = [
    Chem.MolFromSmarts("[CX3;!R](=[O,S,Se,Te])[O,S,Se,Te,N,F,Cl,Br,I;!$(*C#N)]"),
    Chem.MolFromSmarts("[#16X4,#16X3;!R](=O)[O,N,F,Cl,Br,I;!$(*C#N)]"),
]
_CHALCOGEN_KETONE_OK = {
    "acid", "thioic", "peroxoic", "imidic", "sulfonic", "amide", *_CHALCOGEN_AMIDE_CLASSES, "amidine", *_AMIDRAZONE, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES, *_CHALCOGEN_IMIDAMIDE.values(), "hydrazide", *_CHALCOGEN_HYDRAZIDE.values(), *_CHALCOGEN_HYDRAZIDINE.values(), "nitrile", "aldehyde", "ketone",
    *_CHALCOGEN_KETONE_CLASS.values(),
}


def _chalcogen_ketone(mol, atom):
    """A carbon double-bonded to S, Se or Te (thione, selone, tellone)."""
    if atom.GetAtomicNum() != 6:
        return False
    return not _is_isocyanate_carbon(atom) and any(
        n.GetAtomicNum() in (16, 34, 52)
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )


_CHALCOGEN_RANK = {8: 0, 16: 1, 34: 2, 52: 3}
_RING_HETERO_NEIGHBORS = {7, 8, 16, 34, 52}


def _urea_carbon(mol, idx):
    """The carbon of a urea core C(=X)(N)N, X = O, S, Se or Te, the amide of carbonic acid: it ranks below every amide of
    a carboxylic or sulfonic acid (P-66.1.6.1.1.5) and is cited as 'carbamoyl' or 'carbamothioyl' beneath one."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or _acyl_chalcogen_of(mol, atom) is None:
        return False
    return len(_single_neighbors(mol, idx, 7)) == 2


def _acyl_chalcogen_of(mol, atom):
    return next(
        (
            n.GetAtomicNum()
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() in _CHALCOGEN_RANK
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ),
        None,
    )


def _outranks_urea(principal):
    return (
        _is_acid_family(principal)
        or principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES)
    )


def _acyl_chalcogen(mol, idx):
    """Atomic number of the terminal =O/=S/=Se/=Te of an acyl carbon R-C(=X)- (R carbon or H), else None."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return None
    chalcogens = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() in _CHALCOGEN_RANK
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    taken = {n.GetIdx() for n in chalcogens}
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() not in taken]
    if len(chalcogens) != 1 or len(rest) != 2 - (atom.GetTotalNumHs() == 1):
        return None
    if any(n.GetAtomicNum() not in (6, 7) for n in rest):
        return None
    return chalcogens[0].GetAtomicNum()


def _acyl_branch_size(mol, carbon):
    """(ring attached, longest carbon path) of the part of an acyl group beyond its carbon."""
    start = [n.GetIdx() for n in mol.GetAtomWithIdx(carbon).GetNeighbors() if n.GetAtomicNum() == 6]
    if not start:
        return (0, 0)
    ring = int(mol.GetAtomWithIdx(start[0]).IsInRing())
    best, stack = 0, [(start[0], {carbon, start[0]})]
    while stack:
        node, seen = stack.pop()
        best = max(best, len(seen) - 1)
        for n in mol.GetAtomWithIdx(node).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIdx() not in seen:
                stack.append((n.GetIdx(), seen | {n.GetIdx()}))
    return (ring, best)


def _imide_parent(mol, carbon, partner):
    """Whether acyl `carbon` rather than `partner` is the parent amide of an imide nitrogen: the more senior chalcogen
    analogue (O, S, Se, Te), then a ring over a chain and the longer chain (P-66.1.2, P-44.1)."""
    z, other = _acyl_chalcogen(mol, carbon), _acyl_chalcogen(mol, partner)
    if _urea_carbon(mol, carbon) != _urea_carbon(mol, partner):
        return _urea_carbon(mol, partner)
    if z != other:
        return _CHALCOGEN_RANK[z] < _CHALCOGEN_RANK[other]
    mine, theirs = _acyl_branch_size(mol, carbon), _acyl_branch_size(mol, partner)
    return mine > theirs or (mine == theirs and carbon < partner)


def _organyloxy(mol, chalcogen, nitrogen):
    """A chain of chalcogen atoms (-O-, -S-, -Se-, -Te-) on nitrogen that ends in a plain carbon group: an alkoxy,
    sulfanyl or disulfanyl-type prefix (P-63.2.2.1.1, P-63.4.2.2)."""
    previous, atom = nitrogen, chalcogen
    while atom.GetAtomicNum() in _CHALCOGEN_RANK:
        if atom.GetDegree() != 2 or atom.GetFormalCharge() or atom.IsInRing():
            return False
        if mol.GetBondBetweenAtoms(previous.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return False
        previous, atom = atom, next(n for n in atom.GetNeighbors() if n.GetIdx() != previous.GetIdx())
    return (
        atom.GetAtomicNum() == 6
        and mol.GetBondBetweenAtoms(previous.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not _double_oxygens(mol, atom.GetIdx())
        and not is_functional_carbon(mol, atom.GetIdx())
    )


def _carbamimidoyl_on(mol, atom, nitrogen):
    """Whether `atom` is the carbon of an amidine group -C(=NR)-R' bonded to `nitrogen`, R' being hydrogen, a plain carbon
    group or an amino group: the N-substituent is an imidoyl or carbamimidoyl prefix, since the amide outranks the
    amidine (P-41, P-66.4.1.2.1.3, P-66.4.1.3)."""
    if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() or atom.IsInRing():
        return False
    if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False
    ends = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    imino = [n for n in ends if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    rest = [n for n in ends if n not in imino]
    if len(imino) != 1 or len(rest) > 1 or imino[0].GetAtomicNum() != 7 or imino[0].GetFormalCharge() or imino[0].IsInRing():
        return False
    if atom.GetTotalNumHs() != (0 if rest else 1) or not (
        _plain_nitrogen_substituents(mol, imino[0], atom.GetIdx()) or _amino_hydrazone_nitrogen(mol, imino[0], atom.GetIdx())
    ):
        return False
    return not rest or (
        rest[0].GetAtomicNum() == 6 and not is_functional_carbon(mol, rest[0].GetIdx())
        and mol.GetBondBetweenAtoms(atom.GetIdx(), rest[0].GetIdx()).GetBondTypeAsDouble() == 1.0
    ) or (
        rest[0].GetAtomicNum() == 7 and not rest[0].GetFormalCharge() and not rest[0].GetIsAromatic()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), rest[0].GetIdx()).GetBondTypeAsDouble() == 1.0
        and _carbamimidoyl_nitrogen_substituents(mol, rest[0], atom.GetIdx())
    )


def _amino_hydrazone_nitrogen(mol, nitrogen, origin):
    """Whether the only other neighbour of `nitrogen` is a free amino group: the hydrazone nitrogen of a carbamohydrazonoyl
    group."""
    others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != origin]
    return (
        len(others) == 1
        and others[0].GetAtomicNum() == 7
        and others[0].GetDegree() == 1
        and not others[0].GetFormalCharge()
        and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), others[0].GetIdx()).GetBondTypeAsDouble() == 1.0
    )


def _carbamimidoyl_nitrogen_substituents(mol, nitrogen, origin):
    """Whether every neighbour of `nitrogen` other than `origin` is a plain carbon group or a further carbamimidoyl group
    (the N'-carbamimidoylcarbamimidoyl prefix of a condensed guanidine)."""
    return all(
        (
            n.GetAtomicNum() in (6, 7, 8, *HALOGEN_PREFIXES)
            and not is_functional_carbon(mol, n.GetIdx())
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        )
        or _carbamimidoyl_on(mol, n, nitrogen)
        or _plain_acyl_on(mol, n, nitrogen)
        for n in nitrogen.GetNeighbors()
        if n.GetIdx() != origin
    )


def _plain_acyl_on(mol, atom, nitrogen):
    """Whether `atom` is the carbonyl carbon of an acyl group R-CO- (R a plain carbon group) bonded to `nitrogen`: the second
    acyl group of a guanidine is an N-substituent of the carbamimidoyl prefix of the amide parent (P-66.4.1.2.1)."""
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False
    oxo = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() not in (nitrogen.GetIdx(), *(o.GetIdx() for o in oxo))]
    return (
        len(oxo) == 1
        and len(rest) == 1
        and rest[0].GetAtomicNum() == 6
        and mol.GetBondBetweenAtoms(atom.GetIdx(), rest[0].GetIdx()).GetBondTypeAsDouble() == 1.0
        and not is_functional_carbon(mol, rest[0].GetIdx())
    )


def _plain_nitrogen_substituents(mol, nitrogen, origin):
    """Whether every neighbour of `nitrogen` other than `origin` is a plain carbon group joined by a single bond."""
    return all(
        n.GetAtomicNum() == 6 and not is_functional_carbon(mol, n.GetIdx())
        and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in nitrogen.GetNeighbors()
        if n.GetIdx() != origin
    )


def _plain_amide_nitrogen(mol, nitrogen, carbonyl):
    """An amide nitrogen carrying carbon, hydroxy or organyloxy substituents (a hydroxamic acid, P-65.1.3.4, or its
    O-organyl ether) and at most one acyl group -- named with 'N-' prefixes on the amide parent."""
    if (nitrogen.GetFormalCharge() and not (AMINIUM.get() == "amide" and nitrogen.GetFormalCharge() == 1)) or (
        nitrogen.IsInRing() or nitrogen.GetIsAromatic()
    ):
        return False
    others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbonyl]
    hydroxy = [n for n in others if _terminal_heteroatom(mol, n.GetIdx(), 1) and n.GetAtomicNum() in (8, 16, 34, 52)]
    oxy = [n for n in others if _organyloxy(mol, n, nitrogen)]
    acyl = [n for n in others if mol.GetAtomWithIdx(carbonyl).GetAtomicNum() != 6 and _urea_carbon(mol, n.GetIdx())]
    if mol.GetAtomWithIdx(carbonyl).GetAtomicNum() == 6 and _acyl_chalcogen(mol, carbonyl) is not None:
        acyl = [n for n in others if _acyl_chalcogen(mol, n.GetIdx()) is not None]
        if len(acyl) > 1 or (acyl and not _imide_parent(mol, carbonyl, acyl[0].GetIdx())):
            return False
    return bool(others) and len(hydroxy) <= 1 and all(
        n in hydroxy
        or n in oxy
        or n in acyl
        or _acyloxy_on_nitrogen(mol, n, nitrogen)
        or _carbamimidoyl_on(mol, n, nitrogen)
        or _cyano_on(mol, n, nitrogen)
        or n.GetAtomicNum() in HALOGEN_PREFIXES
        or (
            n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and not _double_oxygens(mol, n.GetIdx())
            and not is_functional_carbon(mol, n.GetIdx())
        )
        or _hydride_group_on_nitrogen(mol, n, nitrogen, len(others))
        or _ring_nitrogen_substituent(mol, n, nitrogen)
        or _chalcogen_ylidene_on(mol, n, nitrogen)
        for n in others
    )


def _chalcogen_ylidene_on(mol, atom, nitrogen):
    """The lambda4- or lambda6-chalcogen centre of an N-ylidene group on an amide nitrogen, carrying carbon groups, oxo
    and imino groups only (P-68.4.3.3, P-68.4.3.4)."""
    if (
        atom.GetAtomicNum() not in (16, 34, 52)
        or atom.IsInRing()
        or atom.GetFormalCharge()
        or atom.GetTotalValence() not in (4, 6)
        or mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 2.0
    ):
        return False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == nitrogen.GetIdx():
            continue
        if bond.GetBondTypeAsDouble() == 2.0:
            if not _oxo_or_imino(mol, atom.GetIdx(), other.GetIdx()):
                return False
        elif other.GetAtomicNum() != 6 or bond.GetBondTypeAsDouble() != 1.0 or is_functional_carbon(mol, other.GetIdx()):
            return False
    return True


def _ring_nitrogen_substituent(mol, atom, nitrogen):
    """A neutral three-bonded ring nitrogen joined to `nitrogen` by a single bond whose other neighbors are ring atoms:
    the N-substituent is a 'piperidin-1-yl'- or '1H-pyrazol-1-yl'-type ring group, since the ring nitrogen cannot also
    be the group's own nitrogen (P-62.2.1.3)."""
    if atom.GetAtomicNum() != 7 or not atom.IsInRing() or atom.GetFormalCharge() or atom.GetDegree() != 3:
        return False
    if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False
    return all(n.GetIdx() == nitrogen.GetIdx() or n.IsInRing() for n in atom.GetNeighbors())


def _amidine_ylidene(mol, atom, nitrogen):
    """An acyclic carbon =N- whose other neighbours are amino groups with plain carbon substituents, a plain carbon or
    hydrogen: the N-acyl amidine or guanidine drawn on its imine nitrogen, cited as a (diaminomethylidene)-type prefix on
    the amide (P-41, P-66.4.1.2.1.3)."""
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    if not rest or any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in rest):
        return False
    amino = [n for n in rest if n.GetAtomicNum() == 7]
    carbon = [n for n in rest if n.GetAtomicNum() == 6 and not is_functional_carbon(mol, n.GetIdx())]
    return (
        len(amino) + len(carbon) == len(rest)
        and bool(amino)
        and all(not n.GetFormalCharge() and not n.IsInRing() and _plain_nitrogen_substituents(mol, n, atom.GetIdx()) for n in amino)
    )


def _carbodiimide_ylidene(mol, atom, nitrogen):
    """The central carbon of a carbodiimide drawn on `nitrogen`, cited as an [(organylimino)methylidene] prefix on the
    amide (P-62.3.1.4, P-41)."""
    rest = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    return (
        len(rest) == 1
        and rest[0].GetAtomicNum() == 7
        and mol.GetBondBetweenAtoms(atom.GetIdx(), rest[0].GetIdx()).GetBondTypeAsDouble() == 2.0
        and not rest[0].GetFormalCharge()
        and not rest[0].IsInRing()
        and all(
            mol.GetBondBetweenAtoms(rest[0].GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and (
                n.GetAtomicNum() in (8, 14, 16, *HALOGEN_PREFIXES)
                or n.GetAtomicNum() == 6
            )
            for n in rest[0].GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
        )
    )


def _hydride_group_on_nitrogen(mol, atom, nitrogen, sibling_count):
    """An N-substituent that is a mononuclear hydride group (phosphanyl, silyl, boranyl ...) or, on a nitrogen with no
    other substituent, an ylidene of a hydride atom or of a plain carbon (P-66.1.1.4.3, P-68.1)."""
    from ._hetero_prefixes import MONONUCLEAR_HYDRIDES

    if atom.GetFormalCharge() or atom.GetIsotope():
        return False
    order = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble()
    if atom.IsInRing():
        return order == 2.0 and sibling_count == 1 and atom.GetAtomicNum() == 6
    if order == 1.0:
        return atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES and atom.GetAtomicNum() != 7
    if order == 2.0 and sibling_count == 1:
        if atom.GetAtomicNum() == 6:
            return not _double_oxygens(mol, atom.GetIdx()) and (
                not is_functional_carbon(mol, atom.GetIdx())
                or _amidine_ylidene(mol, atom, nitrogen)
                or _carbodiimide_ylidene(mol, atom, nitrogen)
            )
        return atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES and atom.GetAtomicNum() != 7
    return False


def _ring_nitrogen_acyl(mol, nitrogen, carbonyl):
    """A neutral ring nitrogen without hydrogen whose only acyl group is `carbonyl`, in a saturated or mancude ring:
    the carbonyl is a pseudoketone ('hidden amide', P-66.1.3, P-64.1.2.1(b))."""
    if nitrogen.GetFormalCharge() or not nitrogen.IsInRing() or nitrogen.GetDegree() != 3:
        return False
    for n in nitrogen.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx())
        if n.GetIdx() == carbonyl:
            if bond.GetBondTypeAsDouble() != 1.0:
                return False
        elif bond.GetBondTypeAsDouble() == 1.5:
            if not bond.IsInRing():
                return False
        elif bond.GetBondTypeAsDouble() != 1.0 or not (
            (n.GetAtomicNum() == 6 and (not is_functional_carbon(mol, n.GetIdx()) or _ring_skeleton_carbon(mol, n, nitrogen)))
            or (n.GetAtomicNum() in _RING_HETERO_NEIGHBORS and n.IsInRing() and not n.GetFormalCharge())
        ):
            return False
    return True


def _ring_skeleton_carbon(mol, carbon, nitrogen):
    """A ring carbon whose only multiple bonds to heteroatoms lie inside the ring of `nitrogen` (C2 of a 4,5-dihydro-1H-imidazole):
    it is a member of the ring parent, not a group of its own."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if nitrogen.GetIdx() in r and carbon.GetIdx() in r]
    return bool(rings) and all(
        b.GetBondTypeAsDouble() == 1.0
        or b.GetOtherAtom(carbon).GetAtomicNum() == 6
        or any(b.GetOtherAtomIdx(carbon.GetIdx()) in ring for ring in rings)
        for b in carbon.GetBonds()
    )


def _sulfonyl_group(mol, s_idx, attached):
    """("sulfonic" | "sulfonamide", owned atoms) for an S(=O)(=O)X group whose
    X is OH or an amine nitrogen and whose other neighbor is `attached`."""
    sulfur = mol.GetAtomWithIdx(s_idx)
    group = acid_group_at(mol, s_idx)
    if group is not None and any(n.GetIdx() == attached for n in sulfur.GetNeighbors()):
        if group.spec.kind[0] in ("S", "Se", "Te") and not (group.spec.kind == ("S", 2) and group.spec.plain):
            return group.spec.key, group.owned | {s_idx}
    if sulfur.GetAtomicNum() not in (16, 34, 52) or sulfur.GetFormalCharge():
        return None
    oxygens = [
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 8
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    imides = [
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 7
        and not n.GetFormalCharge()
        and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        and all(m.GetIdx() == s_idx or (m.GetAtomicNum() == 6 and not is_functional_carbon(mol, m.GetIdx())) for m in n.GetNeighbors())
    ]
    replaced = [
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() in _REPLACED_INFIX
        and n.GetDegree() == 1
        and not n.GetFormalCharge()
        and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    hydrazones = {
        n.GetIdx(): end
        for n in sulfur.GetNeighbors()
        if (end := _sulfonyl_hydrazone_end(mol, s_idx, n)) is not None
    }
    rest = [
        n
        for n in sulfur.GetNeighbors()
        if n.GetIdx() not in {*oxygens, *imides, *replaced, *hydrazones} and n.GetIdx() != attached
    ]
    if len(rest) != 1 or sulfur.GetDegree() != len(oxygens) + len(imides) + len(replaced) + len(hydrazones) + 2:
        return None
    other = rest[0]
    if sulfur.GetAtomicNum() == 16 and len(oxygens) == 2 and not imides and not replaced and not hydrazones:
        if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
            return "sulfonic", {s_idx, *oxygens, other.GetIdx()}
    if other.GetAtomicNum() != 7 or not 1 <= len(oxygens) + len(replaced) + len(imides) + len(hydrazones) <= 2:
        return None
    symbols = (
        ["O"] * len(oxygens)
        + [mol.GetAtomWithIdx(i).GetSymbol() for i in replaced]
        + ["NH"] * len(imides)
        + ["NNH2"] * len(hydrazones)
    )
    oxo = make_spec(SULFONYL_CENTRES[sulfur.GetAtomicNum()], symbols, ("O",)).oxo
    owned = {s_idx, *oxygens, *replaced, *imides, *hydrazones, *hydrazones.values(), other.GetIdx()}
    if _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx):
        return SULFONYL_GROUP_NAMES[(sulfur.GetAtomicNum(), oxo, "amide")], owned
    beta = _hydrazide_beta_nitrogen(mol, other, s_idx)
    if beta is not None:
        return SULFONYL_GROUP_NAMES[(sulfur.GetAtomicNum(), oxo, "hydrazide")], {*owned, beta}
    return None


def _sulfonyl_hydrazone_end(mol, s_idx, nitrogen):
    """The terminal nitrogen of an =N-N group double-bonded to the sulfonyl centre `s_idx` (a hydrazono position), else None."""
    if (
        nitrogen.GetAtomicNum() != 7
        or nitrogen.GetFormalCharge()
        or nitrogen.GetDegree() != 2
        or mol.GetBondBetweenAtoms(s_idx, nitrogen.GetIdx()).GetBondTypeAsDouble() != 2.0
    ):
        return None
    end = next((m for m in nitrogen.GetNeighbors() if m.GetIdx() != s_idx), None)
    if (
        end is None
        or end.GetAtomicNum() != 7
        or end.GetFormalCharge()
        or end.IsInRing()
        or end.GetIsAromatic()
        or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), end.GetIdx()).GetBondTypeAsDouble() != 1.0
        or not _plain_nitrogen_arm(mol, end, {nitrogen.GetIdx()}, True)
    ):
        return None
    return end.GetIdx()


def _chalcogen_amide_group(mol, atom, chalcogen, z):
    """(class, owned atoms) of a C(=X)-N carbon, X = S, Se or Te: a thio-, seleno- or telluroamide (P-66.1.4.1), or the
    thione of a 'hidden' amide whose nitrogen is a ring member (P-66.1.4.3)."""
    carbon = atom.GetIdx()
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != chalcogen]
    nitrogens = [n for n in others if n.GetAtomicNum() == 7]
    rest = [n for n in others if n.GetAtomicNum() != 7]
    if any(mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in others):
        return None
    ring_nitrogens = [n for n in nitrogens if n.IsInRing()]
    if len(ring_nitrogens) == 2 and not rest and all(_ring_nitrogen_acyl(mol, n, carbon) for n in ring_nitrogens):
        return _CHALCOGEN_KETONE_CLASS[z], {chalcogen}
    if len(nitrogens) == 2 and len(ring_nitrogens) == 1 and not rest:
        other = next(n for n in nitrogens if not n.IsInRing())
        if _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, carbon):
            return _CHALCOGEN_AMIDE[z], {chalcogen, other.GetIdx()}
        return None
    if len(nitrogens) != 1:
        return None
    if len(rest) > 1 or (rest and rest[0].GetAtomicNum() != 6) or (not rest and atom.GetTotalNumHs() != 1):
        return None
    (nitrogen,) = nitrogens
    if rest and _ring_nitrogen_acyl(mol, nitrogen, carbon):
        return _CHALCOGEN_KETONE_CLASS[z], {chalcogen}
    if _terminal_heteroatom(mol, nitrogen.GetIdx(), 2) or _plain_amide_nitrogen(mol, nitrogen, carbon):
        return _CHALCOGEN_AMIDE[z], {chalcogen, nitrogen.GetIdx()}
    return None


_GROUP_14_ATOMS = {14, 32, 50, 82}
_GROUP_15_ATOMS = {15, 33, 51, 83}


def _diacyl_chalcogen_chain(mol, carbon, first):
    """Whether the acyl `carbon` is joined through `first` to a carbon group by a chain of three or more chalcogen atoms:
    a pseudoketone, not an ester, since no acid has such a chain (P-64.1.2.1, P-65.7.5.1, P-68.4.2.4)."""
    if _homogeneous_chalcogen_chain(mol, carbon, first):
        return True
    chain, previous, atom = [], carbon, first
    while atom.GetAtomicNum() in (8, 16, 34, 52):
        if atom.GetFormalCharge() or atom.IsInRing() or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
            return False
        chain.append(atom.GetIdx())
        onward = [n for n in atom.GetNeighbors() if n.GetIdx() != previous]
        if not onward:
            return len(chain) >= 3 and atom.GetTotalNumHs() == 1
        if len(onward) != 1:
            return False
        previous, atom = atom.GetIdx(), onward[0]
    return len(chain) >= 3 and atom.GetAtomicNum() == 6 and atom.GetIdx() != carbon


def _homogeneous_chalcogen_chain(mol, carbon, first):
    """Three or more identical chalcogen atoms between the acyl `carbon` and a carbon group or hydrogen: the
    compound is a pseudoketone (P-68.4.1.3)."""
    element = first.GetAtomicNum()
    if element not in (8, 16, 34, 52):
        return False
    length, previous, atom = 0, carbon, first
    while atom.GetAtomicNum() == element:
        if atom.GetFormalCharge() or atom.IsInRing() or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
            return False
        length += 1
        onward = [n for n in atom.GetNeighbors() if n.GetIdx() != previous]
        if not onward:
            return length >= 3 and atom.GetTotalNumHs() == 1
        if len(onward) != 1:
            return False
        previous, atom = atom.GetIdx(), onward[0]
    return length >= 3 and atom.GetAtomicNum() == 6


def _is_pseudoketone_heteroatom(mol, atom, carbon):
    """A neutral Group 14 or 15 atom (P-64.1.2.1) that makes its acyl carbon a pseudoketone rather than a member of
    a senior class: no multiple bonds on it, and for Group 15 no oxygen, nitrogen or halogen (phosphinous-type acids)."""
    z = atom.GetAtomicNum()
    if atom.GetFormalCharge() or z not in _GROUP_14_ATOMS | _GROUP_15_ATOMS:
        return False
    if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
        return False
    if any(_double_oxygens(mol, n.GetIdx()) for n in atom.GetNeighbors() if n.GetIdx() != carbon):
        return False
    return z in _GROUP_14_ATOMS or all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors() if n.GetIdx() != carbon)


def _acyl_diazene_nitrogen(mol, atom, carbon):
    """The acylated nitrogen of an acyl diazene R-CO-N=N-R': a ketone with a diazenyl prefix, not an amide
    (P-68.3.1.3.6, P-15.3.3.2.2)."""
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() or atom.GetDegree() != 2 or atom.IsInRing():
        return False
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != carbon]
    return (
        len(others) == 1
        and others[0].GetAtomicNum() == 7
        and not others[0].GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), others[0].GetIdx()).GetBondTypeAsDouble() == 2.0
        and mol.GetBondBetweenAtoms(atom.GetIdx(), carbon).GetBondTypeAsDouble() == 1.0
    )


def _acyloxy_amine_oxygen(mol, oxygen, carbon):
    """The ester oxygen of an acyl-O-N group on an acyclic amine or oxime nitrogen: a pseudoketone (P-65.6.3.4), since only a
    cyclic nitrogen gives a traditional ester."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or oxygen.GetFormalCharge():
        return False
    far = next((n for n in oxygen.GetNeighbors() if n.GetIdx() != carbon), None)
    if far is None or far.GetAtomicNum() != 7 or far.GetFormalCharge() or far.IsInRing() or far.GetIsAromatic():
        return False
    doubles = [b for b in far.GetBonds() if b.GetBondTypeAsDouble() != 1.0]
    if len(doubles) > 1 or any(b.GetBondTypeAsDouble() != 2.0 for b in doubles):
        return False
    return all(n.GetIdx() == oxygen.GetIdx() or n.GetAtomicNum() == 6 for n in far.GetNeighbors())


def _acyloxy_on_nitrogen(mol, oxygen, nitrogen):
    """An acyloxy group -O-CO-R on `nitrogen`, cited as the prefix 'acyloxy' because an acyclic nitrogen gives no
    traditional ester (P-65.6.3.2.3, P-65.6.3.4.1)."""
    acyl = next((n for n in oxygen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()), None)
    return (
        acyl is not None
        and acyl.GetAtomicNum() in (6, 16)
        and bool(_double_oxygens(mol, acyl.GetIdx()))
        and _acyloxy_amine_oxygen(mol, oxygen, acyl.GetIdx())
    )


def _sulfonyl_groups_on(mol, carbon):
    """[(class, owned atoms)] of the sulfonyl-type groups bonded to `carbon`."""
    found = []
    for n in mol.GetAtomWithIdx(carbon).GetNeighbors():
        if n.GetAtomicNum() in (16, 34, 52) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0:
            sulfonyl = _sulfonyl_group(mol, n.GetIdx(), carbon)
            if sulfonyl is not None:
                found.append(sulfonyl)
    return found


def _hydroxy_type_groups(mol, carbon):
    """[(class, terminal atoms)] of the hydroxy, chalcogenol and amino groups on `carbon`, geminal ones together; in an
    anion name only the deprotonated ones form the suffix, a neutral one beside them is a prefix (P-72.2.2.2.2)."""
    found = []
    for z, hydrogens, name in ((8, 1, "alcohol"), (16, 1, "thiol"), (34, 1, "selenol"), (52, 1, "tellurol"), (7, 2, "amine")):
        atoms = [
            n
            for n in _single_neighbors(mol, carbon, z)
            # in an aminium name only the cationic nitrogens are the suffix; a neutral amino group is a prefix
            if _terminal_heteroatom(mol, n, hydrogens)
            and not (AMINIUM.get() is True and z == 7 and not mol.GetAtomWithIdx(n).HasProp("_cationic_amine"))
        ]
        found.append((name, [n for n in atoms if mol.GetAtomWithIdx(n).HasProp(ANION_PROP)] or atoms))
    return found


def _group_weight(mol, carbon, cls):
    """How many groups of class `cls` sit on `carbon` (two sulfonic acid groups on one carbon: methanedisulfonic acid,
    two amino groups: methanediamine)."""
    geminal = next((len(found) for name, found in _hydroxy_type_groups(mol, carbon) if name == cls), 0)
    return sum(1 for c, _ in _sulfonyl_groups_on(mol, carbon) if c == cls) or geminal or 1


def _is_ketene_carbon(mol, atom, oxygen):
    """O=C=C: a ketone-class carbon that ends the chain, its single carbon neighbour joined by the cumulated C=C (P-64.2.2.4)."""
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != oxygen]
    return (
        len(others) == 1
        and others[0].GetAtomicNum() == 6
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), others[0].GetIdx()).GetBondTypeAsDouble() == 2.0
    )


def _group_of(mol, carbon):
    """(class, atoms owned by the group) for a principal-capable group on
    `carbon`, else None. Raises on carbon-bound groups this engine cannot
    name (esters, acid halides, ...)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.HasProp(JUNIOR_GROUP):
        return None
    variant = acid_group_at(mol, carbon)
    if variant is not None and not variant.spec.plain:
        return variant.spec.key, variant.owned
    nitrogens = [
        n.GetIdx()
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 3.0
    ]
    if nitrogens:
        if atom.GetFormalCharge() == -1 and atom.GetDegree() == 1 and mol.GetAtomWithIdx(nitrogens[0]).GetFormalCharge() == 1:
            return None
        others = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogens[0]]
        if (
            len(others) == 1
            and others[0].GetAtomicNum() in (8, 16, 34, 52)
            and others[0].GetDegree() == 2
            and mol.HasSubstructMatch(_CARBOXYLIC_OR_SULFONIC)
        ):
            return None
        if (
            len(others) == 1
            and others[0].GetAtomicNum() in (16, 34, 52)
            and others[0].GetDegree() >= 3
            and mol.HasSubstructMatch(_CARBOXYLIC_OR_SULFONIC)
        ):
            return None
        if len(others) == 1 and others[0].GetAtomicNum() == 7 and not others[0].GetFormalCharge() and not others[0].IsInRing():
            return None
        if (
            len(others) == 1
            and others[0].GetAtomicNum() == 7
            and not others[0].GetFormalCharge()
            and any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_NITRILE)
        ):
            return None
        if (
            len(others) == 1
            and others[0].GetAtomicNum() in (8, 16, 34, 52)
            and (AMINIUM.get() or RING_CENTER.get() or any(a.HasProp(ANION_PROP) for a in mol.GetAtoms()))
        ):
            return None
        if len(others) != 1 or others[0].GetAtomicNum() != 6:
            raise UnsupportedStructure("a cyanide not bonded to carbon is not a nitrile")
        return "nitrile", {nitrogens[0]}
    for z, thione in _CHALCOGEN_KETONE_CLASS.items():
        found = [
            n.GetIdx()
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == z and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if found and all(
            n.GetAtomicNum() == 6
            or _acyl_diazene_nitrogen(mol, n, carbon)
            or (n.GetAtomicNum() in (8, 16, 34, 52) and _diacyl_chalcogen_chain(mol, carbon, n))
            for n in atom.GetNeighbors()
            if n.GetIdx() != found[0]
        ) and atom.GetDegree() == 3:
            return thione, {found[0]}
        if found:
            amide = _chalcogen_amide_group(mol, atom, found[0], z)
            if amide is not None:
                return amide
    if _is_isocyanate_carbon(atom):
        return None
    oxygens = _double_oxygens(mol, carbon)
    if oxygens:
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        hetero = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 6 and n.GetIdx() != oxygens[0]]
        if any(
            mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() != 1.0
            for n in atom.GetNeighbors()
            if n.GetIdx() != oxygens[0]
        ):
            if any(n.IsInRing() for n in carbon_neighbors):
                return None
            if _is_ketene_carbon(mol, atom, oxygens[0]):
                return "ketone", {oxygens[0]}
            raise UnsupportedStructure("a cumulated carbonyl (ketene) is not supported")
        if not hetero:
            if len(carbon_neighbors) == 1 and atom.GetTotalNumHs() == 1:
                return ("ketone" if atom.HasProp(ACYL_ROOT) else "aldehyde"), {oxygens[0]}
            if len(carbon_neighbors) == 2:
                return "ketone", {oxygens[0]}
            raise UnsupportedStructure("a formaldehyde-type carbonyl is not named by the chain engine")
        if len(hetero) == 2:
            hydroxyl = [h for h in hetero if h.GetAtomicNum() == 8 and _terminal_heteroatom(mol, h.GetIdx(), 1)]
            nitro = [h for h in hetero if h.GetAtomicNum() == 7 and is_nitro_nitrogen(mol, h.GetIdx())]
            if hydroxyl and nitro and not carbon_neighbors:
                return "acid", {oxygens[0], hydroxyl[0].GetIdx()}
        if hetero and all(_acyl_diazene_nitrogen(mol, h, carbon) for h in hetero):
            return "ketone", {oxygens[0]}
        ring_nitrogens = [h for h in hetero if h.IsInRing() and h.GetAtomicNum() == 7]
        if len(ring_nitrogens) == 2 and len(hetero) == 2 and not carbon_neighbors:
            if all(_ring_nitrogen_acyl(mol, h, carbon) for h in ring_nitrogens):
                return "ketone", {oxygens[0]}
        ring_centres = [h for h in hetero if h.IsInRing() and h.GetAtomicNum() in _RING_ACID_CENTRES and not h.GetFormalCharge()]
        if len(hetero) == 2 and len(ring_centres) == 1 and not carbon_neighbors:
            other = next(h for h in hetero if h.GetIdx() != ring_centres[0].GetIdx())
            if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
                return "acid", {oxygens[0], other.GetIdx()}
        if len(hetero) == 2 and len(ring_nitrogens) == 1:
            other = next(h for h in hetero if h.GetIdx() != ring_nitrogens[0].GetIdx())
            if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
                return "acid", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and (
                _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, carbon)
            ):
                return "amide", {oxygens[0], other.GetIdx()}
            beta = _hydrazide_beta_nitrogen(mol, other, carbon) if other.GetAtomicNum() == 7 else None
            if beta is not None:
                return "hydrazide", {oxygens[0], other.GetIdx(), beta}
            return None
        if len(hetero) == 1:
            other = hetero[0]
            if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
                return "acid", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and _terminal_heteroatom(mol, other.GetIdx(), 2):
                return "amide", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and _plain_amide_nitrogen(mol, other, carbon):
                return "amide", {oxygens[0], other.GetIdx()}
            if (
                other.GetAtomicNum() == 7
                and AMINIUM.get() == "amide"
                and other.GetDegree() == 1
                and other.GetTotalNumHs() == 3
                and other.GetFormalCharge() == 1
            ):
                return "amide", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and len(carbon_neighbors) <= 1:
                beta = _hydrazide_beta_nitrogen(mol, other, carbon)
                if beta is not None:
                    return "hydrazide", {oxygens[0], other.GetIdx(), beta}
            if carbon_neighbors and other.GetAtomicNum() == 7 and _ring_nitrogen_acyl(mol, other, carbon):
                return "ketone", {oxygens[0]}
            if carbon_neighbors and (
                _is_pseudoketone_heteroatom(mol, other, carbon)
                or _diacyl_chalcogen_chain(mol, carbon, other)
                or _acyloxy_amine_oxygen(mol, other, carbon)
            ):
                return "ketone", {oxygens[0]}
            if not carbon_neighbors and atom.GetTotalNumHs() == 1 and other.GetAtomicNum() == 7 and _ring_nitrogen_acyl(mol, other, carbon):
                return "aldehyde", {oxygens[0]}
        return None
    amidine = _amidine_nitrogens(mol, atom)
    if amidine is not None:
        return "amidine", set(amidine)
    amidrazone = _amidrazone_group(mol, atom)
    if amidrazone is not None:
        return amidrazone
    imine = _acyclic_imine_nitrogen(mol, atom)
    if imine is not None:
        return "imine", {imine}
    sulfonyls = _sulfonyl_groups_on(mol, carbon)
    if sulfonyls:
        same = [owned for cls, owned in sulfonyls if cls == sulfonyls[0][0]]
        return sulfonyls[0][0], set().union(*same)
    for name, found in _hydroxy_type_groups(mol, carbon):
        if found:
            return name, set(found)
    for n in mol.GetAtomWithIdx(carbon).GetNeighbors():
        if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0:
            found = _peroxol_group(mol, n, carbon)
            if found is not None:
                return found
    return None


def _peroxol_group(mol, bridge, anchor):
    """(class, owned atoms) of the -X-YH group of two chalcogens whose first atom `bridge` is bonded to `anchor`: the
    peroxol of Table 6.1 or its chalcogen analogue (P-63.4.1, P-63.4.2.1)."""
    ends = [m for m in bridge.GetNeighbors() if m.GetIdx() != anchor]
    if bridge.GetAtomicNum() not in _PEROXOL_ELEMENTS or bridge.GetDegree() != 2 or bridge.GetTotalNumHs() or len(ends) != 1:
        return None
    end = ends[0]
    if end.GetAtomicNum() not in _PEROXOL_ELEMENTS or not _terminal_heteroatom(mol, end.GetIdx(), 1):
        return None
    if mol.GetBondBetweenAtoms(bridge.GetIdx(), end.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    return peroxol_suffix(bridge.GetAtomicNum(), end.GetAtomicNum()), {bridge.GetIdx(), end.GetIdx()}


def _paths(adj, eligible):
    ends = [a for a in eligible if sum(1 for n in adj[a] if n in eligible) <= 1] or list(eligible)
    paths = []
    for start in ends:
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            onward = [n for n in adj[node] if n in eligible and n not in path]
            if not onward:
                paths.append(path)
            for n in onward:
                stack.append((n, path + [n]))
    return paths


_ACID_OR_CARBOXYLIC_ESTER = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=[O,S,Se,Te])[OX2H1,OX1-,OX2,SX2,SeX2,TeX2]",
        "[#16X4](=O)(=O)([#6])[OX2H1,OX1-]",
        "[#15X4](=O)([#6,#7])[OX2H1,OX1-]",
    )
]


def name_polyfunctional(mol) -> str:
    from ._inorganic_ester_parent import has_inorganic_ester_parent
    from ._isotope_labels import split_isotopes
    from ._linear_phane import has_linear_phane_shape, name_linear_phane
    from ._radical_group import has_radical_group_shape, name_radical_group

    if has_inorganic_ester_parent(mol) and not any(mol.HasSubstructMatch(q) for q in _ACID_OR_CARBOXYLIC_ESTER):
        raise UnsupportedStructure("an ester of a noncarbon acid is the parent and is named by its own namer")
    if has_linear_phane_shape(mol):
        return name_linear_phane(mol)
    if has_radical_group_shape(mol):
        return name_radical_group(mol)
    split = split_isotopes(mol)
    anionic_parent = any(a.HasProp(ANION_PROP) for a in mol.GetAtoms())
    if split is not None:
        ammonium = None if anionic_parent else _aminium_base(split[0])
        if ammonium is not None:
            return _name_isotopic(mol, *split, build=lambda clean, labels: _name_aminium(ammonium, labels))
        center = _ring_center_base(split[0])
        if center is not None:
            return _name_isotopic(mol, *split, build=lambda clean, labels: _name_labelled_ring_center(clean, labels, center))
        return _name_isotopic(mol, *split)
    mixed = None if anionic_parent else _cation_parent_base(mol)
    if mixed is not None:
        with ion_stereo_carried(mol, mixed[0]):
            return _name_aminium(mixed[0], parent=mixed[1], cationic=mixed[2])
    cation = None if anionic_parent else _aminium_base(mol)
    if cation is not None:
        with ion_stereo_carried(mol, cation):
            return _name_aminium(cation)
    amidinium = None if anionic_parent else _amidinium_base(mol)
    if amidinium is not None:
        with ion_stereo_carried(mol, amidinium):
            name = _name_labelled(amidinium, {})
        return _cationic_group_suffix(name, "amidine", sum(a.GetFormalCharge() == 1 for a in mol.GetAtoms()))
    iminium = None if anionic_parent else _iminium_base(mol)
    if iminium is not None:
        with ion_stereo_carried(mol, iminium):
            return _name_aminium(iminium, parent="imine", cationic=sum(1 for a in mol.GetAtoms() if a.GetFormalCharge()))
    for kind in ("amide", "nitrile"):
        acylated = _group_cation_base(mol, kind)
        if acylated is not None:
            with ion_stereo_carried(mol, acylated):
                return _name_aminium(acylated, parent=kind)
    grouped = None if anionic_parent else _ring_cation_with_group(mol)
    if grouped is not None:
        return _name_ring_group_cation(*grouped)
    center = _ring_center_base(mol)
    if center is not None:
        with _oxide_stereo_carried(mol, center[0]):
            return _name_ring_center(*center)
    carbocation = _carbocation_base(mol)
    if carbocation is not None:
        base, atom = carbocation
        return _name_labelled(base, {}, lambda name, placed: _ylium_inserted(name, base, atom, placed))
    return _name_labelled(mol, {})


def _ring_center_base(mol):
    """(neutral analogue, center atom in it, kind, probe) for a single ring-nitrogen cation or N-oxide; the cation
    (P-73.1.1.2) or the N-oxide (P-62.5, P-74.2.1.2) outranks every other group, so the analogue is named with all of
    them as prefixes. `probe` checks the centre's locant when the ring has another nitrogen."""
    from ._diester_ring_diyl import _system_of

    anionic_parent = any(_anionic_group_atom(mol, a) or a.HasProp(ANION_PROP) for a in mol.GetAtoms())
    charged = [
        a
        for a in mol.GetAtoms()
        if a.GetFormalCharge()
        and (not _is_nitro_part(a) or _is_ring_azoxy_part(a))
        and not _anionic_group_atom(mol, a)
        and not _is_chain_amine_oxide_part(a)
        and not (anionic_parent and _cationic_prefix_nitrogen(a))
    ]
    centers = [a for a in charged if a.GetAtomicNum() == 7 and a.IsInRing() and a.GetFormalCharge() == 1]
    if len(centers) > 1:
        return _ring_oxide_centres(mol, centers, charged)
    if len(centers) != 1 or len(charged) - len(centers) > 1:
        return None
    center = centers[0]
    oxides = [a for a in charged if a is not center]
    if oxides and not (
        oxides[0].GetAtomicNum() == 8
        and oxides[0].GetFormalCharge() == -1
        and oxides[0].GetDegree() == 1
        and mol.GetBondBetweenAtoms(oxides[0].GetIdx(), center.GetIdx()) is not None
    ):
        return None
    if mol.GetNumAtoms() > _MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        return None
    system = _system_of(mol, center.GetIdx())
    if system is None:
        return None
    other_nitrogens = [i for i in system[1] if i != center.GetIdx() and mol.GetAtomWithIdx(i).GetAtomicNum() == 7]
    substituted = center.GetDegree() - sum(1 for n in center.GetNeighbors() if n.GetIdx() in system[1]) - len(oxides) > 0
    probe = None
    if other_nitrogens and not substituted:
        probe = _probe_with_marker(mol, center, oxides)
        if probe is None:
            return None
    neutral = Chem.RWMol(mol)
    index = center.GetIdx()
    atom = neutral.GetAtomWithIdx(index)
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(True)
    hydrogens = center.GetTotalNumHs() if oxides else max(center.GetTotalNumHs() - 1, 0)
    atom.SetNumExplicitHs(hydrogens)
    if oxides:
        oxide = oxides[0].GetIdx()
        neutral.RemoveAtom(oxide)
        index -= index > oxide
    base = neutral.GetMol()
    if oxides or center.GetTotalNumHs():
        Chem.SanitizeMol(base)
    else:
        base.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(base)
        if center.GetIsAromatic():
            base.GetAtomWithIdx(index).SetBoolProp("_ring_cation_centre", True)
    return base, index, "oxide" if oxides else "ium", probe


def _oxidized_set_is_fixed(mol, members, centre_atoms):
    """Whether every symmetry of the ring atoms (matched by element) maps the oxidized nitrogens onto themselves, so that
    all permitted numberings give them the same locants."""
    from networkx.algorithms import isomorphism

    graph = nx.Graph()
    graph.add_nodes_from((i, {"z": mol.GetAtomWithIdx(i).GetAtomicNum()}) for i in members)
    graph.add_edges_from((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds() if {b.GetBeginAtomIdx(), b.GetEndAtomIdx()} <= members)
    matcher = isomorphism.GraphMatcher(graph, graph, node_match=lambda a, b: a["z"] == b["z"])
    return all({m[c] for c in centre_atoms} == centre_atoms for m in matcher.isomorphisms_iter())


def _ring_oxide_centres(mol, centers, charged):
    """(neutral analogue, [centre atoms in it], 'oxide', None) for several N-oxides on ring nitrogens whose set the symmetry of
    the ring atoms cannot change, so that the parent numbering alone fixes the locants (P-62.5, P-74.2.1.2); None otherwise."""
    if len(charged) != 2 * len(centers) or mol.GetNumAtoms() > _MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        return None
    oxides = []
    charged_atoms = {a.GetIdx() for a in charged}
    for center in centers:
        found = [
            n
            for n in center.GetNeighbors()
            if n.GetFormalCharge() == -1 and n.GetAtomicNum() == 8 and n.GetDegree() == 1 and n.GetIdx() in charged_atoms
        ]
        if len(found) != 1:
            return None
        oxides.append(found[0].GetIdx())
    members = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    centre_atoms = {c.GetIdx() for c in centers}
    if not centre_atoms <= members or not _oxidized_set_is_fixed(mol, members, centre_atoms):
        return None
    neutral = Chem.RWMol(mol)
    for center in centers:
        atom = neutral.GetAtomWithIdx(center.GetIdx())
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(center.GetTotalNumHs())
    for oxide in sorted(oxides, reverse=True):
        neutral.RemoveAtom(oxide)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    return base, [c.GetIdx() - sum(1 for o in oxides if o < c.GetIdx()) for c in centers], "oxide", None


_CENTRE_MARKERS = ((35, "bromo"), (17, "chloro"), (9, "fluoro"), (53, "iodo"))


def _probe_with_marker(mol, center, oxides):
    """(marked ring cation, its centre, marker prefix): the same ring cation with a halogen on the centre. Its name
    numbers the centre lowest, as the name of an N-oxide or protonated ring (no substituent on the centre) must; the
    marker is the alphabetically first halogen absent from the molecule, so that it is cited first on a tie, and its
    prefix is cut out afterwards."""
    present = {a.GetAtomicNum() for a in mol.GetAtoms()}
    marker = next(((z, word) for z, word in _CENTRE_MARKERS if z not in present), None)
    if marker is None:
        return None
    probe = Chem.RWMol(mol)
    if oxides:
        halogen = probe.GetAtomWithIdx(oxides[0].GetIdx())
        halogen.SetAtomicNum(marker[0])
        halogen.SetFormalCharge(0)
    else:
        added = probe.AddAtom(Chem.Atom(marker[0]))
        probe.AddBond(center.GetIdx(), added, Chem.BondType.SINGLE)
        probe.GetAtomWithIdx(center.GetIdx()).SetNumExplicitHs(max(center.GetTotalNumHs() - 1, 0))
        probe.GetAtomWithIdx(center.GetIdx()).SetNoImplicit(True)
    try:
        probe = probe.GetMol()
        Chem.SanitizeMol(probe)
    except Exception:
        return None
    found = _ring_center_base(probe)
    return None if found is None else (*found[:2], marker[1])


def _without_marker(name, locant, word):
    token = f"{locant}-{word}"
    at = name.find(token)
    if at < 0:
        raise UnsupportedStructure("the centre marker is not cited in the name")
    text = name[:at] + name[at + len(token):]
    if at == 0:
        return text.lstrip("-")
    if text[at - 1:at] == "-" and text[at:at + 1] == "-":
        return text[:at] + text[at + 1:]
    if text[at - 1:at] == "-" and text[at:at + 1].isalpha() and text[at - 2:at - 1].isalpha():
        return text[:at - 1] + text[at:]
    return text


def _attach_oxides(name, placed, centers):
    from ._pin import mark
    from ._substituents import _locant_sort_key

    locants = [placed.get(c) for c in centers]
    if None in locants or " " in name:
        raise UnsupportedStructure("the ring nitrogens are not numbered in the parent hydride")
    word = {2: "dioxide", 3: "trioxide", 4: "tetraoxide"}.get(len(centers))
    if word is None:
        raise UnsupportedStructure("this many N-oxides are not named here")
    text = ",".join(str(x) for x in sorted(locants, key=_locant_sort_key))
    return mark(f"{name} {text}-{word}", "the Blue Book generates the preferred name of two or more amine oxides by the lambda-convention (P-74.2.1.2)")


def _attach(name, placed, center, kind):
    locant = placed.get(center)
    if locant is None or " " in name:
        raise UnsupportedStructure("the ring nitrogen is not numbered in the parent hydride")
    if kind == "oxide":
        return f"{name} {locant}-oxide"
    if not name.endswith("e"):
        raise UnsupportedStructure("the cationic ring parent has no terminal 'e' to elide (P-73.1.1.2)")
    return f"{name[:-1]}-{locant}-ium"


def _name_labelled_ring_center(clean, labels, found):
    base, center, kind, probe = found
    if isinstance(center, list):
        raise UnsupportedStructure("an isotopic modification of several N-oxides is not supported yet")
    shift = clean.GetNumAtoms() - base.GetNumAtoms()
    removed = next((a.GetIdx() for a in clean.GetAtoms() if a.GetFormalCharge() == -1 and a.GetAtomicNum() == 8), None)
    if shift and removed is None:
        raise UnsupportedStructure("the isotopic modification of this ring cation is not supported")
    if STEREO_OF_ISOTOPOLOGUE.get() or removed in labels:
        raise UnsupportedStructure("an isotopic modification of the oxide group or the stereo of the centre is not supported yet")
    moved = labels if probe is not None else {i - (1 if shift and i > removed else 0): e for i, e in labels.items()}
    return _name_ring_center(base, center, kind, probe, moved)


def _mancude_ring(ring):
    """`ring` with the hydrogens of its aromatic atoms assigned afresh so that it is a valid mancude structure: the
    hydrogens a ring carries as oxo or hydro positions leave no aromatic ring as it stands."""
    nitrogens = [a.GetIdx() for a in ring.GetAtoms() if a.GetAtomicNum() == 7 and a.GetIsAromatic()]
    for count in range(len(nitrogens) + 1):
        for with_hydrogen in itertools.combinations(nitrogens, count):
            candidate = Chem.RWMol(ring)
            for atom in candidate.GetAtoms():
                if atom.GetIsAromatic():
                    atom.SetNumExplicitHs(1 if atom.GetIdx() in with_hydrogen else 0)
                    atom.SetNoImplicit(atom.GetAtomicNum() == 7)
            mol = candidate.GetMol()
            try:
                Chem.SanitizeMol(mol)
            except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
                continue
            return mol
    raise UnsupportedStructure("the ring of this cation has no mancude structure")


def _cation_inserted(name, base, center, placed_in_parent):
    """Cite the ring cation as 'ium' on its ring parent hydride, ahead of any suffix or free valence (P-74.1.2)."""
    from ._diester_ring_diyl import _system_of

    system = _system_of(base, center)[1]
    locant = placed_in_parent.get(center)
    if locant is None and any(base.GetAtomWithIdx(i).GetAtomicNum() == 7 and i != center for i in system):
        locant = _prefix_center_locant(base, center, placed_in_parent)
    ring = Chem.RWMol(base)
    for index in system:
        lost = sum(
            int(b.GetBondTypeAsDouble()) for b in base.GetAtomWithIdx(index).GetBonds() if b.GetOtherAtomIdx(index) not in system
        )
        if lost and index != center:
            ring.GetAtomWithIdx(index).SetNumExplicitHs(base.GetAtomWithIdx(index).GetTotalNumHs() + lost)
            ring.GetAtomWithIdx(index).SetNoImplicit(True)
    for index in sorted((a.GetIdx() for a in base.GetAtoms() if a.GetIdx() not in system), reverse=True):
        ring.RemoveAtom(index)
    new_center = sorted(system).index(center)
    ring_mol = ring.GetMol()
    ring_mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(ring_mol)
    lambda_token = CITE_SKELETAL_LAMBDA.set(True)
    try:
        try:
            bare, placed = _bare_ring_name(ring_mol, new_center, hydro_ok=True)
        except ValueError:
            bare, placed = _bare_ring_name(_mancude_ring(ring_mol), new_center, hydro_ok=True)
    finally:
        CITE_SKELETAL_LAMBDA.reset(lambda_token)
    if not bare.endswith("e"):
        raise UnsupportedStructure("the cationic ring parent has no terminal 'e' to elide (P-73.1.1.2)")
    stem = bare[:-1]
    name = _hydro_for_cationic_centre(name, stem, locant or placed[new_center])
    pattern = re.compile(re.escape(stem) + r"(?:e(?![a-z])|(?=-?\d*,?\d*-?yl)|(?=-\d+(?:,\d+)*-[a-z]))")
    found = list(pattern.finditer(name))
    if len(found) != 1:
        raise UnsupportedStructure("the ring parent of the cation is not delimited in the anionic name")
    match = found[0]
    return f"{name[:match.start()]}{stem}-{locant or placed[new_center]}-ium{name[match.end():]}"


_HYDRO_WORDS = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
_HYDRO_COUNTS = {word: count for count, word in _HYDRO_WORDS.items()}


def _hydro_for_cationic_centre(name, stem, locant):
    """A ring group attached through its own nitrogen cites the added hydrogen of the neutral group, '3,4-dihydro
    isoquinolin-2(1H)-yl'; the cationic centre has no added hydrogen, so its hydro prefix takes the locant of the added
    hydrogen and of the centre: '1,2,3,4-tetrahydroisoquinolin-2-ium-2-yl' (P-31.1.4.2.4, P-73.1.1.2)."""
    match = re.search(
        rf"(?P<locants>\d+(?:,\d+)*)-(?P<word>di|tri|tetra|penta|hexa|hepta|octa)hydro{re.escape(stem)}-(?P<centre>{re.escape(str(locant))})\((?P<added>\d+)H\)-yl",
        name,
    )
    if match is None:
        every = re.search(
            rf"(?P<word>tri|tetra|penta|hexa|hepta|octa|nona)hydro{re.escape(stem)}-(?P<centre>{re.escape(str(locant))})\((?P<added>\d+)H\)-yl",
            name,
        )
        word = _HYDRO_WORDS.get(_HYDRO_COUNTS[every.group("word")] + 2) if every else None
        if word is None:
            return name
        return f"{name[:every.start()]}{word}hydro{stem}-{every.group('centre')}-yl{name[every.end():]}"
    locants = sorted({*match.group("locants").split(","), match.group("added"), match.group("centre")}, key=int)
    word = _HYDRO_WORDS.get(len(locants))
    if word is None:
        return name
    return f"{name[:match.start()]}{','.join(locants)}-{word}hydro{stem}-{match.group('centre')}-yl{name[match.end():]}"


def _prefix_center_locant(base, center, placed_in_parent):
    """The locant of the cationic centre in the ring group cited as a prefix: the group's own numbering gives the free
    valence and the heteroatoms their lowest locants, so a second ring nitrogen fixes it (P-31.1.4, P-73.1.1.2)."""
    from ._diester_ring_diyl import _system_of, evaluate_skeleton
    from ._free_valence import suffix_of

    rings, atoms = _system_of(base, center)
    joins = [
        (a, n.GetIdx())
        for a in atoms
        for n in base.GetAtomWithIdx(a).GetNeighbors()
        if n.GetIdx() in placed_in_parent
    ]
    if len(joins) != 1:
        raise UnsupportedStructure("the ring cation is not joined to the parent hydride by a single bond")
    root, parent = joins[0]
    suffix = suffix_of(base.GetBondBetweenAtoms(root, parent).GetBondTypeAsDouble())
    named = None if suffix is None else evaluate_skeleton(base, adjacency(base), "ring", rings, atoms, [root], {parent}, suffix)
    if named is None or center not in named[2]:
        raise UnsupportedStructure("the cationic centre is not numbered in the ring group")
    return named[2][center]


def _carbocation_base(mol):
    """(mol with its one carbocation made neutral, that atom) when an acid anion is the only other charge: the anion
    outranks the cation, which stays on the parent hydride as 'ylium' (P-74.1.2)."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and not _anionic_group_atom(mol, a)]
    if len(charged) != 1 or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    cation = charged[0]
    if cation.GetAtomicNum() != 6 or cation.GetFormalCharge() != 1 or cation.GetIsAromatic():
        return None
    if not any(_anionic_group_atom(mol, a) for a in mol.GetAtoms()):
        return None
    neutral = Chem.RWMol(mol)
    atom = neutral.GetAtomWithIdx(cation.GetIdx())
    atom.SetFormalCharge(0)
    atom.SetNumExplicitHs(cation.GetTotalNumHs() + 1)
    atom.SetNoImplicit(True)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    return base, cation.GetIdx()


def _ylium_inserted(name, base, center, placed):
    from ._numerals import alkane_name
    from ._radical_ion_skeleton import _parent_base, _unsaturation_locants

    locant = placed.get(center)
    if locant is None:
        raise UnsupportedStructure("the carbocation is not on the parent hydride of the anionic name (P-74.1.3)")
    parent = set(placed)
    ring = any(base.GetAtomWithIdx(a).IsInRing() for a in parent)
    size = len(parent)
    stem = ("cyclo" if ring else "") + alkane_name(size)[:-3]
    enes, ynes = _unsaturation_locants(base, parent, placed, size, not ring)
    root = _parent_base(stem, enes, ynes)[:-1]
    if name.count(root) != 1:
        raise UnsupportedStructure("the parent hydride of the anionic name is not delimited")
    start = name.index(root) + len(root)
    tail = name[start + 1 :] if name[start : start + 1] == "e" else name[start:]
    if tail and tail[0] != "-":
        anion_carbons = []
        for a in base.GetAtoms():
            if not _anionic_group_atom(base, a):
                continue
            near = [n for n in a.GetNeighbors() if n.GetAtomicNum() == 6]
            near = near or [m for x in a.GetNeighbors() for m in x.GetNeighbors() if m.GetAtomicNum() == 6]
            anion_carbons += sorted({placed[m.GetIdx()] for m in near if m.GetIdx() in placed})
        anion_carbons.sort()
        if not anion_carbons:
            raise UnsupportedStructure("the anionic suffix locants are not found on the parent hydride")
        tail = "-" + ",".join(map(str, dict.fromkeys(anion_carbons))) + "-" + tail
    return f"{name[:start]}-{locant}-ylium{tail}"


def _name_ring_center(base, center, kind, probe=None, labels=None):
    from ._diester_ring_diyl import _system_of

    if isinstance(center, list):
        if labels:
            raise UnsupportedStructure("an isotopic modification of several N-oxides is not supported yet")
        token = RING_CENTER.set(True)
        cation_token = CATION_PARENT.set(True)
        atom_token = CATION_RING_ATOM.set(center[0])
        try:
            return _name_labelled(base, {}, lambda name, placed: _attach_oxides(name, placed, center))
        finally:
            CATION_RING_ATOM.reset(atom_token)
            CATION_PARENT.reset(cation_token)
            RING_CENTER.reset(token)
    if kind == "ium" and any(_anionic_group_atom(base, a) or a.HasProp(ANION_PROP) for a in base.GetAtoms()):
        marked = Chem.Mol(base)
        marked.GetAtomWithIdx(center).SetBoolProp(CATION_INSERT_PROP, True)
        return _name_labelled(marked, labels or {}, lambda name, placed: _cation_inserted(name, base, center, placed))
    if not labels and len(_system_of(base, center)[1]) == base.GetNumAtoms():
        try:
            name, placed = _bare_ring_name(base, center)
        except UnsupportedStructure:
            if not (_ortho_fused_only(base) or _saturated_rings(base)):
                raise
            name = None
        if name is not None:
            return _attach(name, placed, center, kind)
    token = RING_CENTER.set(True)
    cation_token = CATION_PARENT.set(True)
    atom_token = CATION_RING_ATOM.set(center)
    try:
        if probe is None:
            return _name_labelled(base, labels or {}, lambda name, placed: _attach(name, placed, center, kind))
        marked, marked_center, word = probe

        def attach_marked(name, placed):
            locant = placed.get(marked_center)
            if locant is None:
                raise UnsupportedStructure("the ring nitrogen is not numbered in the parent hydride")
            return _attach(_without_marker(name, locant, word), {marked_center: locant}, marked_center, kind)

        return _name_labelled(marked, labels or {}, attach_marked)
    finally:
        CATION_RING_ATOM.reset(atom_token)
        CATION_PARENT.reset(cation_token)
        RING_CENTER.reset(token)


def _monocycle_locant(base, center):
    """The locant of `center` in a bare heteromonocycle: heteroatoms take the lowest locants, then the senior
    element the lowest of those (P-22.2.2.1.2)."""
    cycle = ring_cycle(
        {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in base.GetAtoms()}, [a.GetIdx() for a in base.GetAtoms()]
    )
    hetero = [a for a in cycle if base.GetAtomWithIdx(a).GetAtomicNum() != 6]
    size = len(cycle)
    best = None
    for start in range(size):
        for step in (1, -1):
            locants = {cycle[(start + step * k) % size]: k + 1 for k in range(size)}
            key = (
                sorted(locants[h] for h in hetero),
                [locants[h] for h in sorted(hetero, key=lambda h: (_HETERO_RANK.get(base.GetAtomWithIdx(h).GetSymbol(), 99), h))],
                locants[center],
            )
            if best is None or key < best[0]:
                best = (key, locants)
    return best[1][center]


def _ortho_fused_only(mol):
    """Whether every two rings of `mol` share at most one bond: a hydrogenated ortho-fused system is named with hydro
    prefixes on its mancude parent, a bridged one needs bridge prefixes."""
    rings = [set(r) for r in Chem.GetSymmSSSR(mol)]
    return all(len(a & b) <= 2 for i, a in enumerate(rings) for b in rings[i + 1 :])


def _saturated_rings(mol):
    """Whether no ring bond is aromatic or multiple, so that a bridged system is a von Baeyer hydride."""
    return all(b.GetBondTypeAsDouble() == 1.0 for b in mol.GetBonds() if b.IsInRing())


def _numbered_monocycle(base, center):
    """(name, {center: locant}) of a bare heteromonocycle numbered for the indicated hydrogen first, then the cationic
    centre, then the hydro prefixes (P-31.1.4.2.4, P-73.1.1.2); None when the ring has no such numbering."""
    from ._ring_diyl_numbering import monocycle_numberings

    graph = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in base.GetAtoms()}
    try:
        numberings = monocycle_numberings(base, ring_cycle(graph, list(graph)), {}, 0)
    except (UnsupportedStructure, ValueError, KeyError, RuntimeError):
        return None
    if not numberings:
        return None
    best = min(numberings, key=lambda n: (n.pre_key, n.position_of[center], n.unsat_key))
    name = best.text((), 0, suffix="")
    return (name, {center: best.position_of[center]}) if name.endswith("e") else None


def _bare_ring_name(base, center, hydro_ok=False):
    """(name, {center: locant}) of an unsubstituted single-heteroatom ring system; with `hydro_ok` a hydrogenated fused
    system gives the name of its mancude skeleton, the stem of a name that cites the hydro prefixes itself."""
    from ._fusion_name import Context, fused_ring_system_name, fusion_name, system_numbering_options
    from ._fused_numbering import _locant_key
    from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic

    if base.GetRingInfo().NumRings() == 1:
        numbered = _numbered_monocycle(base, center)
        if numbered is not None:
            return numbered
        if has_hetero_monocyclic_name(base):
            return name_hetero_monocyclic(base), {center: _monocycle_locant(base, center)}
        from .core import smiles_to_iupac

        try:
            named = smiles_to_iupac(Chem.MolToSmiles(base))
        except UnsupportedStructure:
            raise UnsupportedStructure("this ring has no supported parent name") from None
        return re.sub(r"^(?:\d+[a-z]?(?:,\d+[a-z]?)*H-)+", "", named), {center: _monocycle_locant(base, center)}
    name = fused_ring_system_name(base)
    if name is None and hydro_ok and _ortho_fused_only(base):
        try:
            name = str(fusion_name(base)[0])
        except (UnsupportedStructure, ValueError, KeyError, RuntimeError):
            name = None
    if name is None:
        raise UnsupportedStructure("this ring system has no supported parent name")
    ctx = Context(base)
    fused, root = fusion_name(base)
    options = _with_indicated_hydrogen(base, name, system_numbering_options(ctx, fused, root))
    return name, {center: min((n[center] for n in options), key=_locant_key)}


def _with_indicated_hydrogen(base, name, options):
    """The numberings of a fused ring system that put its NH atoms at the indicated hydrogen locants of `name`; the
    indicated hydrogen is lower than the cationic centre in the order of lowest locants (P-31.1.4.2.4)."""
    prefix = re.match(r"(?:\d+[a-z]?(?:,\d+[a-z]?)*H-)+", name)
    cited = set(re.findall(r"\d+[a-z]?", prefix.group(0))) if prefix else set()
    nh_atoms = [a.GetIdx() for a in base.GetAtoms() if a.GetAtomicNum() == 7 and a.GetTotalNumHs() and a.IsInRing()]
    if not cited or not nh_atoms:
        return options
    kept = [n for n in options if all(str(n[a]) in cited for a in nh_atoms)]
    return kept or options


def _azoxy_charged_atoms(mol):
    atoms = set()
    for atom in mol.GetAtoms():
        found = azoxy_atoms(mol, atom)
        if found is not None:
            atoms |= {atom.GetIdx(), found[0].GetIdx()}
    return atoms


def _plain_oxy_group(mol, atom):
    """A hydroxy or alkoxy group: an oxygen with a hydrogen or one non-functional carbon."""
    if atom.GetAtomicNum() != 8 or atom.GetFormalCharge() or atom.IsInRing():
        return False
    if atom.GetDegree() == 1:
        return atom.GetTotalNumHs() == 1
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return atom.GetDegree() == 2 and len(carbons) == 1 and not (
        is_functional_carbon(mol, carbons[0].GetIdx()) or _double_oxygens(mol, carbons[0].GetIdx())
    )


def _hydride_centre(atom):
    """An acyclic neutral atom of a mononuclear parent hydride other than carbon (silane, germane, ...)."""
    return atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not atom.IsInRing() and not atom.GetFormalCharge()


def _aminium_base(mol, ignore=frozenset()):
    """The mol with its ammonium nitrogen neutralised when the only charge is one N+ bonded to carbon and hydrogen
    (a cation outranks every acid, P-41), else None; the charges of an azoxy group are not ammonium charges."""
    ignore = {*ignore, *_azoxy_charged_atoms(mol)}
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in ignore]
    if not charged or any(a.GetFormalCharge() != 1 or a.GetAtomicNum() != 7 for a in charged):
        return None
    for nitrogen in charged:
        if (
            nitrogen.GetIsAromatic()
            or nitrogen.IsInRing()
            or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 4
            or any(
                (n.GetAtomicNum() != 6 and not _plain_oxy_group(mol, n) and not _hydride_centre(n))
                or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
                or _double_oxygens(mol, n.GetIdx())
                for n in nitrogen.GetNeighbors()
            )
            or not any(n.GetAtomicNum() == 6 or _hydride_centre(n) for n in nitrogen.GetNeighbors())
            or (nitrogen.GetTotalNumHs() == 0 and any(_hydride_centre(n) for n in nitrogen.GetNeighbors()))
        ):
            return None
    if all(n.GetTotalNumHs() == 0 for n in charged):
        return mol
    neutral = Chem.RWMol(mol)
    for nitrogen in charged:
        if not nitrogen.GetTotalNumHs():
            continue
        atom = neutral.GetAtomWithIdx(nitrogen.GetIdx())
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(nitrogen.GetTotalNumHs() - 1)
        atom.SetNoImplicit(True)
        atom.SetBoolProp("_cationic_amine", True)
    Chem.SanitizeMol(neutral)
    return neutral.GetMol()


def _cationic_prefix_nitrogen(atom):
    """An acyclic onium centre that an anionic parent cites as an 'azaniumyl', 'sulfaniumyl' or similar prefix
    (P-74.1.3)."""
    from ._hetero_prefixes import _ONIUM_PREFIX_STEMS, is_diazonium_nitrogen

    return is_diazonium_nitrogen(atom) or (
        (atom.GetAtomicNum() == 7 or atom.GetAtomicNum() in _ONIUM_PREFIX_STEMS)
        and atom.GetFormalCharge() == 1
        and not atom.IsInRing()
        and not atom.GetIsAromatic()
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    )


def _iminium_base(mol, ignore=frozenset()):
    """The mol with its iminium nitrogens neutralised when they still carry a hydrogen, else the cation itself, for the
    =N(+)< centres of acyclic C=N groups whose other substituents are carbon (P-73.1.2.1); else None."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in ignore]
    if not charged or any(a.GetFormalCharge() != 1 or a.GetAtomicNum() != 7 for a in charged):
        return None
    for nitrogen in charged:
        if nitrogen.GetIsAromatic() or nitrogen.IsInRing() or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 3:
            return None
        doubles = [b for b in nitrogen.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
        if len(doubles) != 1 or doubles[0].GetOtherAtom(nitrogen).GetAtomicNum() != 6:
            return None
        if any(n.GetAtomicNum() != 6 for n in nitrogen.GetNeighbors()):
            return None
    if not any(n.GetTotalNumHs() for n in charged):
        return mol
    neutral = Chem.RWMol(mol)
    for nitrogen in charged:
        if nitrogen.GetTotalNumHs():
            atom = neutral.GetAtomWithIdx(nitrogen.GetIdx())
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(nitrogen.GetTotalNumHs() - 1)
            atom.SetNoImplicit(True)
    Chem.SanitizeMol(neutral)
    return neutral.GetMol()


_CATION_KIND_ORDER = ("amide", "amidine", "nitrile", "amine", "imine")


def _cation_kind(mol, atom):
    """The kind of cationic characteristic group an acyclic nitrogen cation stands for (Table 7.4): acylammonium (amide),
    amidinium, nitrilium, ammonium (amine) or iminium; None for any other charged atom."""
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.IsInRing() or atom.GetIsAromatic():
        return None
    if any(b.GetBondTypeAsDouble() == 3.0 for b in atom.GetBonds()):
        return "nitrile"
    doubles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
    valence = atom.GetDegree() + atom.GetTotalNumHs()
    if doubles:
        carbon = doubles[0].GetOtherAtom(atom)
        if len(doubles) != 1 or valence != 3 or carbon.GetAtomicNum() != 6:
            return None
        amino = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != atom.GetIdx() and not n.GetFormalCharge()]
        if atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2 and len(amino) == 1 and carbon.GetDegree() == 3:
            return "amidine"
        return "imine"
    if valence != 4:
        return None
    if any(n.GetAtomicNum() == 6 and _double_oxygens(mol, n.GetIdx()) for n in atom.GetNeighbors()):
        return "amide"
    return "amine"


_NEUTRAL_GROUP_PATTERNS = {
    "amide": Chem.MolFromSmarts("[CX3;!R](=O)[NX3]"),
    "amidine": Chem.MolFromSmarts("[CX3;!R](=[NX2])[NX3]"),
    "nitrile": Chem.MolFromSmarts("[CX2;!R]#[NX1]"),
}


def _demote_neutral_groups(base, kind, cations):
    """`base` with the neutral groups of the cationic kind marked as junior, so that they are cited as prefixes: the
    cation outranks the neutral group of the same class (P-41, P-73.7)."""
    pattern = _NEUTRAL_GROUP_PATTERNS.get(kind)
    if pattern is None:
        return base
    demoted = None
    for match in base.GetSubstructMatches(pattern):
        group_nitrogens = [i for i in match[1:] if base.GetAtomWithIdx(i).GetAtomicNum() == 7]
        if any(i in cations for i in group_nitrogens):
            continue
        demoted = demoted or Chem.Mol(base)
        demoted.GetAtomWithIdx(match[0]).SetBoolProp(JUNIOR_GROUP, True)
    return base if demoted is None else demoted


def _cation_parent_base(mol):
    """(neutral base, parent kind, number of cations of that kind) for the cationic characteristic groups of `mol` when
    the cations are of several kinds, are amidinium groups, or neutral groups of their own class stand beside them: the
    parent carries the most cationic suffixes, then the senior kind (P-73.7); the cations of other kinds stay charged and
    are cited as 'azaniumyl'-type prefixes, neutral groups of the class are prefixes too."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    kinds = {a.GetIdx(): _cation_kind(mol, a) for a in charged}
    if not kinds or None in kinds.values():
        return None
    counts = Counter(kinds.values())
    kind = min(counts, key=lambda k: (-counts[k], _CATION_KIND_ORDER.index(k)))
    ignore = frozenset(i for i, k in kinds.items() if k != kind)
    if kind == "amine":
        base = _aminium_base(mol, ignore=ignore)
    elif kind == "imine":
        base = _iminium_base(mol, ignore=ignore)
    elif kind == "amidine":
        base = _amidinium_base(mol, ignore=ignore)
    else:
        base = _group_cation_base(mol, kind, ignore=ignore, neutral_groups=True)
    if base is None:
        return None
    cations = {i for i, k in kinds.items() if k == kind}
    demoted = _demote_neutral_groups(base, kind, cations)
    if len(counts) == 1 and kind != "amidine" and demoted is base:
        return None
    return demoted, kind, counts[kind]


def _prefix_ammonium(mol, nitrogen):
    return (
        nitrogen.GetDegree() + nitrogen.GetTotalNumHs() == 4
        and nitrogen.GetDegree() > 0
        and all(
            n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and not _double_oxygens(mol, n.GetIdx())
            for n in nitrogen.GetNeighbors()
        )
    )


def _amidinium_base(mol, ignore=frozenset()):
    """The neutral amidine of the charges in `mol` (besides the `ignore`d ring centres), amidinium groups C(NR2)=NH2(+)
    on carbon; each cation takes the suffix 'imidamidium' with the N locants of the neutral amidine (P-73.5.3.1)."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in ignore]
    if not charged:
        return None
    neutral = Chem.RWMol(mol)
    for cation in charged:
        if cation.GetFormalCharge() != 1 or cation.GetAtomicNum() != 7 or cation.GetIsAromatic() or cation.IsInRing():
            return None
        if cation.GetDegree() != 1 or cation.GetTotalNumHs() != 2:
            return None
        carbon = cation.GetNeighbors()[0]
        if carbon.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(cation.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() != 2.0:
            return None
        others = [n for n in carbon.GetNeighbors() if n.GetIdx() != cation.GetIdx()]
        amino = [n for n in others if n.GetAtomicNum() == 7 and not n.GetIsAromatic() and not n.GetFormalCharge()]
        if len(amino) != 1 or len(others) != 2:
            return None
        atom = neutral.GetAtomWithIdx(cation.GetIdx())
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(1)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    return base


def _group_cation_base(mol, kind, ignore=frozenset(), neutral_groups=False):
    """The mol with its hydrogen-bearing cationic nitrogens neutralised, for the acylammonium (amidium, a nitrogen with
    four bonds on an amide carbonyl) or nitrilium (C#N-H) centres of a polycation of one kind (Table 7.4); else None."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in ignore]
    if not charged or any(a.GetFormalCharge() != 1 or a.GetAtomicNum() != 7 or a.GetIsAromatic() or a.IsInRing() for a in charged):
        return None
    if kind == "amide":
        # a cationic amide group outranks an ammonium group, which is then an 'azaniumyl' prefix (P-73.7)
        charged = [a for a in charged if not _prefix_ammonium(mol, a)]
        if not charged:
            return None
    for nitrogen in charged:
        triple = [b for b in nitrogen.GetBonds() if b.GetBondTypeAsDouble() == 3.0]
        if kind == "nitrile":
            if len(triple) != 1 or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 2:
                return None
            if any(
                b.GetBondTypeAsDouble() not in (1.0, 3.0) or b.GetOtherAtom(nitrogen).GetAtomicNum() != 6
                for b in nitrogen.GetBonds()
            ):
                return None
        else:
            acyl = [
                n
                for n in nitrogen.GetNeighbors()
                if n.GetAtomicNum() == 6 and any(
                    b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() == 8 for b in n.GetBonds()
                )
            ]
            if (
                triple
                or len(acyl) != 1
                or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 4
                or any(
                    n.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
                    for n in nitrogen.GetNeighbors()
                )
            ):
                return None
    neutral = Chem.RWMol(mol)
    for nitrogen in charged:
        atom = neutral.GetAtomWithIdx(nitrogen.GetIdx())
        if nitrogen.GetTotalNumHs():
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(nitrogen.GetTotalNumHs() - 1)
    base = neutral.GetMol()
    Chem.SanitizeMol(base)
    groups = Chem.MolFromSmarts("[CX2]#[N]" if kind == "nitrile" else "[CX3](=O)[NX3]")
    if not neutral_groups and len(base.GetSubstructMatches(groups)) > len(charged):
        return None
    return base


_IMIDAMIDE_SUFFIX = re.compile(r"(?P<mult>di|tri|tetra)?(?P<carbo>carbox|carbo)?imidamide$")
_GROUP_SUFFIX = re.compile(r"(?P<mult>di|tri|tetra)?(?P<carbo>carbox|carbo)?(?P<kind>amide|nitrile)$")


def _ring_cation_with_group(mol):
    """(base, ring centre, kind) when one ring nitrogen cation sits beside cationic groups of one kind (aminium,
    amidium, nitrilium) on the same parent: the groups are neutralised in `base`, the ring centre is returned."""
    centres = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() == 1 and a.GetAtomicNum() == 7 and a.IsInRing()]
    others = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in centres]
    if len(centres) != 1 or not others or len(Chem.GetMolFrags(mol)) != 1:
        return None
    for kind, function in (
        ("amine", lambda m: _aminium_base(m, ignore=frozenset(centres))),
        ("amidine", lambda m: _amidinium_base(m, ignore=frozenset(centres))),
        ("amide", lambda m: _group_cation_base(m, "amide", ignore=frozenset(centres))),
        ("nitrile", lambda m: _group_cation_base(m, "nitrile", ignore=frozenset(centres))),
    ):
        base = function(mol)
        if base is not None:
            return base, centres[0], kind, mol
    return None


def _ylium_ring_group_cation(mol, centre, kind):
    """The name of a ring whose centre is a lambda 'ylium' centre (a nitrogen with three ring bonds) together with
    unsubstituted aminium or carboximidamidium groups on its ring atoms, or None for another ring centre."""
    from ._ylium_ring import _is_ylium_centre, name_ylium_ring

    if not _is_ylium_centre(mol.GetAtomWithIdx(centre)):
        return None
    group = {}
    for atom in mol.GetAtoms():
        if not atom.GetFormalCharge() or atom.GetIdx() == centre:
            continue
        if kind == "amine" and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            root, ring_atom = atom, atom.GetNeighbors()[0]
        elif kind == "amidine" and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2:
            root = atom.GetNeighbors()[0]
            ring_atom = next((n for n in root.GetNeighbors() if n.IsInRing()), None)
            amino = [n for n in root.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != atom.GetIdx()]
            if ring_atom is None or len(amino) != 1 or amino[0].GetDegree() != 1 or amino[0].GetTotalNumHs() != 2:
                raise UnsupportedStructure("a substituted amidinium group on a lambda ylium ring is not supported yet")
        else:
            raise UnsupportedStructure("this cationic group on a lambda ylium ring is not supported yet")
        if not ring_atom.IsInRing() or ring_atom.GetIdx() in group:
            raise UnsupportedStructure("this cationic group on a lambda ylium ring is not supported yet")
        group[ring_atom.GetIdx()] = root.GetIdx()
    return name_ylium_ring(mol, group, "aminium" if kind == "amine" else "carboximidamidium")


def _name_ring_group_cation(base, centre, kind, original=None):
    """A ring cation together with a cationic suffix: the skeletal centre is cited as 'ium' after the parent hydride and
    takes the lowest locant before the suffix (P-73.5.3.1, P-73.5.3.2); the suffix names the group."""
    from ._diester_ring_diyl import CATION_CENTRES

    ylium = _ylium_ring_group_cation(original, centre, kind) if original is not None and kind in ("amine", "amidine") else None
    if ylium is not None:
        return ylium
    editable = Chem.RWMol(base)
    atom = editable.GetAtomWithIdx(centre)
    hydrogens = atom.GetTotalNumHs()
    nitrenium = bool(atom.GetNumRadicalElectrons())
    atom.SetFormalCharge(0)
    atom.SetNumRadicalElectrons(0)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(max(hydrogens - 1, 0))
    if atom.GetIsAromatic():
        atom.SetBoolProp("_ring_cation_centre", True)
    ring_base = editable.GetMol()
    ring_base.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(ring_base)
    token = AMINIUM.set(True if kind == "amine" else kind)
    centres_token = CATION_CENTRES.set(((centre, "ium"),))
    try:
        name = _name_labelled(ring_base, {}, lambda text, placed: _cation_inserted(text, ring_base, centre, placed))
    finally:
        CATION_CENTRES.reset(centres_token)
        AMINIUM.reset(token)
    name = _cationic_group_suffix(name, kind)
    return re.sub(r"(\d)-ium(?![a-z])", r"\1-ylium", name, count=1) if nitrenium else name


def _name_aminium(base, labels=None, parent="amine", cationic=None):
    token = AMINIUM.set(True if parent == "amine" else parent)
    cation_token = CATION_PARENT.set(True)
    try:
        name = _name_labelled(base, labels or {})
    finally:
        CATION_PARENT.reset(cation_token)
        AMINIUM.reset(token)
    return _cationic_group_suffix(
        name, parent, cationic or sum(a.HasProp("_cationic_amine") or (a.GetFormalCharge() == 1 and a.GetAtomicNum() == 7) for a in base.GetAtoms())
    )


_MULTIPLIED = {"di": 2, "tri": 3, "tetra": 4, "penta": 5, "hexa": 6}
_UNIT_COUNT = {"di": 2, "bis": 2, "tri": 3, "tris": 3, "tetra": 4, "tetrakis": 4}
_UNIT_WORD = {2: "bis", 3: "tris", 4: "tetrakis"}


def _multiplied_amine_unit(name):
    """(index where the multiplying prefix starts, its prefix, the unit text without its amine ending, the ending) of a
    multiplicative name whose last component is a multiplied amine or aniline unit, else None."""
    ending = next((e for e in ("amine", "aniline") if name.endswith(e) or name.endswith(e + ")")), None)
    if ending is None:
        return None
    if not name.endswith(")"):
        match = re.search(r"(di|tri|tetra)(aniline)$", name)
        return None if match is None or ending != "aniline" else (match.start(), match.group(1), "", "aniline")
    depth = 0
    for index in range(len(name) - 1, -1, -1):
        depth += (name[index] == ")") - (name[index] == "(")
        if depth == 0:
            inner = name[index + 1 : -1]
            before = re.search(r"(di|tri|tetra|bis|tris|tetrakis)$", name[:index])
            if before is None or not inner.endswith(ending):
                return None
            return before.start(), before.group(1), inner[: -len(ending)], ending
    return None


def cationic_multiplied_units(name, cationic):
    """`name` with the identical amine units of a multiplicative name made aminium units, or None when it is not one
    (P-15.3.2, P-73.1.1.2); the locant of an ethane unit is cited because multiplicative names cite all locants (P-14.3.4)."""
    found = _multiplied_amine_unit(name)
    if found is None:
        return None
    start, prefix, unit, ending = found
    count = _UNIT_COUNT[prefix]
    if cationic and count != cationic:
        raise UnsupportedStructure("a neutral amino group beside the cationic units is not cited as a prefix of this parent")
    unit = re.sub(r"(?<!m)ethan$", "ethan-1-", unit)
    return f"{name[:start]}{_UNIT_WORD[count]}({unit}{ending[:-1]}ium)"


def _cationic_group_suffix(name, parent, cationic=None):
    if parent == "amidine":
        match = _IMIDAMIDE_SUFFIX.search(name)
        if match is None:
            raise UnsupportedStructure("the amidinium cation is not named as an amidine parent")
        count = _MULTIPLIED[match.group("mult")] if match.group("mult") else 1
        word = (match.group("carbo") or "") + "imidamidium"
        if count > 1:
            return f"{name[:match.start()]}{_UNIT_WORD[count]}({word})"
        return name[:match.start()] + word
    if parent in ("amide", "nitrile"):
        match = _GROUP_SUFFIX.search(name)
        if match is None or match.group("kind") != parent:
            raise UnsupportedStructure("the cation is not named as an amide or nitrile parent")
        word = (match.group("carbo") or "") + parent[:-1] + "ium"
        if match.group("mult"):
            multiplier = {"di": "bis", "tri": "tris", "tetra": "tetrakis"}[match.group("mult")]
            return f"{name[:match.start()]}{multiplier}({word})"
        return name[:match.start()] + word
    multiplied = cationic_multiplied_units(name, cationic) if parent == "amine" else None
    if multiplied is not None:
        return multiplied
    if not name.endswith(("amine", "aniline", "imine") if parent == "imine" else ("amine", "aniline")):
        raise UnsupportedStructure("the cation is not named as an amine or imine parent")
    imines = re.search(r"(di|tri|tetra)imine$", name) if parent == "imine" else None
    if imines is not None:
        if cationic and _MULTIPLIED[imines.group(1)] != cationic:
            raise UnsupportedStructure("a neutral imine group beside a cationic one is not cited as a prefix of this parent")
        word = {"di": "bis", "tri": "tris", "tetra": "tetrakis"}[imines.group(1)]
        return f"{name[:imines.start()]}{word}(iminium)"
    if parent == "imine" and cationic and cationic > 1:
        raise UnsupportedStructure("several iminium centres need a multiplied imine parent")
    multiple = re.search(r"(di|tri|tetra|penta|hexa)(amine|aniline)$", name)
    if multiple is not None:
        word = {"di": "bis", "tri": "tris", "tetra": "tetrakis", "penta": "pentakis", "hexa": "hexakis"}[multiple.group(1)]
        return f"{name[:multiple.start()]}{word}(aminium)"
    return name[:-1] + "ium"


STEREO_OF_ISOTOPOLOGUE = contextvars.ContextVar("stereo_of_isotopologue", default=None)
STEREO_SUFFIX = contextvars.ContextVar("stereo_suffix", default=None)


def _located_elements(mol):
    located = []
    for kind, idx, code in specified_stereo_elements(mol) or []:
        if kind == "bond":
            bond = mol.GetBondWithIdx(idx)
            located.append(("bond", (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), code))
        else:
            located.append((kind, idx, code))
    return located


@contextlib.contextmanager
def ion_stereo_carried(mol, neutral):
    """The descriptors of an ion carried over to its neutral copy that names it: a charged ligand ranks below its
    neutral twin, and rings that differ only in the charged atom leave a centre unrecognised (P-92.1.4.4, P-92.2.1)."""
    try:
        elements = _located_elements(mol)
    except UnsupportedStructure:
        elements = None
    try:
        carried = elements is not None and elements != _located_elements(neutral)
    except UnsupportedStructure:
        carried = elements is not None
    if not carried:
        yield
        return
    token = STEREO_OF_ISOTOPOLOGUE.set(elements)
    try:
        yield
    finally:
        STEREO_OF_ISOTOPOLOGUE.reset(token)


@contextlib.contextmanager
def _oxide_stereo_carried(mol, base):
    """The descriptors of an N-oxide, whose oxygen is the fourth ligand of the nitrogen (P-93.2.3), moved to the atoms of
    `base`, the amine without its oxide atoms, that names it."""
    oxides = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 8
        and a.GetFormalCharge() == -1
        and a.GetDegree() == 1
        and a.GetNeighbors()[0].GetAtomicNum() == 7
        and a.GetNeighbors()[0].IsInRing()
    }
    if not oxides or mol.GetNumAtoms() - len(oxides) != base.GetNumAtoms():
        yield
        return
    try:
        elements = _located_elements(mol)
    except UnsupportedStructure:
        yield
        return

    def moved(idx):
        return idx - sum(1 for o in oxides if o < idx)

    located = [
        (kind, moved(idx) if kind == "atom" else tuple(moved(i) for i in idx), code)
        for kind, idx, code in elements
        if kind == "bond" or idx not in oxides
    ]
    token = STEREO_OF_ISOTOPOLOGUE.set(located)
    try:
        yield
    finally:
        STEREO_OF_ISOTOPOLOGUE.reset(token)


JUNIOR_GROUP = "_junior_group"
ACYL_ROOT = "_acyl_root"


def _modifications(entry):
    return [*([entry["skeleton"]] if entry["skeleton"] else []), *(n for n, c in entry["H"].items() for _ in range(c))]


def _nuclide_rank(nuclide):
    """Higher atomic number first, then higher mass number (P-82.2.2.2)."""
    table = Chem.GetPeriodicTable()
    symbol = "".join(ch for ch in nuclide if ch.isalpha())
    return -table.GetAtomicNumber(symbol), -int("".join(ch for ch in nuclide if ch.isdigit()))


def _demote_junior_groups(mol, labels):
    """Characteristic groups that are modified in different ways cannot be multiplied: the group with the most
    modifications, then the nuclide of higher atomic number and mass number, stays the suffix and the others are cited
    as prefixes (P-82.2.2.2)."""
    if not labels:
        return mol
    principal, owned_groups = _principal_owned(mol)
    found = {frozenset(g) for g in owned_groups}
    substituted = {frozenset({a.GetIdx()}) for a in mol.GetAtoms() if _substituted_amine_nitrogen(mol, a)}
    if principal in (None, "amine"):
        found |= substituted
        principal = "amine" if found else None
    groups = {g for g in found if not any(g < other for other in found)}
    if principal is None or len(groups) < 2:
        return mol

    def signature(group):
        return tuple(sorted((n for a in group if a in labels for n in _modifications(labels[a])), key=_nuclide_rank))

    senior = min(groups, key=lambda g: (-len(signature(g)), [_nuclide_rank(n) for n in signature(g)]))
    others = [g for g in groups if signature(g) != signature(senior)]
    if not others:
        return mol
    demoted = Chem.Mol(mol)
    for group in others:
        if group in substituted:
            demoted.GetAtomWithIdx(next(iter(group))).SetBoolProp(JUNIOR_GROUP, True)
            continue
        anchors = [a for a in group if mol.GetAtomWithIdx(a).GetAtomicNum() == 6] or [
            n.GetIdx() for a in group for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in group
        ]
        for atom in anchors:
            demoted.GetAtomWithIdx(atom).SetBoolProp(JUNIOR_GROUP, True)
    return demoted


def _hydrogen_isotope_bonds(original, known):
    """E/Z codes of double bonds whose stereo rests on a hydrogen isotope, which RDKit's perception drops."""
    if not any(a.GetAtomicNum() == 1 and a.GetIsotope() for a in original.GetAtoms()):
        return []
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)
    try:
        probe = Chem.Mol(original)
        Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
        assign_cip_labels(probe)
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)
    return [
        ("bond", b.GetIdx(), b.GetProp("_CIPCode"))
        for b in probe.GetBonds()
        if b.HasProp("_CIPCode") and b.GetProp("_CIPCode") in ("E", "Z") and b.GetIdx() not in known
    ]


def _name_isotopic(original, clean, labels, index_map, build=None):
    """Stereodescriptors come from the molecule with its nuclides (a CHD centre is a stereocentre, P-82.4) and are
    carried over to the atoms of the unlabelled molecule."""
    located = []
    elements = list(specified_stereo_elements(original) or [])
    extra = _hydrogen_isotope_bonds(original, {idx for kind, idx, _ in elements if kind == "bond"})
    elements += extra
    for kind, idx, code in elements:
        if kind == "bond":
            bond = original.GetBondWithIdx(idx)
            ends = (index_map.get(bond.GetBeginAtomIdx()), index_map.get(bond.GetEndAtomIdx()))
            if None in ends:
                raise UnsupportedStructure("a double bond to an isotopically labelled hydrogen is not supported")
            located.append(("bond", ends, code))
        elif idx in index_map:
            located.append(("atom", index_map[idx], code))
    specified_bonds = [
        b for b in original.GetBonds()
        if b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
    ]
    if len(specified_bonds) + len(extra) != sum(1 for k, _, _ in located if k == "bond") or any(
        a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED and a.GetAtomicNum() != 1 and a.GetIdx() not in {i for k, i, _ in located if k == "atom"} and False
        for a in original.GetAtoms()
    ):
        raise UnsupportedStructure("a stereo element of this isotopically modified structure is not cited by a supported name")
    token = STEREO_OF_ISOTOPOLOGUE.set(located)
    try:
        return (build or _name_labelled)(clean, labels)
    finally:
        STEREO_OF_ISOTOPOLOGUE.reset(token)


def _with_labels(mol, labels, consumed, name, parts, reselect):
    """`name` with the isotopic descriptor of the labelled atoms of `mol`: parent atoms cite their locants, atoms of
    the principal group take letter locants or go before the suffix (P-82.2, P-82.6)."""
    from ._substituents import ISOTOPE_LABELS

    cited = (ISOTOPE_LABELS.get() or {}).get("cited", ())
    in_parent = {a: e for a, e in labels.items() if a in parts[4] and a not in cited}
    if in_parent and len(parts[4]) > 1 and not _locants_omitted(mol, in_parent, parts):
        force_token = FORCE_LOCANTS.set(True)
        try:
            name, parts = reselect()
        finally:
            FORCE_LOCANTS.reset(force_token)
    # a suffix group cannot also be cited as a prefix, so labels that candidate names consumed there are still open
    parent_groups = [
        group for group in _principal_owned(mol)[1]
        if any(a in parts[4] or any(n.GetIdx() in parts[4] for n in mol.GetAtomWithIdx(a).GetNeighbors()) for a in group)
    ]
    suffix_atoms = set().union(*parent_groups) if parent_groups else set()
    rest = {a: e for a, e in labels.items() if a not in in_parent and a not in cited and (a not in consumed or a in suffix_atoms)}
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in rest):
        name, parts, suffix_text = _carbon_group_label(mol, rest, reselect)
        front = []
    else:
        front, suffix_text = _group_atom_labels(mol, rest, parts[4]) if rest else ([], None)
    if suffix_text and suffix_text[1] == "al" and suffix_text[2] == 1 and name.endswith("aldehyde"):
        # a retained aldehyde name has no 'al' suffix to modify: the nuclide stands in front of the parent as for an amide
        front = [(nuclide, None, 1, False) for nuclide in suffix_text[0].strip("()").split(",")]
        suffix_text = None
    if suffix_text:
        name = _before_suffix(name, *suffix_text)
    if in_parent or front:
        bare = len(parts[4]) == 1 or _locants_omitted(mol, in_parent, parts)
        capacity = {a: mol.GetAtomWithIdx(a).GetTotalNumHs() for a in in_parent}
        name = _label_parent(mol, name, parts, in_parent, front, capacity, bare)
    return name


_HYDRO_PREFIX = re.compile(r"(-?)(?:(\d[0-9a-z',]*)-)?([a-z]*?hydro)(?!xy|per)")


def _label_parent(mol, name, parts, labels, front, capacity, bare):
    """`name` with the isotopic descriptor after the hydro prefixes, which are detachable prefixes at the front of the
    parent hydride (P-82.2.3). A hydrogen added by a hydro prefix beyond the hydrogens of the mancude parent is the
    modified one, cited as 'di[(2H)hydro]' before the unmodified prefixes."""
    from ._isotope_labels import descriptor

    split = len(parts[0]) if parts[0] is not None else (parts[6] if len(parts) > 6 else None)
    if split is None:
        raise UnsupportedStructure("the parent hydride part of this name is not delimited")
    head, rest = name[:split], name[split:]
    match = _HYDRO_PREFIX.match(rest)
    hydro_text = ""
    if match is not None:
        rest = rest[match.end() :]
        labels, capacity, hydro_text = _split_hydro_labels(mol, parts[4], labels, capacity, match)
    text = (
        descriptor(labels, parts[4], bare, front, capacity, _sole_heteroatoms(mol, labels, parts[4]), bare and len(parts[4]) > 1)
        if labels or front
        else ""
    )
    return head + hydro_text + text + ("-" if text and rest[:1].isdigit() else "") + rest


def _split_hydro_labels(mol, position_of, labels, capacity, match):
    """(descriptor labels, capacities, hydro prefix text) once the hydrogens that the hydro prefixes added are taken
    out of the labels."""
    lead, locant_text, _ = match.groups()
    locants = locant_text.split(",") if locant_text else None
    atom_of = {str(locant): atom for atom, locant in position_of.items()}
    if locants is not None and any(x not in atom_of for x in locants):
        raise UnsupportedStructure("a hydro prefix locant is not a position of the parent hydride")
    hydro_atoms = {atom_of[x] for x in locants} if locants is not None else set(position_of)
    kept, kept_capacity, modified = {}, dict(capacity), {}
    for atom, entry in labels.items():
        total = sum(entry["H"].values())
        parent_hydrogens = max(mol.GetAtomWithIdx(atom).GetTotalNumHs() - (atom in hydro_atoms), 0)
        if atom in hydro_atoms and total > parent_hydrogens:
            if len(entry["H"]) != 1 or total - parent_hydrogens != 1:
                raise UnsupportedStructure("the hydrogens added by a hydro prefix are modified in different ways")
            (nuclide,) = entry["H"]
            modified[atom] = nuclide
            entry = {"skeleton": entry["skeleton"], "H": {nuclide: parent_hydrogens} if parent_hydrogens else {}}
        if entry["skeleton"] or entry["H"]:
            kept[atom] = entry
        if atom in hydro_atoms:
            kept_capacity[atom] = parent_hydrogens
    if not modified:
        return labels, kept_capacity, match.group(0)
    if locants is None:
        raise UnsupportedStructure("hydro prefixes without locants cannot be split by their isotopic modification")
    groups = {}
    for x in locants:
        groups.setdefault(modified.get(atom_of[x]), []).append(x)
    if any(len(group) % 2 for group in groups.values()):
        return labels, kept_capacity, match.group(0)
    ordered = sorted(groups, key=lambda nuclide: (nuclide is None, nuclide or ""))
    text = "-".join(
        f"{','.join(groups[n])}-{multiplied_word(len(groups[n]), 'hydro' if n is None else f'[({n})hydro]')}" for n in ordered
    )
    return kept, kept_capacity, lead + text


def _carbon_group_label(mol, rest, reselect):
    """The nuclide on the carbon of a 'carbo' suffix is cited before the suffix of the systematic name, because the
    retained name numbers no such position (P-82.6.3.2: 'benzene(13C)carboxylic acid')."""

    principal, owned_groups = _principal_owned(mol)
    carbon_groups = {frozenset(g) for g in owned_groups if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in g)}
    modified = {g for g in carbon_groups if any(a in rest for a in g)}

    def nuclides(group):
        """The descriptor items of the group: the nuclides of its atoms, those of the hydrogens of its carbon with the
        number of hydrogens modified (P-82.2.1)."""
        items = []
        for a in group:
            if a in rest:
                if rest[a]["skeleton"]:
                    items.append((rest[a]["skeleton"], rest[a]["skeleton"]))
                for nuclide, count in rest[a]["H"].items():
                    on_carbon = mol.GetAtomWithIdx(a).GetAtomicNum() == 6
                    items.extend((nuclide, f"{nuclide}{count}" if on_carbon else nuclide) for _ in range(1 if on_carbon else count))
        return [text for _, text in sorted(items, key=lambda item: _nuclide_sort_key(item[0]))]

    from ._isotope_labels import _nuclide_sort_key

    if (
        principal is None
        or not carbon_groups
        or modified != carbon_groups
        or any(a not in set().union(*carbon_groups) for a in rest)
        or len({tuple(nuclides(g)) for g in carbon_groups}) != 1
        or any(sum(1 for a in g if a in rest and mol.GetAtomWithIdx(a).GetAtomicNum() == 6) != 1 for g in carbon_groups)
    ):
        raise UnsupportedStructure("an isotopically modified carbon of this characteristic group is not supported yet")
    try:
        word = _SUFFIX_WORDS[_RING_SUFFIX[principal]]
    except KeyError:
        raise UnsupportedStructure("this characteristic group has no carbon suffix to modify") from None
    if not word.startswith("carb"):
        raise UnsupportedStructure("the carbon of this characteristic group is not a suffix carbon")
    token = NO_RETAINED_BENZENE.set(True)
    try:
        name, parts = reselect()
    finally:
        NO_RETAINED_BENZENE.reset(token)
    text = f"({','.join(nuclides(next(iter(carbon_groups))))})"
    count = len(carbon_groups)
    if count == 1:
        return name, parts, (text, word)
    tail = multiplied_word(count, word)
    if not name.endswith(tail):
        raise UnsupportedStructure("the suffix of this name is not delimited")
    return name[: -len(tail)] + multiplying_prefix(count) + "[" + text + word + "]", parts, None


def _isotope_entries(labels):
    return [
        ("isotope", atom, "|".join(filter(None, [entry["skeleton"], *(n for n, c in entry["H"].items() for _ in range(c))])))
        for atom, entry in labels.items()
    ]


def _hydride_centre_with_amino_groups(mol):
    """A Group 13/14 atom with two or more amino groups: the hydride holding the most principal groups is the parent,
    named with the suffix 'diamine', not multiplicatively through the atom (P-44.1.1, P-68.2.5)."""
    return any(
        a.GetAtomicNum() in (*_GROUP_14, *_GROUP_13_METALS)
        and not a.IsInRing()
        and sum(
            n.GetAtomicNum() == 7
            and not n.IsInRing()
            and not n.GetFormalCharge()
            and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and not any(
                b.GetBondTypeAsDouble() > 1.0 and b.GetOtherAtom(m).GetAtomicNum() != 6
                for m in n.GetNeighbors()
                if m.GetIdx() != a.GetIdx()
                for b in m.GetBonds()
            )
            for n in a.GetNeighbors()
        )
        >= 2
        for a in mol.GetAtoms()
    )


def _name_labelled(mol, labels, finish=None):
    from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS

    mol = _demote_junior_groups(mol, labels)

    stereo = _check_scope(mol) + _isotope_entries(labels)
    context = {
        "atoms": {where: code for kind, where, code in stereo if kind == "atom"},
        "bonds": {where: code for kind, where, code in stereo if kind == "bond"},
        "used": set(),
    }
    token = BRANCH_STEREO.set(context if stereo else None)
    suffix_token = None
    isotope_context = {"labels": labels, "consumed": set(), "mol": mol}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if labels else None)
    try:
        try:
            multiplicative = (
                None
                if _ionic_linker_centre(mol) or _cationic_amine_parent(mol) or RING_CENTER.get() or _hydride_centre_with_amino_groups(mol)
                else _multiplicative_name(mol, stereo)
            )
        except UnsupportedStructure:
            if not labels:
                raise
            multiplicative = None
        if multiplicative is not None:
            if not set(labels) <= isotope_context["consumed"]:
                raise UnsupportedStructure("isotopic modification of a multiplicative name is not supported yet")
            return multiplicative if finish is None else finish(multiplicative, {})
        isotope_context["consumed"].clear()
        isotope_context.pop("cache", None)
        isotope_context.pop("cited", None)
        owned = set().union(*_principal_owned(mol)[1]) if stereo else frozenset()
        suffix_token = STEREO_SUFFIX.set((mol, owned) if owned else None)
        _, name, parts = _select(mol, stereo=stereo)
        LAST_POSITIONS.set((mol, dict(parts[4])))
        if labels:
            name = _with_labels(
                mol, labels, isotope_context["consumed"], name, parts, lambda: _select(mol, stereo=stereo)[1:]
            )
        if finish is not None:
            name = finish(name, parts[4])
        if AMINIUM.get() is True:
            owned |= {w for kind, w, _ in stereo if kind == "atom" and mol.GetAtomWithIdx(w).GetAtomicNum() == 7}
        prefix = _stereo_prefix(stereo, parts[4], ring_parent=parts[5], used=context["used"], suffix=owned)
        return name if prefix and name.startswith(prefix) else prefix + name
    finally:
        if suffix_token is not None:
            STEREO_SUFFIX.reset(suffix_token)
        BRANCH_STEREO.reset(token)
        ISOTOPE_LABELS.reset(isotope_token)


def _cationic_amine_parent(mol):
    """A single substituted ammonium centre is the parent beside neutral groups (P-41, class 6), so it is never the linker
    of a multiplied name; several centres are the identical parents that a linker joins."""
    centres = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and (a.GetFormalCharge() == 1 or a.HasProp("_cationic_amine"))]
    if AMINIUM.get() is not True:
        return False
    return (len(centres) == 1 and centres[0].GetDegree() > 1) or _neutral_amine_beside_cation(mol, centres)


def _neutral_amine_beside_cation(mol, centres):
    """A neutral amino group beside the cationic nitrogens is a prefix of the aminium parent, not one of its units
    (P-41, class 6), so units that differ in this way are not identical."""
    cationic = {a.GetIdx() for a in centres}
    return bool(cationic) and any(
        a.GetAtomicNum() == 7
        and a.GetIdx() not in cationic
        and not a.GetFormalCharge()
        and not a.IsInRing()
        and not a.GetIsAromatic()
        and all(n.GetAtomicNum() == 6 and not _double_oxygens(mol, n.GetIdx()) for n in a.GetNeighbors())
        and all(b.GetBondTypeAsDouble() == 1.0 for b in a.GetBonds())
        for a in mol.GetAtoms()
    )


def _ionic_linker_centre(mol):
    """An anionic centre with two or more bonds cannot be a repeated unit: it is the parent, senior to every neutral
    group (P-41), unless the terminal anionic groups outnumber the centres in between (P-72.7(a))."""
    marked = [a for a in mol.GetAtoms() if a.HasProp(ANION_PROP)]
    inner = sum(1 for a in marked if a.GetDegree() > 1)
    return bool(inner) and len(marked) - inner <= inner


_SUFFIX_WORD = {
    "alcohol": "ol", "ketone": "one", "aldehyde": "al", "peroxol": "peroxol", "thiol": "thiol", "selenol": "selenol", "tellurol": "tellurol",
    "thione": "thione", "selone": "selone", "tellone": "tellone",
    **{name: name for name in _PEROXOL_ANALOGUES},
}
_SUFFIX_ELEMENT = {"thiol": 16, "selenol": 34, "tellurol": 52, "thione": 16, "selone": 34, "tellone": 52}


def _principal_owned(mol):
    """(principal class, [atoms owned by each principal group])."""
    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], []).append(found[1])
    for cls, _, owned in _ring_occurrences(mol):
        groups.setdefault(cls, []).append(owned)
    principal = _principal_class(set(groups))
    owned = groups.get(principal, [])
    return principal, [group for group in owned if not any(group < other for other in owned)]


_NITROGEN_LABEL_CLASSES = frozenset(
    (
        "amide", "amine", "imine", "nitrile", "amidine", "hydrazide", "sulfonamide",
        *_CHALCOGEN_AMIDE_CLASSES, *_SULFONYL_AMIDE_ORDER, *_SULFONYL_HYDRAZIDE_ORDER,
    )
)


def _group_atom_labels(mol, rest, parent_atoms=None):
    """([(nuclide, locant, count, repeatable)], suffix descriptor text) for nuclides on atoms of the principal
    characteristic group: an amide or amine nitrogen and the oxygens of an acid are cited in front of the parent
    with a letter locant (P-82.2.4, P-82.2.5, P-82.6.2); a hydroxy, oxo or aldehyde oxygen is inserted before the
    suffix (P-82.2.1, P-82.6.1.1). Atoms of the same element in the substituents of the parent do not make the
    letter locant ambiguous; those inside `parent_atoms` do."""
    from ._isotope_labels import _nuclide_sort_key

    principal, owned_groups = _principal_owned(mol)
    if principal in (None, "amine"):
        substituted = [{a.GetIdx()} for a in mol.GetAtoms() if _substituted_amine_nitrogen(mol, a)]
        if substituted:
            principal, owned_groups = "amine", [*owned_groups, *substituted]
    owned = set().union(*owned_groups) if owned_groups else set()
    if principal is None or any(a not in owned for a in rest):
        raise UnsupportedStructure("an isotopically modified atom outside the parent and the principal group is not supported yet")
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in rest}
    before_suffix = principal in _SUFFIX_WORD and elements == {_SUFFIX_ELEMENT.get(principal, 8)}
    if not before_suffix and any(
        a.GetAtomicNum() in elements
        and a.GetIdx() not in owned
        and a.GetIdx() not in rest
        and not (a.IsInRing() and a.GetIdx() in (parent_atoms or ()))
        and (parent_atoms is None or a.GetIdx() in parent_atoms)
        for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("the modified atom of the characteristic group needs a locant that is not defined yet")
    heavy = [a for a, e in rest.items() if e["skeleton"]]
    if not before_suffix and len(elements) == 1 and len(heavy) > 1 and _same_label_on_every_group_atom(mol, rest, owned, next(iter(elements))):
        return [(rest[heavy[0]]["skeleton"], None, len(heavy), False)], None
    if before_suffix:
        if principal == "peroxol" and any(mol.GetAtomWithIdx(a).GetDegree() != 1 for a in rest):
            raise UnsupportedStructure("the oxygen of a peroxol that is bonded to carbon has no defined descriptor")
        per_group = []
        for group in owned_groups:
            nuclides = []
            for a in group:
                if a in rest:
                    nuclides.extend([rest[a]["skeleton"]] if rest[a]["skeleton"] else [])
                    nuclides.extend(n for n, c in rest[a]["H"].items() for _ in range(c))
            per_group.append(tuple(sorted(nuclides, key=_nuclide_sort_key)))
        nuclides = min(per_group, key=lambda found: (-len(found), [_nuclide_rank(n) for n in found]))
        counted = [f"{n}{nuclides.count(n) if nuclides.count(n) > 1 else ''}" for n in dict.fromkeys(nuclides)]
        return [], ("(" + ",".join(counted) + ")", _SUFFIX_WORD[principal], len(owned_groups))
    if principal == "nitrile" and any(e["H"] for e in rest.values()):
        raise UnsupportedStructure("a nitrile nitrogen carries no hydrogen")
    nitrogen_atoms = {a for a in rest if principal in _NITROGEN_LABEL_CLASSES and mol.GetAtomWithIdx(a).GetAtomicNum() == 7}
    chalcogen_classes = ("amide", "hydrazide", "sulfonamide", *_CHALCOGEN_AMIDE_CLASSES, *_SULFONYL_AMIDE_ORDER, *_SULFONYL_HYDRAZIDE_ORDER)
    chalcogen_atoms = {
        a for a in rest
        if (_is_acid_family(principal) or principal in chalcogen_classes) and mol.GetAtomWithIdx(a).GetAtomicNum() in (8, 16, 34, 52)
    }
    if set(rest) - nitrogen_atoms - chalcogen_atoms:
        raise UnsupportedStructure("an isotopically modified atom of this characteristic group is not supported yet")
    items = []
    group_of = {a: k for k, group in enumerate(owned_groups) for a in group}
    hydrazide = principal == "hydrazide" or principal in _SULFONYL_HYDRAZIDE_ORDER
    centres = {
        n.GetIdx()
        for o in owned
        if mol.GetAtomWithIdx(o).GetAtomicNum() in (8, 16, 34, 52) and mol.GetAtomWithIdx(o).GetDegree() == 1
        for n in mol.GetAtomWithIdx(o).GetNeighbors()
    }
    roles = {}
    if nitrogen_atoms and principal in (*_CHALCOGEN_IMIDAMIDE.values(), *_CHALCOGEN_HYDRAZIDINE.values()):
        if len(owned_groups) != 1:
            raise UnsupportedStructure("several imidamide groups cannot be told apart by their N locants")
        roles = _labelled_imidamide_roles(mol, owned_groups[0], set(rest))
    for atom, entry in rest.items():
        heavy_label = entry["skeleton"]
        symbol = mol.GetAtomWithIdx(atom).GetSymbol()
        nitrogen = atom in nitrogen_atoms
        several_atoms = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == symbol) > 1
        only_heavy = not entry["H"]
        locant = symbol
        letter = symbol
        if len(owned_groups) > 1 and not nitrogen:
            locant += "'" * group_of[atom]
            letter += "'" * group_of[atom]
        if atom in roles:
            locant = roles[atom]
            letter = roles[atom]
        elif principal in _GROUP_LOCANTED_N_CLASSES and isinstance(parent_atoms, dict) and sum(
            any(_group_locant(mol, a, parent_atoms) is not None for a in g if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
            for g in owned_groups
        ) > 1:
            number = _group_locant(mol, atom, parent_atoms)
            if number is not None:
                locant += str(number)
                letter += str(number)
        if atom not in roles and (
            (hydrazide and not any(n.GetIdx() in centres for n in mol.GetAtomWithIdx(atom).GetNeighbors()))
            or (nitrogen and _imido_nitrogen(mol, atom, owned))
        ):
            locant += "'"
            letter += "'"
        if heavy_label:
            items.append((heavy_label, letter if several_atoms and only_heavy else None, 1, False))
        for nuclide, count in entry["H"].items():
            items.append((nuclide, locant, count, nitrogen and _capacity(mol, atom) > 1))
    return items, None


_GROUP_LOCANTED_N_CLASSES = frozenset(("amine", "amide", "sulfonamide", *_CHALCOGEN_AMIDE_CLASSES))


def _group_locant(mol, nitrogen, parent_atoms):
    """The locant of the parent atom that bears the group of `nitrogen`, which several groups of one parent tell apart:
    N1, N4 of a diamide or diamine (P-66.1.1.4, P-62.2.4.1.2)."""
    for neighbor in mol.GetAtomWithIdx(nitrogen).GetNeighbors():
        if neighbor.GetIdx() in parent_atoms:
            return parent_atoms[neighbor.GetIdx()]
        found = [n.GetIdx() for n in neighbor.GetNeighbors() if n.GetIdx() in parent_atoms]
        if neighbor.GetAtomicNum() in (6, 16, 34, 52) and len(found) == 1:
            return parent_atoms[found[0]]
    return None


def _labelled_imidamide_roles(mol, owned, modified):
    """{nitrogen: N | N' | N''} of the one sulfonimidamide-type group `owned`; of equivalent imido nitrogens the
    `modified` one takes the lower primes."""
    from ._substituents import ISOTOPE_LABELS

    graph = adjacency(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    token = ISOTOPE_LABELS.set(None)
    try:
        cited = _imidamide_cited(mol, graph, halogen_substituents(mol), aromatic_atoms, owned)
    finally:
        ISOTOPE_LABELS.reset(token)
    return dict(_imidamide_roles(mol, graph, owned, cited, modified))


def _imido_nitrogen(mol, atom, owned):
    """The imino nitrogen of an amidine-type group takes N', the amino nitrogen N (P-66.4.1.1); a second imido
    nitrogen on the same atom (N'') is not supported."""
    nitrogen = mol.GetAtomWithIdx(atom)
    for bond in nitrogen.GetBonds():
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        centre = bond.GetOtherAtom(nitrogen)
        amino = [
            n for n in centre.GetNeighbors()
            if n.GetIdx() in owned and n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(centre.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if not amino:
            return False
        if sum(1 for b in centre.GetBonds() if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(centre).GetAtomicNum() == 7) > 1:
            raise UnsupportedStructure("the locant N'' of a second imido nitrogen is not supported yet")
        return True
    return False


def _same_label_on_every_group_atom(mol, rest, owned, element):
    """Every atom of the element in the characteristic groups carries the same heavy nuclide and nothing else, so one
    descriptor with a count names them all and needs no locant (P-82.6.1.3)."""
    atoms = [a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == element]
    labels = {rest[a]["skeleton"] for a in atoms if a in rest}
    return (
        all(a in rest for a in atoms)
        and len(labels) == 1
        and None not in labels
        and "" not in labels
        and not any(rest[a]["H"] for a in atoms)
        and len(atoms) == len(rest)
    )


def _capacity(mol, atom):
    target = mol.GetAtomWithIdx(atom)
    return target.GetIntProp("_capacity") if target.HasProp("_capacity") else target.GetTotalNumHs()


def _sole_heteroatoms(mol, in_parent, positions):
    """Nuclides of a heteroatom that is the only atom of its element in the parent, so its locant is implied."""
    counts = {}
    for a in positions:
        number = mol.GetAtomWithIdx(a).GetAtomicNum()
        counts[number] = counts.get(number, 0) + 1
    return frozenset(
        e["skeleton"] for a, e in in_parent.items()
        if e["skeleton"] and mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and counts[mol.GetAtomWithIdx(a).GetAtomicNum()] == 1
    )


def _locants_omitted(mol, in_parent, parts):
    """P-82.6.1.3 (every parent position modified in the same way, none keeping a hydrogen; or every hydrogen of the
    compound replaced by one nuclide, acetonitrile) and the one modified atom of a bare hydrocarbon whose positions are
    all equivalent (benzene, ethane; P-82.6.1.1)."""
    positions = set(parts[4])
    if set(in_parent) == positions and len({repr(sorted(e["H"])) + str(e["skeleton"]) for e in in_parent.values()}) == 1:
        if all(sum(in_parent[a]["H"].values()) == mol.GetAtomWithIdx(a).GetTotalNumHs() for a in positions):
            return True
    if (
        in_parent
        and not any(e["skeleton"] for e in in_parent.values())
        and len({tuple(sorted(e["H"])) for e in in_parent.values()}) == 1
        and sum(sum(e["H"].values()) for e in in_parent.values()) == sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    ):
        return True
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    return (
        len(in_parent) == 1
        and sum(bool(e["skeleton"]) + sum(e["H"].values()) for e in in_parent.values()) == 1
        and len(positions) == mol.GetNumAtoms()
        and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in positions)
        and len({ranks[a] for a in positions}) == 1
    )


_GROUP_MULTIPLIER = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}


def _before_suffix(name, text, word, count=1):
    """`name` with the descriptor in front of the suffix `word`; identical modified groups are multiplied with the
    modified suffix enclosed: 'pentane-2,4-di[(18O)one]' (P-82.2.2.2)."""
    if count > 1:
        tail = _GROUP_MULTIPLIER[count] + word
        if not name.endswith(tail):
            raise UnsupportedStructure("the suffix of this name is not delimited")
        return name[: -len(tail)] + _GROUP_MULTIPLIER[count] + "[" + text + word + "]"
    if not name.endswith(word):
        raise UnsupportedStructure("the suffix of this name is not delimited")
    return name[: -len(word)] + text + word


def _check_scope(mol):
    """The specified stereo elements as ("atom", idx, code) or
    ("bond", (a, b), code); [] when the molecule has none."""
    if mol.GetNumAtoms() > _MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("this molecule is out of scope for the polyfunctional chain engine")
    if STEREO_OF_ISOTOPOLOGUE.get() is not None:
        return list(STEREO_OF_ISOTOPOLOGUE.get())
    elements = specified_stereo_elements(mol) or []
    located = []
    for kind, idx, code in elements:
        if kind == "bond":
            bond = mol.GetBondWithIdx(idx)
            located.append(("bond", (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), code))
        else:
            located.append((kind, idx, code))
    return located


def _with_descriptor(name, parts, text):
    """`text` inserted before the parent hydride part of `name`, with a hyphen before a leading locant (P-82.2.1)."""
    split = len(parts[0]) if parts[0] is not None else (parts[6] if len(parts) > 6 else None)
    if split is None:
        raise UnsupportedStructure("the parent hydride part of this name is not delimited")
    return name[:split] + text + ("-" if name[split:split + 1].isdigit() else "") + name[split:]


def _bond_locant(first, second):
    if isinstance(first, PrimedLocant) and isinstance(second, PrimedLocant) and first.primes != second.primes:
        return JunctionBondLocant(*sorted((first, second), key=lambda locant: locant.primes))
    return min(first, second)


def _stereo_entries(stereo, position_of, ring_parent=False, used=frozenset(), suffix=frozenset()):
    """[(locant, code)] for the stereo elements found on the parent, plus
    whether every element was placed (on the parent or inside a substituent). A double bond inside the atoms of the
    principal suffix group, whose carbon the parent does not number, takes no locant (P-93.4.2.1.3)."""
    entries, complete, unlocated = [], True, []
    carried = STEREO_SUFFIX.get()
    if not suffix and carried is not None:
        suffix = carried[1]
    for kind, where, code in stereo or []:
        if kind == "isotope":
            continue
        if kind == "atom":
            if where in position_of:
                entries.append((position_of[where], code))
            elif where in suffix:
                unlocated.append((_suffix_position(carried, position_of, (where,)), code))
            elif ("atom", where) not in used:
                complete = False
        else:
            a, b = where
            if ring_parent or a not in position_of or b not in position_of:
                inside = [x for x in (a, b) if x in position_of]
                if len(inside) == 2:
                    entries.append((_bond_locant(position_of[a], position_of[b]), code))
                elif ("bond", where) in used:
                    continue
                elif len(inside) == 1:
                    # P-93.4.2.1, P-93.5.1.4.2.1: the configuration of a double bond to an ylidene group takes the locant
                    # of the parent atom
                    entries.append((position_of[inside[0]], code))
                elif a in suffix and b in suffix:
                    unlocated.append((_suffix_position(carried, position_of, (a, b)), code))
                else:
                    complete = False
            else:
                entries.append((_bond_locant(position_of[a], position_of[b]), code))
    unlocated.sort(key=lambda found: (found[0] is None, found[0], found[1]))
    return sorted(entries) + [("", code) for _, code in unlocated], complete


def _suffix_position(carried, position_of, atoms):
    """The position in the parent of the atom that bears a principal group whose own atoms the parent does not number,
    which orders the descriptors that take no locant (P-93.4.2.1.2)."""
    if carried is None:
        return None
    mol, owned = carried
    for atom in atoms:
        for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
            if neighbor.GetIdx() in position_of and neighbor.GetIdx() not in owned:
                return position_of[neighbor.GetIdx()]
    return None


def _stereo_rank(stereo, position_of, ring_parent=False):
    entries, _ = _stereo_entries(stereo, position_of, ring_parent)
    marked = [(w, code) for k, w, code in stereo or [] if k == "isotope" and w in position_of]
    return (
        _isotope_counts(stereo, marked),
        *_isotope_locant_rank(stereo, position_of),
        _substituent_isotope_locants(position_of),
        tuple(0 if code in ("R", "Z", "r", "seqCis", "M", "m") else 1 for _, code in entries),
    )


def _isotope_locant_rank(stereo, position_of):
    """P-44.4.1.11.4-6: the lowest locants for the isotopically modified atoms of a parent in its own numbering, then
    for the nuclides of higher atomic number and of higher mass number."""
    marked = [(w, code) for k, w, code in stereo or [] if k == "isotope" and w in position_of]
    return tuple(sorted(position_of[w] for w, code in marked for _ in code.split("|"))), _nuclide_precedence(marked, position_of)


def _substituent_isotope_locants(position_of, mol=None, labels=None):
    """P-45.4: the parent gives the lowest locants to its isotopically modified substituent groups, then to those
    holding the nuclide of higher atomic number, then of higher mass number."""
    from ._isotope_labels import _nuclide_sort_key
    from ._substituents import ISOTOPE_LABELS, _locant_sort_key

    if mol is None:
        context = ISOTOPE_LABELS.get()
        if not context:
            return ()
        mol, labels = context["mol"], context["labels"]
    table = Chem.GetPeriodicTable()
    groups = []
    for parent_atom, locant in position_of.items():
        for neighbor in mol.GetAtomWithIdx(parent_atom).GetNeighbors():
            if neighbor.GetIdx() in position_of:
                continue
            group, stack = {neighbor.GetIdx()}, [neighbor.GetIdx()]
            while stack:
                for onward in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                    if onward.GetIdx() not in group and onward.GetIdx() not in position_of:
                        group.add(onward.GetIdx())
                        stack.append(onward.GetIdx())
            found = {n for a in group if a in labels for n in [labels[a]["skeleton"], *labels[a]["H"]] if n}
            if found:
                groups.append((locant, found))
    order = sorted(
        {n for _, found in groups for n in found},
        key=lambda n: (-table.GetAtomicNumber(_nuclide_sort_key(n)[0]), -_nuclide_sort_key(n)[1]),
    )
    return (
        tuple(sorted(_locant_sort_key(locant) for locant, _ in groups)),
        tuple(tuple(sorted(_locant_sort_key(locant) for locant, found in groups if n in found)) for n in order),
    )


def _parent_configuration_rank(stereo, position_of, ring_parent=False, directional=True):
    """P-44.4.1.11.1-3, P-44.4.1.12: between candidate parents the one with more isotopically modified atoms is senior,
    then the one with more Z double bonds (and lower locants for them), more like descriptor pairs, and R ahead of S.
    A chain numbered in either direction leaves R ahead of S to the stereo rank, after the isotopic locants (P-31.1.4.3.4)."""
    entries, _ = _stereo_entries(stereo, position_of, ring_parent)
    marked = [(w, code) for k, w, code in stereo or [] if k == "isotope" and w in position_of]
    cis = [locant for locant, code in entries if code == "Z" and locant != ""]
    centres = [code for _, code in entries if code in ("R", "S", "r", "s")]
    like = sum(1 for code in centres[1:] if code == centres[0])
    return (
        _isotope_counts(stereo, marked),
        -len(cis),
        tuple(cis),
        -like,
        tuple(0 if code in "Rr" else 1 for code in centres) if directional else (),
    )


def _isotope_counts(stereo, marked):
    """P-44.4.1.11.1-3: the parent with more isotopically modified atoms is senior, then the one with more nuclides of
    higher atomic number, then of higher mass number; smaller keys are senior."""
    from ._isotope_labels import _nuclide_sort_key
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context or not marked:
        return 0, ()
    counts = {}
    for atom, _ in marked:
        entry = context["labels"].get(atom)
        if entry is None:
            continue
        if entry["skeleton"]:
            counts[entry["skeleton"]] = counts.get(entry["skeleton"], 0) + 1
        for nuclide, number in entry["H"].items():
            counts[nuclide] = counts.get(nuclide, 0) + number
    table = Chem.GetPeriodicTable()
    everywhere = {entry["skeleton"] for entry in context["labels"].values() if entry["skeleton"]}
    everywhere |= {n for entry in context["labels"].values() for n in entry["H"]}
    order = sorted(
        everywhere,
        key=lambda n: (-table.GetAtomicNumber(_nuclide_sort_key(n)[0]), -_nuclide_sort_key(n)[1]),
    )
    return -sum(counts.values()), tuple(-counts.get(n, 0) for n in order)


def _nuclide_precedence(marked, position_of):
    """Lowest locants to the nuclide of higher atomic number, then of higher mass number (P-82.5.2)."""
    table = Chem.GetPeriodicTable()
    located = {}
    for atom, code in marked:
        for nuclide in filter(None, code.split("|")):
            symbol = "".join(ch for ch in nuclide if ch.isalpha())
            located.setdefault((-table.GetAtomicNumber(symbol), -int("".join(ch for ch in nuclide if ch.isdigit()))), []).append(position_of[atom])
    return tuple(tuple(sorted(located[key])) for key in sorted(located))


def _stereo_prefix(stereo, position_of, ring_parent=False, used=frozenset(), suffix=frozenset()):
    if not stereo:
        return ""
    entries, _ = _stereo_entries(stereo, position_of, ring_parent, used, suffix)
    if not entries:
        return ""
    if len(entries) == 1 and entries[0][1].startswith("("):
        return f"{entries[0][1]}-"
    if len(position_of) == 1 and len(entries) == 1:
        # P-93.5: the one skeletal atom of a mononuclear parent needs no locant
        return f"({entries[0][1]})-"
    text = ",".join(f"{locant}{code}" for locant, code in entries)
    return f"[{text}]-" if "(" in text else f"({text})-"


def _select(mol, attach=None, n_names=(), stereo=None):
    classes = {_group_of(mol, a.GetIdx())[0] for a in mol.GetAtoms() if _group_of(mol, a.GetIdx())} | {
        c for c, _, _ in _ring_occurrences(mol)
    }
    principal = _principal_class(classes)
    token = EXTENDED_PREFIXES.set(
        bool(AMINIUM.get() or RING_CENTER.get() or any(a.HasProp(ANION_PROP) for a in mol.GetAtoms()))
        or (_is_acid_family(principal) or principal in ("thioic", "peroxoic", "imidic") if principal else False)
    )
    try:
        return _select_with_prefixes(mol, attach, n_names, stereo)
    finally:
        EXTENDED_PREFIXES.reset(token)


def _select_with_prefixes(mol, attach=None, n_names=(), stereo=None):
    """(key, name, parts) of the best parent. `attach`: an atom that carries a
    free valence (a multiplicative unit); it takes the lowest locant after
    the principal groups and multiple bonds."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())

    groups = {}
    for atom in mol.GetAtoms():
        found = _group_of(mol, atom.GetIdx())
        if found is not None:
            groups.setdefault(found[0], {})[atom.GetIdx()] = found[1]
    ring_groups = _ring_occurrences(mol)
    marked = {a.GetIdx() for a in mol.GetAtoms() if a.HasProp(ANION_PROP)}
    ide_extra = {}
    if marked:
        for idx in marked:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7 and any(
                b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds()
            ):
                carbon = next(b.GetOtherAtom(atom).GetIdx() for b in atom.GetBonds())
                groups.setdefault("imine", {})[carbon] = {idx}
            _peroxy_group(mol, atom, groups, ring_groups)
            _carbonyl_variant_group(mol, atom, groups, ring_groups)
            if atom.GetAtomicNum() != 6:
                continue
            if atom.IsInRing():
                ring_groups += [("ide", idx, {idx})] * anion_weight(atom)
            else:
                groups.setdefault("ide", {})[idx] = {idx}
        carbon_marks = {i for i in marked if mol.GetAtomWithIdx(i).GetAtomicNum() == 6}
        has_group = any(cls != "ide" and any(o & marked - carbon_marks for o in members.values()) for cls, members in groups.items()) or any(
            g[0] != "ide" and g[2] & marked - carbon_marks for g in ring_groups
        )
        if has_group and carbon_marks:
            groups.pop("ide", None)
            ring_groups = [g for g in ring_groups if g[0] != "ide"]
            ide_extra = {i: anion_weight(mol.GetAtomWithIdx(i)) for i in carbon_marks}
        groups = {
            cls: kept
            for cls, members in groups.items()
            if (kept := {c: owned for c, owned in members.items() if owned & marked})
        }
        ring_groups = [g for g in ring_groups if g[2] & marked]
        groups, ring_groups = _senior_anionic_centres(mol, groups, ring_groups, marked)
        secondary = any(_substituted_amine_nitrogen(mol, mol.GetAtomWithIdx(i)) for i in marked)
        if not groups and not ring_groups and not secondary:
            raise UnsupportedStructure("this anionic group is not a principal-capable characteristic group")
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = _principal_class(classes)
    token = IDE_EXTRA.set(ide_extra)
    imine_token = IMINE_PARENT.set(principal == "imine")
    try:
        results = [
            _select_with_principal(mol, graph, halogens, aromatic_atoms, kept, kept_rings, principal, attach, n_names, stereo)
            for amidine_kept, amidine_rings in _diamidide_alternatives(mol, groups, ring_groups)
            for kept, kept_rings in _joined_imine_alternatives(mol, graph, amidine_kept, amidine_rings)
        ]
        best = min(results, key=lambda result: result[0])
    finally:
        IMINE_PARENT.reset(imine_token)
        IDE_EXTRA.reset(token)
    if ide_extra:
        anions = _carbanion_parent(mol, graph, halogens, aromatic_atoms, ide_extra, attach, n_names, stereo)
        if anions is not None and _prefer_anion_parent(mol, best, anions):
            return anions
    return best


def _carbanion_parent(mol, graph, halogens, aromatic_atoms, centres, attach, n_names, stereo):
    """The parent named by its carbanide centres alone, the anionic acid groups being prefixes (P-72.7b)."""
    ide_groups = {"ide": {idx: {idx} for idx in centres if not mol.GetAtomWithIdx(idx).IsInRing()}}
    ide_rings = [("ide", idx, {idx}) for idx, weight in centres.items() if mol.GetAtomWithIdx(idx).IsInRing() for _ in range(weight)]
    token = IDE_EXTRA.set({})
    try:
        return _select_with_principal(
            mol, graph, halogens, aromatic_atoms, {k: v for k, v in ide_groups.items() if v}, ide_rings, "ide", attach, n_names, stereo
        )
    except UnsupportedStructure:
        return None
    finally:
        IDE_EXTRA.reset(token)


def _anionic_centres(result):
    return -result[0][0] if isinstance(result[0][0], int) else None


def _ide_centres_in_parent(mol, result):
    placed = result[2][4] if isinstance(result[2], tuple) and len(result[2]) > 4 and isinstance(result[2][4], dict) else {}
    return sum(
        anion_weight(atom) for i in placed if (atom := mol.GetAtomWithIdx(i)).HasProp(ANION_PROP) and atom.GetAtomicNum() == 6
    )


def _prefer_anion_parent(mol, group_parent, ide_parent):
    """P-72.7: the parent with more anionic centres, suffixes included, is senior (a); then the one with more 'ide'
    centres (b)."""
    group_total, ide_total = _anionic_centres(group_parent), _anionic_centres(ide_parent)
    if group_total is None or ide_total is None:
        return False
    if ide_total != group_total:
        return ide_total > group_total
    return _ide_centres_in_parent(mol, ide_parent) > _ide_centres_in_parent(mol, group_parent)


def _diamidide_alternatives(mol, groups, ring_groups):
    """The (group table, ring group list) pairs in which each pair of amidine groups sharing an amino nitrogen keeps
    one group: the other carbon is cited as an N-imidoyl prefix (P-66.4.1.6)."""
    nitrogens = {}
    for c, owned in groups.get("amidine", {}).items():
        nitrogens[c] = {a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7}
    for g in ring_groups:
        if g[0] == "amidine":
            carbon = next((a for a in g[2] if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), g[1])
            nitrogens[carbon] = {a for a in g[2] if mol.GetAtomWithIdx(a).GetAtomicNum() == 7}
    shared = [(x, y) for x, y in itertools.combinations(sorted(nitrogens), 2) if nitrogens[x] & nitrogens[y]]
    if not shared:
        return [(groups, ring_groups)]
    if len(shared) > 1:
        raise UnsupportedStructure("several amidine groups sharing nitrogen atoms are not supported yet")
    alternatives = []
    for dropped in shared[0]:
        kept_groups = dict(groups)
        kept_groups["amidine"] = {c: o for c, o in groups.get("amidine", {}).items() if c != dropped}
        kept_rings = [
            g
            for g in ring_groups
            if not (g[0] == "amidine" and next((a for a in g[2] if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), g[1]) == dropped)
        ]
        alternatives.append((kept_groups, kept_rings))
    return alternatives


def _joined_imine_alternatives(mol, graph, groups, ring_groups):
    """The (group table, ring group list) pairs in which two imine groups joined through the nitrogen of one of them
    keep one group: the other C=N is cited as a (ylideneamino) prefix on the nitrogen it substitutes (P-62.3.1.2)."""
    members = {c: next(iter(owned)) for c, owned in groups.get("imine", {}).items()}
    for cls, ring_atom, owned in ring_groups:
        if cls == "imine":
            members[ring_atom] = next(iter(owned))
    joined = [
        (x, y)
        for x, y in itertools.combinations(sorted(members), 2)
        if any(
            members[b] in _arm_atoms(graph, n, members[a])
            for a, b in ((x, y), (y, x))
            for n in graph[members[a]]
            if n != a
        )
    ]
    if not joined:
        return [(groups, ring_groups)]
    if len(joined) > 1:
        raise UnsupportedStructure("several imine groups joined through nitrogen atoms are not supported yet")
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    if ranks[joined[0][0]] == ranks[joined[0][1]]:
        raise UnsupportedStructure("imine groups joined through nitrogen need a multiplicative name (P-15.3)")
    alternatives = []
    for dropped in joined[0]:
        kept_groups = dict(groups)
        kept_groups["imine"] = {c: o for c, o in groups.get("imine", {}).items() if c != dropped}
        kept_rings = [g for g in ring_groups if not (g[0] == "imine" and g[1] == dropped)]
        alternatives.append((kept_groups, kept_rings))
    return alternatives


_AMINE_STEREO_ELEMENTS = {6, 7, 8, 9, 14, 16, 17, 32, 34, 35, 52, 53}


def _skeletal_heteroatom(atom):
    """An atom of another element that links two groups is a heteroatom of the chain, named by replacement (P-15.4.3)."""
    return atom.GetAtomicNum() not in _AMINE_STEREO_ELEMENTS and sum(n.GetDegree() > 1 for n in atom.GetNeighbors()) > 1


def _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal, attach, n_names, stereo):
    if (
        any(_chalcogen_ketone(mol, a) for a in mol.GetAtoms())
        and principal not in _CHALCOGEN_KETONE_OK
        and not (principal and _is_variant(principal))
        and not (AMINIUM.get() or RING_CENTER.get() or any(a.HasProp(ANION_PROP) for a in mol.GetAtoms()))
    ):
        raise UnsupportedStructure("a thioketone-type group outranks the parents this engine can build here")
    for atom in mol.GetAtoms():
        if atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("isotopes and radicals are not supported by the polyfunctional chain engine")
        if atom.GetFormalCharge() and not _is_nitro_part(atom) and not _anionic_group_atom(mol, atom) and not (
            AMINIUM.get() and atom.GetAtomicNum() == 7
        ) and not (RING_CENTER.get() and _is_chain_amine_oxide_part(atom)) and not (
            _cationic_prefix_nitrogen(atom)
            and any(_anionic_group_atom(mol, a) or a.HasProp(ANION_PROP) for a in mol.GetAtoms())
        ):
            raise UnsupportedStructure("charged atoms are not supported by the polyfunctional chain engine")
        if (
            atom.GetAtomicNum() == 6
            and not atom.IsInRing()
            and not _is_acid_family(principal)
            and not (AMINIUM.get() and principal in (None, "amine", "imine"))
            and not RING_CENTER.get()
            and principal not in ("peroxoic", "thioic", "imidic")
            and not (
                principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES, "hydrazide", *_CHALCOGEN_HYDRAZIDE.values(), *_CHALCOGEN_HYDRAZIDINE.values())
                and atom.GetIdx() in groups.get(principal, {})
            )
            and not (principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "amidine", *_AMIDRAZONE) and atom.GetIdx() in groups.get("hydrazide", {}))
            and _is_ester_like(mol, atom.GetIdx())
            and not (principal in _AMIDE_FAMILY and _nitrogen_acyl(mol, atom.GetIdx()))
            and not (_urea_carbon(mol, atom.GetIdx()) and _outranks_urea(principal))
            and not (FORCED_PRINCIPAL.get() == principal and principal == "nitrile")
            and not (FORCED_PRINCIPAL.get() in ("amine", "imine") and principal in (None, FORCED_PRINCIPAL.get()))
            and principal != "ide"
            and not any(a.HasProp(ANION_PROP) for a in mol.GetAtoms())
        ):
            raise UnsupportedStructure("an ester outranks every parent this engine can build except an acid")
    if principal in (None, "amine") and attach is None and not n_names and not RING_CENTER.get():
        boranyl = _boranyl_amines(mol, graph, halogens, aromatic_atoms)
        if boranyl is not None:
            return boranyl
    if (
        principal in (None, "amine")
        and not SUBSTITUTED_AMINE_PREFIX.get()
        and not (RING_CENTER.get() and not AMINIUM.get())
        and any(_substituted_amine_nitrogen(mol, a) for a in mol.GetAtoms())
    ):
        if attach is not None:
            raise UnsupportedStructure("N-substituted amines inside a unit are not handled by the chain engine")
        if any(kind != "isotope" for kind, _, _ in stereo or ()) and any(_skeletal_heteroatom(a) for a in mol.GetAtoms()):
            raise UnsupportedStructure("stereodescriptors with an N-substituted amine parent are not supported beside other heteroatoms yet")
        return _substituted_amine(mol, graph, halogens, aromatic_atoms, groups, ring_groups)

    if principal is None:
        if attach is not None:
            raise UnsupportedStructure("a unit without a principal group is not supported")
        return _plain_parent(mol, graph, halogens, aromatic_atoms, stereo)

    if _is_variant(principal):
        imidic_ns = _imidic_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if imidic_ns:
            n_names = imidic_ns

    if principal in ("amide", "nitrile", *_CHALCOGEN_AMIDE_CLASSES, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES):
        amide_ns = _amide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if amide_ns:
            n_names = amide_ns

    if principal in _CHALCOGEN_IMIDAMIDE.values() or principal in _CHALCOGEN_HYDRAZIDINE.values():
        imidamide_ns = _imidamide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if imidamide_ns:
            n_names = imidamide_ns

    if principal == "hydrazide" or principal in _CHALCOGEN_HYDRAZIDE.values():
        hydrazide_ns = _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if hydrazide_ns:
            n_names = hydrazide_ns

    if principal in _AMIDRAZONE:
        amidrazone_ns = _amidrazone_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if amidrazone_ns:
            n_names = amidrazone_ns

    if principal == "amidine":
        amidine_ns = _amidine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups)
        if amidine_ns:
            n_names = amidine_ns

    if principal == "imine":
        imine_ns = _imine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups)
        if imine_ns:
            n_names = imine_ns

    anchors = set(groups.get(principal, {}))
    for cls, _, owned in ring_groups:
        if cls == principal:
            anchors.add(min(owned, key=lambda a: (mol.GetAtomWithIdx(a).GetAtomicNum() != 6, a)))
    total_principal = sum(
        _group_weight(mol, a, principal) if mol.GetAtomWithIdx(a).GetAtomicNum() == 6 else 1 for a in anchors
    )
    group_atoms = set().union(*groups.get(principal, {}).values(), *(g[2] for g in ring_groups if g[0] == principal))
    group_atoms |= set(groups.get(principal, {}))

    def _finish(result):
        carried = -result[0][0] if isinstance(result[0][0], int) else 0
        if attach is None and not n_names and carried < total_principal and _identical_group_units(
            mol, graph, group_atoms
        ):
            raise UnsupportedStructure(
                "identical parents joined through a linking group need a multiplicative name (P-15.3)"
            )
        return result

    chain_best = None
    chain_error = None
    stereo_context = BRANCH_STEREO.get()
    used_before = set(stereo_context["used"]) if stereo_context else set()
    chain_used = set()
    from ._substituents import ISOTOPE_LABELS

    isotope_context = ISOTOPE_LABELS.get()
    consumed_before = set(isotope_context["consumed"]) if isotope_context else set()
    chain_consumed = set()
    principal_atoms = groups.get(principal, {})
    if principal_atoms:
        owned = set().union(*principal_atoms.values())
        eligible = {
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and not a.IsInRing()
            and (
                a.GetIdx() in principal_atoms
                or not is_functional_carbon(mol, a.GetIdx())
                or (
                    (
                        _is_acid_family(principal)
                        or principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "hydrazide", "amidine", *_AMIDRAZONE)
                        or ((AMINIUM.get() or FORCED_PRINCIPAL.get() in ("amine", "imine")) and principal in ("amine", "imine"))
                    )
                    and _junior_end_group(mol, a.GetIdx())
                )
            )
        }
        try:
            for path in _paths(graph, eligible):
                for chain in (path, path[::-1]):
                    if stereo_context:
                        stereo_context["used"] &= used_before
                    if isotope_context:
                        isotope_context["consumed"] &= consumed_before
                    candidate = _evaluate(
                        mol, graph, halogens, aromatic_atoms, chain, principal, principal_atoms, owned, attach, n_names, stereo
                    )
                    if chain_best is None or candidate[0] < chain_best[0]:
                        chain_best = candidate
                        chain_used = set(stereo_context["used"]) - used_before if stereo_context else set()
                        chain_consumed = set(isotope_context["consumed"]) - consumed_before if isotope_context else set()
        except UnsupportedStructure as error:
            chain_error, chain_best = error, None
        if isotope_context:
            isotope_context["consumed"] &= consumed_before
        if chain_best is not None and chain_best[0][0] == 1:
            chain_best = None
    chain_count = -chain_best[0][0] if chain_best else 0
    if stereo_context:
        stereo_context["used"] &= used_before

    assembly = _assembly_parent(
        mol, graph, halogens, aromatic_atoms, principal, [g for g in ring_groups if g[0] == principal], stereo, n_names
    )
    if assembly is not None and assembly[0] >= chain_count:
        if attach is not None:
            raise UnsupportedStructure("a ring assembly inside a unit is not supported yet")
        return _finish(assembly[1])

    ring_best = _best_ring(
        mol, graph, halogens, aromatic_atoms, principal, [g for g in ring_groups if g[0] == principal], n_names, stereo
    )
    ring_count = ring_best[0] if ring_best else 0

    ring_needed = total_principal
    if chain_error is not None and ring_count:
        chain_capacity = max(
            (
                sum(
                    anion_weight(mol.GetAtomWithIdx(a)) if principal == "ide" else _group_weight(mol, a, principal)
                    for a in path
                    if a in principal_atoms
                )
                for path in _paths(graph, eligible)
            ),
            default=0,
        )
        ring_needed = max(min(total_principal, chain_capacity), 1)
    if chain_error is not None and not (ring_count and ring_count >= ring_needed):
        raise chain_error
    if ring_count and ring_count >= chain_count:
        if attach is not None:
            raise UnsupportedStructure("a ring parent inside a multiplicative unit is not supported yet")
        return _finish(ring_best[1])
    if chain_best is None:
        raise UnsupportedStructure("no parent carries the principal group")
    if stereo_context:
        stereo_context["used"] |= chain_used
    if isotope_context:
        isotope_context["consumed"] |= chain_consumed
    if _is_terminal(principal) and chain_count < len(principal_atoms):
        used_with_chain = set(stereo_context["used"]) if stereo_context else set()
        carbo = _carbo_best(mol, graph, halogens, aromatic_atoms, principal, principal_atoms, stereo, n_names)
        if carbo is not None and -carbo[0][0] > chain_count:
            return _finish(carbo)
        if stereo_context:
            stereo_context["used"] &= used_with_chain
    if (
        _is_terminal(principal)
        and not (principal in ("amidine", *_AMIDRAZONE) and attach is None)
        and len(chain_best[2][4]) == 1
        and not (_is_variant(principal) and attach is None)
        and not chain_best[1].endswith("formic acid")
        and not (principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "hydrazide") and attach is None)
    ):
        raise UnsupportedStructure("one-carbon acid, amide, nitrile and aldehyde parents use retained names")
    return _finish(chain_best)


def _isolated_candidates(evaluations):
    """[(result, consumed nuclide labels)] of the evaluations that succeed. Each one starts from the same consumed set,
    so the outcome never depends on the order in which equally senior parents are tried."""
    from ._substituents import _labels_consumed

    consumed = _labels_consumed()
    base = set(consumed)
    found = []
    for evaluate in evaluations:
        consumed.clear()
        consumed.update(base)
        try:
            found.append((evaluate(), set(consumed)))
        except UnsupportedStructure:
            continue
    consumed.clear()
    consumed.update(base)
    return found


def _choose_isolated(candidates, key):
    from ._substituents import _labels_consumed

    result, gained = min(candidates, key=lambda c: key(c[0]))
    consumed = _labels_consumed()
    consumed.clear()
    consumed.update(gained)
    return result


def _plain_ring_parent(mol, graph, halogens, aromatic_atoms, ring, stereo):
    """(sort key, result) for `ring` as the parent of a molecule without a
    principal group; the key ranks equally senior rings by their lowest
    substituent locants, then by name."""
    ring_set = set(ring)
    roots = [
        (r, n.GetIdx())
        for r in ring
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring_set
    ]
    from ._substituents import ISOTOPE_LABELS

    labels = (ISOTOPE_LABELS.get() or {}).get("labels", {})
    modified = _isotope_counts(stereo, [(a, None) for a in ring_set if a in labels])
    spec = monocycle_spec(mol, ring)
    if spec is None:
        from ._diester_ring_diyl import CATION_CENTRES, evaluate_skeleton

        cation = CATION_RING_ATOM.get()
        centres_token = CATION_CENTRES.set(((cation, "ium"),)) if cation in ring_set else None
        try:
            found = evaluate_skeleton(mol, graph, "ring", [ring], ring_set, [], set(), "")
        finally:
            if centres_token is not None:
                CATION_CENTRES.reset(centres_token)
        if found is None:
            raise UnsupportedStructure("this ring has no supported name")
        from ._diester_ring_diyl import PARENT_START

        placed = found[2]
        name = _without_stereo(found[1])
        return (modified, _isotope_locant_rank(stereo, placed), _parent_configuration_rank(stereo, placed), -len(roots), tuple(sorted(placed[r] for r, _ in roots)), alphanumerical_name_key(name), name), (
            (0,), name, (None, None, None, 0, placed, True, PARENT_START.get())
        )
    from ._substituents import ISOTOPE_LABELS

    isotope_context = ISOTOPE_LABELS.get()
    ring_labelled = isotope_context is not None and any(a in ring_set for a in isotope_context["labels"])
    if not roots:
        if not ring_labelled:
            raise UnsupportedStructure("an unsubstituted ring is not a polyfunctional case")
        locants = min(numberings(spec), key=lambda option: _stereo_rank(stereo, option))
        return (modified, _isotope_locant_rank(stereo, locants), _parent_configuration_rank(stereo, locants), 0, (), spec.parent), ((0,), spec.parent, (None, None, None, 0, locants, True, 0))
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    best = None
    for locants in numberings(spec):
        key = (
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
            _stereo_rank(stereo, locants),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    if len(entries) == 1 and spec.hetero is None and not ring_labelled:
        _, only_name, only_compound = entries[0]
        prefix_text = format_substituent_prefixes(
            {only_name: {"locants": [1], "compound": only_compound}}, omit_locants=True
        )
    else:
        prefix_text = _prefix_text(entries, best[1])
    name = _join(prefix_text, spec.parent)
    return (modified, _isotope_locant_rank(stereo, best[1]), _parent_configuration_rank(stereo, best[1]), -len(roots), best[0][0], alphanumerical_name_key(name), name), ((0,), name, (None, None, None, 0, best[1], True, len(name) - len(spec.parent)))


_GROUP_14 = (14, 32, 50, 82)
_GROUP_13_METALS = (13, 31, 49, 81)
_CHALCOGENOL_WORDS = {8: "ol", 16: "thiol", 34: "selenol", 52: "tellurol"}


_PARENT_HYDRIDE_ORDER = (7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81)


def _junior_hydride_atom(mol, atom, center_z):
    """Whether `atom` is a mononuclear hydride atom of an element junior to the parent hydride `center_z` (P-68.1.5.2.3):
    a boranyl group on silicon, not the other way round."""
    z = mol.GetAtomWithIdx(atom).GetAtomicNum()
    return (
        z in _PARENT_HYDRIDE_ORDER
        and center_z in _PARENT_HYDRIDE_ORDER
        and _PARENT_HYDRIDE_ORDER.index(z) > _PARENT_HYDRIDE_ORDER.index(center_z)
    )


_NITROGEN_GROUP_PREFIXES = {
    "nitro", "nitroso", "azido", "isocyano", "isocyanato", "isothiocyanato", "isoselenocyanato", "isotellurocyanato",
}


def _nitrogen_group_prefix(mol, graph, atom, parent, halogens, aromatic_atoms):
    """Whether the nitrogen `atom` bonded to the hydride atom `parent` is a nitro, nitroso, azido or isocyano group."""
    if mol.GetAtomWithIdx(atom).GetAtomicNum() != 7:
        return False
    try:
        name, _ = name_branch(graph, atom, parent, halogens, aromatic_atoms, mol=mol)
    except UnsupportedStructure:
        return False
    return name in _NITROGEN_GROUP_PREFIXES


def _mononuclear_parent(mol, graph, halogens, aromatic_atoms, center):
    """A single Si, Ge, P, B, ... atom is the senior parent hydride when there is no principal group (P-44.1.2):
    'trimethyl(phenyl)silane', 'methoxy(trimethyl)silane'. On a Group 14 atom a hydroxy or amino group is the
    suffix (P-68.2): 'trimethylsilanol', 'trimethylsilanamine' ('1,1,1-trimethyl-N-(trimethylsilyl)silanamine' once N-locants occur)."""
    from ._substituents import centre_nuclides, format_mononuclear_prefixes

    z = center.GetAtomicNum()
    stem, _, valence = MONONUCLEAR_HYDRIDES[z]
    index = center.GetIdx()
    nuclide = centre_nuclides(index, consume=False)
    neighbors = [n for n in graph[index]]
    if (
        center.GetFormalCharge()
        or center.GetIsotope()
        or sum(mol.GetBondBetweenAtoms(index, n).GetBondTypeAsDouble() for n in neighbors) > valence
        or any(
            mol.GetBondBetweenAtoms(index, n).GetBondTypeAsDouble() != 1.0
            for n in neighbors
            if mol.GetAtomWithIdx(n).GetAtomicNum() != 6 and not _junior_hydride_atom(mol, n, z)
        )
    ):
        return None
    chalcogenols = {
        word: [n for n in neighbors if _terminal_heteroatom(mol, n, 1) and mol.GetAtomWithIdx(n).GetAtomicNum() == z]
        for z, word in _CHALCOGENOL_WORDS.items()
    }
    principal_word = next((w for w, atoms in chalcogenols.items() if atoms), None)
    junior_chalcogenols = {n for w, atoms in chalcogenols.items() if w != principal_word for n in atoms}
    amines = [
        n
        for n in neighbors
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 7
        and not mol.GetAtomWithIdx(n).GetFormalCharge()
        and not mol.GetAtomWithIdx(n).IsInRing()
        and all(
            mol.GetAtomWithIdx(m).GetAtomicNum() in (6, 8, 16, 34, 52, *MONONUCLEAR_HYDRIDES) or _terminal_heteroatom(mol, m, 1)
            for m in graph[n]
            if m != index
        )
        and not _nitrogen_group_prefix(mol, graph, n, index, halogens, aromatic_atoms)
    ]
    suffix_atoms = chalcogenols[principal_word] if principal_word else []
    amino_prefixes = amines if suffix_atoms else []
    if suffix_atoms:
        amines = []
    others = [n for n in neighbors if n not in suffix_atoms and n not in amines]
    if any(
        mol.GetAtomWithIdx(n).GetAtomicNum() not in (6, 8, 16, 34, 52, *HALOGEN_PREFIXES)
        and n not in junior_chalcogenols
        and n not in amino_prefixes
        and not _nitrogen_group_prefix(mol, graph, n, index, halogens, aromatic_atoms)
        and not _junior_hydride_atom(mol, n, z)
        for n in others
    ):
        return None
    if (suffix_atoms or amines) and (
        z not in (*_GROUP_14, *_GROUP_13_METALS) and not (amines and not suffix_atoms and z == 5)
    ):
        return None
    entries = [name_branch(graph, n, index, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in others]
    if suffix_atoms:
        suffix = (
            stem[:-1] + "ol" if principal_word == "ol" and len(suffix_atoms) == 1 else stem + multiplied_word(len(suffix_atoms), principal_word)
        )
        return format_mononuclear_prefixes(entries) + nuclide + suffix
    if amines:
        grouped = group_substituents({1: entries} if entries else {})
        n_entries = [
            (*name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True), "N" + "'" * k)
            for k, nitrogen in enumerate(amines)
            for n in graph[nitrogen]
            if n != index
        ]
        merged = _with_n_names(grouped, n_entries)
        suffix = stem[:-1] + "amine" if len(amines) == 1 else stem + multiplied_word(len(amines), "amine")
        return format_substituent_prefixes(merged, omit_locants=len(amines) > 1 and not n_entries) + nuclide + suffix
    return format_mononuclear_prefixes(entries) + nuclide + stem


def _hydroxy_suffix_count(name):
    """The number of hydroxy or chalcogenol suffixes of a hydride name: stannanol 1, silanetriol 3 (P-41, P-44.1.1)."""
    match = re.search(r"(di|tri|tetra)?(?:ol|thiol|selenol|tellurol)$", name)
    return {"di": 2, "tri": 3, "tetra": 4}.get(match.group(1), 1) if match else 0


def _amine_count(name):
    """The number of amine suffixes of a hydride name: boranamine 1, boranediamine 2 (P-44.1.1)."""
    match = re.search(r"(di|tri|tetra)amine$", name)
    return {"di": 2, "tri": 3, "tetra": 4}[match.group(1)] if match else int(name.endswith("amine"))


def _azane_ylidyne_parent(mol, graph, halogens, aromatic_atoms):
    """A nitrogen joined by a triple bond to a mononuclear hydride atom is the senior hydride, azane (P-44.1.2): N#P is
    'phosphanylidyneazane'."""
    from ._substituents import format_mononuclear_prefixes

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 1:
            continue
        other = atom.GetNeighbors()[0]
        if (
            other.GetAtomicNum() in MONONUCLEAR_HYDRIDES
            and not other.IsInRing()
            and mol.GetBondBetweenAtoms(atom.GetIdx(), other.GetIdx()).GetBondTypeAsDouble() == 3.0
        ):
            entry = name_branch(graph, other.GetIdx(), atom.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
            return format_mononuclear_prefixes([entry]) + "azane"
    return None


_HYDROXYLAMINE_PREFIX = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}
_CHALCOGEN_SYMBOL_LOCANT = {8: "O", 16: "S", 34: "Se", 52: "Te"}


def _hydroxylamine_parent(mol, graph, halogens, aromatic_atoms):
    """P-68.3.1.1.1.2-3: a nitrogen bonded to chalcogens and hydrogen only is the centre of hydroxylamine, a functional
    parent with full substitution: 'N-methoxy-O-methylhydroxylamine', the other chalcogen groups N prefixes."""
    from ._substituents import format_substituent_prefixes

    nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 1 or mol.GetRingInfo().NumRings():
        return None
    nitrogen = nitrogens[0]
    n_idx = nitrogen.GetIdx()
    neighbors = [n.GetIdx() for n in nitrogen.GetNeighbors()]
    if (
        nitrogen.GetFormalCharge()
        or len(neighbors) < 2
        or any(mol.GetAtomWithIdx(n).GetAtomicNum() not in _HYDROXYLAMINE_PREFIX for n in neighbors)
        or any(mol.GetBondBetweenAtoms(n_idx, n).GetBondTypeAsDouble() != 1.0 for n in neighbors)
        or any(mol.GetAtomWithIdx(n).GetFormalCharge() for n in neighbors)
    ):
        return None
    senior = min((mol.GetAtomWithIdx(n).GetAtomicNum() for n in neighbors), key=list(_HYDROXYLAMINE_PREFIX).index)
    candidates = []
    for centre in (n for n in neighbors if mol.GetAtomWithIdx(n).GetAtomicNum() == senior):
        entries = []
        for n in neighbors:
            if n != centre:
                name, compound = name_branch(graph, n, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.append(("N", name, compound))
        for r in graph[centre]:
            if r != n_idx:
                name, compound = name_branch(graph, r, centre, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.append((_CHALCOGEN_SYMBOL_LOCANT[senior], name, compound))
        entries.sort(key=lambda e: alpha_sort_key(e[1]))
        grouped = {}
        for locant, name, compound in entries:
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
        parent = f"{_HYDROXYLAMINE_PREFIX[senior]}hydroxylamine"
        if senior != 8 and grouped:
            parent = f"({parent})"
        text = format_substituent_prefixes(grouped) + parent
        candidates.append(([e[0] for e in entries], text))
    return min(candidates)[1]


def _hydrazine_parent(mol, aromatic_atoms):
    """The hydrazine as the parent of a molecule without a principal group and without a ring that contains nitrogen
    (P-44.1.2, P-68.3.1.2.1): one plain acyclic N-N bond whose nitrogens carry only carbon groups."""
    from ._hydrazine import hydrazine_with_substituents

    if any(a.GetAtomicNum() == 7 and a.IsInRing() for a in mol.GetAtoms()) or not _only_plain_hydrazines(mol):
        return None
    pairs = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in mol.GetBonds()
        if b.GetBeginAtom().GetAtomicNum() == 7 and b.GetEndAtom().GetAtomicNum() == 7 and not b.IsInRing()
    ]
    if len(pairs) != 1:
        return None
    return hydrazine_with_substituents(mol, *pairs[0], aromatic_atoms, unsaturated=True)


def _joins_two_rings(bond):
    """A single bond between ring atoms of two rings: an N-N bond of this kind joins ring systems, it is no hydrazine."""
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return bond.GetBondTypeAsDouble() == 1.0 and a.IsInRing() and b.IsInRing() and not bond.IsInRing()


def _plain_diazene(bond):
    """An acyclic N=N bond whose nitrogens carry only a carbon group or hydrogen, or one of them an amino nitrogen that
    carries only carbon groups (a triaz-1-ene): a ring that contains nitrogen is the parent beside it, the chain being
    the prefix 'phenyldiazenyl' or '(dimethylamino)diazenyl' (P-44.1.2.2, P-68.3.1.3)."""
    atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
    inner = {x.GetIdx() for x in atoms}
    if bond.GetBondTypeAsDouble() != 2.0:
        return False
    amino = 0
    for a in atoms:
        if a.IsInRing() or a.GetFormalCharge() or a.GetDegree() > 2:
            return False
        for n in a.GetNeighbors():
            if n.GetIdx() in inner or n.GetAtomicNum() == 6 or (n.GetAtomicNum() == 7 and n.IsInRing()):
                continue
            if not (
                n.GetAtomicNum() == 7
                and not n.IsInRing()
                and not n.GetFormalCharge()
                and all(m.GetAtomicNum() == 6 for m in n.GetNeighbors() if m.GetIdx() != a.GetIdx())
            ):
                return False
            amino += 1
    return amino <= 1


def _only_plain_hydrazines(mol):
    """Every acyclic N-N bond is a single bond between chain nitrogens that bear only hydrogen, nitrogen and carbon
    without a double bond, so the hydrazine is a prefix beside a ring that contains nitrogen (P-44.1.2.2); one C=N on
    a hydrazine nitrogen is its ylidene substituent (P-68.3.1.2.2), two make an azine."""
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if (
            a.GetAtomicNum() != 7
            or b.GetAtomicNum() != 7
            or bond.IsInRing()
            or _joins_two_rings(bond)
            or a.IsInRing() != b.IsInRing()
            or _plain_diazene(bond)
        ):
            continue
        if bond.GetBondTypeAsDouble() != 1.0 or a.IsInRing() or b.IsInRing() or _is_azide_part(a):
            return False
        ylidene_nitrogens = 0
        for atom in (a, b):
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() not in (6, 7):
                    return False
                if neighbor.GetAtomicNum() != 6:
                    continue
                for x in neighbor.GetBonds():
                    other = x.GetOtherAtom(neighbor)
                    if x.GetBondTypeAsDouble() != 2.0 or other.GetAtomicNum() == 6 or x.IsInRing():
                        continue
                    if other.GetIdx() != atom.GetIdx():
                        return False
                    ylidene_nitrogens += 1
        if ylidene_nitrogens > 1:
            return False
    return True


_SENIOR_ATOM_ORDER = (*_PARENT_HYDRIDE_ORDER, 8, 16, 34, 52, 6)


def _ring_outranks_hydrides(mol, centers):
    """P-44.1.2: the parent holds the most senior atom of N > P > As > Sb > Bi > Si > Ge > Sn > Pb > B > Al > Ga > In > Tl > O > S >
    Se > Te > C, and a ring wins over a chain of the same element."""
    rank = {z: i for i, z in enumerate(_SENIOR_ATOM_ORDER)}
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    ring_rank = min((rank.get(mol.GetAtomWithIdx(a).GetAtomicNum(), len(rank)) for a in ring_atoms), default=len(rank))
    return ring_rank <= min(rank[c.GetAtomicNum()] for c in centers)


def _plain_parent(mol, graph, halogens, aromatic_atoms, stereo=None):
    """Parent without a principal group: the monocycle when there is one
    (P-44.1.2.2), else the longest chain, with every substituent a prefix."""
    from ._hydride_chain import name_hydride_chain

    chain_name = name_hydride_chain(mol, graph, halogens, aromatic_atoms)
    if chain_name is not None:
        return ((0,), chain_name, (None, None, None, 0, {}, False))
    ring_cation = RING_CENTER.get()
    if not ring_cation:
        hydroxylamine = _hydroxylamine_parent(mol, graph, halogens, aromatic_atoms)
        if hydroxylamine is not None:
            return ((0,), hydroxylamine, (None, None, None, 0, {}, False))
        azane = _azane_ylidyne_parent(mol, graph, halogens, aromatic_atoms)
        if azane is not None:
            return ((0,), azane, (None, None, None, 0, {}, False))
    centers = [] if ring_cation else [a for a in mol.GetAtoms() if a.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not a.IsInRing()]
    if centers and _ring_outranks_hydrides(mol, centers):
        centers = []
    if centers:
        named = [
            (
                -_hydroxy_suffix_count(n),
                -_amine_count(n),
                _PARENT_HYDRIDE_ORDER.index(c.GetAtomicNum()) if c.GetAtomicNum() in _PARENT_HYDRIDE_ORDER else 99,
                -len(graph[c.GetIdx()]),
                n,
                c.GetIdx(),
            )
            for c in centers
            if (n := _mononuclear_parent(mol, graph, halogens, aromatic_atoms, c))
        ]
        if named:
            *_, name, center = min(named)
            from ._substituents import centre_nuclides

            cited = bool(centre_nuclides(center))
            return ((0,), name, (None, None, None, 0, {} if cited else {center: 1}, False))
    if centers or (
        not ring_cation
        and not (any(a.GetAtomicNum() == 7 and a.IsInRing() for a in mol.GetAtoms()) and _only_plain_hydrazines(mol))
        and any(
            b.GetBeginAtom().GetAtomicNum() == 7
            and b.GetEndAtom().GetAtomicNum() == 7
            and not b.IsInRing()
            and not _joins_two_rings(b)
            and not _is_azide_part(b.GetBeginAtom())
            and not _is_azide_bond_end(b.GetBeginAtom(), b.GetEndAtom())
            for b in mol.GetBonds()
        )
    ):
        hydrazine = None if centers else _hydrazine_parent(mol, aromatic_atoms)
        if hydrazine is not None:
            return ((0,), hydrazine, (None, None, None, 0, {}, False))
        raise UnsupportedStructure("a heteroatom hydride is the senior parent when there is no principal group (P-44.1.2.2)")
    if not ring_cation and chalcogen_link_is_parent(mol):
        raise UnsupportedStructure("a heteroatom hydride is the senior parent when there is no principal group (P-44.1.2.2)")
    ring_info = mol.GetRingInfo()
    rings = [r for r in ring_info.AtomRings()]
    if len(rings) >= 2:
        assembly = _assembly_parent(mol, graph, halogens, aromatic_atoms, None, [], stereo) or _assembly_beside_systems(
            mol, graph, halogens, aromatic_atoms, stereo
        )
        if assembly is not None:
            return assembly[1]
    if rings and any(ring_info.NumAtomRings(a) != 1 for r in rings for a in r):
        return _fused_plain_parent(mol, graph, rings)
    if rings:
        if any(ring_info.NumAtomRings(a) != 1 for r in rings for a in r):
            raise UnsupportedStructure("several rings without a principal group are not named by the chain engine")
        cationic = CATION_RING_ATOM.get()
        if cationic is not None and any(cationic in r for r in rings):
            rings = [r for r in rings if cationic in r]  # the cationic ring outranks every neutral ring (P-41)
        ranked = sorted(rings, key=lambda r: _ring_rank(mol, r[0]))
        top = _ring_rank(mol, ranked[0][0])
        tied = [r for r in ranked if _ring_rank(mol, r[0]) == top]
        for i, first in enumerate(tied):
            for second in tied[:i]:
                bonded = any(mol.GetBondBetweenAtoms(a, b) is not None for a in first for b in second)
                if bonded and _bare_key(mol, set(first)) == _bare_key(mol, set(second)):
                    raise UnsupportedStructure("identical rings joined directly form a ring assembly (P-28)")
        candidates = _isolated_candidates(
            [lambda ring=ring: _plain_ring_parent(mol, graph, halogens, aromatic_atoms, ring, stereo) for ring in tied]
        )
        if not candidates:
            raise UnsupportedStructure("none of the equally senior rings has a supported name")
        return _choose_isolated(candidates, key=lambda c: c[0])[1]

    eligible = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not is_functional_carbon(mol, a.GetIdx())
    }
    from ._substituents import ISOTOPE_LABELS

    best = None
    stereo_context = BRANCH_STEREO.get()
    used_before = set(stereo_context["used"]) if stereo_context else set()
    best_used = set()
    isotope_context = ISOTOPE_LABELS.get()
    consumed_before = set(isotope_context["consumed"]) if isotope_context else set()
    best_consumed = set()
    paths = _paths(graph, eligible)
    longest = max((len(p) for p in paths), default=0)
    for path in (p for p in paths if len(p) == longest):
        for chain in (path, path[::-1]):
            if stereo_context:
                stereo_context["used"] &= used_before
            if isotope_context:
                isotope_context["consumed"] &= consumed_before
            candidate = _evaluate_plain(mol, graph, halogens, aromatic_atoms, chain, stereo)
            if best is None or candidate[0] < best[0]:
                best = candidate
                best_used = set(stereo_context["used"]) - used_before if stereo_context else set()
                best_consumed = set(isotope_context["consumed"]) - consumed_before if isotope_context else set()
    if stereo_context:
        stereo_context["used"] &= used_before
        stereo_context["used"] |= best_used
    if isotope_context:
        isotope_context["consumed"] &= consumed_before
        isotope_context["consumed"] |= best_consumed
    if best is None:
        raise UnsupportedStructure("no chain carbons")
    return best


def _without_stereo(name):
    return _ALL_DESCRIPTORS.sub("", name, count=1) if _ALL_DESCRIPTORS.match(name) else name


def _fused_plain_parent(mol, graph, rings):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    systems = []
    for ring in rings:
        system_rings, system_atoms = _system_of(mol, ring[0])
        if not any(set(system_atoms) == set(seen[1]) for seen in systems):
            systems.append((system_rings, system_atoms))
    cationic = CATION_RING_ATOM.get()
    if cationic is not None and any(cationic in atoms for _, atoms in systems):
        systems = [system for system in systems if cationic in system[1]]
    from ._substituents import ISOTOPE_LABELS

    labels = (ISOTOPE_LABELS.get() or {}).get("labels", {})

    def rank(system):
        count, nuclides = _isotope_counts(None, [(a, None) for a in system[1] if a in labels])
        return ring_seniority_key(mol, system[1]), count, nuclides

    ranked = sorted(systems, key=rank)
    tied = [system for system in ranked if rank(system) == rank(ranked[0])]
    found = None
    if len(tied) > 1:
        found = _configuration_senior_system(mol, graph, tied)
    else:
        system_rings, system_atoms = ranked[0]
        if len(system_rings) > 1:
            _require_mancude_system(mol, system_atoms)
        found = evaluate_skeleton(mol, graph, "ring", system_rings, system_atoms, [], set(), "")
    if found is None:
        raise UnsupportedStructure("this fused ring system has no supported numbering")
    from ._diester_ring_diyl import PARENT_START

    return ((0,), _without_stereo(found[1]), (None, None, None, 0, found[2], True, PARENT_START.get()))


_DESCRIPTOR = (
    r"(?:\d+[a-z]*['\u2032]*(?:\(\d+[a-z]*['\u2032]*\))?)?"
    r"(?:[RSEZrsMPmp]|seqCis|seqTrans|\((?:TBPY-5|OC-6|SPY-5|SS-4|TS-3)-\d+['\u2032]*(?:-[AC])?\))"
)
_ALL_DESCRIPTORS = re.compile(rf"[(\[]{_DESCRIPTOR}(?:,{_DESCRIPTOR})*[)\]]-")


def _configuration_senior_system(mol, graph, tied):
    """The numbered parent among equally senior ring systems that differ only in isotopic modification or configuration:
    the lowest locants for the modified atoms (P-44.4.1.11.4-6), then Z before E and R before S (P-44.4.1.12.1,
    P-44.4.1.12.2, P-45.6.2); identical ones need a multiplicative or assembly name."""
    from ._diester_ring_diyl import PARENT_START, evaluate_skeleton

    candidates, starts = [], []
    for rings, atoms in tied:
        if len(rings) > 1:
            _require_mancude_system(mol, atoms)
        found = evaluate_skeleton(mol, graph, "ring", rings, atoms, [], set(), "")
        if found is None:
            raise UnsupportedStructure("this fused ring system has no supported numbering")
        candidates.append(found)
        starts.append(PARENT_START.get())
    from ._substituents import ISOTOPE_LABELS

    entries = _isotope_entries((ISOTOPE_LABELS.get() or {}).get("labels", {}))
    ranks = [_isotope_locant_rank(entries, c[2]) for c in candidates]
    keep = [i for i, rank in enumerate(ranks) if rank == min(ranks)]
    candidates, starts = [candidates[i] for i in keep], [starts[i] for i in keep]
    if len(candidates) == 1:
        return candidates[0]
    if len({_ALL_DESCRIPTORS.sub("", c[1]) for c in candidates}) > 1 or len({c[0] for c in candidates}) == 1:
        return _by_substituents(mol, graph, [t for i, t in enumerate(tied) if i in keep], candidates, starts)
    return min(candidates, key=lambda c: c[0])


def _by_substituents(mol, graph, tied, candidates, starts):
    """Identical ring systems that the units of a multiplicative name do not match (different junction atoms or different
    substituents) are told apart as parent by the number of substituents, their locants and the alphanumerical order of
    the names (P-45.2.1, P-45.2.2, P-45.2.3, P-45.5)."""
    from ._diester_ring_diyl import PARENT_START
    from ._multiplicative import name_if_multiplicative

    if name_if_multiplicative(mol) is not None:
        raise UnsupportedStructure("several equally senior ring systems need a multiplicative or assembly name")

    def key(index):
        found = candidates[index]
        roots = sum(1 for a in tied[index][1] for n in graph[a] if n not in tied[index][1])
        return (-roots, found[0], alphanumerical_name_key(found[1]), found[1])

    best = min(range(len(candidates)), key=key)
    PARENT_START.set(starts[best])
    return candidates[best]


def _evaluate_plain(mol, graph, halogens, aromatic_atoms, chain, stereo=None):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    entries = {}
    for atom in chain:
        for neighbor in graph[atom]:
            if neighbor in chain_set:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    from ._substituents import ISOTOPE_LABELS

    if not entries and not ISOTOPE_LABELS.get():
        raise UnsupportedStructure("an unsubstituted chain is not a polyfunctional case")
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    length = len(chain)
    single = total_count == 1
    prefix = format_substituent_prefixes(
        grouped, omit_locants=(length == 1 or (length == 2 and single)) and not FORCE_LOCANTS.get()
    )
    if length == 2 and (ene or yne):
        body = "ethene" if ene else "ethyne"
    else:
        body = name_from_substituents(length, ene, yne, "e")
    name = prefix + body
    key = (
        -length,
        -(len(ene) + len(yne)),
        -len(ene),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        _parent_configuration_rank(stereo, position_of, directional=False),
        -total_count,
        locant_set,
        citation,
        alphanumerical_name_key(name),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (None, None, None, 0, position_of, False, len(prefix))


_CARBO_WORDS = {
    "acid": "carboxylic acid",
    "amide": "carboxamide",
    **{name: f"carbo{name}" for name in _CHALCOGEN_AMIDE_CLASSES},
    "amidine": "carboximidamide",
    "hydrazonamide": "carbohydrazonamide",
    "imidohydrazide": "carboximidohydrazide",
    "hydrazonohydrazide": "carbohydrazonohydrazide",
    "hydrazide": "carbohydrazide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
}


def _carbo_best(mol, graph, halogens, aromatic_atoms, principal, principal_atoms, stereo, n_names=()):
    """The chain whose attached group carbons are cited as 'carbo' suffixes
    (propane-1,2,3-tricarboxylic acid), or None."""
    owned = set().union(*principal_atoms.values())
    group_carbons = set(principal_atoms)
    eligible = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6
        and not a.IsInRing()
        and a.GetIdx() not in group_carbons
        and not is_functional_carbon(mol, a.GetIdx())
    }
    best = None
    for path in _paths(graph, eligible):
        for chain in (path, path[::-1]):
            candidate = _evaluate_carbo(
                mol, graph, halogens, aromatic_atoms, chain, principal, group_carbons, owned, stereo, n_names
            )
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    return best


def _evaluate_carbo(mol, graph, halogens, aromatic_atoms, chain, principal, group_carbons, owned, stereo, n_names=()):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    attached = {}
    for g in group_carbons:
        host = [n for n in graph[g] if n in chain_set]
        if host:
            attached[g] = host[0]
    if not attached:
        return None
    suffix_locants = sorted(position_of[a] for a in attached.values())
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    entries = {}
    for atom in chain:
        for neighbor in graph[atom]:
            if neighbor in chain_set or neighbor in attached:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    count = len(attached)
    length = len(chain)
    word = carbo_suffix(spec_from_key(principal), 1) if _is_variant(principal) else _CARBO_WORDS[principal]
    prefix = format_substituent_prefixes(_with_n_names(grouped, n_names, position_of, count))
    body = name_from_substituents(length, ene, yne, multiplied_word(count, word), suffix_locants)
    name = prefix + body
    key = (
        -count,
        -length,
        tuple(suffix_locants),
        -(len(ene) + len(yne)),
        lowest_locant_set(ene + yne),
        _parent_configuration_rank(stereo, position_of, directional=False),
        -total_count,
        locant_set,
        citation,
        alphanumerical_name_key(name),
        _n_group_positions(n_names, position_of),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (prefix, body, "", 0, position_of, False)


def _locant_order(locant):
    """Primed locants follow the unprimed one of the same number: 1 < 1' < 2."""
    return locant[1], locant[0]


def _replacement_locants(mol, locants):
    """(sort key, text) of the skeletal replacement prefixes of a cycloalkane assembly numbering (P-28.4.2)."""
    from ._ring_assembly import _replacement_text

    return _replacement_text(mol, {atom: f"{number}{chr(39) * prime}" for atom, (prime, number) in locants.items()})


def _assembly_replacement(mol, specs, locants, bracketed=True):
    """The skeletal replacement prefixes written ahead of the assembly name: '3,3'-dioxa-1,1'-bi(...)', or directly
    before the bracket of '3,3'-dioxa[1,1'-bi(...)]-2-ol' (P-28.4.1, P-28.4.2)."""
    if not any(sp.replaced for sp in specs):
        return ""
    text = _replacement_locants(mol, locants)[1]
    return text.rstrip("-") if bracketed else text


def _replacement_junction_ok(mol, specs, join):
    """A replacement prefix names a skeletal atom of the hydrocarbon assembly, which the junction atoms must be part of
    as carbon (P-28.4.1)."""
    return not any(sp.replaced for sp in specs) or all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in join)


def _assembly_stereo_prefix(locants, stereo):
    """'(1R,1'S,2S)-' for the CIP descriptors `stereo` ({atom: code}) of ring atoms of an assembly numbered by `locants`."""
    from ._substituents import _branch_stereo_prefix

    return _branch_stereo_prefix(sorted((PrimedLocant(*locants[a]), code) for a, code in stereo.items() if a in locants))


def _assembly_hydro_positions(mol, hydro_rings, locants):
    """The (prime count, position) of the atoms of the hydro rings of a two-ring assembly that take hydro prefixes: every
    position of a six-membered ring, every one but the heteroatom of a five-membered ring (P-54.3, P-31.1.4.2.4)."""
    from ._ring_assembly_chain import hydrogenated_atoms

    found = []
    for ring in hydro_rings:
        for atom in hydrogenated_atoms(mol, ring):
            if len(ring) == 5 and locants[atom][1] == 1:
                continue
            found.append(locants[atom])
    return sorted(found, key=_locant_order)


def _assembly_numbering(graph, rings, join, marked, entries, specs=None, cite_marked=None, mol=None, stereo=None, hydro_rings=()):
    """Atom -> (prime count, locant) for the best numbering of two directly
    joined identical rings. Benzene and cycloalkane rings number from the
    junction (1 and 1'); heteroaromatic rings keep their fixed numbering and
    the junction takes the lowest locants, then the replacement heteroatoms of
    a large ring, the marked atoms (principal groups or a free valence), then
    the substituents."""
    from ._diester_ring_diyl import _isotope_key

    hetero = bool(specs) and specs[0].hetero is not None
    replaced = bool(specs) and mol is not None and any(sp.replaced for sp in specs)
    cycles = [ring_cycle(graph, rings[0]), ring_cycle(graph, rings[1])]

    def orientations(index):
        if hetero:
            return [numbering for numbering in numberings(specs[index])]
        cycle = cycles[index]
        start = cycle.index(join[index])
        rotated = cycle[start:] + cycle[:start]
        return [
            {atom: i + 1 for i, atom in enumerate(order)} for order in (rotated, [rotated[0]] + rotated[:0:-1])
        ]

    best = None
    unsaturated = bool(specs) and any(sp.kind == "cycloalkene" for sp in specs) and not hydro_rings
    for unprimed in (0, 1):
        for first in orientations(unprimed):
            for second in orientations(1 - unprimed):
                locants = {atom: (0, number) for atom, number in first.items()}
                locants.update({atom: (1, number) for atom, number in second.items()})
                ene, closing = (), 0
                if unsaturated:
                    found = sum(_assembly_multiple_locants(specs, locants), [])
                    size = len(specs[0].cycle)
                    closing = sum(1 for loc in found if loc[1] == size)
                    ene = tuple(sorted(_locant_order((loc[0], 1 if loc[1] == size else loc[1])) for loc in found))
                key = (
                    (locants[join[unprimed]][1], locants[join[1 - unprimed]][1]),
                    _replacement_locants(mol, locants)[0] if replaced else (),
                    tuple(sorted(_locant_order(locants[a]) for a in marked)),
                    tuple(_locant_order(locants[a]) for a in cite_marked) if cite_marked else (),
                    closing,
                    ene,
                    tuple(_locant_order(loc) for loc in _assembly_hydro_positions(mol, hydro_rings, locants)) if hydro_rings else (),
                    tuple(sorted(_locant_order(locants[r]) for r, _, _ in entries)),
                    _citation_key([(_locant_order(locants[r]), name) for r, name, _ in entries]),
                    _isotope_key({a: _locant_order(loc) for a, loc in locants.items()}, locants),
                    _branch_stereo_rank(sorted((_locant_order(locants[a]), code) for a, code in (stereo or {}).items() if a in locants)),
                )
                if best is None or key < best[0]:
                    best = (key, locants, ene)
    return best[1]


def _assembly_multiple_locants(specs, locants):
    """([(prime count, locant)] of the ring double bonds, [...] of the triple bonds) of an assembly numbering."""
    ene, yne = [], []
    for spec in specs:
        if spec.kind != "cycloalkene":
            continue
        prime = locants[spec.cycle[0]][0]
        found_ene, found_yne = multiple_locants(spec, {a: locants[a][1] for a in spec.cycle})
        ene += [(prime, n) for n in found_ene]
        yne += [(prime, n) for n in found_yne]
    return sorted(ene, key=_locant_order), sorted(yne, key=_locant_order)


def _assembly_unsaturation(specs, locants):
    """The ending of an assembly of saturated components with ring double or triple bonds, cited after the bracket
    (P-31.1.7.1): '1,2'-diene'; '' when no ring is unsaturated."""
    ene, yne = _assembly_multiple_locants(specs, locants)
    if not ene and not yne:
        return ""

    size = len(specs[0].cycle)

    def spot(prime, n):
        mark = chr(39) * prime
        return f"1{mark}({n}{mark})" if n == size else f"{n}{mark}"

    def cite(found):
        # P-31.1.4.2.4: a bond closing the ring is cited from its lower end with the other end in parentheses and
        # takes its place in the order by that lower end
        ordered = sorted(found, key=lambda item: (1 if item[1] == size else item[1], item[0]))
        return ",".join(spot(prime, n) for prime, n in ordered)

    ene_word, yne_word = multiplied_word(len(ene), "ene"), multiplied_word(len(yne), "yne")
    if ene and yne:
        return f"{cite(ene)}-{ene_word[:-1]}-{cite(yne)}-{yne_word}"
    return f"{cite(ene)}-{ene_word}" if ene else f"{cite(yne)}-{yne_word}"


def _compatible_assembly_rings(mol, rings, specs):
    """Whether two rings are the same component of a ring assembly: unsaturation of cycloalkane components is
    cited as endings, so a cycloalkene counts as its cycloalkane (P-31.1.7.1)."""
    if any(sp is None or sp.kind == "pyrrole" for sp in specs) or len(rings[0]) != len(rings[1]):
        return False
    if {sp.kind for sp in specs} <= {"cycloalkane", "cycloalkene"}:
        return True
    if specs[0].name or specs[1].name:
        return specs[0].name == specs[1].name
    return specs[0].kind == specs[1].kind and _bare_key(mol, set(rings[0])) == _bare_key(mol, set(rings[1]))


def _junction_is_ylidene(mol, join, specs):
    """False for a single-bond junction, True for a double bond between two saturated rings, else None."""
    order = mol.GetBondBetweenAtoms(*join).GetBondTypeAsDouble()
    if order == 1.0:
        return False
    if order == 2.0 and all(sp.kind == "cycloalkane" for sp in specs):
        return True
    return None


def _assembly_base(specs, locants, join, elide, ylidene=False):
    """'1,1'-biphenyl', '1,1'-bi(cyclohexane)', '2,2'-bipyridine', '1,1'-bi(cyclohexylidene)'; the final
    'e' goes before a vowel-initial suffix."""
    spec = specs[0]
    if spec.kind == "benzene":
        return "1,1'-biphenyl"
    unprimed, primed = sorted(join, key=lambda atom: locants[atom][0])
    spots = f"{locants[unprimed][1]},{locants[primed][1]}'"
    if ylidene:
        return f"{spots}-bi({spec.parent[:-3]}ylidene)"
    parent = spec.parent
    stem = parent[:-1] if elide and not KEEP_ASSEMBLY_ENDING.get() and parent.endswith("e") else parent
    marked = spec.hetero is None or re.search(r"[\d\[-]", stem)
    return f"{spots}-bi({stem})" if marked else f"{spots}-bi{stem}"


def _assembly_spec(mol, graph, ring, replaceable=False):
    """The component of a ring assembly: a benzene, cycloalkane or heteroaromatic ring, or a saturated heteromonocycle
    named as its parent hydride (P-28.4.1) with the numberings of that name; with `replaceable`, a ring of more than ten
    members whose heteroatoms are cited as skeletal replacement prefixes of the cycloalkane assembly (P-28.4.1)."""
    spec = spec_of(mol, ring)
    if spec is not None:
        return spec
    from ._ring_assembly import _replacement_ring_kind, _saturated_heteroring_kind
    from ._ring_diyl_numbering import _hetero_monocycle

    found = _saturated_heteroring_kind(mol, graph, ring)
    if found is None:
        if replaceable and _replacement_ring_kind(mol, ring) is not None:
            return RingSpec(ring_cycle(graph, list(ring)), "cycloalkane", None, replaced=True)
        return None
    cycle = ring_cycle(graph, list(ring))
    options = tuple(position_of for position_of, *_ in _hetero_monocycle(mol, cycle, None)[1])
    hetero = next(a for a in cycle if mol.GetAtomWithIdx(a).GetAtomicNum() != 6)
    return RingSpec(cycle, found[0], hetero, name=found[1], fixed=options)


def assembly_substituent(mol, graph, root, coming_from, halogens, aromatic_atoms, place=None, place_monocycles=False):
    """(name, True) of a substituent group made of two identical directly
    joined rings ([1,1'-biphenyl]-4-yl), entered at `root`; None otherwise. For a group of fused, chained or von
    Baeyer components `place`, a dict, receives the locant of every ring atom and the offset of the assembly name in the
    group name, so that an isotopic descriptor can be cited in front of it. A double or triple bond to
    `coming_from` makes the group an ylidene or ylidyne (P-29.2)."""
    from ._free_valence import suffix_of

    named = _assembly_substituent_yl(mol, graph, root, coming_from, halogens, aromatic_atoms, place, place_monocycles)
    if named is None:
        return None
    suffix = suffix_of(mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble())
    if suffix == "yl" or named[0].endswith(("ylidene", "ylidyne")):
        return named
    return (named[0][: -len("yl")] + suffix, named[1]) if suffix and named[0].endswith("yl") else None


def _assembly_substituent_yl(mol, graph, root, coming_from, halogens, aromatic_atoms, place, place_monocycles):
    from ._system_assembly import system_assembly

    fused = system_assembly(
        mol, graph, halogens, aromatic_atoms, None, [], None, free=(root, coming_from), within=_arm_atoms(graph, root, coming_from),
        place=place,
    )
    if fused is not None:
        return fused
    from ._chain_assembly import chain_assembly

    chained = chain_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=(root, coming_from), place=place)
    if chained is not None:
        return chained
    from ._vb_ring_assembly import vb_assembly

    von_baeyer = vb_assembly(
        mol, graph, halogens, aromatic_atoms, None, [], None, free=(root, coming_from), within=_arm_atoms(graph, root, coming_from),
        place=place,
    )
    if von_baeyer is not None:
        return von_baeyer
    ring_info = mol.GetRingInfo()
    own = next((list(r) for r in ring_info.AtomRings() if root in r), None)
    if own is None:
        return None
    branch = _arm_atoms(graph, root, coming_from)
    inside = [list(r) for r in ring_info.AtomRings() if set(r) <= branch]
    if len(inside) != 2 or any(ring_info.NumAtomRings(a) != 1 for r in inside for a in r):
        return None
    other = next(r for r in inside if r != own)
    joins = [(a, b) for a in own for b in other if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1:
        return None
    specs = [_assembly_spec(mol, graph, r, True) for r in (own, other)]
    if not _compatible_assembly_rings(mol, (own, other), specs) or not _replacement_junction_ok(mol, specs, joins[0]):
        return None
    ylidene = _junction_is_ylidene(mol, joins[0], specs)
    if ylidene is None:
        return None
    ring_atoms = set(own) | set(other)
    roots = [
        (r, n.GetIdx())
        for r in ring_atoms
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring_atoms and not (r == root and n.GetIdx() == coming_from)
    ]
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    rings = [own, other]
    context = BRANCH_STEREO.get()
    ring_stereo = {a: code for a, code in (context["atoms"] if context else {}).items() if a in ring_atoms}
    locants = _assembly_numbering(graph, rings, joins[0], [root], entries, specs, mol=mol, stereo=ring_stereo)
    grouped = {}
    for r, name, compound in entries:
        spot = locants[r]
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(f"{spot[1]}{chr(39) * spot[0]}")
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(text.rstrip(chr(39))), text.count(chr(39))))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    ending = _assembly_unsaturation(specs, locants)
    base = _assembly_base(specs, locants, joins[0], elide=not ending, ylidene=ylidene)
    spot = locants[root]
    enes = f"{ending[:-1]}-" if ending else ""
    core = _assembly_replacement(mol, specs, locants) + f"[{base}]-{enes}{spot[1]}{chr(39) * spot[0]}-yl"
    core = _assembly_group_label(mol, graph, ring_atoms, root, locants, core)
    from ._substituents import _branch_stereo_entries

    _branch_stereo_entries({a: PrimedLocant(*loc) for a, loc in locants.items()}, ring=True, record=True)
    name = _assembly_stereo_prefix(locants, ring_stereo) + assembly_join(prefix, core)
    if place is not None and place_monocycles:
        place.update(locants={a: PrimedLocant(*loc) for a, loc in locants.items()}, start=len(name) - len(core), atoms=ring_atoms)
    return name, True


def _assembly_group_label(mol, graph, ring_atoms, root, locants, core):
    """`core` with the isotopic descriptor of the modified ring atoms in front of the assembly name (P-82.2.1, P-82.5.2)."""
    from ._isotope_labels import descriptor
    from ._multiplicative_text import PrimedLocant
    from ._substituents import ISOTOPE_LABELS, _fully_modified, _hydrogen_capacity

    context = ISOTOPE_LABELS.get()
    labelled = {a: e for a, e in (context["labels"] if context else {}).items() if a in ring_atoms and a not in context["consumed"]}
    if not labelled:
        return core
    positions = {a: PrimedLocant(*locants[a]) for a in ring_atoms}
    capacity = _hydrogen_capacity(context, labelled)
    expected = {a for a in ring_atoms if a != root and (mol.GetAtomWithIdx(a).GetTotalNumHs() or any(n not in ring_atoms for n in graph[a]))}
    omit = _fully_modified(labelled, positions, capacity, expected)
    context["consumed"].update(labelled)
    return descriptor(labelled, positions, omit, capacity=capacity) + core


def _assembly_parent(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, n_names=()):
    """P-28: two identical benzene or cycloalkane rings joined by a single
    bond -- [1,1'-biphenyl]-4-ol, 4-nitro-1,1'-biphenyl. Returns
    (principal count, (key, name, parts)) or None when `mol` is not one."""
    if principal is not None and not occurrences:
        return None
    from ._system_assembly import system_assembly

    fused = system_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, n_names=n_names)
    if fused is not None:
        return fused
    from ._chain_assembly import chain_assembly

    chained = chain_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, n_names=n_names)
    if chained is not None:
        return chained
    from ._vb_ring_assembly import vb_assembly

    von_baeyer = vb_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo)
    if von_baeyer is not None:
        return _without_n_prefixes(von_baeyer, n_names)
    if occurrences and len(_joined_system_sets(mol, 2)) > 1:
        fused = _joined_pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, n_names)
        if fused is not None:
            return fused
    ring_info = mol.GetRingInfo()
    all_rings = [list(r) for r in ring_info.AtomRings()]
    if occurrences:
        # P-44.1.1: the assembly holding the principal group is the parent, so every other ring system is a substituent
        all_rings = [r for r in all_rings if all(ring_info.NumAtomRings(a) == 1 for a in r)]
    if len(all_rings) < 2 or any(ring_info.NumAtomRings(a) != 1 for r in all_rings for a in r):
        return None
    ring_specs = [_assembly_spec(mol, graph, r, True) for r in all_rings]
    if len(all_rings) > 2 and any(sp is not None and sp.kind == "cycloalkene" for sp in ring_specs):
        return None
    if not occurrences and not all(
        _compatible_assembly_rings(mol, (all_rings[0], r), (ring_specs[0], sp)) or _hydro_pair_rings(mol, graph, [all_rings[0], r])
        for r, sp in zip(all_rings, ring_specs)
    ):
        return None
    best = None
    for i, first in enumerate(all_rings):
        for second in all_rings[i + 1 :]:
            found = _pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, [first, second], n_names)
            if found is None:
                continue
            ylidene = _junction_is_ylidene(mol, found[2], [_assembly_spec(mol, graph, first, True), _assembly_spec(mol, graph, second, True)])
            # P-28.2.2: a double-bond junction is a two-ring assembly only; the pair holding the double bond has the
            # parent's multiple bond, so it ranks first
            key = (-found[0], not ylidene, ring_seniority_key(mol, set(first) | set(second), assembly=True), found[1][1])
            if best is None or key < best[0]:
                best = (key, (found[0], found[1]))
    return best[1] if best else None


def _hydro_pair_rings(mol, graph, rings):
    """The hydro ring of two directly joined rings that are the saturated and the mancude form of one ring, pyridine and
    piperidine (P-54.3); empty for benzene beside cyclohexane, which keeps its substitutive name (cyclohexylbenzene)."""
    from ._ring_assembly_chain import _has_ring_multiple_bond, ring_kinds_with_hydro

    kinds, hydro = ring_kinds_with_hydro(mol, graph, [tuple(r) for r in rings])
    if None in kinds or len(set(kinds)) != 1 or len(hydro) != 1:
        return []
    if kinds[0] == ("aromatic", 6) and not _has_ring_multiple_bond(mol, next(iter(hydro))):
        return []
    return hydro


def _joined_pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, n_names):
    """The assembly of two identical ring systems that hold every principal group, named by the system assembly with the
    other ring systems of the molecule cited as substituents (P-44.1.1, P-28.2.3); None when no pair is one."""
    from ._system_assembly import system_assembly

    wanted = {o[1] for o in occurrences}
    best = None
    for atoms, members in _joined_system_sets(mol, 2):
        if len(members) != 2 or not wanted <= atoms:
            continue
        try:
            found = system_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo, within=atoms, n_names=n_names)
        except UnsupportedStructure:
            continue
        if found is None:
            continue
        key = (-found[0], ring_seniority_key(mol, atoms, assembly=True), found[1][1])
        if best is None or key < best[0]:
            best = (key, found)
    return best[1] if best else None


def _without_n_prefixes(found, n_names):
    if n_names:
        raise UnsupportedStructure("N-prefixes on a bridged ring assembly are not supported yet")
    return found


def assembly_n_names(grouped, n_names, position_of, group_count):
    """`grouped` plus the N-prefixes of the amine groups of a ring assembly, ordered by the position of their group
    (N4,N4'-dimethyl, P-62.2.4.1.2); `position_of` maps the atoms of the assembly to their locants."""
    if not n_names:
        return grouped

    def position(entry):
        found = _n_position(entry[2], position_of) if len(entry) > 2 and isinstance(entry[2], tuple) else None
        return found if found is not None else min(position_of.values())

    return _with_n_names(grouped, sorted(n_names, key=position), position_of, group_count)


def _assembly_beside_systems(mol, graph, halogens, aromatic_atoms, stereo):
    """The ring assembly of two identical rings as the parent of a molecule that also holds other ring systems, when
    no other ring system outranks it (P-44.2.2.2.7); None otherwise."""
    from ._multiplicative import _ring_systems

    ring_info = mol.GetRingInfo()
    monocycles = [list(r) for r in ring_info.AtomRings() if all(ring_info.NumAtomRings(a) == 1 for a in r)]
    best, named = None, set()
    for i, first in enumerate(monocycles):
        for second in monocycles[i + 1 :]:
            found = _pair_assembly(mol, graph, halogens, aromatic_atoms, None, [], [first, second])
            if found is None:
                continue
            ylidene = _junction_is_ylidene(mol, found[2], [_assembly_spec(mol, graph, first, True), _assembly_spec(mol, graph, second, True)])
            atoms = set(first) | set(second)
            named.add(frozenset(atoms))
            key = (not ylidene, ring_seniority_key(mol, atoms, assembly=True), found[1][1])
            if best is None or key < best[0]:
                best = (key, atoms, (found[0], found[1]))
    for atoms, members in _joined_system_sets(mol):
        found = _system_set_assembly(mol, graph, halogens, aromatic_atoms, stereo, atoms, members)
        if found is None:
            continue
        joined = frozenset(found[1][2][4])
        if joined in named:
            continue
        named.add(joined)
        ylidene = any(
            b.GetBondTypeAsDouble() == 2.0
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in joined and b.GetEndAtomIdx() in joined and not b.IsInRing()
        )
        key = (not ylidene, ring_seniority_key(mol, joined, assembly=True), found[1][1])
        if best is None or key < best[0]:
            best = (key, joined, found)
    if best is None:
        return None
    if any(ring_seniority_key(mol, system) < best[0][1] for system in _ring_systems(mol) if not system <= best[1]):
        return None
    return best[2]


def _joined_system_sets(mol, largest=6):
    """(atoms, member systems) of every set of two to `largest` ring systems linked through direct bonds between them."""
    from ._multiplicative import _ring_systems

    systems = [frozenset(s) for s in _ring_systems(mol)]
    owner = {a: i for i, system in enumerate(systems) for a in system}
    links = {i: set() for i in range(len(systems))}
    for bond in mol.GetBonds():
        i, j = owner.get(bond.GetBeginAtomIdx()), owner.get(bond.GetEndAtomIdx())
        if i is not None and j is not None and i != j:
            links[i].add(j)
            links[j].add(i)
    found, grown = set(), {frozenset([i]) for i in links}
    for _ in range(largest - 1):
        grown = {chosen | {n} for chosen in grown for i in chosen for n in links[i] if n not in chosen}
        found |= grown
    ordered = sorted(found, key=sorted)
    return [(frozenset().union(*(systems[i] for i in chosen)), [systems[i] for i in chosen]) for chosen in ordered]


def _system_set_assembly(mol, graph, halogens, aromatic_atoms, stereo, atoms, members):
    """The assembly name of ring systems that are the same component joined by direct bonds, as
    (count, (key, name, parts)); None when they are not one or cannot be named."""
    from ._chain_assembly import chain_assembly
    from ._system_assembly import _skeleton_key, system_assembly

    if len({_skeleton_key(mol, m) for m in members}) != 1:
        return None
    from ._vb_ring_assembly import vb_assembly

    build = system_assembly if len(members) == 2 else chain_assembly
    for builder in (build, vb_assembly):
        try:
            found = builder(mol, graph, halogens, aromatic_atoms, None, [], stereo, within=atoms)
        except UnsupportedStructure:
            continue
        if found is not None:
            return found
    return None


def _pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, rings, n_names=()):
    """The assembly name of two directly joined identical rings with every other part of `mol` cited as a
    substituent; None when the pair is not such an assembly or leaves a principal group outside it."""
    if set(rings[0]) & set(rings[1]):
        return None
    if occurrences and any(o[1] not in set(rings[0]) | set(rings[1]) for o in occurrences):
        return None
    specs = [_assembly_spec(mol, graph, r, True) for r in rings]
    hydro_rings = _hydro_pair_rings(mol, graph, rings) if all(sp is not None for sp in specs) else []
    name_specs = specs
    if hydro_rings:
        mancude = next(sp for sp, r in zip(specs, rings) if frozenset(r) not in hydro_rings)
        name_specs = [mancude, mancude]
    elif not _compatible_assembly_rings(mol, rings, specs):
        return None
    joins = [(a, b) for a in rings[0] for b in rings[1] if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1 or not _replacement_junction_ok(mol, name_specs, joins[0]):
        return None
    ylidene = _junction_is_ylidene(mol, joins[0], name_specs)
    if ylidene is None:
        return None
    owned = set().union(*(o[2] for o in occurrences)) if occurrences else set()
    ring_atoms = set(rings[0]) | set(rings[1])
    roots = [
        (r, n.GetIdx())
        for r in ring_atoms
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring_atoms and n.GetIdx() not in owned
    ]
    entries = [(r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots]
    locants = _assembly_numbering(
        graph, rings, joins[0], [o[1] for o in occurrences], entries, specs, mol=mol, hydro_rings=hydro_rings
    )

    def cite(locant):
        return f"{locant[1]}{chr(39) * locant[0]}"

    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(text.rstrip(chr(39))), text.count(chr(39))))
    count = len(occurrences)
    grouped = assembly_n_names(grouped, n_names, {a: PrimedLocant(*loc) for a, loc in locants.items()}, count)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    ending = _assembly_unsaturation(name_specs, locants)
    replacement = _assembly_replacement(mol, name_specs, locants, bool(principal) or bool(ending))
    hydro = ""
    if hydro_rings:
        spots = [cite(loc) for loc in _assembly_hydro_positions(mol, hydro_rings, locants)]
        hydro = f"{','.join(spots)}-{multiplied_word(len(spots), 'hydro')}"
    if principal is None:
        core = _assembly_base(name_specs, locants, joins[0], elide=False, ylidene=ylidene)
        core = f"{replacement}[{core}]-{ending}" if ending else replacement + core
        core = f"{hydro}-{core}" if hydro else core
    else:
        word = multiplied_word(count, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        base = _assembly_base(name_specs, locants, joins[0], elide=word[0] in "aeiouy" and not ending, ylidene=ylidene)
        spots = ",".join(cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: _locant_order(locants[o[1]])))
        enes = f"{ending[:-1] if word[0] in 'aeiouy' else ending}-" if ending else ""
        core = f"{hydro}{replacement}[{base}]-{enes}{spots}-{word}"
    name = assembly_join(prefix, core)
    return count, ((-count,), name, (None, None, None, 0, {a: PrimedLocant(*loc) for a, loc in locants.items()}, True, len(name) - len(core))), joins[0]


_CYCLOPENTA_A_PHENANTHRENE = Chem.MolFromSmarts(
    "[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2"
)


def _require_mancude_system(mol, atoms):
    """Only fully aromatic fused systems (arenes, mancude heterocycles): partly
    hydrogenated, bridged and spiro systems need hydro/von Baeyer names. Beyond three rings only
    all-six-membered systems and peri-fused ones with a retained numbering are verified."""
    if not any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in atoms) and any(
        set(match) == set(atoms) for match in mol.GetSubstructMatches(_CYCLOPENTA_A_PHENANTHRENE)
    ):
        raise UnsupportedStructure("a saturated cyclopenta[a]phenanthrene skeleton is a steroid parent hydride (P-101), not a hydro fusion name")
    member_rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= set(atoms)]
    only_spiro_contacts = all(len(set(a) & set(b)) <= 1 for i, a in enumerate(member_rings) for b in member_rings[i + 1 :])
    if len(member_rings) > 3 and not only_spiro_contacts and any(len(r) not in (3, 4, 5, 6, 7) for r in member_rings):
        if not any(sum(a in r for r in member_rings) > 2 for a in atoms):
            raise UnsupportedStructure("the numbering of this larger fused system is not verified here")


class _FusedSuffix(dict):
    def __missing__(self, key):
        if _is_variant(key):
            return carbo_suffix(spec_from_key(key), 1)
        raise KeyError(key)

    def __contains__(self, key):
        return dict.__contains__(self, key) or _is_variant(key)


_FUSED_SUFFIX = _FusedSuffix({
    "ide": "ide",
    "acid": "carboxylic acid",
    "sulfonic": "sulfonic acid",
    "amide": "carboxamide",
    **{name: f"carbo{name}" for name in _CHALCOGEN_AMIDE_CLASSES},
    "amidine": "carboximidamide",
    "hydrazonamide": "carbohydrazonamide",
    "imidohydrazide": "carboximidohydrazide",
    "hydrazonohydrazide": "carbohydrazonohydrazide",
    "sulfonamide": "sulfonamide",
    **{name: name for name in _CHALCOGEN_SULFONAMIDE_CLASSES},
    **{name: name for name in _CHALCOGEN_HYDRAZIDE.values()},
    **{name: name for name in _CHALCOGEN_HYDRAZIDINE.values()},
    **{name: name for name in _CHALCOGEN_IMIDAMIDE.values()},
    "hydrazide": "carbohydrazide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
    "alcohol": "ol",
    "peroxol": "peroxol",
    **{name: name for name in _PEROXOL_ANALOGUES},
    "ketone": "one",
    "imine": "imine",
    "thione": "thione",
    "selone": "selone",
    "tellone": "tellone",
    "thiol": "thiol",
    "selenol": "selenol",
    "tellurol": "tellurol",
    "amine": "amine",
})


def _fused_parent(mol, graph, principal, occurrences, here, n_names, stereo):
    """A fused ring system bearing the principal groups, named through the
    ring-system numbering and parent names used for diyl groups."""
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    if principal not in _FUSED_SUFFIX:
        raise UnsupportedStructure("this fused-ring parent is not supported by the chain engine yet")
    rings, atoms = _system_of(mol, here[0][1])
    if len(rings) > 1:
        _require_mancude_system(mol, atoms)
    on_system = [o for o in occurrences if o[1] in atoms]
    attach = [o[1] for o in on_system]
    blocked = set().union(*(o[2] for o in on_system))
    token = OMIT_SUFFIX_LOCANT.set(omits_group_locants(mol, atoms, len(attach)))
    try:
        found = evaluate_skeleton(mol, graph, "ring", rings, atoms, attach, blocked, _FUSED_SUFFIX[principal], n_names=n_names)
    finally:
        OMIT_SUFFIX_LOCANT.reset(token)
    if found is None:
        raise UnsupportedStructure("this fused ring system has no supported numbering")
    from ._diester_ring_diyl import PARENT_START

    count = len(on_system)
    return count, ((-count,), _without_stereo(found[1]), (None, None, None, 0, found[2], True, PARENT_START.get()))


class _RingSuffix(dict):
    def __missing__(self, key):
        if _is_variant(key):
            return key
        raise KeyError(key)


_RING_SUFFIX = _RingSuffix({
    "ide": "ide",
    "peroxoic": "peroxoic",
    "peroxol": "peroxol",
    **{name: name for name in _PEROXOL_ANALOGUES},
    "thioic": "thioic",
    "imidic": "imidic",
    "acid": "carboxylic_acid",
    "sulfonic": "sulfonic_acid",
    "amide": "amide",
    **{name: name for name in _CHALCOGEN_AMIDE_CLASSES},
    "amidine": "amidine",
    "hydrazonamide": "hydrazonamide",
    "imidohydrazide": "imidohydrazide",
    "hydrazonohydrazide": "hydrazonohydrazide",
    "sulfonamide": "sulfonamide",
    **{name: name for name in _CHALCOGEN_SULFONAMIDE_CLASSES},
    **{name: name for name in _CHALCOGEN_HYDRAZIDE.values()},
    **{name: name for name in _CHALCOGEN_HYDRAZIDINE.values()},
    **{name: name for name in _CHALCOGEN_IMIDAMIDE.values()},
    "hydrazide": "hydrazide",
    "nitrile": "nitrile",
    "aldehyde": "aldehyde",
    "ketone": "ketone",
    "thione": "thione",
    "selone": "selone",
    "tellone": "tellone",
    "imine": "imine",
    "alcohol": "alcohol",
    "thiol": "thiol",
    "selenol": "selenol",
    "tellurol": "tellurol",
    "amine": "amine",
})


# P-66.4.1.2.1: an acylated guanidine is named with an amide as the parent and the guanidine as a prefix, not multiplicatively
_DIACYLGUANIDINE = Chem.MolFromSmarts("[CX3](=O)[NX3][CX3](=[NX2])[NX3][CX3]=O")


def _best_ring(mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names=(), stereo=None):
    """(principal group count, (key, name, parts)) of the monocycle bearing
    the most principal groups, or None."""
    if not occurrences:
        return None
    ring_info = mol.GetRingInfo()
    candidates = []
    seen_systems = []
    for ring in ring_info.AtomRings():
        here = [o for o in occurrences if o[1] in ring]
        if not here:
            continue
        if any(ring_info.NumAtomRings(a) != 1 for a in ring):
            from ._diester_ring_diyl import _system_of

            _, system_atoms = _system_of(mol, here[0][1])
            if any(set(system_atoms) == seen for seen in seen_systems):
                continue
            seen_systems.append(set(system_atoms))
            here = [o for o in occurrences if o[1] in system_atoms]
        if any(ring_info.NumAtomRings(a) != 1 for a in ring):
            candidates.append((len(here), ring, None, here))
            continue
        spec = monocycle_spec(mol, ring)
        if spec is not None and _exocyclic_oxo(mol, set(ring)) and not (principal == "imine" and spec.kind == "cycloalkane"):
            spec = None
        candidates.append((len(here), ring, spec, here))
    if not candidates:
        return None
    top = max(c[0] for c in candidates)
    leading = [c for c in candidates if c[0] == top]
    if len(leading) > 1:
        # P-44.1.2.2, P-59.2.1.5: equal numbers of principal groups leave the senior ring system as the parent
        from ._diester_ring_diyl import _system_of

        def system_rank(candidate):
            return ring_seniority_key(mol, _system_of(mol, candidate[1][0])[1])

        senior = min(system_rank(c) for c in leading)
        leading = [c for c in leading if system_rank(c) == senior]
    if len(leading) == 1:
        return _ring_parent(mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names, stereo, leading[0])[:2]
    ranks = CanonicalRankAtoms(mol, breakTies=False)
    if len({ranks[next(iter(c[3]))[1]] for c in leading}) == 1 and not mol.HasSubstructMatch(_DIACYLGUANIDINE):
        raise UnsupportedStructure("several rings bear the principal group; a multiplicative name is needed")
    ranked = _isolated_candidates(
        [
            lambda candidate=candidate: _ring_parent(
                mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names, stereo, candidate
            )
            for candidate in leading
        ]
    )
    if not ranked:
        raise UnsupportedStructure("several rings bear the principal group and none has a supported name")
    count, found, _ = _choose_isolated(ranked, key=lambda r: r[2])
    return count, found


def _ring_parent(mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names, stereo, candidate):
    count, ring, spec, here = candidate
    if spec is None:
        found = _fused_parent(mol, graph, principal, occurrences, here, n_names, stereo)
        placed = found[1][2][4]
        taken = set().union(*(o[2] for o in here))
        attached = [
            placed[r]
            for r in ring
            for n in mol.GetAtomWithIdx(r).GetNeighbors()
            if r in placed and n.GetIdx() not in ring and n.GetIdx() not in taken
        ]
        locant_set = tuple(sorted(attached))
        return found[0], found[1], (
            -len(attached), tuple(_prefix_locant_key(loc) for loc in locant_set), locant_set, re.sub(r"[^a-z]", "", found[1][1]),
            tuple(int(n) for n in re.findall(r"\d+", found[1][1])),
        )
    owned = set().union(*(o[2] for o in here))
    ring_set = set(ring)
    roots = [
        (r, n.GetIdx())
        for r in ring
        for n in mol.GetAtomWithIdx(r).GetNeighbors()
        if n.GetIdx() not in ring_set and n.GetIdx() not in owned
    ]
    from ._diester_ring_diyl import _system_of

    for _, n in roots:
        neighbor_system = _system_of(mol, n)
        if neighbor_system is None or neighbor_system[1] & ring_set:
            continue
        if _bare_key(mol, neighbor_system[1]) == _bare_key(mol, ring_set):
            raise UnsupportedStructure("identical rings joined directly form a ring assembly (P-28)")
    entries = [
        (r, *name_branch(graph, n, r, halogens, aromatic_atoms, mol=mol, unsaturated=True)) for r, n in roots
    ]
    principal_atoms = [o[1] for o in here]
    ide_atoms = [a for a, w in IDE_EXTRA.get().items() if a in ring_set for _ in range(w)]
    best = None
    for locants in numberings(spec):
        key = (
            tuple(sorted(locants[a] for a in ide_atoms)),
            tuple(sorted(locants[a] for a in principal_atoms)),
            tuple(sorted(locants[r] for r, _, _ in entries)),
            _citation_key([(locants[r], name) for r, name, _ in entries]),
            _n_group_positions(n_names, locants),
            _stereo_rank(stereo, locants, True),
        )
        if best is None or key < best[0]:
            best = (key, locants)
    locants = best[1]
    suffix_name = _RING_SUFFIX[principal]
    suffix_locants = [locants[a] for a in principal_atoms]
    prefix_text = _ring_prefix_text(entries, locants, n_names, len(here))
    if ide_atoms:
        core = _ring_compound_core(spec, suffix_name, locants, ide_atoms, principal_atoms)
    elif (
        count == 1
        and not entries
        and not FORCE_LOCANTS.get()
        and (
            spec.kind == "cycloalkane"
            or (spec.kind == "benzene" and (suffix_name not in _RETAINED_BENZENE or NO_RETAINED_BENZENE.get()))
            or (spec.kind != "benzene" and omits_group_locants(mol, ring, 1))
        )
    ):
        word = _SUFFIX_WORDS[suffix_name]
        stem = spec.parent[:-1] if word[0] in "aeiouy" else spec.parent
        core = stem + word
    else:
        every_position_modified = (
            count > 1 and not entries and not n_names and not any(mol.GetAtomWithIdx(r).GetTotalNumHs() for r in ring)
        )
        core, _ = _suffix_text(spec.parent, suffix_name, suffix_locants, spec, every_position_modified)
    name = _join(prefix_text, core)
    count += len(ide_atoms)
    letters = re.sub(r"[^a-z]", "", name)
    prefix_locants = sorted(
        [
            *(_prefix_locant_key(locants[r]) for r, _, _ in entries),
            *(
                _prefix_locant_key("N")
                for entry in n_names
                if len(entry) < 3 or not isinstance(entry[2], tuple) or _n_position(entry[2], locants) is not None
            ),
        ]
    )
    rank = (-len(prefix_locants), tuple(prefix_locants), best[0][3], letters, tuple(int(n) for n in re.findall(r"\d+", name)))
    return count, ((-count,), name, (None, None, None, 0, locants, True, len(name) - len(core))), rank


def _ring_compound_core(spec, suffix_name, locants, ide_atoms, principal_atoms):
    if spec.kind not in ("benzene", "cycloalkane"):
        raise UnsupportedStructure("an anionic carbon beside a group on this ring is not supported yet")
    word = _SUFFIX_WORDS[suffix_name]
    ide_locants = sorted(locants[a] for a in ide_atoms)
    ide_word = multiplied_word(len(ide_locants), "ide")
    if word[0] in "aeiouy":
        ide_word = ide_word[:-1]
    group_locants = sorted(locants[a] for a in principal_atoms)
    stem = spec.parent[:-1]
    return (
        f"{stem}-{','.join(str(x) for x in ide_locants)}-{ide_word}-"
        f"{','.join(str(x) for x in group_locants)}-{multiplied_word(len(group_locants), word)}"
    )


def amine_unit_name(unit, attach):
    """(name, locant of the atom `attach`) of a unit that is one N-substituted amine, when that atom lies in the parent
    chosen for it; None otherwise (P-15.3: identical parent structures are multiplied)."""
    token = PREFERRED_AMINE_ATOM.set(attach)
    try:
        _, name, parts = _select(unit)
    except UnsupportedStructure:
        return None
    finally:
        PREFERRED_AMINE_ATOM.reset(token)
    locant = parts[4].get(attach) if isinstance(parts[4], dict) else None
    return None if locant is None or not name.endswith("amine") else (name, locant)


def _substituted_amine(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """A single secondary or tertiary amine: the ring or chain bound to
    nitrogen is the parent (ring first, P-44.1.2.2) and every other group on
    nitrogen is cited as an 'N-' prefix (P-62.2.2)."""
    amine_nitrogens = [a for a in mol.GetAtoms() if _substituted_amine_nitrogen(mol, a)]
    primaries = [a for a in groups.get("amine", {})] + [r for c, r, _ in ring_groups if c == "amine"]
    from ._substituents import ISOTOPE_LABELS

    if (
        len(amine_nitrogens) > 1
        and not primaries
        and not ISOTOPE_LABELS.get()
        and not SUBSTITUTED_AMINE_PREFIX.get()
        and not any(a.GetFormalCharge() for a in amine_nitrogens)
        and not _share_chain_backbone(mol, graph, amine_nitrogens)
    ):
        named = _backbone_polyamine(mol, graph, halogens, aromatic_atoms, amine_nitrogens)
        if named is not None:
            return named
        return _one_of_substituted_amines(mol, graph, halogens, aromatic_atoms, groups, ring_groups, amine_nitrogens)
    if (
        len(amine_nitrogens) > 2
        and not primaries
        and not ISOTOPE_LABELS.get()
        and not SUBSTITUTED_AMINE_PREFIX.get()
        and not any(a.GetFormalCharge() for a in amine_nitrogens)
    ):
        named = _largest_backbone_polyamine(mol, graph, halogens, aromatic_atoms, amine_nitrogens)
        if named is not None:
            return named
    if len(amine_nitrogens) > 1 and not SUBSTITUTED_AMINE_PREFIX.get() and (
        _labels_on_backbone_only(mol, graph, amine_nitrogens) or _labels_on_alike_nitrogens(mol, graph, amine_nitrogens)
    ):
        named = _backbone_polyamine(mol, graph, halogens, aromatic_atoms, amine_nitrogens)
        if named is not None:
            return named
    if (
        len(amine_nitrogens) > 1
        and AMINIUM.get() is True
        and not primaries
        and not ISOTOPE_LABELS.get()
        and not SUBSTITUTED_AMINE_PREFIX.get()
        and all(a.GetFormalCharge() == 1 or a.HasProp("_cationic_amine") for a in amine_nitrogens)
    ):
        named = _backbone_polyamine(mol, graph, halogens, aromatic_atoms, amine_nitrogens)
        if named is not None:
            return named
    primary_nitrogens = [
        mol.GetAtomWithIdx(n) for owned in (*groups.get("amine", {}).values(), *(o for c, _, o in ring_groups if c == "amine"))
        for n in owned if mol.GetAtomWithIdx(n).GetAtomicNum() == 7
    ]
    if len(amine_nitrogens) != 1 or (
        primaries
        and ISOTOPE_LABELS.get()
        and not _labels_on_backbone_only(mol, graph, amine_nitrogens)
        and not _labels_on_alike_nitrogens(mol, graph, [*amine_nitrogens, *primary_nitrogens])
    ):
        raise UnsupportedStructure("several amine groups with N-substitution are not handled by the chain engine")
    nitrogen = amine_nitrogens[0]
    n_idx = nitrogen.GetIdx()
    if nitrogen.IsInRing() or nitrogen.GetTotalNumHs() > 1:
        raise UnsupportedStructure("a ring nitrogen is not an acyclic amine parent")
    neighbors = [n.GetIdx() for n in nitrogen.GetNeighbors()]
    if any(
        mol.GetBondBetweenAtoms(n_idx, c).GetBondTypeAsDouble() != 1.0
        or (
            mol.GetAtomWithIdx(c).GetAtomicNum() == 6
            and (_double_oxygens(mol, c) or (is_functional_carbon(mol, c) and not mol.GetAtomWithIdx(c).IsInRing()))
        )
        or (mol.GetAtomWithIdx(c).GetAtomicNum() in (16, 34, 52) and _double_oxygens(mol, c))
        for c in neighbors
    ):
        raise UnsupportedStructure("this nitrogen is not a plain amine nitrogen")
    carbon_neighbors = [
        c for c in neighbors if mol.GetAtomWithIdx(c).GetAtomicNum() == 6 or _ring_nitrogen_parent(mol.GetAtomWithIdx(c))
    ]
    arms = {c: _arm_atoms(graph, c, n_idx) for c in carbon_neighbors}
    if sum(len(a) for a in arms.values()) != len(set().union(*arms.values())) or n_idx in set().union(*arms.values()):
        raise UnsupportedStructure("a nitrogen closing a ring is not an acyclic amine parent")
    context = ISOTOPE_LABELS.get()
    settled = set(context["consumed"]) if context else set()
    stereo_context = BRANCH_STEREO.get()
    used_before = set(stereo_context["used"]) if stereo_context else set()
    used_by = {}
    results = []
    for c in carbon_neighbors:
        others = [o for o in neighbors if o != c]
        if context:
            context["consumed"] = set(settled)
        if stereo_context:
            stereo_context["used"] &= used_before
        n_names = [name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True) for o in others]
        parent, mapped = _amine_parent_molecule(mol, arms[c], c, n_idx)
        if primaries:
            n_names = [(name, compound, ("N", mapped)) for name, compound in n_names]
        unlabelled = ISOTOPE_LABELS.set(None)
        try:
            with _renumbered_branch_stereo(sorted(set(arms[c]) | {n_idx})):
                result = _select(parent, None, n_names)
        finally:
            ISOTOPE_LABELS.reset(unlabelled)
        if stereo_context:
            used_by[c] = set(stereo_context["used"]) - used_before
        if primaries and (mapped not in result[2][4] or any(name not in result[1] for name, *_ in n_names)):
            continue
        ring = parent.GetAtomWithIdx(mapped).IsInRing()
        count = -result[0][0] if primaries else 0
        results.append((-count, ring, _ring_rank(parent, mapped), _chain_size(parent, mapped), result, c))
        if context:
            keep = sorted(set(arms[c]) | {n_idx})
            placed = {keep[new]: locant for new, locant in result[2][4].items() if new < len(keep)}
            placed[n_idx] = "N"
            results[-1] += (_substituent_isotope_locants(placed, mol, context["labels"]),)
    prefix_used = set()
    if stereo_context:
        stereo_context["used"] &= used_before
    if primaries:
        token = SUBSTITUTED_AMINE_PREFIX.set(True)
        try:
            as_prefix = _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, "amine", None, (), None)
        except UnsupportedStructure:
            as_prefix = None
        finally:
            SUBSTITUTED_AMINE_PREFIX.reset(token)
        if stereo_context:
            prefix_used = set(stereo_context["used"]) - used_before
            stereo_context["used"] &= used_before
        if as_prefix is not None:
            ring = as_prefix[2][5]
            anchor = next(iter(as_prefix[2][4]))
            results.append(
                (as_prefix[0][0], ring, _ring_rank(mol, anchor), len(as_prefix[2][4]) if not ring else 0, as_prefix, None)
            )
        if not results:
            raise UnsupportedStructure("no parent carries the amine groups of this substituted amine")
    preferred = PREFERRED_AMINE_ATOM.get()
    results.sort(
        key=lambda r: (
            r[0], not r[1], r[2], not _is_assembly_result(r[4]), -r[3], r[6] if len(r) > 6 else (),
            preferred is None or r[5] is None or preferred not in arms[r[5]],
            r[4][1],
        )
    )
    best = results[0]
    if stereo_context:
        stereo_context["used"] |= prefix_used if best[5] is None else used_by.get(best[5], set())
    if best[5] is None:
        return best[4]
    best = (best[1], best[2], best[3], best[4], best[5])
    if context:
        context["consumed"] = set(settled)
        others = [o for o in neighbors if o != best[4]]
        n_names = [name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True) for o in others]
        keep = sorted(set(arms[best[4]]) | {n_idx})
        new_index = {old: new for new, old in enumerate(keep)}
        parent_labels = {new_index[a]: e for a, e in context["labels"].items() if a in new_index}
        if parent_labels:
            parent, mapped = _amine_parent_molecule(mol, arms[best[4]], best[4], n_idx)
            if primaries:
                n_names = [(name, compound, ("N", mapped)) for name, compound in n_names]
            sub_context = {"labels": parent_labels, "consumed": set(), "mol": parent}
            sub_token = ISOTOPE_LABELS.set(sub_context)
            try:
                with _renumbered_branch_stereo(keep):
                    selected = _select(parent, None, n_names)
                    name = _with_labels(
                        parent, parent_labels, sub_context["consumed"], selected[1], selected[2],
                        lambda: _select(parent, None, n_names)[1:],
                    )
            finally:
                ISOTOPE_LABELS.reset(sub_token)
            context["consumed"].update(a for a in context["labels"] if a in new_index)
            context.setdefault("cited", set()).update(a for a in context["labels"] if a in new_index)
            placed = {}
            if stereo_context and (stereo_context["atoms"] or stereo_context["bonds"]):
                placed = {keep[new]: locant for new, locant in selected[2][4].items() if new < len(keep)}
            return selected[0], name, (*selected[2][:4], placed, *selected[2][5:])
    kept = sorted(set(arms[best[4]]) | {n_idx})
    positions = best[3][2][4]
    if isinstance(positions, dict) and positions:
        remapped = {kept[new]: locant for new, locant in positions.items() if new < len(kept)}
        return best[3][0], best[3][1], (*best[3][2][:4], remapped, *best[3][2][5:])
    return best[3]


def _hydroxy_or_alkoxy_on_nitrogen(atom):
    """The oxygen of an N-hydroxy or N-alkoxy group: an N-prefix of the amine, not a second parent (P-68.3.1.1.1.1)."""
    return (
        atom.GetAtomicNum() == 8
        and not atom.GetFormalCharge()
        and not atom.GetIsotope()
        and all(n.GetAtomicNum() in (6, 7) for n in atom.GetNeighbors())
        and sum(n.GetAtomicNum() == 7 for n in atom.GetNeighbors()) == 1
    )


def _is_assembly_result(result):
    """Whether the parent of `result` is a ring assembly, whose locants are primed (P-28.2.1)."""
    return any(isinstance(locant, PrimedLocant) for locant in result[2][4].values())


def _labels_on_alike_nitrogens(mol, graph, nitrogens):
    """Every amine nitrogen carries the same nuclides and no other atom off the backbone is modified, so the amine
    groups are modified alike and are multiplied (P-82.2.2.2)."""
    from ._substituents import ISOTOPE_LABELS

    labels = ISOTOPE_LABELS.get()["labels"]
    indexes = {n.GetIdx() for n in nitrogens}
    on_nitrogen = {tuple(sorted(_modifications(labels[i]))) for i in indexes if i in labels}
    if len(on_nitrogen) != 1 or not all(i in labels for i in indexes):
        return False
    stripped = {a: e for a, e in labels.items() if a not in indexes}
    token = ISOTOPE_LABELS.set({**ISOTOPE_LABELS.get(), "labels": stripped})
    try:
        return _labels_on_backbone_only(mol, graph, nitrogens)
    finally:
        ISOTOPE_LABELS.reset(token)


def _labels_on_backbone_only(mol, graph, nitrogens):
    """No nuclide, or only nuclides on the carbon chain or ring that joins the nitrogens, not on a nitrogen or a group of it."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if not context:
        return True
    indexes = {n.GetIdx() for n in nitrogens}
    carbons = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
    chain = carbons | {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (8, 16, 34, 52) and a.GetDegree() == 2}
    reached, stack = set(), [n for i in indexes for n in graph[i] if n in carbons]
    while stack:
        atom = stack.pop()
        if atom not in reached:
            reached.add(atom)
            stack.extend(n for n in graph[atom] if n in chain and n not in reached)
    arm_roots = {x for i in indexes for x in graph[i] if x not in reached}
    return all(a in reached and not any(a in _arm_atoms(graph, root, i) for i in indexes for root in arm_roots if root in graph[i]) for a in context["labels"])


def _backbone_polyamine(mol, graph, halogens, aromatic_atoms, nitrogens):
    """Several N-substituted amines on one chain: the chain with an amine suffix at each of those nitrogens is the parent
    and the other groups on the nitrogens are N-prefixes (P-62.2.2, P-62.2.4.1.2); None when the nitrogens do not hang on a single backbone."""
    from ._substituents import ISOTOPE_LABELS

    if any(
        n.IsInRing() or (n.GetFormalCharge() and not (AMINIUM.get() is True and n.GetFormalCharge() == 1))
        or n.GetTotalNumHs() > 1 or any(b.GetBondTypeAsDouble() != 1.0 for b in n.GetBonds())
        or any(
            not _hydroxy_or_alkoxy_on_nitrogen(x)
            and not _ring_nitrogen_parent(x)
            and (x.GetAtomicNum() != 6 or is_functional_carbon(mol, x.GetIdx()) or _double_oxygens(mol, x.GetIdx()))
            for x in n.GetNeighbors()
        )
        for n in nitrogens
    ):
        return None
    indexes = {n.GetIdx() for n in nitrogens}
    seen, components = set(indexes), []
    for start in range(mol.GetNumAtoms()):
        if start in seen:
            continue
        component, stack = {start}, [start]
        seen.add(start)
        while stack:
            for nxt in graph[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    component.add(nxt)
                    stack.append(nxt)
        components.append(component)
    joined = [c for c in components if sum(any(n in c for n in graph[i]) for i in indexes) > 1]
    if len(joined) != 1:
        return None
    backbone = joined[0]
    arms, carbons = {}, {}
    for i in indexes:
        attached = [x for x in graph[i] if x in backbone]
        if len(attached) != 1:
            return None
        carbons[i] = attached[0]
        arms[i] = [x for x in graph[i] if x not in backbone]
    def attempt(kept):
        cut = indexes - kept
        kept_atoms = set(kept)
        for i in kept:
            for arm in arms[i]:
                kept_atoms |= _arm_atoms(graph, arm, i)
        keep = sorted(backbone | cut | kept_atoms)
        editable = Chem.RWMol(mol)
        position_of_carbon = {c: position for position, c in enumerate(sorted({carbons[i] for i in cut}), start=1)}
        for i in cut:
            editable.GetAtomWithIdx(carbons[i]).SetAtomMapNum(position_of_carbon[carbons[i]])
            editable.GetAtomWithIdx(i).SetBoolProp("_cut_amine_n", True)
        for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
            editable.RemoveAtom(idx)
        parent = editable.GetMol()
        for a in parent.GetAtoms():
            if a.HasProp("_cut_amine_n"):
                a.SetFormalCharge(0)
                a.SetIntProp("_capacity", 2)
                a.SetNumExplicitHs(2)
                a.SetNoImplicit(True)
                a.SetBoolProp("_cationic_amine", True)
        Chem.SanitizeMol(parent)
        mapped = {a.GetAtomMapNum(): a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum()}
        for a in parent.GetAtoms():
            a.SetAtomMapNum(0)
        n_names = []
        for i in sorted(cut):
            for arm in arms[i]:
                name, compound = name_branch(graph, arm, i, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                n_names.append((name, compound, ("N", mapped[position_of_carbon[carbons[i]]], i)))
        token = SUBSTITUTED_AMINE_PREFIX.set(True)
        unlabelled = ISOTOPE_LABELS.set(None)
        try:
            with _renumbered_branch_stereo(keep):
                result = _select(parent, None, n_names)
        except UnsupportedStructure:
            return None, set()
        finally:
            ISOTOPE_LABELS.reset(unlabelled)
            SUBSTITUTED_AMINE_PREFIX.reset(token)
        missing = {i for i in cut if mapped[position_of_carbon[carbons[i]]] not in result[2][4]}
        if missing or any(name not in result[1] for name, *_ in n_names):
            return None, missing
        positions = {keep[new]: locant for new, locant in result[2][4].items() if new < len(keep)}
        return (result[0], result[1], (*result[2][:4], positions, *result[2][5:])), set()

    named, missing = attempt(frozenset())
    if named is None and missing and len(indexes - missing) > 1:
        # P-62.2.4.1.3: an amine beyond the largest chain of amine suffixes is an amino prefix that keeps its substituents
        named, _ = attempt(frozenset(missing))
    return named


def _chain_backbone_groups(mol, graph, nitrogens):
    """`nitrogens` split into the sets joined through acyclic carbons only: each set holds the amino groups of one chain
    (P-62.2.4.1.2)."""
    chain = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    by_index = {n.GetIdx(): n for n in nitrogens}
    joined = {first: set() for first in by_index}
    for first in by_index:
        reached, stack = {first}, [first]
        while stack:
            for nxt in graph[stack.pop()]:
                if nxt in chain and nxt not in reached:
                    reached.add(nxt)
                    stack.append(nxt)
        for atom in reached:
            for other in graph[atom]:
                if other in by_index and other != first:
                    joined[first].add(other)
                    joined[other].add(first)
    groups, seen = [], set()
    for start in by_index:
        if start in seen:
            continue
        group, stack = set(), [start]
        while stack:
            i = stack.pop()
            if i not in group:
                group.add(i)
                stack.extend(joined[i])
        seen |= group
        groups.append([by_index[i] for i in by_index if i in group])
    return groups


def _share_chain_backbone(mol, graph, nitrogens):
    return any(len(group) > 1 for group in _chain_backbone_groups(mol, graph, nitrogens))


def _largest_backbone_polyamine(mol, graph, halogens, aromatic_atoms, amine_nitrogens):
    """Amines on several backbones: the backbone with the most amine groups gives the polyamine parent (P-44.1.1) and the
    amines on the others are amino prefixes of the N-substituents (P-62.2.2); of equally large backbones the senior
    parent is chosen (P-62.2.4.1.3). None unless a single backbone holds a part of the amines."""
    groups = _chain_backbone_groups(mol, graph, amine_nitrogens)
    most = max(len(group) for group in groups)
    if most < 2 or most == len(amine_nitrogens):
        return None
    named = [
        _backbone_polyamine(mol, graph, halogens, aromatic_atoms, group) for group in groups if len(group) == most
    ]
    named = [result for result in named if result is not None]
    return min(named, key=_by_parent_then_prefix_count) if named else None


def _by_parent_then_prefix_count(result):
    """Rank of a chain parent among others: the chain key puts locants before the number of prefixes, which suits the
    numbering of one parent, but between parents the substituent count comes first (P-45.2.1, then P-45.2.2)."""
    key = result[0]
    return (*key[:7], key[13], *key[7:13], *key[14:], result[1])


def _one_of_substituted_amines(mol, graph, halogens, aromatic_atoms, groups, ring_groups, amine_nitrogens):
    """Several N-substituted amines, none primary: each in turn is the parent amine and the others are amino prefixes;
    the oxidized nitrogen of an N-oxide, when marked, is the only candidate (P-62.5)."""
    marked = [a for a in amine_nitrogens if a.HasProp("_oxidized_amine")]
    candidates = []
    for nitrogen in marked or amine_nitrogens:
        n_idx = nitrogen.GetIdx()
        if nitrogen.IsInRing() or nitrogen.GetTotalNumHs() > 1:
            continue
        neighbors = [n.GetIdx() for n in nitrogen.GetNeighbors()]
        carbons = [c for c in neighbors if mol.GetAtomWithIdx(c).GetAtomicNum() == 6]
        if not carbons or any(is_functional_carbon(mol, c) or _double_oxygens(mol, c) for c in carbons):
            continue
        for c in carbons:
            arm = _arm_atoms(graph, c, n_idx)
            if n_idx in arm or any(o in arm for o in neighbors if o != c):
                continue
            n_names = [
                name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                for o in neighbors
                if o != c
            ]
            parent, mapped = _amine_parent_molecule(mol, arm, c, n_idx)
            token = SUBSTITUTED_AMINE_PREFIX.set(True)
            try:
                with _renumbered_branch_stereo(sorted(arm | {n_idx})):
                    result = _select(parent, None, n_names)
            except UnsupportedStructure:
                continue
            finally:
                SUBSTITUTED_AMINE_PREFIX.reset(token)
            if mapped not in result[2][4]:
                continue
            ring = parent.GetAtomWithIdx(mapped).IsInRing()
            kept = sorted(set(arm) | {n_idx})
            positions = {kept[new]: locant for new, locant in result[2][4].items() if new < len(kept)}
            ranked = (not ring, _ring_rank(parent, mapped), -_chain_size(parent, mapped), result[1])
            candidates.append((ranked, (result[0], result[1], (*result[2][:4], positions, *result[2][5:]))))
    if not candidates:
        raise UnsupportedStructure("no parent carries the amine nitrogen of these N-substituted amines")
    return min(candidates, key=lambda candidate: candidate[0])[1]


def _boranyl_nitrogens(mol):
    """[(nitrogen, parent carbon, [hydride atoms])] for nitrogens bonded to one carbon and to mononuclear hydride groups
    such as boranyl, when together with the primary amines they number at least two, else []."""
    nitrogens = []
    primaries = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() or atom.IsInRing() or atom.GetIsAromatic():
            continue
        carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        hydrides = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in MONONUCLEAR_HYDRIDES]
        if hydrides and len(carbons) == 1 and _amidine_nitrogens(mol, carbons[0]) is not None:
            continue
        if hydrides and len(carbons) == 1 and len(carbons) + len(hydrides) == atom.GetDegree():
            nitrogens.append((atom.GetIdx(), carbons[0].GetIdx(), [h.GetIdx() for h in hydrides]))
        elif atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2 and len(carbons) == 1:
            primaries += 1
    return nitrogens if nitrogens and len(nitrogens) + primaries >= 2 else []


def boranyl_amine_parent_applies(mol):
    return bool(_boranyl_nitrogens(mol))


def _boranyl_amines(mol, graph, halogens, aromatic_atoms):
    """Nitrogens bonded to one carbon of a shared parent and carrying mononuclear hydride groups such as boranyl: the
    carbon parent whose amine groups outnumber those of any hydride parent is senior (P-44.1.1, P-68.1.5.2.1), and the
    hydride groups become N<locant> prefixes. None when the structure is not of this kind."""
    nitrogens = _boranyl_nitrogens(mol)
    if not nitrogens:
        return None
    removed = set()
    for n_idx, _, hydrides in nitrogens:
        for h in hydrides:
            arm = _arm_atoms(graph, h, n_idx)
            if arm & removed or any(a in arm for _, c, _ in nitrogens for a in (n_idx, c)):
                return None
            removed |= arm
    n_names = []
    editable = Chem.RWMol(mol)
    numbers = {}
    for n_idx, carbon, hydrides in nitrogens:
        number = numbers.setdefault(carbon, len(numbers) + 1)
        editable.GetAtomWithIdx(carbon).SetAtomMapNum(number)
        for h in hydrides:
            n_names.append((*name_branch(graph, h, n_idx, halogens, aromatic_atoms, mol=mol), number, n_idx))
    kept = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    parent = editable.GetMol()
    try:
        Chem.SanitizeMol(parent)
    except Exception:
        return None
    if len(Chem.GetMolFrags(parent)) != 1:
        return None
    mapped = {a.GetAtomMapNum(): a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum()}
    for a in parent.GetAtoms():
        a.SetAtomMapNum(0)
    tagged = [(name, compound, ("N", mapped[number], n_idx)) for name, compound, number, n_idx in n_names]
    with _renumbered_branch_stereo(kept):
        result = _select(parent, None, tagged)
    positions = result[2][4]
    if isinstance(positions, dict) and any(carbon not in positions for _, _, (_, carbon, _) in tagged):
        return None
    if isinstance(positions, dict) and positions:
        positions = {kept[new]: locant for new, locant in positions.items() if new < len(kept)}
        return result[0], result[1], (*result[2][:4], positions, *result[2][5:])
    return result


@contextlib.contextmanager
def _renumbered_branch_stereo(kept):
    """The stereo elements cited inside substituents, renumbered for the molecule made of the sorted atoms `kept`
    (as `_amine_parent_molecule` renumbers them); what is cited there counts as cited in the whole molecule."""
    from ._substituents import BRANCH_STEREO

    outer = BRANCH_STEREO.get()
    if not outer:
        yield
        return
    new_of = {old: new for new, old in enumerate(kept)}
    old_of = {new: old for old, new in new_of.items()}
    inner = {
        **outer,
        "atoms": {new_of[a]: code for a, code in outer["atoms"].items() if a in new_of},
        "bonds": {(new_of[a], new_of[b]): code for (a, b), code in outer["bonds"].items() if a in new_of and b in new_of},
        "used": set(),
    }
    token = BRANCH_STEREO.set(inner)
    try:
        yield
    finally:
        BRANCH_STEREO.reset(token)
        for kind, where in inner["used"]:
            outer["used"].add((kind, old_of[where] if kind == "atom" else tuple(old_of[w] for w in where)))


def _amine_parent_molecule(mol, atoms, carbon, n_idx):
    editable = Chem.RWMol(mol)
    editable.GetAtomWithIdx(carbon).SetAtomMapNum(1)
    editable.GetAtomWithIdx(n_idx).SetBoolProp("_cut_amine_n", True)
    keep = set(atoms) | {n_idx}
    for idx in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(idx)
    parent = editable.GetMol()
    for a in parent.GetAtoms():
        if a.HasProp("_cut_amine_n") and not a.IsInRing():
            a.SetIntProp("_capacity", mol.GetAtomWithIdx(n_idx).GetTotalNumHs())
            a.SetFormalCharge(0)
            a.SetNumExplicitHs(2)
            a.SetNoImplicit(True)
            a.SetBoolProp("_cationic_amine", True)
    if any(a.HasProp("_ring_cation_centre") or a.HasProp("_anion_word") for a in parent.GetAtoms()):
        parent.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(parent)
    else:
        Chem.SanitizeMol(parent)
    mapped = next(a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum() == 1)
    for a in parent.GetAtoms():
        if not is_prefix_placeholder(a):
            a.SetAtomMapNum(0)
    return parent, mapped


def _ring_rank(mol, idx):
    """P-44.2 key of the ring system holding atom `idx` (smaller is senior); `()` for an acyclic atom."""
    if not mol.GetAtomWithIdx(idx).IsInRing():
        return ()
    from ._diester_ring_diyl import _system_of

    return ring_seniority_key(mol, _system_of(mol, idx)[1])


def _chain_size(mol, idx):
    if mol.GetAtomWithIdx(idx).IsInRing():
        return 0
    graph = adjacency(mol)
    eligible = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    return max((len(p) for p in _paths(graph, eligible) if idx in p), default=0)


def _unit_parts(unit, attach, arm):
    """Parts of the name of the multiplied unit cut from `arm` (atoms of the full molecule); its isotopic descriptor
    is inserted in front of the parent hydride name (P-82.2.1)."""
    from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    index = {atom: i for i, atom in enumerate(sorted(arm))}
    labels = {index[a]: entry for a, entry in context["labels"].items() if a in index} if context else {}
    token = ISOTOPE_LABELS.set(None)
    stereo_token = BRANCH_STEREO.set(None)
    try:
        if not labels:
            return _select(unit, attach)[2]
        entries = _isotope_entries(labels)
        inner = {"labels": labels, "consumed": set(), "mol": unit}
        ISOTOPE_LABELS.set(inner)
        BRANCH_STEREO.set({"atoms": {}, "bonds": {}, "used": set()})
        _, name, parts = _select(unit, attach, stereo=entries)
        inner["consumed"].clear()
        inner.pop("cache", None)
        _, name, parts = _select(unit, attach, stereo=entries)
        name = _with_labels(unit, labels, inner["consumed"], name, parts, lambda: _select(unit, attach, stereo=entries)[1:])
    finally:
        BRANCH_STEREO.reset(stereo_token)
        ISOTOPE_LABELS.reset(token)
    prefix, body, tail, locant = parts[:4]
    stem = name[: len(name) - len(tail)] if tail else name
    return (stem[: len(prefix)] if stem.startswith(prefix) else "", stem[len(prefix):] if stem.startswith(prefix) else stem, tail, locant)


def _constitution_smiles(unit):
    """Canonical SMILES of a cut unit without its stereo, since a configuration is cited after the units are matched."""
    bare = Chem.Mol(unit)
    Chem.RemoveStereochemistry(bare)
    return Chem.MolToSmiles(bare)


def _modified_unit_smiles(unit, atoms):
    """Canonical SMILES of a cut unit, its isotopically modified atoms marked, so that units modified differently are
    not identical (P-82.2.2)."""
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    labels = {a: entry for a, entry in context["labels"].items() if a in atoms} if context else {}
    if not labels:
        return Chem.MolToSmiles(unit)
    marked, codes = Chem.Mol(unit), {}
    for index, atom in zip(sorted(atoms), marked.GetAtoms()):
        if index in labels:
            entry = labels[index]
            atom.SetAtomMapNum(10 + codes.setdefault((entry["skeleton"], tuple(sorted(entry["H"].items()))), len(codes)))
    return Chem.MolToSmiles(marked)


def _identical_group_units(mol, graph, group_atoms):
    """True when two disjoint, identical, non-trivial fragments hanging off
    different acyclic bonds each carry a principal group -- the shape that
    P-15.3 names multiplicatively."""
    sides = []
    for bond in mol.GetBonds():
        if bond.IsInRing() or bond.GetBondTypeAsDouble() != 1.0:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for root, other in ((b, a), (a, b)):
            neighbor = mol.GetAtomWithIdx(other)
            if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                continue
            if neighbor.GetFormalCharge() > 0:
                # P-74.1.3: a cationic centre is cited as a prefix of the anionic parent, never as a linking group
                continue
            atoms = _arm_atoms(graph, root, other)
            if other in atoms or not atoms & group_atoms or not atoms - group_atoms:
                continue
            canonical, _ = _unit_molecule(mol, atoms, root)
            sides.append((bond.GetIdx(), atoms, _modified_unit_smiles(canonical, atoms), root))
    plain = Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=False)
    chiral = Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True)
    for i, (bi, ai, si, ri) in enumerate(sides):
        for bj, aj, sj, rj in sides[:i]:
            if bi != bj and si == sj and not ai & aj:
                if _unsymmetric_chalcogen_linker(mol, graph, ai, aj):
                    continue
                if plain[ri] == plain[rj] and chiral[ri] != chiral[rj]:
                    continue
                # a multiplied parent must express every principal group (P-15.6.1.5)
                extra = set().union(*(ak for bk, ak, sk, _ in sides if sk == si and not ak & (ai | aj)))
                units = ai | aj | extra
                if plain[ri] != plain[rj] and not extra and _unsymmetrical_branches(mol, set(range(mol.GetNumAtoms())) - units, ai, aj):
                    continue
                if group_atoms <= units:
                    return True
    return False


def _unsymmetrical_branches(mol, linker, first_arm, second_arm):
    """Whether the components of the linking group (ring system, carbon chain, run of one heteroatom) differ on the two
    sides of its centre, which no multiplicative group can express (P-15.3.1.2.2.3, P-15.3.3.1, P-15.3.4.1.1)."""
    parent = {a: a for a in linker}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def same_kind(bond):
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ends = (mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(b))
        return bond.IsInRing() or (ends[0].GetAtomicNum() == ends[1].GetAtomicNum() and not any(e.IsInRing() for e in ends))

    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in linker and bond.GetEndAtomIdx() in linker and same_kind(bond):
            parent[find(bond.GetBeginAtomIdx())] = find(bond.GetEndAtomIdx())
    members = {}
    for a in linker:
        members.setdefault(find(a), set()).add(a)
    neighbours = {c: set() for c in members}
    for c, atoms in members.items():
        for a in atoms:
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                if n.GetIdx() in linker and find(n.GetIdx()) != c:
                    neighbours[c].add(find(n.GetIdx()))
    touching = [[c for c, atoms in members.items() if any(n.GetIdx() in arm for a in atoms for n in mol.GetAtomWithIdx(a).GetNeighbors())] for arm in (first_arm, second_arm)]
    if len(members) < 2 or any(len(t) != 1 for t in touching) or any(len(n) > 2 for n in neighbours.values()):
        return False
    order, seen = [touching[0][0]], {touching[0][0]}
    while order[-1] != touching[1][0]:
        onward = [n for n in neighbours[order[-1]] if n not in seen]
        if len(onward) != 1:
            return False
        order.append(onward[0])
        seen.add(onward[0])
    if len(order) != len(members):
        return False

    def signature(component, outside):
        editable = Chem.RWMol()
        index = {}
        for a in sorted(members[component]):
            index[a] = editable.AddAtom(Chem.Atom(mol.GetAtomWithIdx(a).GetAtomicNum()))
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a in index and b in index:
                editable.AddBond(index[a], index[b], bond.GetBondType())
            elif a in index or b in index:
                inner, other = (a, b) if a in index else (b, a)
                dummy = Chem.Atom(0)
                dummy.SetAtomMapNum(1 if other in outside else 2)
                editable.AddBond(index[inner], editable.AddAtom(dummy), bond.GetBondType())
        return Chem.MolToSmiles(editable.GetMol(), isomericSmiles=False)

    outside = [set(first_arm), set(second_arm)]
    low, high = 0, len(order) - 1
    while high - low >= 1:
        if signature(order[low], outside[0]) != signature(order[high], outside[1]):
            return True
        outside[0] |= members[order[low]]
        outside[1] |= members[order[high]]
        low, high = low + 1, high - 1
    return False


def _unsymmetric_chalcogen_linker(mol, graph, first_arm, second_arm):
    """Whether the rest of the molecule is one chain of chalcogen atoms whose element sequence reads differently from
    each end, so no divalent multiplying group names it (P-65.7.5.1: 1-[(acetylperoxy)sulfanyl]ethan-1-one)."""
    rest = set(range(mol.GetNumAtoms())) - first_arm - second_arm
    if len(rest) < 2 or any(mol.GetAtomWithIdx(a).GetAtomicNum() not in (8, 16, 34, 52) for a in rest):
        return False
    start = next((a for a in rest if first_arm & set(graph[a])), None)
    if start is None:
        return False
    order, previous = [start], None
    while True:
        onward = [n for n in graph[order[-1]] if n in rest and n != previous]
        if not onward:
            break
        previous = order[-1]
        order.append(onward[0])
    elements = [mol.GetAtomWithIdx(a).GetAtomicNum() for a in order]
    return len(order) == len(rest) and elements != elements[::-1]


_LINKER_WORDS = {8: "oxy", 16: "sulfanediyl", 34: "selanediyl"}
_DICHALCOGEN_WORDS = {8: "peroxy", 16: "disulfanediyl", 34: "diselanediyl"}


_CHALCOGEN_CHAIN_STEMS = {8: "oxidane", 16: "sulfane", 34: "selane", 52: "tellane"}
_THIO_PREFIXES = {16: "thio", 34: "seleno", 52: "telluro"}


def _chalcogen_linker_word(elements):
    """'trioxidanediyl', 'tetrasulfanediyl', 'dithioxanediyl' for a symmetric chain of chalcogen atoms, else None."""
    if elements != elements[::-1]:
        return None
    if len(set(elements)) == 1:
        return multiplying_prefix(len(elements)) + _CHALCOGEN_CHAIN_STEMS[elements[0]] + "diyl"
    if len(elements) == 3 and elements[1] == 8 and elements[0] in _THIO_PREFIXES:
        return "di" + _THIO_PREFIXES[elements[0]] + "xanediyl"
    half = len(elements) // 2
    centre, arms = (elements[half : half + 1], elements[:half]) if len(elements) % 2 else (elements[half - 1 : half + 1], elements[: half - 1])
    return _chalcogen_run_word(centre, True) + "bis(" + _chalcogen_run_word(arms[::-1], False) + ")"


def _chalcogen_run_word(elements, central):
    """'sulfanediyl', 'disulfanediyl', 'oxy', 'peroxy', 'oxyselanediyl': the divalent group of a run of chalcogen atoms,
    each stretch of identical atoms one word (P-29.3.2.2, P-15.3.1.2)."""
    words = []
    for z, group in itertools.groupby(elements):
        count = len(list(group))
        if z == 8:
            words.append({1: "oxy", 2: "peroxy"}.get(count) or multiplying_prefix(count) + "oxidanediyl")
        else:
            words.append((multiplying_prefix(count) if count > 1 else "") + _CHALCOGEN_CHAIN_STEMS[z][:-1] + "ediyl")
    return "".join(words)


def _diacyl_chalcogen_chain_name(mol, graph, stereo):
    """P-65.7.5.1: two identical acyl groups joined by a symmetric chain of three or more chalcogen atoms are
    pseudoketones named multiplicatively (1,1'-trioxidanediyldi(ethan-1-one))."""
    chain = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() in _CHALCOGEN_CHAIN_STEMS
        and a.GetDegree() == 2
        and not a.IsInRing()
        and not a.GetFormalCharge()
        and not a.GetTotalNumHs()
    ]
    if len(chain) < 3:
        return None
    chain_set = set(chain)
    ends = [a for a in chain if sum(1 for n in graph[a] if n in chain_set) == 1]
    if len(ends) != 2 or any(mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble() != 1.0 for a in chain for b in graph[a]):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(chain):
        onward = [n for n in graph[order[-1]] if n in chain_set and n != previous]
        if len(onward) != 1:
            return None
        previous = order[-1]
        order.append(onward[0])
    word = _chalcogen_linker_word([mol.GetAtomWithIdx(a).GetAtomicNum() for a in order])
    roots = [next(n for n in graph[end] if n not in chain_set) for end in (order[0], order[-1])]
    if word is None or any(
        mol.GetAtomWithIdx(r).GetAtomicNum() != 6 or not _double_oxygens(mol, r) for r in roots
    ):
        return None
    arms = [_arm_atoms(graph, r, end) for r, end in zip(roots, (order[0], order[-1]))]
    if arms[0] & arms[1] or len(arms[0]) + len(arms[1]) + len(chain) != mol.GetNumAtoms():
        return None
    found = _multiplied_units(mol, arms, roots, chalcogen_cut=True)
    if found is None:
        return None
    units, parts = found
    stereo_text = ""
    plain_stereo = [element for element in stereo if element[0] != "isotope"]
    if plain_stereo:
        arms_named = _stereo_arms(mol, plain_stereo, arms, units)
        if arms_named is None:
            return None
        parts, stereo_text = arms_named
    prefix, body, tail, locant = parts[:4]
    if locant is None:
        return None
    marked = _linker_descriptor(order)
    _consume_labels(order)
    if marked:
        word = enclose(marked + word)
    elif "bis(" in word:
        word = f"[{word}]"
    text = prefix + body
    if prefix:
        return f"{stereo_text}{locant},{locant}'-{word}{_bis_unit(text, tail)}"
    return f"{stereo_text}{locant},{locant}'-{word}di{unit_phrase(text, tail)}"


def _arm_atoms(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _multiplicative_name(mol, stereo=None):
    """P-15.3: identical chain parents joined through one oxygen, chalcogen,
    NH or N atom are named multiplicatively (2,2'-oxydi(ethan-1-ol))."""
    from ._chain_multiplicative import _principal_atoms

    graph = adjacency(mol)
    for z in mol.GetAtoms():
        if z.IsInRing() or z.GetFormalCharge() or z.GetIsotope():
            continue
        neighbors = [n.GetIdx() for n in z.GetNeighbors()]
        hydrogens = z.GetTotalNumHs()
        number = z.GetAtomicNum()
        blockers = [z.GetIdx()]
        extra = 1
        partner = next(
            (
                n
                for n in z.GetNeighbors()
                if n.GetAtomicNum() == number and number in _DICHALCOGEN_WORDS and n.GetDegree() == 2 and not n.GetTotalNumHs()
            ),
            None,
        )
        if partner is not None and len(neighbors) == 2 and hydrogens == 0:
            if partner.GetIdx() < z.GetIdx() or partner.IsInRing():
                continue
            linker, arms = _DICHALCOGEN_WORDS[number], 2
            neighbors = [n for n in neighbors if n != partner.GetIdx()] + [
                n.GetIdx() for n in partner.GetNeighbors() if n.GetIdx() != z.GetIdx()
            ]
            blockers = [z.GetIdx(), partner.GetIdx()]
            extra = 2
        elif number in _LINKER_WORDS and len(neighbors) == 2 and hydrogens == 0:
            linker, arms = _LINKER_WORDS[number], 2
        elif number == 7 and len(neighbors) == 2 and hydrogens == 1:
            linker, arms = "azanediyl", 2
        elif number == 7 and len(neighbors) == 3 and hydrogens == 0:
            linker, arms = "nitrilo", 3
        else:
            continue
        if len(neighbors) != arms or any(
            mol.GetAtomWithIdx(n).GetAtomicNum() != 6
            or _double_oxygens(mol, n)
            or is_functional_carbon(mol, n)
            or mol.GetAtomWithIdx(n).GetIsAromatic()
            for n in neighbors
        ):
            continue
        parts_atoms = [_arm_atoms(graph, n, blockers[0] if i == 0 or len(blockers) == 1 else blockers[1]) for i, n in enumerate(neighbors)]
        if sum(len(p) for p in parts_atoms) != mol.GetNumAtoms() - extra or any(
            parts_atoms[i] & parts_atoms[j] for i in range(arms) for j in range(i)
        ) or len({len(p) for p in parts_atoms}) != 1:
            continue
        from ._substituents import ISOTOPE_LABELS

        isotope_context = ISOTOPE_LABELS.get()
        labelled = bool(isotope_context)
        units = [_unit_molecule(mol, atoms, n, refill=labelled) for atoms, n in zip(parts_atoms, neighbors)]
        if labelled:
            keys = {_modified_unit_smiles(unit[0], atoms) for unit, atoms in zip(units, parts_atoms)}
        else:
            keys = {Chem.MolToSmiles(unit[0], isomericSmiles=False) for unit in units}
        if len(keys) != 1:
            continue
        if number == 7 and (_principal_atoms(units[0][0]) or (None,))[0] == "amine":
            continue
        unit, attach = units[0]
        try:
            if labelled:
                parts = _unit_parts(unit, attach, parts_atoms[0])
            else:
                _, _, parts = _select(unit, attach)
        except UnsupportedStructure:
            continue
        stereo_text = ""
        plain_stereo = [element for element in stereo or () if element[0] != "isotope"]
        if plain_stereo:
            arms_named = _stereo_arms(mol, plain_stereo, parts_atoms, units)
            if arms_named is None:
                continue
            parts, stereo_text = arms_named
        if labelled:
            linker_atoms = [z.GetIdx()] if extra == 1 else [z.GetIdx(), partner.GetIdx()]
            _consume_labels(linker_atoms + [a for atoms in parts_atoms for a in atoms])
            marked = _linker_descriptor(linker_atoms, {z.GetIdx(): hydrogens})
            if marked:
                linker = enclose(marked + linker)
        prefix, body, tail, locant = parts[:4]
        lead = ",".join(str(locant) + "'" * i for i in range(arms)) + "-" if locant is not None else ""
        text = prefix + body
        if prefix:
            word = {2: "bis", 3: "tris"}[arms]
            if tail.startswith(" "):
                return f"{stereo_text}{lead}{linker}{word}{enclose(text + tail)}"
            return f"{stereo_text}{lead}{linker}{word}{enclose(text)}{tail}"
        word = {2: "di", 3: "tri"}[arms]
        return f"{stereo_text}{lead}{linker}{word}{unit_phrase(text, tail)}"
    from ._chain_multiplicative import chain_multiplicative_name

    return (
        _diacyl_chalcogen_chain_name(mol, graph, stereo)
        or _natural_product_linker_name(mol, graph, stereo)
        or _ring_linker_name(mol, graph, stereo)
        or _composite_linker_name(mol, graph, stereo)
        or chain_multiplicative_name(mol, stereo)
    )


_CENTER_WORDS = {8: "oxy", 16: "sulfanediyl", 34: "selanediyl"}


def _composite_linker_name(mol, graph, stereo):
    """Two identical chain parents on benzene rings joined through one atom:
    2,2'-[oxybis(4,1-phenylene)]di(ethan-1-ol), 2,2'-[methylenebis(4,1-phenylene)]diethanoic acid."""
    for center in mol.GetAtoms():
        z = center.GetAtomicNum()
        if center.IsInRing() or center.GetDegree() != 2 or center.GetFormalCharge() or center.GetIsotope():
            continue
        if z in _CENTER_WORDS and center.GetTotalNumHs() == 0:
            word = _CENTER_WORDS[z]
        elif z == 7 and center.GetTotalNumHs() == 1:
            word = "azanediyl"
        elif z == 6 and center.GetTotalNumHs() == 2:
            word = "methylene"
        else:
            continue
        index = center.GetIdx()
        sides = []
        for ring_atom in graph[index]:
            ring = next((r for r in mol.GetRingInfo().AtomRings() if ring_atom in r), None)
            if ring is None or any(mol.GetRingInfo().NumAtomRings(a) != 1 for a in ring):
                sides = []
                break
            spec = monocycle_spec(mol, ring)
            if spec is None or spec.kind != "benzene":
                sides = []
                break
            ring_set = set(ring)
            others = [
                (r, n.GetIdx())
                for r in ring
                for n in mol.GetAtomWithIdx(r).GetNeighbors()
                if n.GetIdx() not in ring_set and n.GetIdx() != index
            ]
            if len(others) != 1:
                sides = []
                break
            sides.append((ring, ring_atom, others[0]))
        if len(sides) != 2:
            continue
        arms = [_arm_atoms(graph, root, r) for _, _, (r, root) in sides]
        ring_total = sum(len(ring) for ring, _, _ in sides)
        if arms[0] & arms[1] or len(arms[0]) + len(arms[1]) + ring_total + 1 != mol.GetNumAtoms():
            continue
        if any(
            mol.GetAtomWithIdx(root).GetAtomicNum() != 6
            or mol.GetAtomWithIdx(root).GetIsAromatic()
            or mol.GetBondBetweenAtoms(r, root).GetBondTypeAsDouble() != 1.0
            for _, _, (r, root) in sides
        ):
            continue
        found = _multiplied_units(mol, arms, [root for _, _, (_, root) in sides])
        if found is None:
            continue
        units, parts = found
        components = []
        for ring, ring_atom, (r, root) in sides:
            ring_isotopes = [element for element in stereo if element[0] == "isotope" and element[1] in ring]
            named = name_ring_component(
                mol, ring, [(r, root), (ring_atom, index)], [], None, directed=(r, ring_atom), stereo=ring_isotopes
            )
            components.append(named[0] if named else None)
            _consume_labels(ring)
        if components[0] is None or components[0] != components[1]:
            continue
        stereo_text = ""
        plain_stereo = [element for element in stereo if element[0] != "isotope"]
        if plain_stereo:
            arms_named = _stereo_arms(mol, plain_stereo, arms, units)
            if arms_named is None:
                continue
            parts, stereo_text = arms_named
        _consume_labels([index])
        prefix, body, tail, locant = parts[:4]
        lead = f"{locant},{locant}'-" if locant is not None else ""
        linker = enclose(f"{_linker_descriptor([index], {index: center.GetTotalNumHs()})}{word}bis{enclose(components[0])}")
        text = prefix + body
        if prefix:
            return f"{stereo_text}{lead}{linker}{_bis_unit(text, tail)}"
        return f"{stereo_text}{lead}{linker}di{unit_phrase(text, tail)}"
    return None


def _is_phane_candidate(mol, ring_set):
    return len(ring_set) > 12 and any(len(r) == 6 and set(r) <= ring_set for r in mol.GetRingInfo().AtomRings())


def _natural_product_linker_name(mol, graph, stereo):
    """Two identical chain parents on an Appendix 3 parent: 3,3'-(yohimban-14,18-diyl)dipropanoic acid."""
    from ._appendix3_skeletons import appendix3_central_group

    found = appendix3_central_group(mol, graph)
    if found is None or len(found[1]) != 2:
        return None
    linker, attachments, central = found
    arms = [_arm_atoms(graph, root, r) for r, root in attachments]
    if arms[0] & arms[1] or len(arms[0]) + len(arms[1]) + len(central) != mol.GetNumAtoms():
        return None
    if any(
        mol.GetAtomWithIdx(root).GetAtomicNum() != 6
        or mol.GetAtomWithIdx(root).GetIsAromatic()
        or mol.GetBondBetweenAtoms(r, root).GetBondTypeAsDouble() != 1.0
        for r, root in attachments
    ):
        return None
    found = _multiplied_units(mol, arms, [root for _, root in attachments])
    if found is None:
        return None
    units, parts = found
    arm_atoms = arms[0] | arms[1]
    arm_stereo = [element for element in stereo if element[0] != "isotope" and _stereo_within(element, arm_atoms)]
    stereo_text = ""
    if arm_stereo:
        arms_named = _stereo_arms(mol, arm_stereo, arms, units)
        if arms_named is None:
            return None
        parts, stereo_text = arms_named
    prefix, body, tail, locant = parts[:4]
    lead = f"{locant},{locant}'-" if locant is not None else ""
    text = prefix + body
    if prefix:
        return f"{stereo_text}{lead}{enclose(linker)}{_bis_unit(text, tail)}"
    return f"{stereo_text}{lead}{enclose(linker)}di{unit_phrase(text, tail)}"


def _stereo_within(element, atoms):
    kind, where, _ = element
    return where in atoms if kind in ("atom", "isotope") else all(a in atoms for a in where)


def _ring_linker_name(mol, graph, stereo):
    """Two identical chain parents on one monocycle: 2,2'-(1,4-phenylene)di(ethan-1-ol)."""
    from ._multiplicative import _ring_systems

    ring_info = mol.GetRingInfo()
    for system in _ring_systems(mol):
        ring_set = set(system)
        member_rings = [list(r) for r in ring_info.AtomRings() if set(r) <= ring_set]
        fused = len(member_rings) > 1
        if fused:
            try:
                _require_mancude_system(mol, ring_set)
            except UnsupportedStructure:
                if not _is_phane_candidate(mol, ring_set):
                    continue
        ring = member_rings[0]
        attachments = [
            (r, n.GetIdx())
            for r in ring_set
            for n in mol.GetAtomWithIdx(r).GetNeighbors()
            if n.GetIdx() not in ring_set
        ]
        if len(attachments) != 2:
            continue
        spec = None if fused else monocycle_spec(mol, ring)
        arms = [_arm_atoms(graph, root, r) for r, root in attachments]
        if arms[0] & arms[1] or len(arms[0]) + len(arms[1]) + len(ring_set) != mol.GetNumAtoms():
            continue
        if any(
            mol.GetAtomWithIdx(root).GetAtomicNum() != 6
            or mol.GetAtomWithIdx(root).GetIsAromatic()
            or mol.GetBondBetweenAtoms(r, root).GetBondTypeAsDouble() != 1.0
            for r, root in attachments
        ):
            continue
        from ._substituents import ISOTOPE_LABELS

        labelled = bool(ISOTOPE_LABELS.get())
        units = [_unit_molecule(mol, atoms, root, refill=labelled) for atoms, (_, root) in zip(arms, attachments)]
        if len({_modified_unit_smiles(unit[0], atoms) for (unit, atoms) in zip(units, arms)}) != 1:
            continue
        unit, attach = units[0]
        try:
            parts = _unit_parts(unit, attach, arms[0])
        except UnsupportedStructure:
            continue
        arm_atoms = arms[0] | arms[1]
        ring_stereo = [element for element in stereo if _stereo_within(element, ring_set)]
        arm_isotopes = [element for element in stereo if element[0] == "isotope" and element[1] in arm_atoms]
        arm_stereo = [element for element in stereo if element[0] != "isotope" and _stereo_within(element, arm_atoms)]
        stereo_text = ""
        if arm_stereo:
            arms_named = _stereo_arms(mol, arm_stereo, arms, units)
            if arms_named is None:
                continue
            parts, stereo_text = arms_named
        if spec is not None:
            linker = name_ring_component(mol, ring, attachments, [], None, stereo=ring_stereo)
        else:
            from ._diester_ring_diyl import evaluate_skeleton

            try:
                found = evaluate_skeleton(
                    mol, graph, "ring", member_rings, ring_set, [r for r, _ in attachments], arms[0] | arms[1], "yl"
                )
            except UnsupportedStructure:
                found = None
            linker = (found[1], False) if found else None
            if linker is None:
                from ._phane_general import phane_diyl_name

                phane = phane_diyl_name(mol, ring_set, [r for r, _ in attachments], arms[0] | arms[1])
                linker = (phane, True) if phane else None
        if linker is None:
            continue
        prefix, body, tail, locant = parts[:4]
        lead = f"{locant},{locant}'-" if locant is not None else ""
        text = prefix + body
        isotope_context = ISOTOPE_LABELS.get()
        if isotope_context is not None:
            isotope_context["consumed"].update(atom for kind, atom, _ in ring_stereo + arm_isotopes if kind == "isotope")
        linker_text = f"({linker[0]})" if linker[1] and not ring_stereo else enclose(linker[0])
        if prefix:
            return f"{stereo_text}{lead}{linker_text}{_bis_unit(text, tail)}"
        return f"{stereo_text}{lead}{linker_text}di{unit_phrase(text, tail)}"
    return None


def _stereo_arms(mol, stereo, parts_atoms, units):
    """(parts of the first unit, '(2R,2'S)-' text) for identical units whose stereo differs; None when the units name
    differently once their configurations are cited, so that the substitutive name is used (P-93.6)."""
    from ._substituents import BRANCH_STEREO

    per_arm = []
    token = BRANCH_STEREO.set(None)
    try:
        for atoms, (unit, attach) in zip(parts_atoms, units):
            index = {orig: i for i, orig in enumerate(sorted(atoms))}
            arm = []
            for kind, where, code in stereo:
                if kind == "atom" and where in index:
                    arm.append((kind, index[where], code))
                elif kind == "bond" and all(w in index for w in where):
                    arm.append((kind, tuple(index[w] for w in where), code))
            context = {
                "atoms": {where: code for kind, where, code in arm if kind == "atom"},
                "bonds": {where: code for kind, where, code in arm if kind == "bond"},
                "used": set(),
            }
            BRANCH_STEREO.set(context)
            _, _, parts = _select(unit, attach, stereo=arm)
            entries, complete = _stereo_entries(arm, parts[4], ring_parent=parts[5], used=context["used"])
            if not complete:
                return None
            per_arm.append((parts, entries))
    finally:
        BRANCH_STEREO.reset(token)
    if len({tuple(parts[:4]) for parts, _ in per_arm}) != 1:
        return None
    ordered = sorted(per_arm, key=lambda arm: tuple(0 if code in ("R", "Z", "r", "seqCis", "M", "m") else 1 for _, code in arm[1]))
    labels = [
        f"{locant}{chr(39) * i}{code}"
        for i, (_, entries) in enumerate(ordered)
        for locant, code in entries
    ]
    return per_arm[0][0], f"({','.join(labels)})-" if labels else ""


def _linker_descriptor(atoms, capacity=None):
    """The isotopic descriptor of the nuclides on the atoms of a multiplying linker, in the numbering that gives the
    modified atoms the lower locants; locants are omitted when every position of the linker is equivalent (P-82.5.2,
    P-82.6.1)."""
    from ._isotope_labels import descriptor, modification_key
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    labels = {a: context["labels"][a] for a in atoms if a in context["labels"]} if context else {}
    if not labels:
        return ""
    numberings = ({a: i + 1 for i, a in enumerate(order)} for order in (atoms, atoms[::-1]))
    position_of = min(numberings, key=lambda candidate: modification_key(labels, candidate))
    kinds = {(entry["skeleton"], tuple(sorted(entry["H"].items()))) for entry in labels.values()}
    unlocated = len(kinds) == 1 and (len(atoms) <= 2 or len(labels) == len(atoms))
    return descriptor(labels, position_of, unlocated, capacity=capacity)


def _consume_labels(atoms):
    from ._substituents import ISOTOPE_LABELS

    context = ISOTOPE_LABELS.get()
    if context:
        context["consumed"].update(a for a in atoms if a in context["labels"])


def _multiplied_units(mol, arms, roots, **cut):
    """(units, parts of the first unit) of identical multiplied arms, their nuclides cited inside the parts, or None
    when the arms differ in constitution or modification or have no name (P-82.2.2)."""
    from ._substituents import ISOTOPE_LABELS

    labelled = bool(ISOTOPE_LABELS.get())
    units = [_unit_molecule(mol, atoms, root, refill=labelled, **cut) for atoms, root in zip(arms, roots)]
    if labelled:
        keys = {_modified_unit_smiles(_unit_without_stereo(unit[0]), atoms) for unit, atoms in zip(units, arms)}
    else:
        keys = {_constitution_smiles(unit[0]) for unit in units}
    if len(keys) != 1:
        return None
    try:
        parts = _unit_parts(*units[0], arms[0]) if labelled else _select(*units[0])[2]
    except UnsupportedStructure:
        return None
    _consume_labels(set().union(*arms))
    return units, parts


def _unit_without_stereo(unit):
    bare = Chem.Mol(unit)
    Chem.RemoveStereochemistry(bare)
    return bare


def _bis_unit(text, tail):
    return "bis" + enclose(text + tail) if tail.startswith(" ") else "bis" + enclose(text) + tail


def _is_acyl_root(mol, atoms, attach, chalcogen_cut):
    """The carbonyl carbon of an acyl arm joined to its core through a carbon or a ring nitrogen: once cut it is a
    ketone, not an aldehyde (P-64.2.2.1, P-64.1.2.1(b))."""
    root = mol.GetAtomWithIdx(attach)
    cut = [n for n in root.GetNeighbors() if n.GetIdx() not in atoms]
    return (
        root.GetAtomicNum() == 6
        and not root.IsInRing()
        and not root.GetTotalNumHs()
        and len(cut) == 1
        and (
            cut[0].GetAtomicNum() == 6
            or (cut[0].GetAtomicNum() == 7 and cut[0].IsInRing())
            or (chalcogen_cut and cut[0].GetAtomicNum() in (8, 16, 34, 52))
        )
        and len(_double_oxygens(mol, attach)) == 1
        and root.GetDegree() == 3
    )


def _unit_molecule(mol, atoms, attach, refill=False, chalcogen_cut=False):
    editable = Chem.RWMol(mol)
    root = editable.GetAtomWithIdx(attach)
    root.SetAtomMapNum(1)
    cut = sum(1 for n in root.GetNeighbors() if n.GetIdx() not in atoms)
    if (root.GetIsAromatic() or root.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED or refill) and cut:
        root.SetNumExplicitHs(root.GetNumExplicitHs() + cut)
    for idx in sorted(set(range(mol.GetNumAtoms())) - atoms, reverse=True):
        editable.RemoveAtom(idx)
    unit = editable.GetMol()
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(unit)
    except Chem.rdchem.AtomValenceException:
        unit.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(unit)
    except Chem.rdchem.KekulizeException as error:
        raise UnsupportedStructure("the fragment cut from an aromatic ring cannot be kekulized") from error
    attach_idx = next(a.GetIdx() for a in unit.GetAtoms() if a.GetAtomMapNum() == 1)
    if _is_acyl_root(mol, atoms, attach, chalcogen_cut):
        unit.GetAtomWithIdx(attach_idx).SetBoolProp(ACYL_ROOT, True)
    canonical = Chem.Mol(unit)
    for a in canonical.GetAtoms():
        a.SetAtomMapNum(1 if a.GetIdx() == attach_idx else 0)
    return canonical, attach_idx


def _is_azide_part(atom):
    """The charged atoms of an azido group -N=N(+)=N(-)."""
    if atom.GetAtomicNum() != 7 or atom.GetDegree() > 2:
        return False
    if atom.GetFormalCharge() == 1:
        ends = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetFormalCharge() == -1 and n.GetDegree() == 1]
        return atom.GetDegree() == 2 and len(ends) == 1
    if atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (middle,) = atom.GetNeighbors()
        return middle.GetAtomicNum() == 7 and middle.GetFormalCharge() == 1 and _is_azide_part(middle)
    return False


def _is_azide_bond_end(first, second):
    """The bond joining a nitrogen to the charged middle nitrogen of an azido group."""
    return any(_is_azide_part(a) and a.GetFormalCharge() == 1 for a in (first, second))


def _is_isocyanate_carbon(atom):
    """The central carbon of -N=C=X (X = O, S, Se, Te), a prefix of its own rather than a carbonyl."""
    if atom.GetAtomicNum() != 6 or atom.GetDegree() != 2 or atom.GetFormalCharge():
        return False
    neighbors = sorted(atom.GetNeighbors(), key=lambda n: n.GetAtomicNum())
    nitrogen = neighbors[0]
    return (
        nitrogen.GetAtomicNum() == 7
        and not nitrogen.GetFormalCharge()
        and neighbors[1].GetAtomicNum() in (8, 16, 34, 52)
        and neighbors[1].GetDegree() == 1
        and all(b.GetBondTypeAsDouble() == 2.0 for b in atom.GetBonds())
    )


def _is_ring_azoxy_part(atom):
    """The oxide nitrogen of an azoxy group whose N=N bond lies in a ring, or its oxygen: the ring N-oxide of a heteroatom
    oxide (P-74.2.2.1.4), not a group cited beside another parent."""
    if atom.GetAtomicNum() == 8:
        return atom.GetDegree() == 1 and any(_is_ring_azoxy_part(n) for n in atom.GetNeighbors())
    found = azoxy_atoms(atom.GetOwningMol(), atom) if atom.GetAtomicNum() == 7 else None
    if found is None:
        return False
    bond = atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), found[1].GetIdx())
    return bond.IsInRing()


def _is_chain_amine_oxide_part(atom):
    """The charged atoms of an amine oxide on a non-ring nitrogen, cited as an '(oxo-λ5-azanyl)' prefix beside the
    ring N-oxide that is the parent (P-62.5, method 2)."""
    if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (nitrogen,) = atom.GetNeighbors()
        return nitrogen.GetAtomicNum() == 7 and _is_chain_amine_oxide_part(nitrogen)
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1 or atom.IsInRing() or atom.GetIsAromatic():
        return False
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 and n.GetDegree() == 1]
    return (
        len(oxygens) == 1
        and atom.GetDegree() + atom.GetTotalNumHs() == 4
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors() if n.GetIdx() != oxygens[0].GetIdx())
    )


def _is_nitro_part(atom):
    """The charged atoms of a nitro, azido, isocyano or azoxy group: N+ bonded to two oxygens (one O-)."""
    if _is_azide_part(atom) or is_halogen_oxo_part(atom.GetOwningMol(), atom):
        return True
    if atom.GetAtomicNum() == 7 and azoxy_atoms(atom.GetOwningMol(), atom) is not None:
        return True
    if atom.GetAtomicNum() in (6, 7) and atom.GetFormalCharge() in (-1, 1):
        mol = atom.GetOwningMol()
        if any(
            {n.GetAtomicNum(), atom.GetAtomicNum()} == {6, 7}
            and n.GetFormalCharge() == -atom.GetFormalCharge()
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 3.0
            and (atom if atom.GetAtomicNum() == 6 else n).GetDegree() == 1
            for n in atom.GetNeighbors()
        ):
            return True
    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 and atom.GetDegree() == 4:
        mol = atom.GetOwningMol()
        return azono_ligands(mol, atom.GetIdx(), [n.GetIdx() for n in atom.GetNeighbors()])
    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        return (
            len(oxygens) in (2, 3)
            and not any(o.IsInRing() for o in oxygens)
            and sum(o.GetFormalCharge() for o in oxygens) == -1
            and atom.GetDegree() == 3
            and not atom.GetTotalNumHs()
        )
    if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (n,) = atom.GetNeighbors()
        return n.GetAtomicNum() == 7 and _is_nitro_part(n)
    return False


_CHALCOGEN_SYMBOL = {8: "O", 16: "S", 34: "Se"}


def _peroxy_group(mol, atom, groups, ring_groups):
    if atom.GetAtomicNum() not in _CHALCOGEN_SYMBOL or atom.GetDegree() != 1:
        return
    (y,) = atom.GetNeighbors()
    if y.GetAtomicNum() not in _CHALCOGEN_SYMBOL or y.GetDegree() != 2 or y.IsInRing():
        return
    carbon = next((n for n in y.GetNeighbors() if n.GetIdx() != atom.GetIdx()), None)
    if carbon is None or carbon.GetAtomicNum() != 6:
        return
    oxo = [
        n.GetIdx()
        for n in carbon.GetNeighbors()
        if n.GetAtomicNum() in (8, 16, 34, 52)
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    owned = {atom.GetIdx(), y.GetIdx()}
    c = carbon.GetIdx()
    if carbon.IsInRing():
        if not oxo and ("peroxol", c, owned) not in ring_groups:
            ring_groups.append(("peroxol", c, owned))
        return
    if oxo:
        owned |= {oxo[0]}
        groups.setdefault("peroxoic", {})[c] = owned
        ring_groups.extend(("peroxoic", n.GetIdx(), {c} | owned) for n in carbon.GetNeighbors() if n.IsInRing())
    else:
        groups.setdefault("peroxol", {})[c] = owned


def _carbonyl_variant_group(mol, atom, groups, ring_groups):
    """Thio/seleno/telluro acid and imidic acid anions: a carbon with =Y and -X(-)."""
    if atom.GetAtomicNum() not in (8, 16, 34, 52) or atom.GetDegree() != 1:
        return
    (carbon,) = atom.GetNeighbors()
    if carbon.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return
    double = [
        n
        for n in carbon.GetNeighbors()
        if n.GetIdx() != atom.GetIdx()
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(double) != 1:
        return
    other = double[0]
    if other.GetAtomicNum() == 7 and other.GetTotalNumHs() == 1 and atom.GetAtomicNum() == 8:
        cls = "imidic"
    elif other.GetAtomicNum() in (8, 16, 34, 52) and not (atom.GetAtomicNum() == 8 and other.GetAtomicNum() == 8):
        cls = "thioic"
    else:
        return
    c = carbon.GetIdx()
    owned = {atom.GetIdx(), other.GetIdx()}
    groups.setdefault(cls, {})[c] = owned
    ring_groups.extend((cls, n.GetIdx(), {c} | owned) for n in carbon.GetNeighbors() if n.IsInRing())


def chalcogen_acid_variant(mol):
    """'thio' / 'dithio' / 'seleno' / ... for the marked thio-acid groups, or None."""
    variants = set()
    for atom in mol.GetAtoms():
        if not atom.HasProp(ANION_PROP) or atom.GetAtomicNum() not in (8, 16, 34, 52) or atom.GetDegree() != 1:
            continue
        (carbon,) = atom.GetNeighbors()
        others = [
            n
            for n in carbon.GetNeighbors()
            if n.GetIdx() != atom.GetIdx()
            and n.GetDegree() == 1
            and n.GetAtomicNum() in (8, 16, 34, 52)
            and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(others) != 1:
            continue
        pair = sorted((atom.GetAtomicNum(), others[0].GetAtomicNum()))
        if pair == [8, 8]:
            continue
        variants.add(tuple(pair))
    if len(variants) > 1:
        raise UnsupportedStructure("different chalcogen acid patterns in one name are not supported yet")
    if not variants:
        return None
    pair = variants.pop()
    words = {16: "thio", 34: "seleno", 52: "telluro"}
    chalcogens = [z for z in pair if z != 8]
    if len(chalcogens) == 2 and chalcogens[0] == chalcogens[1]:
        return "di" + words[chalcogens[0]]
    return "".join(words[z] for z in sorted(chalcogens, key=lambda z: {34: 0, 16: 1, 52: 2}[z]))


def peroxy_carbonyl_thio(mol):
    """'thio' / 'seleno' / 'telluro' when the carbonyl chalcogen of a peroxoic acid anion is not oxygen."""
    words = {16: "thio", 34: "seleno", 52: "telluro"}
    found = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in words and atom.GetDegree() == 1:
            (carbon,) = atom.GetNeighbors()
            if mol.GetBondBetweenAtoms(atom.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() == 2.0 and any(
                n.GetAtomicNum() in (8, 16, 34, 52) and n.GetDegree() == 2 and any(
                    m.GetDegree() == 1 and m.HasProp(ANION_PROP) for m in n.GetNeighbors()
                )
                for n in carbon.GetNeighbors()
            ):
                found.add(words[atom.GetAtomicNum()])
    return found.pop() if len(found) == 1 else ""


def peroxy_variant(mol):
    """('' | 'dithio' | 'diseleno' | '(XY-thio' ...) prefix of the 'peroxo' word for the marked groups."""
    variants = set()
    for atom in mol.GetAtoms():
        if not atom.HasProp(ANION_PROP) or atom.GetAtomicNum() not in _CHALCOGEN_SYMBOL or atom.GetDegree() != 1:
            continue
        (y,) = atom.GetNeighbors()
        if y.GetAtomicNum() not in _CHALCOGEN_SYMBOL or y.GetDegree() != 2:
            continue
        pair = (y.GetAtomicNum(), atom.GetAtomicNum())
        if pair[0] == pair[1]:
            variants.add({8: "", 16: "dithio", 34: "diseleno"}[pair[0]])
        else:
            thio = "thio" if 16 in pair else "seleno"
            variants.add(f"({_CHALCOGEN_SYMBOL[pair[0]]}{_CHALCOGEN_SYMBOL[pair[1]]}-{thio}")
    if len(variants) > 1:
        raise UnsupportedStructure("different peroxy chalcogen patterns in one name are not supported yet")
    return variants.pop() if variants else None


_CHAIN_GROUP_WORDS = {
    "acid": ("oic", " acid"),
    "sulfonic": ("sulfonic acid", ""),
    "alcohol": ("ol", ""),
    "thiol": ("thiol", ""),
    "selenol": ("selenol", ""),
    "tellurol": ("tellurol", ""),
    "amine": ("amine", ""),
    "peroxoic": ("peroxoic", " acid"),
    "peroxol": ("peroxol", ""),
    **{name: (name, "") for name in _PEROXOL_ANALOGUES},
    "thioic": ("thioic", " acid"),
    "imidic": ("imidic", " acid"),
}


def _chain_group_word(principal, count):
    if principal not in _CHAIN_GROUP_WORDS:
        raise UnsupportedStructure("an anionic carbon beside this group is not supported yet")
    word, tail = _CHAIN_GROUP_WORDS[principal]
    return multiplied_word(count, word), tail


def _compound_word(ide_locants, group_locants, group_word):
    ide_word = multiplied_word(len(ide_locants), "ide")
    if group_word[0] in "aeiouy":
        ide_word = ide_word[:-1]
    return f"{ide_word}-{','.join(str(x) for x in group_locants)}-{group_word}"


def _is_carboxylic_ester_carbon(mol, carbon):
    """-CO-OR or -CO-NR2 on a carbon chain, kept in the chain as 'alkoxy...oxo' / 'amino...oxo' (P-65.6.3.2.3)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetDegree() != 3 or atom.GetIsAromatic():
        return False
    oxo = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    if len(oxo) != 1 or mol.GetBondBetweenAtoms(carbon, oxo[0].GetIdx()).GetBondTypeAsDouble() != 2.0:
        return False
    ether = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 2]
    if len(ether) == 1:
        return all(n.GetAtomicNum() == 6 for n in ether[0].GetNeighbors() if n.GetIdx() != carbon)
    amine = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7 and not n.IsInRing() and not n.GetFormalCharge()]
    if len(amine) != 1 or not any(a.HasProp(ANION_PROP) for a in mol.GetAtoms()):
        return False
    return all(n.GetAtomicNum() == 6 and not n.IsInRing() or n.GetIdx() == carbon for n in amine[0].GetNeighbors())


def _acyloxy_lambda_centre_oxygen(oxygen, carbon):
    """The oxygen of an acyl-O-X group on a halogen or chalcogen of nonstandard bonding number: iodine and the other
    centres are not pseudoester elements, so the acyloxy group is a prefix of the lambda-n substituent (P-65.6.3.1.2)."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or oxygen.GetFormalCharge():
        return False
    far = next((n for n in oxygen.GetNeighbors() if n.GetIdx() != carbon), None)
    if far is None or far.GetAtomicNum() not in LAMBDA_CENTRE_STEMS or far.GetFormalCharge() or far.IsInRing():
        return False
    return far.GetTotalValence() > LAMBDA_CENTRE_STEMS[far.GetAtomicNum()][1]


_AMIDE_FAMILY = frozenset(
    {"amide", *_CHALCOGEN_AMIDE_CLASSES, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES, "hydrazide", *_CHALCOGEN_HYDRAZIDE.values(), *_CHALCOGEN_HYDRAZIDINE.values()}
)


def _nitrogen_acyl(mol, carbon):
    """An acyl carbon whose only hetero neighbours besides the carbonyl atom are nitrogens: an amide, hydrazide, urea or imide group,
    which is junior to an amide parent and so a prefix beside it (P-41)."""
    atom = mol.GetAtomWithIdx(carbon)
    return not any(
        (n.GetAtomicNum() not in (6, 7) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0)
        or (n.GetAtomicNum() == 7 and any(b.GetBondTypeAsDouble() != 1.0 for b in n.GetBonds()))
        for n in atom.GetNeighbors()
    )


def _ring_nitrogen_pair(mol, carbon):
    """Whether the carbonyl `carbon` is bonded to two ring nitrogens and nothing else: a pseudoketone with two hidden
    amide nitrogens (P-64.3.2, P-66.1.3)."""
    atom = mol.GetAtomWithIdx(carbon)
    others = [n for n in atom.GetNeighbors() if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0]
    return len(others) == 2 and all(n.GetAtomicNum() == 7 and _ring_nitrogen_acyl(mol, n, carbon) for n in others)


def _is_ester_like(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    if not _double_oxygens(mol, carbon):
        return False
    if any(_acyloxy_amine_oxygen(mol, n, carbon) for n in atom.GetNeighbors()):
        return False
    return any(
        n.GetAtomicNum() in (8, 7, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not _diacyl_chalcogen_chain(mol, carbon, n)
        and not _acyloxy_lambda_centre_oxygen(n, carbon)
        and not (
            n.GetAtomicNum() == 7
            and (atom.GetTotalNumHs() == 1 or any(c.GetAtomicNum() == 6 for c in atom.GetNeighbors()) or _ring_nitrogen_pair(mol, carbon))
            and _ring_nitrogen_acyl(mol, n, carbon)
        )
        and not (n.GetAtomicNum() == 8 and _terminal_heteroatom(mol, n.GetIdx(), 1))
        and not (n.GetAtomicNum() == 7 and _terminal_heteroatom(mol, n.GetIdx(), 2))
        and not (n.GetAtomicNum() == 7 and _acyl_diazene_nitrogen(mol, n, carbon))
        for n in atom.GetNeighbors()
    )


def _terminal_chalcogen_hydride(atom):
    """The sulfur, selenium or tellurium of an -SH, -SeH or -TeH group on nitrogen (thiohydroxylamine analogues, P-68.3.1.1.1.6)."""
    return atom.GetAtomicNum() in (16, 34, 52) and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1 and not atom.GetFormalCharge()


def _substituted_amine_nitrogen(mol, atom):
    if atom.HasProp(JUNIOR_GROUP):
        return False
    if FORCED_PRINCIPAL.get() in ("amine", "imine") and any(
        n.GetAtomicNum() == 6 and _is_ester_like(mol, n.GetIdx()) for n in atom.GetNeighbors()
    ):
        return False
    if any(n.GetAtomicNum() in MONONUCLEAR_HYDRIDES for n in atom.GetNeighbors()):
        return False
    if atom.GetAtomicNum() != 7 or (atom.GetFormalCharge() and not (AMINIUM.get() and atom.GetFormalCharge() == 1)) or atom.GetIsAromatic() or atom.IsInRing():
        return False
    if AMINIUM.get() is True and not atom.GetFormalCharge() and not atom.HasProp("_cationic_amine"):
        return False
    if any(
        n.GetAtomicNum() not in (6, 8, 16, 34, 52) and not _ring_nitrogen_parent(n) and not is_oxo_nitrogen(mol, n) and not _terminal_chalcogen_hydride(n)
        for n in atom.GetNeighbors()
    ):
        return False
    if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
        return False
    parents = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 or _ring_nitrogen_parent(n)]
    parent_ids = {n.GetIdx() for n in parents}
    if any(_ring_nitrogen_parent(n) for n in parents) and any(
        n.GetIdx() not in parent_ids and n.GetAtomicNum() not in (8, 16, 34, 52) for n in atom.GetNeighbors()
    ):
        return False
    return len(parents) >= 2 or (len(parents) == 1 and atom.GetDegree() >= 2)


def _ring_nitrogen_parent(atom):
    """A neutral ring nitrogen with all valences on single bonds, which bears the amine nitrogen as a substituent (P-62.2.1.3)."""
    return (
        atom.GetAtomicNum() == 7
        and atom.IsInRing()
        and not atom.GetFormalCharge()
        and not atom.GetTotalNumHs()
        and all(b.GetBondTypeAsDouble() == 1.0 or b.GetIsAromatic() for b in atom.GetBonds())
    )


def _cyano_on(mol, atom, nitrogen):
    """Whether `atom` is the carbon of a cyano group on `nitrogen`: cited as an N-prefix on a senior parent (P-66.1.6.2)."""
    return (
        atom.GetAtomicNum() == 6
        and atom.GetDegree() == 2
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() == 1.0
        and any(
            n.GetAtomicNum() == 7 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 3.0
            for n in atom.GetNeighbors()
        )
    )


def _carbamoyl_on(mol, atom, nitrogen):
    """Whether `atom` is the carbon of a carbamoyl or carbamothioyl group -C(=X)-NR2 (X = O, S, Se, Te) joined to `nitrogen`,
    cited as a prefix on the hydrazide nitrogen since the hydrazide outranks the carbonic acid derivative (P-41, P-65.2.1.5)."""
    if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3 or atom.GetFormalCharge() or atom.IsInRing():
        return False
    if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False
    ends = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    chalcogens = [n for n in ends if n.GetAtomicNum() in _CHALCOGEN_RANK]
    amino = [n for n in ends if n.GetAtomicNum() == 7]
    return (
        len(chalcogens) == 1
        and len(amino) == 1
        and chalcogens[0].GetDegree() == 1
        and not chalcogens[0].GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), chalcogens[0].GetIdx()).GetBondTypeAsDouble() == 2.0
        and not amino[0].GetFormalCharge()
        and not amino[0].IsInRing()
        and mol.GetBondBetweenAtoms(atom.GetIdx(), amino[0].GetIdx()).GetBondTypeAsDouble() == 1.0
        and _plain_nitrogen_substituents(mol, amino[0], atom.GetIdx())
    )


def _hydrazide_beta_nitrogen(mol, alpha, carbonyl):
    """The terminal nitrogen N' of an acyl-NR-NR'R'' group (P-66.3), whose substituents are carbon groups or a
    ylidene (a hydrazone), else None."""
    if alpha.GetFormalCharge() or alpha.IsInRing() or alpha.GetIsAromatic():
        return None
    betas = [n for n in alpha.GetNeighbors() if n.GetAtomicNum() == 7]
    if len(betas) != 1 or mol.GetBondBetweenAtoms(alpha.GetIdx(), betas[0].GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    beta = betas[0]
    if beta.GetFormalCharge() or beta.IsInRing() or beta.GetIsAromatic():
        return None
    for nitrogen, partner in ((alpha, beta), (beta, alpha)):
        orders = []
        for n in nitrogen.GetNeighbors():
            if n.GetIdx() in (partner.GetIdx(), carbonyl):
                continue
            if (nitrogen is beta and _carbamoyl_on(mol, n, nitrogen)) or _cyano_on(mol, n, nitrogen):
                continue
            if n.GetAtomicNum() != 6 or n.GetFormalCharge() or is_functional_carbon(mol, n.GetIdx()):
                return None
            orders.append(mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble())
        if (orders and 2.0 in orders and (nitrogen is alpha or len(orders) > 1)) or 3.0 in orders:
            return None
    return beta.GetIdx()


def _plain_nitrogen_arm(mol, nitrogen, partners, allow_ylidene):
    """Whether every other neighbour of `nitrogen` is a plain carbon group (a ylidene on a terminal nitrogen)."""
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() in partners:
            continue
        if n.GetAtomicNum() != 6 or n.GetFormalCharge() or is_functional_carbon(mol, n.GetIdx()):
            return False
        order = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble()
        if order != 1.0 and not (allow_ylidene and order == 2.0 and nitrogen.GetDegree() == 2):
            return False
    return True


def _carbonic_ether_oxygen(mol, oxygen, carbon):
    """Whether `oxygen` joins `carbon` to a second carbonic carbon -O-C(=N)-N, a dicarbonic derivative (P-65.2.3.1)."""
    ends = [n for n in oxygen.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
    return (
        oxygen.GetAtomicNum() == 8
        and oxygen.GetDegree() == 2
        and len(ends) == 1
        and ends[0].GetAtomicNum() == 6
        and sum(n.GetAtomicNum() == 7 for n in ends[0].GetNeighbors()) >= 1
        and any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(ends[0]).GetAtomicNum() in (7, 8, 16)
            for b in ends[0].GetBonds()
        )
    )


def _amidrazone_group(mol, atom):
    """("hydrazonamide" | "imidohydrazide", owned nitrogens) of a -C(-N)=N-N or -C(=N)-N-N group on a carbon with at
    most one carbon neighbour (P-66.4.2.1); the owned atoms are the three nitrogens."""
    carbon = atom.GetIdx()
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    others = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 7]
    ring_nitrogens = [n for n in nitrogens if n.IsInRing()]
    if len(nitrogens) == 3 and len(ring_nitrogens) == 1 and not others and not ring_nitrogens[0].GetFormalCharge():
        nitrogens = [n for n in nitrogens if not n.IsInRing()]
    forced = FORCED_PRINCIPAL.get() in _AMIDRAZONE
    if len(nitrogens) != 2 or len(others) > 1 or any(
        n.GetAtomicNum() != 6 and not (forced and _carbonic_ether_oxygen(mol, n, atom)) for n in others
    ):
        return None
    double = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    single = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0]
    if len(double) != 1 or len(single) != 1:
        return None
    imino, amino = double[0], single[0]
    if any(n.GetFormalCharge() or n.IsInRing() or n.GetIsAromatic() for n in (imino, amino)):
        return None
    imino_far = [m for m in imino.GetNeighbors() if m.GetAtomicNum() == 7 and m.GetIdx() != carbon]
    amino_far = [m for m in amino.GetNeighbors() if m.GetAtomicNum() == 7 and m.GetIdx() != carbon]
    if len(imino_far) == 1 and len(amino_far) == 1:
        hydrazone_end, hydrazide_end = imino_far[0], amino_far[0]
        if any(
            mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx()).GetBondTypeAsDouble() != 1.0
            for a, b in ((imino, hydrazone_end), (amino, hydrazide_end))
        ) or any(e.GetFormalCharge() or e.IsInRing() or e.GetIsAromatic() for e in (hydrazone_end, hydrazide_end)):
            return None
        if not (
            _plain_nitrogen_arm(mol, imino, {carbon, hydrazone_end.GetIdx()}, False)
            and _plain_nitrogen_arm(mol, amino, {carbon, hydrazide_end.GetIdx()}, False)
            and _plain_nitrogen_arm(mol, hydrazone_end, {imino.GetIdx()}, True)
            and _plain_nitrogen_arm(mol, hydrazide_end, {amino.GetIdx()}, True)
        ):
            return None
        if imino.GetDegree() != 2 or amino.GetDegree() > 3:
            return None
        return "hydrazonohydrazide", {imino.GetIdx(), amino.GetIdx(), hydrazone_end.GetIdx(), hydrazide_end.GetIdx()}
    if len(imino_far) == 1 and not amino_far:
        terminal, cls = imino_far[0], "hydrazonamide"
        arms = ((imino, {carbon, terminal.GetIdx()}, False), (amino, {carbon}, False))
    elif len(amino_far) == 1 and not imino_far:
        terminal, cls = amino_far[0], "imidohydrazide"
        arms = ((imino, {carbon}, False), (amino, {carbon, terminal.GetIdx()}, False))
    else:
        return None
    if terminal.GetFormalCharge() or terminal.IsInRing() or terminal.GetIsAromatic():
        return None
    if mol.GetBondBetweenAtoms(imino.GetIdx() if cls == "hydrazonamide" else amino.GetIdx(), terminal.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    center = imino if cls == "hydrazonamide" else amino
    if not all(_plain_nitrogen_arm(mol, n, partners, allow) for n, partners, allow in arms):
        return None
    if not _plain_nitrogen_arm(mol, terminal, {center.GetIdx()}, True):
        return None
    if imino.GetDegree() > (2 if cls == "hydrazonamide" else 2) or amino.GetDegree() > 3:
        return None
    return cls, {imino.GetIdx(), amino.GetIdx(), terminal.GetIdx()}


def _imidoyl_carbon(mol, carbon, nitrogen):
    """Whether `carbon` is the imidoyl carbon R-C(=NR')- of a diamidide (P-66.4.1.6), joined to `nitrogen`."""
    if carbon.IsInRing() or carbon.GetFormalCharge():
        return False
    others = [n for n in carbon.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    imino = [n for n in others if n.GetAtomicNum() == 7]
    rest = [n for n in others if n.GetAtomicNum() != 7]
    return (
        len(imino) == 1
        and len(rest) <= 1
        and all(n.GetAtomicNum() == 6 and not is_functional_carbon(mol, n.GetIdx()) for n in rest)
        and mol.GetBondBetweenAtoms(carbon.GetIdx(), imino[0].GetIdx()).GetBondTypeAsDouble() == 2.0
        and not imino[0].GetFormalCharge()
        and all(m.GetIdx() == carbon.GetIdx() or m.GetAtomicNum() == 6 for m in imino[0].GetNeighbors())
    )


def _cyano_on(mol, carbon, nitrogen):
    """Whether `carbon` is the carbon of a cyano group on `nitrogen`: a nitrile ranks below the amidine (P-41)."""
    return (
        carbon.GetAtomicNum() == 6
        and carbon.GetDegree() == 2
        and mol.GetBondBetweenAtoms(carbon.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() == 1.0
        and any(
            n.GetAtomicNum() == 7 and n.GetIdx() != nitrogen.GetIdx() and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 3.0
            for n in carbon.GetNeighbors()
        )
    )


def _amidine_nitrogens(mol, atom):
    """(imino N, amino N) of a -C(=NR)-NR2 group (P-66.4.1.1) whose carbon has at most one carbon neighbour; the
    imino nitrogen carries hydrogen, a carbon group or an oxygen (an amidoxime, P-66.4.4)."""
    carbon = atom.GetIdx()
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    others = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 7]
    ring_nitrogens = [n for n in nitrogens if n.IsInRing()]
    if len(nitrogens) == 3 and len(ring_nitrogens) == 1 and not others and not ring_nitrogens[0].GetFormalCharge():
        nitrogens = [n for n in nitrogens if not n.IsInRing()]
    if len(nitrogens) != 2 or len(others) > 1 or any(n.GetAtomicNum() != 6 for n in others):
        return None
    imino = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    amino = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0]
    if len(imino) != 1 or len(amino) != 1:
        return None
    imino, amino = imino[0], amino[0]
    if any(n.GetFormalCharge() or n.IsInRing() or n.GetIsAromatic() for n in (imino, amino)):
        return None
    allowed = (6, 8, 9, 17, 35, 53, *(z for z in MONONUCLEAR_HYDRIDES if z not in (7, 8, 16, 34, 52)))
    for nitrogen in (imino, amino):
        for n in nitrogen.GetNeighbors():
            if n.GetIdx() == carbon:
                continue
            if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
                return None
            if _ring_nitrogen_substituent(mol, n, nitrogen):
                continue
            if n.GetAtomicNum() not in allowed or n.GetFormalCharge():
                return None
            if is_functional_carbon(mol, n.GetIdx()) and not (
                _imidoyl_carbon(mol, n, nitrogen) or _carbamimidoyl_on(mol, n, nitrogen) or _cyano_on(mol, n, nitrogen)
            ):
                return None
    if imino.GetDegree() > 2 or amino.GetDegree() > 3:
        return None
    return imino.GetIdx(), amino.GetIdx()


def _acyclic_imine_nitrogen(mol, atom):
    """The =N- atom of a C=N group (P-62.3) whose carbon has only carbon neighbours, or a ring nitrogen that makes the
    group the imine of a pseudoketone (P-66.4.1.1), and whose nitrogen carries hydrogen, a carbon group or an oxygen
    (an oxime or its ether, P-68.3.1.1.2); None otherwise."""
    carbon = atom.GetIdx()
    imines = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if len(imines) != 1:
        return None
    nitrogen = imines[0]
    cationic = AMINIUM.get() == "imine" and nitrogen.GetFormalCharge() == 1
    if (nitrogen.GetFormalCharge() and not cationic) or nitrogen.IsInRing() or nitrogen.GetIsAromatic():
        return None
    if nitrogen.GetDegree() > (3 if cationic else 2):
        return None
    if any(
        n.GetAtomicNum() != 6
        and not _nitro_or_nitroso(mol, n)
        and not _ring_nitrogen_acyl(mol, n, carbon)
        and not _acyl_diazene_nitrogen(mol, n, carbon)
        for n in atom.GetNeighbors()
        if n.GetIdx() != nitrogen.GetIdx()
    ):
        return None
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == carbon:
            continue
        if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
        if n.GetFormalCharge() or (
            n.GetAtomicNum() not in ((6,) if cationic else (6, 8, 16, 34, 52))
            and not (_ring_nitrogen_parent(n) and not cationic)
        ):
            return None
    return nitrogen.GetIdx()


def _nitro_or_nitroso(mol, atom):
    """The nitrogen of a nitro or nitroso group, which makes an oxime a nitrolic or nitrosolic acid (P-68.3.1.1.3)."""
    return atom.GetAtomicNum() == 7 and (
        is_nitro_nitrogen(mol, atom.GetIdx())
        or (atom.GetDegree() == 2 and any(n.GetAtomicNum() == 8 and n.GetDegree() == 1 for n in atom.GetNeighbors()))
    )


def _ring_imine_nitrogen(mol, nitrogen, ring_atom):
    """A neutral =N-H, =N-R, =N-OH, =N-OR, =N-halogen or =N-(ring nitrogen) on a saturated ring carbon (P-62.3, P-68.3.1.1)."""
    if nitrogen.GetFormalCharge() or nitrogen.GetDegree() > 2:
        return False
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == ring_atom:
            continue
        if n.GetAtomicNum() == 8:
            if n.GetFormalCharge() or n.GetDegree() > 2:
                return False
        elif n.GetAtomicNum() in (16, 34, 52):
            if n.GetFormalCharge() or n.GetDegree() > 2:
                return False
        elif n.GetAtomicNum() not in (6, *HALOGEN_PREFIXES) and not _ring_nitrogen_parent(n):
            return False
    return True


def _hydroperoxy_oxygen(mol, oxygen, ring_atom):
    """The terminal OH oxygen of an -O-OH group whose first oxygen is bonded to `ring_atom`, else None."""
    if oxygen.GetDegree() != 2 or oxygen.GetFormalCharge():
        return None
    other = next((n for n in oxygen.GetNeighbors() if n.GetIdx() != ring_atom), None)
    if other is None or other.GetAtomicNum() != 8 or other.GetDegree() != 1 or other.GetTotalNumHs() != 1 or other.GetFormalCharge():
        return None
    return other.GetIdx()


def _ring_occurrences(mol):
    """[(class, ring_atom, owned atoms)] for every principal-capable group
    sitting directly on a ring atom (or, for a ketone, the ring carbonyl)."""
    found = []
    for atom in mol.GetAtoms():
        if not atom.IsInRing() or atom.HasProp(JUNIOR_GROUP):
            continue
        r = atom.GetIdx()
        senior_ring = atom.GetAtomicNum() in (16, 34, 52) and _ring_outranks_centre(mol, r)
        for n in atom.GetNeighbors():
            if n.IsInRing():
                continue
            z, i = n.GetAtomicNum(), n.GetIdx()
            order = mol.GetBondBetweenAtoms(r, i).GetBondTypeAsDouble()
            if senior_ring and order == 2.0 and z in (7, 8):
                continue
            if z == 8 and order == 2.0 and n.GetDegree() == 1:
                found.append(("ketone", r, {i}))
            elif z in _CHALCOGEN_KETONE_CLASS and order == 2.0 and n.GetDegree() == 1 and not n.GetFormalCharge():
                found.append((_CHALCOGEN_KETONE_CLASS[z], r, {i}))
            elif z == 8 and _terminal_heteroatom(mol, i, 1):
                found.append(("alcohol", r, {i}))
            elif order == 1.0 and (peroxol := _peroxol_group(mol, n, r)) is not None:
                found.append((peroxol[0], r, peroxol[1]))
            elif z in (16, 34, 52) and _terminal_heteroatom(mol, i, 1):
                found.append(({16: "thiol", 34: "selenol", 52: "tellurol"}[z], r, {i}))
            elif z in (16, 34, 52) and order == 1.0:
                sulfonyl = _sulfonyl_group(mol, i, r)
                if sulfonyl is not None:
                    found.append((sulfonyl[0], r, sulfonyl[1]))
            elif z == 7 and _terminal_heteroatom(mol, i, 2):
                if AMINIUM.get() is True and not n.HasProp("_cationic_amine"):
                    continue
                found.append(("amine", r, {i}))
            elif z == 7 and order == 2.0 and _ring_imine_nitrogen(mol, n, r):
                found.append(("imine", r, {i}))
            elif z == 6:
                group = _group_of(mol, i)
                if group is not None and _is_terminal(group[0]):
                    found.append((group[0], r, {i} | group[1]))
    return found


def _n_position(located, position_of):
    for atom in (located[1], *(located[3] if len(located) > 3 else ())):
        if atom in position_of:
            return position_of[atom]
    return None


def _prefix_locant_key(locant):
    """P-45.2.2: the N locants of the prefixes on the principal groups count with the numbered ones, letters before
    numerals; the number of the group that carries them follows the numbering of the chain and is not compared."""
    if isinstance(locant, str):
        return 0, 0, locant.count("'")
    return 1, locant, 0


def _with_n_names(grouped, n_names, position_of=None, group_count=2):
    """`grouped` plus the N-prefixes; a locant ("N", atom) of one of several groups reads N<position of atom>
    (P-66.1.1.4.2: 'N1,N5-dimethylpentanediamide'); groups on one atom are told apart by primes, the one with the
    most substituents unprimed (P-16.9.3)."""
    if not n_names:
        return grouped
    merged = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    groups = {}
    for name, compound, *locant in n_names:
        if locant and isinstance(locant[0], tuple) and len(locant[0]) > 2:
            entry = groups.setdefault(locant[0][2], [0, name, locant[0]])
            entry[0] += 1
            entry[1] = min(entry[1], name, key=lambda text: (alphanumerical_name_key(text), text))
    primes = {}
    if group_count > 1:
        by_position = {}
        for group, (count, first, located) in groups.items():
            by_position.setdefault(_n_position(located, position_of), []).append(
                (-count, alphanumerical_name_key(first), first, group)
            )
        for members in by_position.values():
            for prime, (*_, group) in enumerate(sorted(members)):
                primes[group] = prime
    stereo_context = BRANCH_STEREO.get()
    for name, compound, *locant in n_names:
        located = locant[0] if locant else "N"
        if isinstance(located, tuple):
            position = _n_position(located, position_of)
            if position is None:
                continue
            if stereo_context:
                stereo_context["used"] |= stereo_context.get("parked", {}).get(located, set())
            if group_count == 1:
                located = located[0]
            else:
                mark = "'" * primes.get(located[2], 0) if len(located) > 2 else ""
                if len(located) > 5:
                    mark *= located[5]
                elif len(located) > 4:
                    mark *= 2
                located = f"{located[0]}{mark}{'' if len(position_of) == 1 else position}"
        merged.setdefault(name, {"locants": [], "compound": compound})["locants"].append(located)
    return merged


def _n_group_positions(n_names, position_of):
    """Sorted numbers of the groups that carry N-prefixes, for choosing the numbering."""
    located = [entry[2] for entry in n_names if len(entry) > 2 and isinstance(entry[2], tuple)]
    return tuple(
        sorted((p, e[0].count("'")) for e in located if (p := _n_position(e, position_of)) is not None)
    )


def _ring_prefix_text(entries, locants, n_names, group_count=2):
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locants[r])
    return format_substituent_prefixes(_with_n_names(grouped, n_names, locants, group_count)) if (grouped or n_names) else ""


def _imidic_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls):
    """N-prefix names of the imidic or hydrazonic acid groups that carry substituents on their =N atom; with several
    groups each is located by the position of its attachment (P-16.9.2: 'N1,N4-dihydroxybutanediimidic acid')."""
    carbon_centre = spec_from_key(cls).center == "C"
    members = {}
    for carbon, owned in groups.get(cls, {}).items():
        centre = carbon if carbon_centre else next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)), carbon)
        members[centre] = (carbon, owned)
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            centre = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in ((6,) if carbon_centre else (6, 16, 34, 52)))
            members[centre] = (ring_atom, owned)
    entries = []
    for centre, (anchor, owned) in members.items():
        group = acid_group_at(mol, centre)
        if group is None or not group.n_substituents:
            continue
        located = ("N", anchor, centre, tuple(a for a in graph[centre] if a not in owned)) if len(members) > 1 else "N"
        for n in group.n_substituents:
            name, compound = _n_branch(located, graph, n, group.n_atom, halogens, aromatic_atoms, mol)
            entries.append((name, compound, located))
    return entries


def _n_branch(located, graph, root, nitrogen, halogens, aromatic_atoms, mol):
    """`name_branch` of an N-substituent of one of several groups; the stereo elements it cites count as cited only once
    that group proves to be the parent (`_with_n_names`), not while another group may still be."""
    context = BRANCH_STEREO.get()
    before = set(context["used"]) if context else set()
    named = name_branch(graph, root, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
    if context and isinstance(located, tuple):
        added = context["used"] - before
        context["used"] -= added
        context.setdefault("parked", {}).setdefault(located, set()).update(added)
    return named


def _imidamide_cited(mol, graph, halogens, aromatic_atoms, owned):
    """{nitrogen: [(name, compound)]} of the substituents on the nitrogens of a sulfonimidamide-type group."""
    centre = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52))
    nitrogens = [a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    return {
        a: [
            name_branch(graph, n, a, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            for n in graph[a]
            if n != centre and n not in nitrogens
        ]
        for a in nitrogens
    }


def _imidamide_roles(mol, graph, owned, cited, preferred=frozenset()):
    """[(nitrogen, N | N' | N'' ...)] of a sulfonimidamide-type group: N on the nitrogen bonded to the centre by a single
    bond, then primes on the terminal nitrogen of a hydrazide, the imido nitrogens (more substituents, then the
    alphanumerically first, then the `preferred` modified one) and the terminal nitrogens of hydrazono groups
    (P-66.4.1.1, P-14.3.5)."""
    centre = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52))
    nitrogens = [a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    bonded = [a for a in nitrogens if mol.GetBondBetweenAtoms(centre, a) is not None]
    amino = next(a for a in bonded if mol.GetBondBetweenAtoms(centre, a).GetBondTypeAsDouble() == 1.0)
    hydrazone_nitrogens = [a for a in bonded if a != amino and any(n in nitrogens for n in graph[a])]
    remote = [a for a in nitrogens if a not in bonded]

    def ranked(atoms):
        return sorted(
            atoms, key=lambda a: (-len(cited[a]), sorted(alpha_sort_key(name) for name, _ in cited[a]), a not in preferred, a)
        )

    hydrazide_end = [a for a in remote if amino in graph[a]]
    hydrazone_end = [a for a in remote if a not in hydrazide_end]
    imido = [a for a in bonded if a != amino and a not in hydrazone_nitrogens]
    sequence = [*hydrazide_end, *ranked(imido), *ranked(hydrazone_end)]
    return [(amino, "N"), *((a, "N" + chr(39) * k) for k, a in enumerate(sequence, 1))]


def _imidamide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls):
    """N, N', N'' ... prefix entries of a sulfonimidamide-, sulfonimidohydrazide- or sulfonohydrazonamide-type group:
    N on the nitrogen bonded to the centre by a single bond, then the primes on the terminal nitrogen of a hydrazide, on
    the imido nitrogens and on the terminal nitrogens of hydrazono groups (P-66.4.1.1, P-66.4.2.1); with several groups
    the locant carries the position of its group."""
    members = _group_members(mol, groups, ring_groups, cls)
    entries = []
    for carbon, identity, owned, anchor in members:
        cited = _imidamide_cited(mol, graph, halogens, aromatic_atoms, owned)
        rest = tuple(a for a in graph[carbon] if a not in owned)
        for nitrogen, locant in _imidamide_roles(mol, graph, owned, cited):
            located = (locant, anchor, identity, rest, "paired") if len(members) > 1 else locant
            entries.extend((name, compound, located) for name, compound in cited[nitrogen])
    return entries


def _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls="hydrazide"):
    """N and N' prefix entries (name, compound, locant) of the hydrazide groups: N on the acylated nitrogen, N' on
    the terminal nitrogen; with several groups the locant carries the position of its group (P-66.3.3.2)."""
    members = _group_members(mol, groups, ring_groups, cls)
    entries = []
    for carbon, identity, owned, anchor in members:
        acyl = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)), carbon)
        alpha = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and acyl in graph[a])
        beta = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and a != alpha)
        rest = tuple(a for a in graph[carbon] if a not in owned)
        for nitrogen, partner, base in ((alpha, beta, "N"), (beta, alpha, "N'")):
            located = (base, anchor, identity, rest, "paired") if len(members) > 1 else base
            for n in graph[nitrogen]:
                if n in (partner, acyl):
                    continue
                name, compound = _n_branch(located, graph, n, nitrogen, halogens, aromatic_atoms, mol)
                entries.append((name, compound, located))
    return entries


def _amidrazone_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls):
    """N-prefix entries of an amidrazone group (P-66.4.2.1): N on the amino nitrogen, N' on the terminal nitrogen of
    the hydrazine, and N'' on the imino nitrogen of an imidohydrazide."""
    members = {c: (owned, c) for c, owned in groups.get(cls, {}).items()}
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            members[next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)] = (owned, ring_atom)
    step = 2 if cls == "hydrazonamide" else 3
    entries = []
    for carbon_key, (owned, anchor) in members.items():
        nitrogens = [a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
        centre = next(
            c
            for c in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(c).GetAtomicNum() == 6 and sum(1 for n in graph[c] if n in nitrogens) == 2
        )
        order = {n: mol.GetBondBetweenAtoms(centre, n).GetBondTypeAsDouble() for n in nitrogens if n in graph[centre]}
        imino = next(n for n, o in order.items() if o == 2.0)
        amino = next(n for n, o in order.items() if o == 1.0)
        if cls == "hydrazonohydrazide":
            hydrazone_end = next(n for n in nitrogens if n in graph[imino] and n != centre)
            hydrazide_end = next(n for n in nitrogens if n in graph[amino] and n != centre)
            roles = ((amino, "N"), (hydrazide_end, "N'"), (hydrazone_end, "N''"))
        elif cls == "hydrazonamide":
            terminal = next(n for n in nitrogens if n not in (imino, amino))
            roles = ((amino, "N"), (terminal, "N'"))
        else:
            terminal = next(n for n in nitrogens if n not in (imino, amino))
            roles = ((amino, "N"), (terminal, "N'"), (imino, "N''"))
        rest = tuple(a for a in graph[centre] if a not in owned)
        for nitrogen, base in roles:
            located = (base, anchor, centre, rest, "amidrazone", step) if len(members) > 1 else base
            for n in graph[nitrogen]:
                if n in nitrogens or n == centre:
                    continue
                name, compound = _n_branch(located, graph, n, nitrogen, halogens, aromatic_atoms, mol)
                entries.append((name, compound, located))
    return entries


def _amidine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """N and N' prefix entries (name, compound, locant) of the amidine groups: N on the amino nitrogen, N' on the
    imino nitrogen; with several groups the locant carries the position of its group (P-66.4.1.4.1) and geminal groups
    continue the primes (P-66.4.1.4.2)."""
    members = {c: (owned, c) for c, owned in groups.get("amidine", {}).items()}
    for cls, ring_atom, owned in ring_groups:
        if cls == "amidine":
            key = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)
            members[key] = (owned, ring_atom)
    entries = []
    for carbon, (owned, anchor) in members.items():
        rest = tuple(a for a in graph[carbon] if a not in owned)
        for nitrogen in (a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7):
            atom = mol.GetAtomWithIdx(nitrogen)
            imino = any(
                b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds()
            )
            base = "N'" if imino else "N"
            located = (base, anchor, carbon, rest, "amidine") if len(members) > 1 else base
            for n in graph[nitrogen]:
                if n in members:
                    continue
                name, compound = _n_branch(located, graph, n, nitrogen, halogens, aromatic_atoms, mol)
                if compound and re.fullmatch(r"[a-z]+imidoyl", name) and _imidoyl_carbon(mol, mol.GetAtomWithIdx(n), atom):
                    compound = False
                entries.append((name, compound, located))
    return entries


def _imine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """N-prefix names when an imine group is substituted on its =N atom, else []; with several groups each entry
    carries the position of its group (N2,N3-dihydroxybutane-2,3-diimine)."""
    members = {c: (owned, c) for c, owned in groups.get("imine", {}).items()}
    for cls, ring_atom, owned in ring_groups:
        if cls == "imine":
            members[ring_atom] = (owned, ring_atom)
    substituted = []
    for carbon, (owned, anchor) in members.items():
        nitrogen = next(iter(owned))
        subs = [n for n in graph[nitrogen] if mol.GetBondBetweenAtoms(nitrogen, n).GetBondTypeAsDouble() == 1.0]
        if subs:
            substituted.append((carbon, anchor, owned, nitrogen, subs))
    if not substituted:
        return []
    if len(members) == 1:
        _, _, _, nitrogen, subs = substituted[0]
        return [name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in subs]
    entries = []
    for carbon, anchor, owned, nitrogen, subs in substituted:
        for n in subs:
            if any(
                next(iter(other_owned)) in _arm_atoms(graph, n, nitrogen)
                for other_carbon, (other_owned, _) in members.items()
                if other_carbon != carbon
            ):
                raise UnsupportedStructure("imine groups joined through nitrogen need a multiplicative name (P-15.3)")
        # an imine has no carbo suffix: a group beside the parent is a prefix, so its neighbours never locate its N
        located = ("N", anchor, carbon, ())
        entries += [
            (*name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True), located)
            for n in subs
        ]
    return entries


def _group_members(mol, groups, ring_groups, cls):
    """[(key atom, group identity, owned atoms, anchor)] of the groups of class `cls`. Sulfonyl-type groups on one carbon
    are merged in `groups` and several groups on one ring atom share their key, so each is told apart by its own atoms."""
    def hetero(owned):
        return frozenset(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() != 6)

    members = {}
    for carbon, owned in groups.get(cls, {}).items():
        split = [part for kind, part in _sulfonyl_groups_on(mol, carbon) if kind == cls]
        for part in split if len(split) > 1 else [owned]:
            members[carbon, hetero(part)] = (part, carbon)
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            key = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)
            members[key, hetero(owned)] = (owned, ring_atom)
    shared = Counter(key for key, _ in members)
    return [(key, key if shared[key] == 1 else min(owned), owned, anchor) for (key, _), (owned, anchor) in members.items()]


def _amide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls="amide"):
    """N-prefix names when the single amide group is N-substituted, else []."""
    members = _group_members(mol, groups, ring_groups, cls)
    entries = []
    for carbon, identity, owned, anchor in members:
        nitrogen = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
        center = next((a for a in graph[nitrogen] if a in owned), carbon)
        located = ("N", anchor, identity, tuple(a for a in graph[carbon] if a not in owned)) if len(members) > 1 else "N"
        for n in graph[nitrogen]:
            if n != center:
                name, compound = _n_branch(located, graph, n, nitrogen, halogens, aromatic_atoms, mol)
                entries.append((name, compound, located))
    return entries


_BONDING_LEVELS = (7, 6, 5, 4, 3)


def _bonding_rank(mol, graph, parent_atoms, position_of):
    """P-45.3: among equally ranked parents, the one with more substituents of the higher bonding number attached
    directly to it, then with lower locants for them."""
    found = [
        (nonstandard_bonding(mol.GetAtomWithIdx(n)), position_of[a])
        for a in parent_atoms
        for n in graph[a]
        if n not in parent_atoms and nonstandard_bonding(mol.GetAtomWithIdx(n))
    ]
    return (
        tuple(-sum(1 for b, _ in found if b == level) for level in _BONDING_LEVELS),
        tuple(tuple(sorted(p for b, p in found if b == level)) for level in _BONDING_LEVELS),
    )


def _evaluate(
    mol, graph, halogens, aromatic_atoms, chain, principal, principal_atoms, owned, attach=None, n_names=(), stereo=None
):
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    chain_set = set(chain)
    if attach is not None and attach not in chain_set:
        return ((1,), "", None)
    on_chain = [a for a in principal_atoms if a in chain_set]
    owned = set().union(*(principal_atoms[a] for a in on_chain)) if on_chain else owned
    if principal == "ide":
        on_chain = [a for a in on_chain for _ in range(anion_weight(mol.GetAtomWithIdx(a)))]
    else:
        on_chain = [a for a in on_chain for _ in range(_group_weight(mol, a, principal))]
    if not on_chain or (_is_terminal(principal) and any(position_of[a] not in (1, len(chain)) for a in on_chain)):
        return ((1,), "", None)
    ene, yne = [], []
    for a, b in zip(chain, chain[1:]):
        order = mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        if order == 2.0:
            ene.append(position_of[a])
        elif order == 3.0:
            yne.append(position_of[a])
    entries = {}
    for atom in chain:
        for neighbor in graph[atom]:
            if neighbor in chain_set or neighbor in owned:
                continue
            name, compound = name_branch(graph, neighbor, atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            entries.setdefault(position_of[atom], []).append((name, compound))
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    suffix_locants = sorted(position_of[a] for a in on_chain)
    count = len(on_chain)
    ide_locants = sorted(position_of[a] for a, w in IDE_EXTRA.get().items() if a in chain_set for _ in range(w))
    length = len(chain)
    force = (attach is not None and length != 1) or (principal == "ide" and bool(grouped) and length > 1) or FORCE_LOCANTS.get()
    completely_substituted = (
        length == 2
        and _is_terminal(principal)
        and len(on_chain) == 1
        and not ene
        and not yne
        and not (
            principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "hydrazide", "amidine", *_AMIDRAZONE)
            and any(mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and mol.GetAtomWithIdx(a).GetTotalNumHs() for a in owned)
        )
        and (
            mol.GetAtomWithIdx(next(a for a in chain if a != on_chain[0])).GetTotalNumHs() == 0
            # P-14.3.4.3: the hydrogen of a formyl group and the carbon of a cyano group are not substitutable
            or principal in ("nitrile", "aldehyde")
        )
    )
    with_n = _with_n_names(grouped, n_names, position_of, len(on_chain))
    nitrogens = {a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7}
    elsewhere = set().union(*(principal_atoms[a] | {a} for a in principal_atoms if a not in chain_set))
    # a parent whose N-substituent holds another principal group yields to the one citing that group as an amido prefix
    holds_group = any(
        _arm_atoms(graph, neighbor, n) & elsewhere
        for n in nitrogens
        for neighbor in graph[n]
        if neighbor not in owned and neighbor not in chain_set
    )
    prefix_locants = [loc for info in with_n.values() for loc in info["locants"]]
    uncited_locants = length > 1 and not force and omits_all_locants(mol, chain_set | nitrogens, with_n, owned - nitrogens)
    # the printed PINs 3-amino-2,3-dioxopropanethioic acid and 4-bromo-3,4-dioxobutanoic acid cite locants on a carbon chain
    carbon_chain_beyond_two = length > 2 and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in chain_set)
    single_arrangement = (
        length > 1 and not force and not ene and not yne and not carbon_chain_beyond_two
        and has_single_arrangement(mol, chain_set | nitrogens, formyl_oxygens=owned if principal == "aldehyde" else frozenset())
    )
    uncited_locants = uncited_locants or single_arrangement
    if completely_substituted and not force and not n_names:
        prefix = single_site_prefixes(grouped)
    else:
        keeps_nitrogen_hydrogen = principal == "amine" and count == 1 and any(mol.GetAtomWithIdx(a).GetTotalNumHs() for a in nitrogens)
        prefix = format_substituent_prefixes(with_n, omit_locants=length == 1 and not force and not n_names and not keeps_nitrogen_hydrogen)
    tail = ""
    if principal == "ide" and attach is None and length == 2 and not grouped and (ene or yne) and (count == 1 or yne):
        word = multiplied_word(count, "ide")
        stem = "ethyn" if yne else "ethen"
        body = stem + word if word[0] in "aeiouy" else stem + "e" + word
    elif ide_locants:
        word, tail = _chain_group_word(principal, count)
        body = name_from_substituents(
            length, ene, yne, _compound_word(ide_locants, suffix_locants, word), ide_locants, force_own_locant=True
        )
    elif _is_variant(principal):
        spec = spec_from_key(principal)
        retained = None
        if attach is None and spec.anion and spec.oxo == ("O",) and spec.y == ("O",) and spec.center == "C":
            retained = retained_chain_acid(grouped, length, ene, yne, count, "anion")
        if retained is None and attach is None and spec.anion and spec.center == "C" and length == 1 and set(grouped) == {"carboxy"}:
            retained = "hydrogen oxalate"
        if retained is not None:
            (prefix, body), tail = split_retained_prefix(grouped, retained), ""
        elif spec.center == "C":
            word, tail = chain_suffix(spec, count)
            body = name_from_substituents(length, ene, yne, word)
        else:
            body = name_from_substituents(
                length, ene, yne, carbo_suffix(spec, count), suffix_locants, force_own_locant=force,
                substituted=bool(grouped),
            )
    elif principal == "acid":
        retained = None if attach is not None else retained_chain_acid(grouped, length, ene, yne, count, "acid")
        if retained is not None:
            (prefix, body), tail = split_retained_prefix(grouped, retained), ""
        else:
            body, tail = name_from_substituents(length, ene, yne, multiplied_word(count, "oic")), " acid"
            if attach is not None and length == 2 and count == 1 and not ene and not yne:
                body = "acetic"
    elif principal in ("peroxoic", "thioic", "imidic"):
        body, tail = name_from_substituents(length, ene, yne, multiplied_word(count, principal)), " acid"
    elif principal == "peroxol" or principal in _PEROXOL_ANALOGUES:
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, principal), suffix_locants, substituted=bool(grouped)
        )
    elif principal == "amide" or principal in _CHALCOGEN_AMIDE_CLASSES:
        body = name_from_substituents(length, ene, yne, multiplied_word(count, principal))
    elif principal == "amidine":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "imidamide"))
    elif principal in _AMIDRAZONE:
        body = name_from_substituents(length, ene, yne, multiplied_word(count, principal))
    elif principal == "hydrazide":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "hydrazide"))
    elif principal == "nitrile":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "nitrile"))
    elif principal == "aldehyde":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "al"))
    else:
        word = CHAIN_SUFFIX_WORD[principal]
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, word), None if single_arrangement else suffix_locants,
            force_own_locant=force, substituted=bool(grouped) and not single_arrangement,
        )
    if prefix and uncited_locants and not any(ch.isdigit() for ch in body):
        prefix = format_substituent_prefixes(with_n, omit_all=True)
    name = prefix + body + tail
    attach_locant = position_of[attach] if attach is not None else 0
    reported_attach = None if (attach is not None and length == 1) else attach_locant
    key = (
        -(count + len(ide_locants)),
        -len(ide_locants),
        -count,
        tuple(ide_locants),
        -length,
        -(len(ene) + len(yne)),
        -len(ene),
        tuple(suffix_locants),
        lowest_locant_set(ene + yne),
        lowest_locant_set(ene),
        attach_locant,
        _parent_configuration_rank(stereo, position_of, directional=False),
        holds_group,
        -len(prefix_locants),
        tuple(sorted(_prefix_locant_key(loc) for loc in prefix_locants)),
        tuple(
            _prefix_locant_key(loc)
            for cited in sorted(with_n, key=citation_order_key)
            for loc in sorted(with_n[cited]["locants"], key=_prefix_locant_key)
        ),
        _bonding_rank(mol, graph, chain_set, position_of),
        alphanumerical_name_key(name),
        _n_group_positions(n_names, position_of),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (prefix, body, tail, reported_attach, position_of, False)


CHAIN_SUFFIX_WORD = {
    "ide": "ide",
    "ketone": "one",
    "thione": "thione",
    "selone": "selone",
    "tellone": "tellone",
    "alcohol": "ol",
    "thiol": "thiol",
    "selenol": "selenol",
    "tellurol": "tellurol",
    "amine": "amine",
    "imine": "imine",
    "sulfonic": "sulfonic acid",
    "sulfonamide": "sulfonamide",
    **{name: name for name in _CHALCOGEN_SULFONAMIDE_CLASSES},
    **{name: name for name in _CHALCOGEN_HYDRAZIDE.values()},
    **{name: name for name in _CHALCOGEN_HYDRAZIDINE.values()},
    **{name: name for name in _CHALCOGEN_IMIDAMIDE.values()},
}


def _stereo_free(mol):
    return not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) and not any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
