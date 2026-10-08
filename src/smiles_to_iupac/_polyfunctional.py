"""Acyclic chain parents carrying any mix of supported groups (P-41, P-44.1.1,
P-44.3, P-45): the most senior class present becomes the suffix and every
other group, ring or branch is cited as a substituent prefix through
`name_branch` (P-29.3.3, P-29.4). Rings may only be substituents; a principal
group on a ring is not handled here.
"""

import contextvars
import itertools
from types import MappingProxyType
import re

from rdkit import Chem, rdBase
from rdkit.Chem import CanonicalRankAtoms, rdCIPLabeler

from ._common import (
    alphanumerical_name_key,
    assembly_join,
    alpha_sort_key,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    nonstandard_bonding,
    name_from_substituents,
    ring_cycle,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._anion import ANION_PROP, anion_weight
from ._numerals import multiplying_prefix
from ._chalcogenourea import is_oxo_nitrogen
from ._functional_prefixes import is_nitro_nitrogen
from ._hetero_prefixes import (
    CATION_PARENT,
    EXTENDED_PREFIXES,
    LAMBDA_CENTRE_STEMS,
    MONONUCLEAR_HYDRIDES,
    is_functional_carbon,
    is_halogen_oxo_part,
)
from ._multiplicative import _bare_key
from ._multiplicative_text import PrimedLocant, enclose, unit_phrase
from ._multiplicative_ring import (
    NO_RETAINED_BENZENE,
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
    parent_text,
    spec_of,
)
from ._fused_numbering import HETERO_RANK as _HETERO_RANK
from ._ring_diyl_numbering import _exocyclic_oxo, is_hydro_fusion_system
from ._ring_system_seniority import ring_seniority_key
from ._acid_groups import acid_group_at
from ._acid_lexicon import carbo_suffix, chain_suffix, make_spec, rank_key, spec_from_key
from ._retained_acids import retained_chain_acid, single_site_prefixes
from ._locant_omission import omits_all_locants
from ._substituents import format_substituent_prefixes, name_branch

_CHALCOGEN_HYDRAZIDE = {
    (16, 2): "sulfonohydrazide",
    (16, 1): "sulfinohydrazide",
    (34, 2): "selenonohydrazide",
    (34, 1): "seleninohydrazide",
    (52, 2): "telluronohydrazide",
    (52, 1): "tellurinohydrazide",
}
_CHALCOGEN_HYDRAZIDINE = {
    (16, 1): "sulfonohydrazonohydrazide",
    (16, 0): "sulfinohydrazonohydrazide",
    (34, 1): "selenonohydrazonohydrazide",
    (34, 0): "seleninohydrazonohydrazide",
    (52, 1): "telluronohydrazonohydrazide",
    (52, 0): "tellurinohydrazonohydrazide",
    (16, 1, "amide"): "sulfonohydrazonamide",
    (16, 0, "amide"): "sulfinohydrazonamide",
    (34, 1, "amide"): "selenonohydrazonamide",
    (34, 0, "amide"): "seleninohydrazonamide",
    (52, 1, "amide"): "telluronohydrazonamide",
    (52, 0, "amide"): "tellurinohydrazonamide",
}
_AMIDRAZONE = ("hydrazonamide", "imidohydrazide", "hydrazonohydrazide")
_CHALCOGEN_IMIDAMIDE = {}
for _z, (_on, _in) in ((16, ("sulfon", "sulfin")), (34, ("selenon", "selenin")), (52, ("telluron", "tellurin"))):
    _CHALCOGEN_IMIDAMIDE[(_z, 1, 1)] = f"{_on}imidamide"
    _CHALCOGEN_IMIDAMIDE[(_z, 0, 2)] = f"{_on}odiimidamide"
    _CHALCOGEN_IMIDAMIDE[(_z, 0, 1)] = f"{_in}imidamide"
_CHALCOGEN_AMIDE = {16: "thioamide", 34: "selenoamide", 52: "telluroamide"}
_CHALCOGEN_AMIDE_CLASSES = tuple(_CHALCOGEN_AMIDE.values())
_REPLACED_INFIX = {16: "thio", 34: "seleno", 52: "telluro"}


def _chalcogen_sulfonamide_name(slots, replaced):
    """'sulfonothioamide', 'sulfonodithioamide', 'sulfinothioamide', ...: the amide of a sulfonic (two =X) or sulfinic
    (one =X) acid whose =O are replaced by S, Se or Te, infixes in alphanumerical order (P-66.1.4.1.1, Table 4.4)."""
    counts = {}
    for z in replaced:
        counts[_REPLACED_INFIX[z]] = counts.get(_REPLACED_INFIX[z], 0) + 1
    infix = "".join(multiplied_word(counts[word], word) for word in sorted(counts))
    return ("sulfon" if slots == 2 else "sulfin") + "o" + infix + "amide"


_CHALCOGEN_SULFONAMIDE = {}
for _slots, _size in ((2, 1), (2, 2)):
    for _combo in itertools.combinations_with_replacement((16, 34, 52), _size):
        _CHALCOGEN_SULFONAMIDE[(_slots, _combo)] = _chalcogen_sulfonamide_name(_slots, _combo)
_CHALCOGEN_SULFONAMIDE[(1, ())] = "sulfinamide"
for _combo in itertools.combinations_with_replacement((16, 34, 52), 1):
    _CHALCOGEN_SULFONAMIDE[(1, _combo)] = _chalcogen_sulfonamide_name(1, _combo)
_ACID_AMIDE = {}
for _z, _stem in ((34, "selen"), (52, "tellur")):
    for _slots, _ending in ((2, "on"), (1, "in")):
        _ACID_AMIDE[(_z, _slots)] = _CHALCOGEN_SULFONAMIDE[(_z, _slots)] = f"{_stem}{_ending}amide"
_ACID_AMIDE[(16, 1)] = "sulfinamide"
_CHALCOGEN_SULFONAMIDE_CLASSES = tuple(_CHALCOGEN_SULFONAMIDE.values())
_SENIORITY = [
    "ide", "acid", "thioic", "peroxoic", "imidic", "sulfonic", "amide", *_CHALCOGEN_AMIDE_CLASSES, "amidine", *_AMIDRAZONE, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES, *_CHALCOGEN_IMIDAMIDE.values(), "hydrazide", *_CHALCOGEN_HYDRAZIDE.values(), *_CHALCOGEN_HYDRAZIDINE.values(), "nitrile", "aldehyde", "ketone", "thione", "selone", "tellone", "alcohol", "peroxol",
    "thiol", "selenol", "tellurol", "amine", "imine",
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
FORCE_LOCANTS = contextvars.ContextVar("force_locants", default=False)
LAST_POSITIONS = contextvars.ContextVar("last_positions", default=None)
FORCED_PRINCIPAL = contextvars.ContextVar("forced_principal", default=None)
SUBSTITUTED_AMINE_PREFIX = contextvars.ContextVar("substituted_amine_prefix", default=False)


def _principal_class(classes):
    if FORCED_PRINCIPAL.get():
        return FORCED_PRINCIPAL.get() if FORCED_PRINCIPAL.get() in classes else None
    if RING_CENTER.get():
        return None
    if AMINIUM.get():
        parent = AMINIUM.get() if AMINIUM.get() in ("imine", "amide", "nitrile") else "amine"
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


def _urea_carbon(mol, idx):
    """The carbonyl carbon of a urea core C(=O)(N)N, the amide of carbonic acid: it ranks below every amide of a
    carboxylic or sulfonic acid (P-66.1.6.1.1.5) and is cited as 'carbamoyl' beneath one."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing() or atom.GetDegree() != 3 or not _double_oxygens(mol, idx):
        return False
    return len(_single_neighbors(mol, idx, 7)) == 2


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


def _organyloxy(mol, oxygen, nitrogen):
    """An -O-R group on nitrogen whose R is a plain carbon group, cited as an alkoxy or aryloxy prefix (P-63.2.2.1.1)."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or oxygen.GetFormalCharge():
        return False
    (carbon,) = [n for n in oxygen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    return (
        carbon.GetAtomicNum() == 6
        and mol.GetBondBetweenAtoms(oxygen.GetIdx(), carbon.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not _double_oxygens(mol, carbon.GetIdx())
        and not is_functional_carbon(mol, carbon.GetIdx())
    )


def _carbamimidoyl_on(mol, atom, nitrogen):
    """Whether `atom` is the carbon of an unsubstituted carbamimidoyl group -C(=NH)NH2 bonded to `nitrogen`."""
    if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3 or atom.GetFormalCharge() or atom.IsInRing():
        return False
    if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return False
    ends = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
    orders = sorted(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in ends)
    return (
        orders == [1.0, 2.0]
        and all(n.GetAtomicNum() == 7 and n.GetDegree() == 1 and not n.GetFormalCharge() for n in ends)
        and sum(n.GetTotalNumHs() for n in ends) == 3
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
        or _carbamimidoyl_on(mol, n, nitrogen)
        or n.GetAtomicNum() in HALOGEN_PREFIXES
        or (
            n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and not _double_oxygens(mol, n.GetIdx())
            and not is_functional_carbon(mol, n.GetIdx())
        )
        or _hydride_group_on_nitrogen(mol, n, nitrogen, len(others))
        for n in others
    )


def _hydride_group_on_nitrogen(mol, atom, nitrogen, sibling_count):
    """An N-substituent that is a mononuclear hydride group (phosphanyl, silyl, boranyl ...) or, on a nitrogen with no
    other substituent, an ylidene of a hydride atom or of a plain carbon (P-66.1.1.4.3, P-68.1)."""
    from ._hetero_prefixes import MONONUCLEAR_HYDRIDES

    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
        return False
    order = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), atom.GetIdx()).GetBondTypeAsDouble()
    if order == 1.0:
        return atom.GetAtomicNum() in MONONUCLEAR_HYDRIDES and atom.GetAtomicNum() != 7
    if order == 2.0 and sibling_count == 1:
        if atom.GetAtomicNum() == 6:
            return not _double_oxygens(mol, atom.GetIdx()) and not is_functional_carbon(mol, atom.GetIdx())
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
            n.GetAtomicNum() == 6 and not is_functional_carbon(mol, n.GetIdx())
        ):
            return False
    return True


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
    hydrazones = [
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 7
        and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        and n.GetDegree() == 2
        and any(m.GetAtomicNum() == 7 and _terminal_heteroatom(mol, m.GetIdx(), 2) for m in n.GetNeighbors())
    ]
    rest = [
        n for n in sulfur.GetNeighbors() if n.GetIdx() not in oxygens + imides + replaced + hydrazones and n.GetIdx() != attached
    ]
    if len(rest) != 1 or sulfur.GetDegree() != len(oxygens) + len(imides) + len(replaced) + len(hydrazones) + 2:
        return None
    other = rest[0]
    if (
        len(hydrazones) == 1
        and not imides
        and not replaced
        and other.GetAtomicNum() == 7
        and _terminal_heteroatom(mol, other.GetIdx(), 2)
        and (name := _CHALCOGEN_HYDRAZIDINE.get((sulfur.GetAtomicNum(), len(oxygens), "amide"))) is not None
    ):
        terminal = next(m.GetIdx() for m in mol.GetAtomWithIdx(hydrazones[0]).GetNeighbors() if m.GetIdx() != s_idx)
        return name, {s_idx, *oxygens, hydrazones[0], terminal, other.GetIdx()}
    if len(hydrazones) == 1 and not imides and not replaced and other.GetAtomicNum() == 7 and (
        (name := _CHALCOGEN_HYDRAZIDINE.get((sulfur.GetAtomicNum(), len(oxygens)))) is not None
    ):
        beta = _hydrazide_beta_nitrogen(mol, other, s_idx)
        if beta is not None and _terminal_heteroatom(mol, beta, 2) and other.GetTotalNumHs() == 1:
            terminal = next(m.GetIdx() for m in mol.GetAtomWithIdx(hydrazones[0]).GetNeighbors() if m.GetIdx() != s_idx)
            return name, {s_idx, *oxygens, hydrazones[0], terminal, other.GetIdx(), beta}
        return None
    if hydrazones:
        return None
    if sulfur.GetAtomicNum() == 16 and replaced and not imides and len(oxygens) + len(replaced) in (1, 2):
        if other.GetAtomicNum() == 7 and (
            _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx)
        ):
            zs = tuple(sorted(mol.GetAtomWithIdx(i).GetAtomicNum() for i in replaced))
            name = _CHALCOGEN_SULFONAMIDE[(len(oxygens) + len(replaced), zs)]
            return name, {s_idx, *oxygens, *replaced, other.GetIdx()}
        return None
    plain_amide = _ACID_AMIDE.get((sulfur.GetAtomicNum(), len(oxygens)))
    if plain_amide is not None and not imides and not replaced and other.GetAtomicNum() == 7 and (
        _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx)
    ):
        return plain_amide, {s_idx, *oxygens, other.GetIdx()}
    if sulfur.GetAtomicNum() == 16 and len(oxygens) == 2 and not imides:
        if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
            return "sulfonic", {s_idx, *oxygens, other.GetIdx()}
        if other.GetAtomicNum() == 7 and (
            _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx)
        ):
            return "sulfonamide", {s_idx, *oxygens, other.GetIdx()}
    imidamide = _CHALCOGEN_IMIDAMIDE.get((sulfur.GetAtomicNum(), len(oxygens), len(imides)))
    if imidamide is not None and other.GetAtomicNum() == 7 and (
        _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx)
    ):
        return imidamide, {s_idx, *oxygens, *imides, other.GetIdx()}
    hydrazide = _CHALCOGEN_HYDRAZIDE.get((sulfur.GetAtomicNum(), len(oxygens)))
    if hydrazide is not None and not imides and other.GetAtomicNum() == 7:
        beta = _hydrazide_beta_nitrogen(mol, other, s_idx)
        if beta is not None:
            return hydrazide, {s_idx, *oxygens, other.GetIdx(), beta}
    return None


def _chalcogen_amide_group(mol, atom, chalcogen, z):
    """(class, owned atoms) of a C(=X)-N carbon, X = S, Se or Te: a thio-, seleno- or telluroamide (P-66.1.4.1), or the
    thione of a 'hidden' amide whose nitrogen is a ring member (P-66.1.4.3)."""
    carbon = atom.GetIdx()
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != chalcogen]
    nitrogens = [n for n in others if n.GetAtomicNum() == 7]
    rest = [n for n in others if n.GetAtomicNum() != 7]
    if len(nitrogens) != 1 or any(mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in others):
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
    """Whether the acyl `carbon` is joined through `first` to another acyl carbon by a chain of three or more
    chalcogen atoms (P-65.7.5.1: a diacyl trioxidane or tetrasulfane is a pseudoketone, not an ester)."""
    if _homogeneous_chalcogen_chain(mol, carbon, first):
        return True
    chain, previous, atom = [], carbon, first
    while atom.GetAtomicNum() in (8, 16, 34, 52):
        if atom.GetFormalCharge() or atom.GetDegree() != 2 or atom.IsInRing():
            return False
        chain.append(atom.GetIdx())
        onward = next(n for n in atom.GetNeighbors() if n.GetIdx() != previous)
        previous, atom = atom.GetIdx(), onward
    if len(chain) < 3 or atom.GetAtomicNum() != 6 or atom.GetIdx() == carbon:
        return False
    return bool(_double_oxygens(mol, atom.GetIdx())) and any(
        n.GetAtomicNum() == 6 for n in atom.GetNeighbors()
    )


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
    """The ester oxygen of an acyl-O-N group on an acyclic amine nitrogen: a pseudoketone (P-65.6.3.4.1), since only a
    cyclic nitrogen gives a traditional ester."""
    if oxygen.GetAtomicNum() != 8 or oxygen.GetDegree() != 2 or oxygen.GetFormalCharge():
        return False
    far = next((n for n in oxygen.GetNeighbors() if n.GetIdx() != carbon), None)
    if far is None or far.GetAtomicNum() != 7 or far.GetFormalCharge() or far.IsInRing() or far.GetIsAromatic():
        return False
    if any(b.GetBondTypeAsDouble() != 1.0 for b in far.GetBonds()):
        return False
    return all(
        n.GetIdx() == oxygen.GetIdx() or (n.GetAtomicNum() == 6 and not _double_oxygens(mol, n.GetIdx()))
        for n in far.GetNeighbors()
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


def _group_weight(mol, carbon, cls):
    """How many groups of class `cls` sit on `carbon` (two sulfonic acid groups on one carbon: methanedisulfonic acid)."""
    return sum(1 for c, _ in _sulfonyl_groups_on(mol, carbon) if c == cls) or 1


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
            n.GetAtomicNum() == 6 or _acyl_diazene_nitrogen(mol, n, carbon) for n in atom.GetNeighbors() if n.GetIdx() != found[0]
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
            raise UnsupportedStructure("a cumulated carbonyl (ketene) is not supported")
        if not hetero:
            if len(carbon_neighbors) == 1 and atom.GetTotalNumHs() == 1:
                return "aldehyde", {oxygens[0]}
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
        if len(hetero) == 2 and len(ring_nitrogens) == 1:
            other = next(h for h in hetero if h.GetIdx() != ring_nitrogens[0].GetIdx())
            if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
                return "acid", {oxygens[0], other.GetIdx()}
            if other.GetAtomicNum() == 7 and (
                _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, carbon)
            ):
                return "amide", {oxygens[0], other.GetIdx()}
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
    for z, hydrogens, name in ((8, 1, "alcohol"), (16, 1, "thiol"), (34, 1, "selenol"), (52, 1, "tellurol"), (7, 2, "amine")):
        for n in _single_neighbors(mol, carbon, z):
            # in an aminium name only the cationic nitrogens are the suffix; a neutral amino group is a prefix
            if AMINIUM.get() is True and z == 7 and not mol.GetAtomWithIdx(n).HasProp("_cationic_amine"):
                continue
            if _terminal_heteroatom(mol, n, hydrogens):
                return name, {n}
    for n in _single_neighbors(mol, carbon, 8):
        bridge = mol.GetAtomWithIdx(n)
        ends = [m for m in bridge.GetNeighbors() if m.GetIdx() != carbon]
        if bridge.GetDegree() == 2 and len(ends) == 1 and ends[0].GetAtomicNum() == 8 and _terminal_heteroatom(mol, ends[0].GetIdx(), 1):
            return "peroxol", {n, ends[0].GetIdx()}
    return None


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


def name_polyfunctional(mol) -> str:
    from ._isotope_labels import split_isotopes
    from ._linear_phane import has_linear_phane_shape, name_linear_phane
    from ._radical_group import has_radical_group_shape, name_radical_group

    if has_linear_phane_shape(mol):
        return name_linear_phane(mol)
    if has_radical_group_shape(mol):
        return name_radical_group(mol)
    split = split_isotopes(mol)
    if split is not None:
        ammonium = _aminium_base(split[0])
        if ammonium is not None:
            return _name_isotopic(mol, *split, build=lambda clean, labels: _name_aminium(ammonium, labels))
        center = _ring_center_base(split[0])
        if center is not None:
            return _name_isotopic(mol, *split, build=lambda clean, labels: _name_labelled_ring_center(clean, labels, center))
        return _name_isotopic(mol, *split)
    anionic_parent = any(a.HasProp(ANION_PROP) for a in mol.GetAtoms())
    cation = None if anionic_parent else _aminium_base(mol)
    if cation is not None:
        return _name_aminium(cation)
    iminium = None if anionic_parent else _iminium_base(mol)
    if iminium is not None:
        return _name_aminium(iminium, parent="imine")
    for kind in ("amide", "nitrile"):
        acylated = _group_cation_base(mol, kind)
        if acylated is not None:
            return _name_aminium(acylated, parent=kind)
    grouped = _ring_cation_with_group(mol)
    if grouped is not None:
        return _name_ring_group_cation(*grouped)
    center = _ring_center_base(mol)
    if center is not None:
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

    charged = [
        a for a in mol.GetAtoms() if a.GetFormalCharge() and not _is_nitro_part(a) and not _anionic_group_atom(mol, a)
    ]
    centers = [a for a in charged if a.GetAtomicNum() == 7 and a.IsInRing() and a.GetFormalCharge() == 1]
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
        probe = _probe_with_iodine(mol, center, oxides)
        if probe is None:
            return None
    neutral = Chem.RWMol(mol)
    index = center.GetIdx()
    atom = neutral.GetAtomWithIdx(index)
    atom.SetFormalCharge(0)
    atom.SetNoImplicit(True)
    hydrogens = 0 if oxides else max(center.GetTotalNumHs() - 1, 0)
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


def _probe_with_iodine(mol, center, oxides):
    """The same ring cation with an iodine on the centre: its name numbers the centre lowest, as the name of an
    N-oxide or protonated ring (no substituent on the centre) must; the iodo prefix is cut out afterwards."""
    if any(a.GetAtomicNum() == 53 for a in mol.GetAtoms()):
        return None
    probe = Chem.RWMol(mol)
    if oxides:
        iodine = probe.GetAtomWithIdx(oxides[0].GetIdx())
        iodine.SetAtomicNum(53)
        iodine.SetFormalCharge(0)
    else:
        added = probe.AddAtom(Chem.Atom(53))
        probe.AddBond(center.GetIdx(), added, Chem.BondType.SINGLE)
        probe.GetAtomWithIdx(center.GetIdx()).SetNumExplicitHs(0)
        probe.GetAtomWithIdx(center.GetIdx()).SetNoImplicit(True)
    try:
        probe = probe.GetMol()
        Chem.SanitizeMol(probe)
    except Exception:
        return None
    found = _ring_center_base(probe)
    return None if found is None else found[:2]


def _without_iodo(name, locant):
    token = f"{locant}-iodo"
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
    shift = clean.GetNumAtoms() - base.GetNumAtoms()
    removed = next((a.GetIdx() for a in clean.GetAtoms() if a.GetFormalCharge() == -1 and a.GetAtomicNum() == 8), None)
    if shift and removed is None:
        raise UnsupportedStructure("the isotopic modification of this ring cation is not supported")
    centre_in_clean = center + (1 if shift and center >= removed else 0)
    if STEREO_OF_ISOTOPOLOGUE.get() or any(i in (removed, centre_in_clean) for i in labels):
        raise UnsupportedStructure("an isotopic modification of the cationic centre or its stereo is not supported yet")
    moved = {i - (1 if shift and i > removed else 0): e for i, e in labels.items()}
    return _name_ring_center(base, center, kind, probe, moved)


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
    bare, placed = _bare_ring_name(ring_mol, new_center)
    if not bare.endswith("e"):
        raise UnsupportedStructure("the cationic ring parent has no terminal 'e' to elide (P-73.1.1.2)")
    stem = bare[:-1]
    pattern = re.compile(re.escape(stem) + r"(?:e(?![a-z])|(?=-?\d*,?\d*-?yl)|(?=-\d+(?:,\d+)*-[a-z]))")
    found = list(pattern.finditer(name))
    if len(found) != 1:
        raise UnsupportedStructure("the ring parent of the cation is not delimited in the anionic name")
    match = found[0]
    return f"{name[:match.start()]}{stem}-{locant or placed[new_center]}-ium{name[match.end():]}"


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

    if kind == "ium" and any(_anionic_group_atom(base, a) or a.HasProp(ANION_PROP) for a in base.GetAtoms()):
        return _name_labelled(base, labels or {}, lambda name, placed: _cation_inserted(name, base, center, placed))
    if not labels and len(_system_of(base, center)[1]) == base.GetNumAtoms():
        name, placed = _bare_ring_name(base, center)
        return _attach(name, placed, center, kind)
    token = RING_CENTER.set(True)
    cation_token = CATION_PARENT.set(True)
    try:
        if probe is None:
            return _name_labelled(base, labels or {}, lambda name, placed: _attach(name, placed, center, kind))
        marked, marked_center = probe

        def attach_marked(name, placed):
            locant = placed.get(marked_center)
            if locant is None:
                raise UnsupportedStructure("the ring nitrogen is not numbered in the parent hydride")
            return _attach(_without_iodo(name, locant), {marked_center: locant}, marked_center, kind)

        return _name_labelled(marked, labels or {}, attach_marked)
    finally:
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


def _bare_ring_name(base, center):
    """(name, {center: locant}) of an unsubstituted single-heteroatom ring system."""
    from ._fusion_name import Context, fused_ring_system_name, fusion_name, system_numbering_options
    from ._fused_numbering import _locant_key
    from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic

    if base.GetRingInfo().NumRings() == 1:
        if has_hetero_monocyclic_name(base):
            return name_hetero_monocyclic(base), {center: _monocycle_locant(base, center)}
        from .core import smiles_to_iupac

        try:
            named = smiles_to_iupac(Chem.MolToSmiles(base))
        except UnsupportedStructure:
            raise UnsupportedStructure("this ring has no supported parent name") from None
        return re.sub(r"^(?:\d+[a-z]?(?:,\d+[a-z]?)*H-)+", "", named), {center: _monocycle_locant(base, center)}
    name = fused_ring_system_name(base)
    if name is None:
        raise UnsupportedStructure("this ring system has no supported parent name")
    ctx = Context(base)
    fused, root = fusion_name(base)
    options = system_numbering_options(ctx, fused, root)
    return name, {center: min((n[center] for n in options), key=_locant_key)}


def _aminium_base(mol, ignore=frozenset()):
    """The mol with its ammonium nitrogen neutralised when the only charge is one N+ bonded to carbon and hydrogen
    (a cation outranks every acid, P-41), else None."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() and a.GetIdx() not in ignore]
    if not charged or any(a.GetFormalCharge() != 1 or a.GetAtomicNum() != 7 for a in charged):
        return None
    for nitrogen in charged:
        if (
            nitrogen.GetIsAromatic()
            or nitrogen.IsInRing()
            or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 4
            or any(
                n.GetAtomicNum() != 6
                or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
                or _double_oxygens(mol, n.GetIdx())
                for n in nitrogen.GetNeighbors()
            )
        ):
            return None
    if all(n.GetTotalNumHs() == 0 for n in charged):
        return mol
    neutral = Chem.RWMol(mol)
    for nitrogen in charged:
        atom = neutral.GetAtomWithIdx(nitrogen.GetIdx())
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(nitrogen.GetTotalNumHs() - 1 if nitrogen.GetTotalNumHs() else 0)
        atom.SetNoImplicit(True)
        atom.SetBoolProp("_cationic_amine", True)
    Chem.SanitizeMol(neutral)
    return neutral.GetMol()


def _cationic_prefix_nitrogen(atom):
    """An acyclic N+ that an anionic parent cites as an 'azaniumyl' or 'azaniumylidene' prefix (P-74.1.3)."""
    return (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == 1
        and not atom.IsInRing()
        and not atom.GetIsAromatic()
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    )


def _iminium_base(mol):
    """The mol with its iminium nitrogen neutralised when it still carries a hydrogen, else the cation itself, for the
    single =N(+)< centre of an acyclic C=N group whose other substituents are carbon (P-73.1.2.1); else None."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() != 1 or charged[0].GetAtomicNum() != 7:
        return None
    nitrogen = charged[0]
    if nitrogen.GetIsAromatic() or nitrogen.IsInRing() or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 3:
        return None
    doubles = [b for b in nitrogen.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
    if len(doubles) != 1 or doubles[0].GetOtherAtom(nitrogen).GetAtomicNum() != 6:
        return None
    if any(n.GetAtomicNum() != 6 for n in nitrogen.GetNeighbors()):
        return None
    if nitrogen.GetTotalNumHs() == 0:
        return mol
    neutral = Chem.RWMol(mol)
    atom = neutral.GetAtomWithIdx(nitrogen.GetIdx())
    atom.SetFormalCharge(0)
    atom.SetNumExplicitHs(nitrogen.GetTotalNumHs() - 1)
    atom.SetNoImplicit(True)
    Chem.SanitizeMol(neutral)
    return neutral.GetMol()


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


def _group_cation_base(mol, kind, ignore=frozenset()):
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
            if len(triple) != 1 or nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 1:
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
    return base


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
        ("amide", lambda m: _group_cation_base(m, "amide", ignore=frozenset(centres))),
        ("nitrile", lambda m: _group_cation_base(m, "nitrile", ignore=frozenset(centres))),
    ):
        base = function(mol)
        if base is not None:
            return base, centres[0], kind
    return None


def _name_ring_group_cation(base, centre, kind):
    """A ring cation together with a cationic suffix: the skeletal centre is cited as 'ium' after the parent hydride and
    takes the lowest locant before the suffix (P-73.5.3.1, P-73.5.3.2); the suffix names the group."""
    from ._diester_ring_diyl import CATION_CENTRES

    editable = Chem.RWMol(base)
    atom = editable.GetAtomWithIdx(centre)
    hydrogens = atom.GetTotalNumHs()
    atom.SetFormalCharge(0)
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
    return _cationic_group_suffix(name, kind)


def _name_aminium(base, labels=None, parent="amine"):
    token = AMINIUM.set(True if parent == "amine" else parent)
    try:
        name = _name_labelled(base, labels or {})
    finally:
        AMINIUM.reset(token)
    return _cationic_group_suffix(name, parent)


def _cationic_group_suffix(name, parent):
    if parent in ("amide", "nitrile"):
        match = _GROUP_SUFFIX.search(name)
        if match is None or match.group("kind") != parent:
            raise UnsupportedStructure("the cation is not named as an amide or nitrile parent")
        word = (match.group("carbo") or "") + parent[:-1] + "ium"
        if match.group("mult"):
            multiplier = {"di": "bis", "tri": "tris", "tetra": "tetrakis"}[match.group("mult")]
            return f"{name[:match.start()]}{multiplier}({word})"
        return name[:match.start()] + word
    if not name.endswith(("amine", "aniline", "imine") if parent == "imine" else ("amine", "aniline")):
        raise UnsupportedStructure("the cation is not named as an amine or imine parent")
    multiple = re.search(r"(di|tri|tetra|penta|hexa)(amine|aniline)$", name)
    if multiple is not None:
        word = {"di": "bis", "tri": "tris", "tetra": "tetrakis", "penta": "pentakis", "hexa": "hexakis"}[multiple.group(1)]
        return f"{name[:multiple.start()]}{word}(aminium)"
    return name[:-1] + "ium"


STEREO_OF_ISOTOPOLOGUE = contextvars.ContextVar("stereo_of_isotopologue", default=None)
JUNIOR_GROUP = "_junior_group"


def _modifications(entry):
    return [*([entry["skeleton"]] if entry["skeleton"] else []), *(n for n, c in entry["H"].items() for _ in range(c))]


def _nuclide_rank(nuclide):
    """Higher atomic number first, then higher mass number (P-82.2.2.2)."""
    table = Chem.GetPeriodicTable()
    symbol = "".join(ch for ch in nuclide if ch.isalpha())
    return -table.GetAtomicNumber(symbol), -int("".join(ch for ch in nuclide if ch.isdigit()))


def _demote_junior_groups(mol, labels):
    """Carbon-suffix groups on a ring that are modified in different ways cannot be multiplied: the group with the
    most modifications, then the nuclide of higher atomic number and mass number, stays the suffix and the others are
    cited as prefixes (P-82.2.2.2)."""
    if not labels:
        return mol
    principal, owned_groups = _principal_owned(mol)
    groups = {frozenset(g) for g in owned_groups if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in g)}
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
        for atom in group:
            if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6:
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
        rdCIPLabeler.AssignCIPLabels(probe)
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
    from ._isotope_labels import descriptor

    in_parent = {a: e for a, e in labels.items() if a in parts[4]}
    if in_parent and len(parts[4]) > 1 and not _locants_omitted(mol, in_parent, parts):
        force_token = FORCE_LOCANTS.set(True)
        try:
            name, parts = reselect()
        finally:
            FORCE_LOCANTS.reset(force_token)
    # a suffix group cannot also be cited as a prefix, so labels that candidate names consumed there are still open
    suffix_atoms = set().union(*_principal_owned(mol)[1])
    rest = {a: e for a, e in labels.items() if a not in in_parent and (a not in consumed or a in suffix_atoms)}
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in rest):
        name, parts, suffix_text = _carbon_group_label(mol, rest, reselect)
        front = []
    else:
        front, suffix_text = _group_atom_labels(mol, rest) if rest else ([], None)
    if suffix_text:
        name = _before_suffix(name, *suffix_text)
    if in_parent or front:
        bare = len(parts[4]) == 1 or _locants_omitted(mol, in_parent, parts)
        capacity = {a: mol.GetAtomWithIdx(a).GetTotalNumHs() for a in in_parent}
        name = _with_descriptor(name, parts, descriptor(in_parent, parts[4], bare, front, capacity, _sole_heteroatoms(mol, in_parent, parts)))
    return name


def _carbon_group_label(mol, rest, reselect):
    """The nuclide on the carbon of a 'carbo' suffix is cited before the suffix of the systematic name, because the
    retained name numbers no such position (P-82.6.3.2: 'benzene(13C)carboxylic acid')."""

    principal, owned_groups = _principal_owned(mol)
    carbon_groups = {frozenset(g) for g in owned_groups if any(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in g)}
    modified = {g for g in carbon_groups if any(a in rest for a in g)}

    def nuclides(group):
        return sorted(
            (n for a in group if a in rest for n in _modifications(rest[a])), key=_nuclide_sort_key
        )

    from ._isotope_labels import _nuclide_sort_key

    if (
        principal is None
        or not carbon_groups
        or modified != carbon_groups
        or any(a not in set().union(*carbon_groups) for a in rest)
        or len({tuple(nuclides(g)) for g in carbon_groups}) != 1
        or any(e["H"] for a, e in rest.items() if mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
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


def _name_labelled(mol, labels, finish=None):
    from ._isotope_labels import descriptor
    from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS

    mol = _demote_junior_groups(mol, labels)

    stereo = _check_scope(mol) + [
        ("isotope", atom, "|".join(filter(None, [entry["skeleton"], *(n for n, c in entry["H"].items() for _ in range(c))])))
        for atom, entry in labels.items()
    ]
    context = {
        "atoms": {where: code for kind, where, code in stereo if kind == "atom"},
        "bonds": {where: code for kind, where, code in stereo if kind == "bond"},
        "used": set(),
    }
    token = BRANCH_STEREO.set(context if stereo else None)
    isotope_context = {"labels": labels, "consumed": set(), "mol": mol}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if labels else None)
    try:
        try:
            multiplicative = _multiplicative_name(mol, stereo)
        except UnsupportedStructure:
            if not labels:
                raise
            multiplicative = None
        if multiplicative is not None:
            if labels:
                raise UnsupportedStructure("isotopic modification of a multiplicative name is not supported yet")
            return multiplicative if finish is None else finish(multiplicative, {})
        isotope_context["consumed"].clear()
        isotope_context.pop("cache", None)
        _, name, parts = _select(mol, stereo=stereo)
        LAST_POSITIONS.set((mol, dict(parts[4])))
        if labels:
            name = _with_labels(
                mol, labels, isotope_context["consumed"], name, parts, lambda: _select(mol, stereo=stereo)[1:]
            )
        if finish is not None:
            name = finish(name, parts[4])
        return _stereo_prefix(stereo, parts[4], ring_parent=parts[5], used=context["used"]) + name
    finally:
        BRANCH_STEREO.reset(token)
        ISOTOPE_LABELS.reset(isotope_token)


_SUFFIX_WORD = {"alcohol": "ol", "ketone": "one", "aldehyde": "al"}


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
    return principal, groups.get(principal, [])


def _group_atom_labels(mol, rest):
    """([(nuclide, locant, count, repeatable)], suffix descriptor text) for nuclides on atoms of the principal
    characteristic group: an amide or amine nitrogen and the oxygens of an acid are cited in front of the parent
    with a letter locant (P-82.2.4, P-82.2.5, P-82.6.2); a hydroxy, oxo or aldehyde oxygen is inserted before the
    suffix (P-82.2.1, P-82.6.1.1)."""
    from ._isotope_labels import _nuclide_sort_key

    principal, owned_groups = _principal_owned(mol)
    owned = set().union(*owned_groups) if owned_groups else set()
    if principal is None or any(a not in owned for a in rest):
        raise UnsupportedStructure("an isotopically modified atom outside the parent and the principal group is not supported yet")
    elements = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in rest}
    if len(elements) != 1:
        raise UnsupportedStructure("isotopic modification of several elements of a characteristic group is not supported yet")
    (element,) = elements
    if any(
        a.GetAtomicNum() == element and a.GetIdx() not in owned and a.GetIdx() not in rest for a in mol.GetAtoms()
    ):
        raise UnsupportedStructure("the modified atom of the characteristic group needs a locant that is not defined yet")
    symbol = mol.GetAtomWithIdx(next(iter(rest))).GetSymbol()
    heavy = [a for a, e in rest.items() if e["skeleton"]]
    if len(heavy) > 1:
        raise UnsupportedStructure("several isotopically modified atoms in one characteristic group are not supported yet")
    if principal in _SUFFIX_WORD and element == 8:
        if len(owned_groups) != 1:
            raise UnsupportedStructure("this isotopically modified oxygen has no defined suffix descriptor")
        nuclides = []
        for entry in rest.values():
            if entry["skeleton"]:
                nuclides.append(entry["skeleton"])
            nuclides.extend(n for n, c in entry["H"].items() for _ in range(c))
        if len(nuclides) != len(set(nuclides)) or len(rest) != 1:
            raise UnsupportedStructure("several isotopically modified atoms at the suffix are not supported yet")
        return [], ("(" + ",".join(sorted(nuclides, key=_nuclide_sort_key)) + ")", _SUFFIX_WORD[principal])
    nitrogen = principal in ("amide", "amine", "nitrile") and element == 7
    if principal == "nitrile" and any(e["H"] for e in rest.values()):
        raise UnsupportedStructure("a nitrile nitrogen carries no hydrogen")
    acid = _is_acid_family(principal) and element == 8
    if not (nitrogen or acid):
        raise UnsupportedStructure("an isotopically modified atom of this characteristic group is not supported yet")
    items = []
    for atom, entry in rest.items():
        heavy_label = entry["skeleton"]
        locant = f"{heavy_label[:-len(symbol)]}{symbol}" if heavy_label else symbol
        if heavy_label:
            items.append((heavy_label, None, 1, False))
        for nuclide, count in entry["H"].items():
            items.append((nuclide, locant, count, nitrogen and _capacity(mol, atom) > 1))
    return items, None


def _capacity(mol, atom):
    target = mol.GetAtomWithIdx(atom)
    return target.GetIntProp("_capacity") if target.HasProp("_capacity") else target.GetTotalNumHs()


def _sole_heteroatoms(mol, in_parent, parts):
    """Nuclides of a heteroatom that is the only atom of its element in the parent, so its locant is implied."""
    counts = {}
    for a in parts[4]:
        number = mol.GetAtomWithIdx(a).GetAtomicNum()
        counts[number] = counts.get(number, 0) + 1
    return frozenset(
        e["skeleton"] for a, e in in_parent.items()
        if e["skeleton"] and mol.GetAtomWithIdx(a).GetAtomicNum() != 6 and counts[mol.GetAtomWithIdx(a).GetAtomicNum()] == 1
    )


def _locants_omitted(mol, in_parent, parts):
    """P-82.6.1.3 (every parent position modified in the same way, none keeping a hydrogen) and the one modified atom
    of a bare hydrocarbon whose positions are all equivalent (benzene, ethane; P-82.6.1.1)."""
    positions = set(parts[4])
    if set(in_parent) == positions and len({repr(sorted(e["H"])) + str(e["skeleton"]) for e in in_parent.values()}) == 1:
        if all(sum(in_parent[a]["H"].values()) == mol.GetAtomWithIdx(a).GetTotalNumHs() for a in positions):
            return True
    ranks = Chem.CanonicalRankAtoms(mol, breakTies=False)
    return (
        len(in_parent) == 1
        and sum(bool(e["skeleton"]) + sum(e["H"].values()) for e in in_parent.values()) == 1
        and len(positions) == mol.GetNumAtoms()
        and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in positions)
        and len({ranks[a] for a in positions}) == 1
    )


def _before_suffix(name, text, word):
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


def _stereo_entries(stereo, position_of, ring_parent=False, used=frozenset()):
    """[(locant, code)] for the stereo elements found on the parent, plus
    whether every element was placed (on the parent or inside a substituent)."""
    entries, complete = [], True
    for kind, where, code in stereo or []:
        if kind == "isotope":
            continue
        if kind == "atom":
            if where in position_of:
                entries.append((position_of[where], code))
            elif ("atom", where) not in used:
                complete = False
        else:
            a, b = where
            if ring_parent or a not in position_of or b not in position_of:
                inside = [x for x in (a, b) if x in position_of]
                if len(inside) == 2:
                    entries.append((min(position_of[a], position_of[b]), code))
                elif ("bond", where) in used:
                    continue
                elif len(inside) == 1:
                    # P-93.4.2.1, P-93.5.1.4.2.1: the configuration of a double bond to an ylidene group takes the locant
                    # of the parent atom
                    entries.append((position_of[inside[0]], code))
                else:
                    complete = False
            else:
                entries.append((min(position_of[a], position_of[b]), code))
    return sorted(entries), complete


def _stereo_rank(stereo, position_of, ring_parent=False):
    entries, _ = _stereo_entries(stereo, position_of, ring_parent)
    marked = [(w, code) for k, w, code in stereo or [] if k == "isotope" and w in position_of]
    isotopic = tuple(sorted(position_of[w] for w, code in marked for _ in code.split("|")))
    return (
        _isotope_counts(stereo, marked),
        isotopic,
        _nuclide_precedence(marked, position_of),
        _substituent_isotope_locants(position_of),
        tuple(0 if code in "RZr" else 1 for _, code in entries),
    )


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


def _parent_configuration_rank(stereo, position_of, ring_parent=False):
    """P-44.4.1.11.1-3, P-44.4.1.12: between candidate parents the one with more isotopically modified atoms is senior,
    then the one with more Z double bonds (and lower locants for them), more like descriptor pairs, and R ahead of S."""
    entries, _ = _stereo_entries(stereo, position_of, ring_parent)
    marked = [(w, code) for k, w, code in stereo or [] if k == "isotope" and w in position_of]
    cis = [locant for locant, code in entries if code == "Z"]
    centres = [code for _, code in entries if code in ("R", "S", "r", "s")]
    like = sum(1 for code in centres[1:] if code == centres[0])
    return (
        _isotope_counts(stereo, marked),
        -len(cis),
        tuple(cis),
        -like,
        tuple(0 if code in "Rr" else 1 for code in centres),
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


def _stereo_prefix(stereo, position_of, ring_parent=False, used=frozenset()):
    if not stereo:
        return ""
    entries, complete = _stereo_entries(stereo, position_of, ring_parent, used)
    if not complete:
        raise UnsupportedStructure("stereodescriptors outside the parent are not supported by the chain engine yet")
    if not entries:
        return ""
    if len(position_of) == 1 and len(entries) == 1:
        # P-93.5: the one skeletal atom of a mononuclear parent needs no locant
        return f"({entries[0][1]})-"
    return "(" + ",".join(f"{locant}{code}" for locant, code in entries) + ")-"


def _select(mol, attach=None, n_names=(), stereo=None):
    classes = {_group_of(mol, a.GetIdx())[0] for a in mol.GetAtoms() if _group_of(mol, a.GetIdx())} | {
        c for c, _, _ in _ring_occurrences(mol)
    }
    principal = _principal_class(classes)
    token = EXTENDED_PREFIXES.set((_is_acid_family(principal) or principal in ("thioic", "peroxoic", "imidic")) if principal else False)
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
        secondary = any(_substituted_amine_nitrogen(mol, mol.GetAtomWithIdx(i)) for i in marked)
        if not groups and not ring_groups and not secondary:
            raise UnsupportedStructure("this anionic group is not a principal-capable characteristic group")
    classes = set(groups) | {c for c, _, _ in ring_groups}
    principal = _principal_class(classes)
    token = IDE_EXTRA.set(ide_extra)
    try:
        results = [
            _select_with_principal(mol, graph, halogens, aromatic_atoms, kept, kept_rings, principal, attach, n_names, stereo)
            for kept, kept_rings in _diamidide_alternatives(mol, groups, ring_groups)
        ]
        best = min(results, key=lambda result: result[0])
    finally:
        IDE_EXTRA.reset(token)
    if len(ide_extra) > 1 and not ring_groups:
        anions = _carbanion_parent(mol, graph, halogens, aromatic_atoms, ide_extra, groups, attach, n_names, stereo)
        if anions is not None and _prefer_anion_parent(best[1], anions[1]):
            return anions
    return best


def _carbanion_parent(mol, graph, halogens, aromatic_atoms, centres, groups, attach, n_names, stereo):
    """The parent chain named by its carbanide centres alone, the anionic acid groups being prefixes (P-72.7b)."""
    ide_groups = {"ide": {idx: {idx} for idx in centres}}
    token = IDE_EXTRA.set({})
    try:
        return _select_with_principal(mol, graph, halogens, aromatic_atoms, ide_groups, [], "ide", attach, n_names, stereo)
    except UnsupportedStructure:
        return None
    finally:
        IDE_EXTRA.reset(token)


_IDE_COUNT = re.compile(r"-(\d+(?:,\d+)*)-(?:di|tri|tetra)?ide$")
_ID_COUNT = re.compile(r"-(\d+(?:,\d+)*)-(?:di|tri|tetra)?id-")


def _prefer_anion_parent(acid_name, anion_name):
    """P-72.7: with as many anionic centres in the parent, the parent with more 'ide' centres is senior."""
    anion = _IDE_COUNT.search(anion_name)
    centre = _ID_COUNT.search(acid_name)
    if anion is None or centre is None:
        return False
    return len(anion.group(1).split(",")) > len(centre.group(1).split(","))


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


def _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal, attach, n_names, stereo):
    if any(_chalcogen_ketone(mol, a) for a in mol.GetAtoms()) and principal not in _CHALCOGEN_KETONE_OK and not (principal and _is_variant(principal)):
        raise UnsupportedStructure("a thioketone-type group outranks the parents this engine can build here")
    for atom in mol.GetAtoms():
        if atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("isotopes and radicals are not supported by the polyfunctional chain engine")
        if atom.GetFormalCharge() and not _is_nitro_part(atom) and not _anionic_group_atom(mol, atom) and not (
            AMINIUM.get() and atom.GetAtomicNum() == 7
        ) and not (
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
            and not (principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES) and atom.GetIdx() in groups.get("hydrazide", {}))
            and _is_ester_like(mol, atom.GetIdx())
            and not (_urea_carbon(mol, atom.GetIdx()) and _outranks_urea(principal))
            and not (FORCED_PRINCIPAL.get() == principal and principal == "nitrile")
            and principal != "ide"
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
        if attach is not None or n_names:
            raise UnsupportedStructure("N-substituted amines inside a unit are not handled by the chain engine")
        if any(kind != "isotope" for kind, _, _ in stereo or ()):
            raise UnsupportedStructure("stereodescriptors with an N-substituted amine parent are not supported yet")
        return _substituted_amine(mol, graph, halogens, aromatic_atoms, groups, ring_groups)

    if principal is None:
        if attach is not None:
            raise UnsupportedStructure("a unit without a principal group is not supported")
        return _plain_parent(mol, graph, halogens, aromatic_atoms, stereo)

    if _is_variant(principal):
        imidic_ns = _imidic_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if imidic_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted imidic acid inside a unit is not handled by the chain engine")
            n_names = imidic_ns

    if principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "sulfonamide", *_CHALCOGEN_SULFONAMIDE_CLASSES):
        amide_ns = _amide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if amide_ns:
            if n_names:
                raise UnsupportedStructure("an N-substituted amide inside a unit is not handled by the chain engine")
            n_names = amide_ns

    if principal in _CHALCOGEN_IMIDAMIDE.values():
        imidamide_ns = _imidamide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if imidamide_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted imidamide inside a unit is not handled by the chain engine")
            n_names = imidamide_ns

    if principal == "hydrazide" or principal in _CHALCOGEN_HYDRAZIDE.values():
        hydrazide_ns = _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if hydrazide_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted hydrazide inside a unit is not handled by the chain engine")
            n_names = hydrazide_ns

    if principal in _AMIDRAZONE:
        amidrazone_ns = _amidrazone_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if amidrazone_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted amidrazone inside a unit is not handled by the chain engine")
            n_names = amidrazone_ns

    if principal == "amidine":
        amidine_ns = _amidine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups)
        if amidine_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted amidine inside a unit is not handled by the chain engine")
            n_names = amidine_ns

    if principal == "imine":
        imine_ns = _imine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups)
        if imine_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted imine inside a unit is not handled by the chain engine")
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
                or ((_is_acid_family(principal) or principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES, "hydrazide", "amidine", *_AMIDRAZONE)) and _junior_end_group(mol, a.GetIdx()))
            )
        }
        try:
            for path in _paths(graph, eligible):
                for chain in (path, path[::-1]):
                    candidate = _evaluate(
                        mol, graph, halogens, aromatic_atoms, chain, principal, principal_atoms, owned, attach, n_names, stereo
                    )
                    if chain_best is None or candidate[0] < chain_best[0]:
                        chain_best = candidate
        except UnsupportedStructure as error:
            chain_error, chain_best = error, None
        if chain_best is not None and chain_best[0][0] == 1:
            chain_best = None
    chain_count = -chain_best[0][0] if chain_best else 0

    assembly = _assembly_parent(
        mol, graph, halogens, aromatic_atoms, principal, [g for g in ring_groups if g[0] == principal], stereo
    )
    if assembly is not None and assembly[0] >= chain_count:
        if attach is not None or n_names:
            raise UnsupportedStructure("a ring assembly inside a unit is not supported yet")
        return _finish(assembly[1])

    ring_best = _best_ring(
        mol, graph, halogens, aromatic_atoms, principal, [g for g in ring_groups if g[0] == principal], n_names, stereo
    )
    ring_count = ring_best[0] if ring_best else 0

    if chain_error is not None and not (ring_count and ring_count >= total_principal):
        raise chain_error
    if ring_count and ring_count >= chain_count:
        if attach is not None:
            raise UnsupportedStructure("a ring parent inside a multiplicative unit is not supported yet")
        return _finish(ring_best[1])
    if chain_best is None:
        raise UnsupportedStructure("no parent carries the principal group")
    if _is_terminal(principal) and chain_count < len(principal_atoms):
        carbo = _carbo_best(mol, graph, halogens, aromatic_atoms, principal, principal_atoms, stereo, n_names)
        if carbo is not None and -carbo[0][0] > chain_count:
            return _finish(carbo)
    if (
        _is_terminal(principal)
        and not (principal in ("amidine", *_AMIDRAZONE) and attach is None and not mol.GetRingInfo().NumRings())
        and len(chain_best[2][4]) == 1
        and not (_is_variant(principal) and attach is None)
        and not chain_best[1].endswith("formic acid")
        and not (principal in ("amide", *_CHALCOGEN_AMIDE_CLASSES) and attach is None)
    ):
        raise UnsupportedStructure("one-carbon acid, amide, nitrile and aldehyde parents use retained names")
    return _finish(chain_best)


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
        from ._diester_ring_diyl import evaluate_skeleton

        found = evaluate_skeleton(mol, graph, "ring", [ring], ring_set, [], set(), "")
        if found is None:
            raise UnsupportedStructure("this ring has no supported name")
        from ._diester_ring_diyl import PARENT_START

        placed = found[2]
        name = _without_stereo(found[1])
        return (modified, _parent_configuration_rank(stereo, placed), -len(roots), tuple(sorted(placed[r] for r, _ in roots)), alphanumerical_name_key(name), name), (
            (0,), name, (None, None, None, 0, placed, True, PARENT_START.get())
        )
    from ._substituents import ISOTOPE_LABELS

    isotope_context = ISOTOPE_LABELS.get()
    ring_labelled = isotope_context is not None and any(a in ring_set for a in isotope_context["labels"])
    if not roots:
        if not ring_labelled:
            raise UnsupportedStructure("an unsubstituted ring is not a polyfunctional case")
        locants = min(numberings(spec), key=lambda option: _stereo_rank(stereo, option))
        return (modified, _parent_configuration_rank(stereo, locants), 0, (), spec.parent), ((0,), spec.parent, (None, None, None, 0, locants, True, 0))
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
    return (modified, _parent_configuration_rank(stereo, best[1]), -len(roots), best[0][0], alphanumerical_name_key(name), name), ((0,), name, (None, None, None, 0, best[1], True, len(name) - len(spec.parent)))


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
    from ._substituents import format_mononuclear_prefixes

    z = center.GetAtomicNum()
    stem, _, valence = MONONUCLEAR_HYDRIDES[z]
    index = center.GetIdx()
    neighbors = [n for n in graph[index]]
    if (
        center.GetFormalCharge()
        or center.GetIsotope()
        or sum(mol.GetBondBetweenAtoms(index, n).GetBondTypeAsDouble() for n in neighbors) > valence
        or any(
            mol.GetBondBetweenAtoms(index, n).GetBondTypeAsDouble() != 1.0
            for n in neighbors
            if mol.GetAtomWithIdx(n).GetAtomicNum() != 6
        )
    ):
        return None
    chalcogenols = {
        word: [n for n in neighbors if _terminal_heteroatom(mol, n, 1) and mol.GetAtomWithIdx(n).GetAtomicNum() == z]
        for z, word in _CHALCOGENOL_WORDS.items()
    }
    principal_word = next((w for w, atoms in chalcogenols.items() if atoms), None)
    hydroxyls = chalcogenols["ol"] if principal_word == "ol" else []
    junior_chalcogenols = {n for w, atoms in chalcogenols.items() if w != principal_word for n in atoms}
    amines = [
        n
        for n in neighbors
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 7
        and not mol.GetAtomWithIdx(n).GetFormalCharge()
        and not mol.GetAtomWithIdx(n).IsInRing()
        and all(
            mol.GetAtomWithIdx(m).GetAtomicNum() in (6, *MONONUCLEAR_HYDRIDES) or _terminal_heteroatom(mol, m, 1)
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
        return format_mononuclear_prefixes(entries) + suffix
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
        return format_substituent_prefixes(merged, omit_locants=not n_entries) + suffix
    return format_mononuclear_prefixes(entries) + stem


def _hydroxy_suffix_count(name):
    """The number of hydroxy or chalcogenol suffixes of a hydride name: stannanol 1, silanetriol 3 (P-41, P-44.1.1)."""
    match = re.search(r"(di|tri|tetra)?(?:ol|thiol|selenol|tellurol)$", name)
    return {"di": 2, "tri": 3, "tetra": 4}.get(match.group(1), 1) if match else 0


def _amine_count(name):
    """The number of amine suffixes of a hydride name: boranamine 1, boranediamine 2 (P-44.1.1)."""
    match = re.search(r"(di|tri|tetra)amine$", name)
    return {"di": 2, "tri": 3, "tetra": 4}[match.group(1)] if match else int(name.endswith("amine"))


def _plain_parent(mol, graph, halogens, aromatic_atoms, stereo=None):
    """Parent without a principal group: the monocycle when there is one
    (P-44.1.2.2), else the longest chain, with every substituent a prefix."""
    from ._hydride_chain import name_hydride_chain

    chain_name = name_hydride_chain(mol, graph, halogens, aromatic_atoms)
    if chain_name is not None:
        return ((0,), chain_name, (None, None, None, 0, {}, False))
    ring_cation = RING_CENTER.get()
    centers = [] if ring_cation else [a for a in mol.GetAtoms() if a.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not a.IsInRing()]
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
            return ((0,), name, (None, None, None, 0, {center: 1}, False))
    if centers or (
        not ring_cation
        and any(
            b.GetBeginAtom().GetAtomicNum() == 7
            and b.GetEndAtom().GetAtomicNum() == 7
            and not b.IsInRing()
            and not _is_azide_part(b.GetBeginAtom())
            and not _is_azide_bond_end(b.GetBeginAtom(), b.GetEndAtom())
            for b in mol.GetBonds()
        )
    ):
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
        if len(rings) >= 4:
            raise UnsupportedStructure("four or more rings need a phane or ring-assembly name")
        ranked = sorted(rings, key=lambda r: _ring_rank(mol, r[0]))
        top = _ring_rank(mol, ranked[0][0])
        tied = [r for r in ranked if _ring_rank(mol, r[0]) == top]
        for i, first in enumerate(tied):
            for second in tied[:i]:
                bonded = any(mol.GetBondBetweenAtoms(a, b) is not None for a in first for b in second)
                if bonded and _bare_key(mol, set(first)) == _bare_key(mol, set(second)):
                    raise UnsupportedStructure("identical rings joined directly form a ring assembly (P-28)")
        candidates = [_plain_ring_parent(mol, graph, halogens, aromatic_atoms, ring, stereo) for ring in tied]
        return min(candidates, key=lambda c: c[0])[1]

    eligible = {
        a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not is_functional_carbon(mol, a.GetIdx())
    }
    best = None
    for path in _paths(graph, eligible):
        for chain in (path, path[::-1]):
            candidate = _evaluate_plain(mol, graph, halogens, aromatic_atoms, chain, stereo)
            if best is None or candidate[0] < best[0]:
                best = candidate
    if best is None:
        raise UnsupportedStructure("no chain carbons")
    return best


def _without_stereo(name):
    return re.sub(r"^\(\d+[a-z]*[RSEZrs](?:,\d+[a-z]*[RSEZrs])*\)-", "", name)


def _fused_plain_parent(mol, graph, rings):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    systems = []
    for ring in rings:
        system_rings, system_atoms = _system_of(mol, ring[0])
        if not any(set(system_atoms) == set(seen[1]) for seen in systems):
            systems.append((system_rings, system_atoms))
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


_ALL_DESCRIPTORS = re.compile(r"\((?:\d+[a-z]*['\u2032]*)?[RSEZrs](?:,(?:\d+[a-z]*['\u2032]*)?[RSEZrs])*\)-")


def _configuration_senior_system(mol, graph, tied):
    """The numbered parent among equally senior ring systems that differ only in configuration: Z before E and R before
    S (P-44.4.1.12.1, P-44.4.1.12.2, P-45.6.2); identical configurations need a multiplicative or assembly name."""
    from ._diester_ring_diyl import evaluate_skeleton

    candidates = []
    for rings, atoms in tied:
        if len(rings) > 1:
            _require_mancude_system(mol, atoms)
        found = evaluate_skeleton(mol, graph, "ring", rings, atoms, [], set(), "")
        if found is None:
            raise UnsupportedStructure("this fused ring system has no supported numbering")
        candidates.append(found)
    if len({_ALL_DESCRIPTORS.sub("", c[1]) for c in candidates}) > 1 or len({c[0] for c in candidates}) == 1:
        raise UnsupportedStructure("several equally senior ring systems need a multiplicative or assembly name")
    return min(candidates, key=lambda c: c[0])


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
        _parent_configuration_rank(stereo, position_of),
        -total_count,
        locant_set,
        citation,
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
        _parent_configuration_rank(stereo, position_of),
        -total_count,
        locant_set,
        citation,
        _n_group_positions(n_names, position_of),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (prefix, body, "", 0, position_of, False)


def _locant_order(locant):
    """Primed locants follow the unprimed one of the same number: 1 < 1' < 2."""
    return locant[1], locant[0]


def _assembly_numbering(graph, rings, join, marked, entries, specs=None, cite_marked=None):
    """Atom -> (prime count, locant) for the best numbering of two directly
    joined identical rings. Benzene and cycloalkane rings number from the
    junction (1 and 1'); heteroaromatic rings keep their fixed numbering and
    the junction takes the lowest locants, then the marked atoms (principal
    groups or a free valence), then the substituents."""
    hetero = bool(specs) and specs[0].hetero is not None
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
    unsaturated = bool(specs) and any(sp.kind == "cycloalkene" for sp in specs)
    for unprimed in (0, 1):
        for first in orientations(unprimed):
            for second in orientations(1 - unprimed):
                locants = {atom: (0, number) for atom, number in first.items()}
                locants.update({atom: (1, number) for atom, number in second.items()})
                ene = ()
                if unsaturated:
                    ene = tuple(sorted(_locant_order(loc) for loc in _assembly_multiple_locants(specs, locants)[0] + _assembly_multiple_locants(specs, locants)[1]))
                key = (
                    (locants[join[unprimed]][1], locants[join[1 - unprimed]][1]),
                    tuple(sorted(_locant_order(locants[a]) for a in marked)),
                    tuple(_locant_order(locants[a]) for a in cite_marked) if cite_marked else (),
                    ene,
                    tuple(sorted(_locant_order(locants[r]) for r, _, _ in entries)),
                    _citation_key([(_locant_order(locants[r]), name) for r, name, _ in entries]),
                )
                if best is None or key < best[0]:
                    best = (key, locants, ene)
    if unsaturated and any(loc[0] == len(rings[0]) for loc in best[2]):
        raise UnsupportedStructure("a ring double bond closing the numbering (1(n) locant) is not supported in an assembly")
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

    def cite(found):
        return ",".join(f"{n}{chr(39) * prime}" for prime, n in found)

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
    stem = parent[:-1] if elide and parent.endswith("e") else parent
    return f"{spots}-bi({stem})" if spec.hetero is None else f"{spots}-bi{stem}"


def assembly_substituent(mol, graph, root, coming_from, halogens, aromatic_atoms):
    """(name, True) of a substituent group made of two identical directly
    joined rings ([1,1'-biphenyl]-4-yl), entered at `root`; None otherwise."""
    from ._system_assembly import system_assembly

    fused = system_assembly(
        mol, graph, halogens, aromatic_atoms, None, [], None, free=(root, coming_from), within=_arm_atoms(graph, root, coming_from)
    )
    if fused is not None:
        return fused
    from ._chain_assembly import chain_assembly

    chained = chain_assembly(mol, graph, halogens, aromatic_atoms, None, [], None, free=(root, coming_from))
    if chained is not None:
        return chained
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
    specs = [spec_of(mol, r) for r in (own, other)]
    if not _compatible_assembly_rings(mol, (own, other), specs):
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
    locants = _assembly_numbering(graph, rings, joins[0], [root], entries, specs)
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
    core = f"[{base}]-{enes}{spot[1]}{chr(39) * spot[0]}-yl" if not ending else f"[{base}]-{enes}{spot[1]}{chr(39) * spot[0]}-yl"
    return (assembly_join(prefix, core)), True


def _assembly_parent(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo):
    """P-28: two identical benzene or cycloalkane rings joined by a single
    bond -- [1,1'-biphenyl]-4-ol, 4-nitro-1,1'-biphenyl. Returns
    (principal count, (key, name, parts)) or None when `mol` is not one."""
    if principal is not None and not occurrences:
        return None
    from ._system_assembly import system_assembly

    fused = system_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo)
    if fused is not None:
        return fused
    from ._chain_assembly import chain_assembly

    chained = chain_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, stereo)
    if chained is not None:
        return chained
    ring_info = mol.GetRingInfo()
    all_rings = [list(r) for r in ring_info.AtomRings()]
    if len(all_rings) < 2 or any(ring_info.NumAtomRings(a) != 1 for r in all_rings for a in r):
        return None
    ring_specs = [spec_of(mol, r) for r in all_rings]
    if len(all_rings) > 2 and any(sp is not None and sp.kind == "cycloalkene" for sp in ring_specs):
        return None
    if not all(_compatible_assembly_rings(mol, (all_rings[0], r), (ring_specs[0], sp)) for r, sp in zip(all_rings, ring_specs)):
        return None
    best = None
    for i, first in enumerate(all_rings):
        for second in all_rings[i + 1 :]:
            found = _pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, [first, second])
            if found is None:
                continue
            ylidene = _junction_is_ylidene(mol, found[2], [spec_of(mol, first), spec_of(mol, second)])
            # P-28.2.2: a double-bond junction is a two-ring assembly only; the pair holding the double bond has the
            # parent's multiple bond, so it ranks first
            key = (-found[0], not ylidene, ring_seniority_key(mol, set(first) | set(second), assembly=True), found[1][1])
            if best is None or key < best[0]:
                best = (key, (found[0], found[1]))
    return best[1] if best else None


def _assembly_beside_systems(mol, graph, halogens, aromatic_atoms, stereo):
    """The ring assembly of two identical rings as the parent of a molecule that also holds other ring systems, when
    no other ring system outranks it (P-44.2.2.2.7); None otherwise."""
    from ._multiplicative import _ring_systems

    ring_info = mol.GetRingInfo()
    monocycles = [list(r) for r in ring_info.AtomRings() if all(ring_info.NumAtomRings(a) == 1 for a in r)]
    best = None
    for i, first in enumerate(monocycles):
        for second in monocycles[i + 1 :]:
            found = _pair_assembly(mol, graph, halogens, aromatic_atoms, None, [], [first, second])
            if found is None:
                continue
            ylidene = _junction_is_ylidene(mol, found[2], [spec_of(mol, first), spec_of(mol, second)])
            atoms = set(first) | set(second)
            key = (not ylidene, ring_seniority_key(mol, atoms, assembly=True), found[1][1])
            if best is None or key < best[0]:
                best = (key, atoms, (found[0], found[1]))
    if best is None:
        return None
    if any(ring_seniority_key(mol, system) < best[0][1] for system in _ring_systems(mol) if not system <= best[1]):
        return None
    return best[2]


def _pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, rings):
    """The assembly name of two directly joined identical rings with every other part of `mol` cited as a
    substituent; None when the pair is not such an assembly or leaves a principal group outside it."""
    if set(rings[0]) & set(rings[1]):
        return None
    if occurrences and any(o[1] not in set(rings[0]) | set(rings[1]) for o in occurrences):
        return None
    specs = [spec_of(mol, r) for r in rings]
    if not _compatible_assembly_rings(mol, rings, specs):
        return None
    joins = [(a, b) for a in rings[0] for b in rings[1] if mol.GetBondBetweenAtoms(a, b) is not None]
    if len(joins) != 1:
        return None
    ylidene = _junction_is_ylidene(mol, joins[0], specs)
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
    locants = _assembly_numbering(graph, rings, joins[0], [o[1] for o in occurrences], entries, specs)

    def cite(locant):
        return f"{locant[1]}{chr(39) * locant[0]}"

    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(cite(locants[r]))
    for info in grouped.values():
        info["locants"].sort(key=lambda text: (int(text.rstrip(chr(39))), text.count(chr(39))))
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    count = len(occurrences)
    ending = _assembly_unsaturation(specs, locants)
    if principal is None:
        core = _assembly_base(specs, locants, joins[0], elide=False, ylidene=ylidene)
        if ending:
            core = f"[{core}]-{ending}"
    else:
        word = multiplied_word(count, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        base = _assembly_base(specs, locants, joins[0], elide=word[0] in "aeiouy" and not ending, ylidene=ylidene)
        spots = ",".join(cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: _locant_order(locants[o[1]])))
        enes = f"{ending[:-1] if word[0] in 'aeiouy' else ending}-" if ending else ""
        core = f"[{base}]-{enes}{spots}-{word}"
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
    if len(member_rings) > 3 and not only_spiro_contacts and any(len(r) not in (5, 6, 7) for r in member_rings):
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
    "ketone": "one",
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
    found = evaluate_skeleton(mol, graph, "ring", rings, atoms, attach, blocked, _FUSED_SUFFIX[principal], n_names=n_names)
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
    if len({ranks[next(iter(c[3]))[1]] for c in leading}) == 1:
        raise UnsupportedStructure("several rings bear the principal group; a multiplicative name is needed")
    ranked = []
    for candidate in leading:
        try:
            ranked.append(
                _ring_parent(mol, graph, halogens, aromatic_atoms, principal, occurrences, n_names, stereo, candidate)
            )
        except UnsupportedStructure:
            continue
    if not ranked:
        raise UnsupportedStructure("several rings bear the principal group and none has a supported name")
    count, found, _ = min(ranked, key=lambda r: r[2])
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
        return found[0], found[1], (-len(attached), locant_set, locant_set, re.sub(r"[^a-z]", "", found[1][1]))
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
        and (spec.kind == "cycloalkane" or (spec.kind == "benzene" and (suffix_name not in _RETAINED_BENZENE or NO_RETAINED_BENZENE.get())))
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
    rank = (-len(entries), best[0][2], best[0][3], letters)
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
        return _one_of_substituted_amines(mol, graph, halogens, aromatic_atoms, groups, ring_groups, amine_nitrogens)
    if len(amine_nitrogens) != 1 or (primaries and ISOTOPE_LABELS.get()):
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
            and (_double_oxygens(mol, c) or is_functional_carbon(mol, c))
        )
        for c in neighbors
    ):
        raise UnsupportedStructure("this nitrogen is not a plain amine nitrogen")
    carbon_neighbors = [c for c in neighbors if mol.GetAtomWithIdx(c).GetAtomicNum() == 6]
    arms = {c: _arm_atoms(graph, c, n_idx) for c in carbon_neighbors}
    if sum(len(a) for a in arms.values()) != len(set().union(*arms.values())) or n_idx in set().union(*arms.values()):
        raise UnsupportedStructure("a nitrogen closing a ring is not an acyclic amine parent")
    context = ISOTOPE_LABELS.get()
    settled = set(context["consumed"]) if context else set()
    results = []
    for c in carbon_neighbors:
        others = [o for o in neighbors if o != c]
        if context:
            context["consumed"] = set(settled)
        n_names = [name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True) for o in others]
        parent, mapped = _amine_parent_molecule(mol, arms[c], c, n_idx)
        if primaries:
            n_names = [(name, compound, ("N", mapped)) for name, compound in n_names]
        unlabelled = ISOTOPE_LABELS.set(None)
        try:
            result = _select(parent, None, n_names)
        finally:
            ISOTOPE_LABELS.reset(unlabelled)
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
    if primaries:
        token = SUBSTITUTED_AMINE_PREFIX.set(True)
        try:
            as_prefix = _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, "amine", None, (), None)
        except UnsupportedStructure:
            as_prefix = None
        finally:
            SUBSTITUTED_AMINE_PREFIX.reset(token)
        if as_prefix is not None:
            ring = as_prefix[2][5]
            anchor = next(iter(as_prefix[2][4]))
            results.append(
                (as_prefix[0][0], ring, _ring_rank(mol, anchor), len(as_prefix[2][4]) if not ring else 0, as_prefix, None)
            )
        if not results:
            raise UnsupportedStructure("no parent carries the amine groups of this substituted amine")
    results.sort(key=lambda r: (r[0], not r[1], r[2], -r[3], r[6] if len(r) > 6 else (), r[4][1]))
    best = results[0]
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
            parent, _ = _amine_parent_molecule(mol, arms[best[4]], best[4], n_idx)
            sub_context = {"labels": parent_labels, "consumed": set(), "mol": parent}
            sub_token = ISOTOPE_LABELS.set(sub_context)
            try:
                selected = _select(parent, None, n_names)
                name = _with_labels(
                    parent, parent_labels, sub_context["consumed"], selected[1], selected[2],
                    lambda: _select(parent, None, n_names)[1:],
                )
            finally:
                ISOTOPE_LABELS.reset(sub_token)
            context["consumed"].update(a for a in context["labels"] if a in new_index)
            return selected[0], name, (*selected[2][:4], {}, *selected[2][5:])
    kept = sorted(set(arms[best[4]]) | {n_idx})
    positions = best[3][2][4]
    if isinstance(positions, dict) and positions:
        remapped = {kept[new]: locant for new, locant in positions.items() if new < len(kept)}
        return best[3][0], best[3][1], (*best[3][2][:4], remapped, *best[3][2][5:])
    return best[3]


def _share_chain_backbone(mol, graph, nitrogens):
    """Two of `nitrogens` joined through acyclic carbons only: they are amino groups of one chain (P-62.2.4.1.2)."""
    chain = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and not a.IsInRing()}
    for i, first in enumerate(nitrogens):
        reached, stack = {first.GetIdx()}, [first.GetIdx()]
        while stack:
            for nxt in graph[stack.pop()]:
                if nxt in chain and nxt not in reached:
                    reached.add(nxt)
                    stack.append(nxt)
        reached_nitrogens = {n for atom in reached for n in graph[atom]}
        if any(second.GetIdx() in reached_nitrogens for second in nitrogens[i + 1 :]):
            return True
    return False


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
        if len(carbons) != len(neighbors) or any(is_functional_carbon(mol, c) or _double_oxygens(mol, c) for c in carbons):
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
    for number, (n_idx, carbon, hydrides) in enumerate(nitrogens, start=1):
        editable.GetAtomWithIdx(carbon).SetAtomMapNum(number)
        for h in hydrides:
            n_names.append((*name_branch(graph, h, n_idx, halogens, aromatic_atoms, mol=mol), number))
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
    tagged = [(name, compound, ("N", mapped[number])) for name, compound, number in n_names]
    result = _select(parent, None, tagged)
    positions = result[2][4]
    if isinstance(positions, dict) and positions:
        positions = {kept[new]: locant for new, locant in positions.items() if new < len(kept)}
        return result[0], result[1], (*result[2][:4], positions, *result[2][5:])
    return result


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
    if any(a.HasProp("_ring_cation_centre") for a in parent.GetAtoms()):
        parent.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(parent)
    else:
        Chem.SanitizeMol(parent)
    mapped = next(a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum() == 1)
    for a in parent.GetAtoms():
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
            atoms = _arm_atoms(graph, root, other)
            if other in atoms or not atoms & group_atoms or not atoms - group_atoms:
                continue
            canonical, _ = _unit_molecule(mol, atoms, root)
            sides.append((bond.GetIdx(), atoms, Chem.MolToSmiles(canonical)))
    for i, (bi, ai, si) in enumerate(sides):
        for bj, aj, sj in sides[:i]:
            if bi != bj and si == sj and not ai & aj:
                if _unsymmetric_chalcogen_linker(mol, graph, ai, aj):
                    continue
                # a multiplied parent must express every principal group (P-15.6.1.5)
                units = ai | aj | set().union(*(ak for bk, ak, sk in sides if sk == si and not ak & (ai | aj)))
                if group_atoms <= units:
                    return True
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
    return None


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
    if stereo:
        raise UnsupportedStructure("stereodescriptors in a multiplicative name are not supported yet")
    from .core import smiles_to_iupac

    names = []
    for atoms, root in zip(arms, roots):
        editable = Chem.RWMol(mol)
        editable.GetAtomWithIdx(root).SetNumExplicitHs(1)
        editable.GetAtomWithIdx(root).SetNoImplicit(True)
        for idx in sorted(set(range(mol.GetNumAtoms())) - atoms, reverse=True):
            editable.RemoveAtom(idx)
        unit = editable.GetMol()
        Chem.SanitizeMol(unit)
        names.append(smiles_to_iupac(Chem.MolToSmiles(unit)))
    if names[0] != names[1]:
        return None
    if names[0] == "acetaldehyde":
        unit = "ethan-1-one"
    elif names[0].endswith(("anal", "enal", "ynal")):
        unit = names[0][:-2] + "-1-one"
    else:
        return None
    if unit[0].isdigit():
        return f"1,1'-{word}bis({unit})"
    return f"1,1'-{word}di({unit})"


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
        ):
            continue
        units = [_unit_molecule(mol, atoms, n) for atoms, n in zip(parts_atoms, neighbors)]
        keys = {Chem.MolToSmiles(unit[0], isomericSmiles=False) for unit in units}
        if len(keys) != 1:
            continue
        if number == 7 and (_principal_atoms(units[0][0]) or (None,))[0] == "amine":
            continue
        stereo_text = ""
        if stereo:
            parts, stereo_text = _stereo_arms(mol, stereo, parts_atoms, units)
        else:
            unit, attach = units[0]
            try:
                _, _, parts = _select(unit, attach)
            except UnsupportedStructure:
                continue
        prefix, body, tail, locant = parts[:4]
        lead = ",".join(str(locant) + "'" * i for i in range(arms)) + "-" if locant is not None else ""
        text = prefix + body
        if prefix:
            word = {2: "bis", 3: "tris"}[arms]
            if tail.startswith(" "):
                return f"{stereo_text}{lead}{linker}{word}({text}{tail})"
            return f"{stereo_text}{lead}{linker}{word}({text}){tail}"
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
        units = [_unit_molecule(mol, atoms, root) for atoms, (_, _, (_, root)) in zip(arms, sides)]
        if len({Chem.MolToSmiles(unit[0]) for unit in units}) != 1:
            continue
        components = []
        for ring, ring_atom, (r, root) in sides:
            named = name_ring_component(
                mol, ring, [(r, root), (ring_atom, index)], [], None, directed=(r, ring_atom)
            )
            components.append(named[0] if named else None)
        if components[0] is None or components[0] != components[1]:
            continue
        unit, attach = units[0]
        try:
            _, _, parts = _select(unit, attach)
        except UnsupportedStructure:
            continue
        if stereo:
            raise UnsupportedStructure("stereodescriptors in a multiplicative name are not supported yet")
        prefix, body, tail, locant = parts[:4]
        lead = f"{locant},{locant}'-" if locant is not None else ""
        linker = enclose(f"{word}bis({components[0]})")
        text = prefix + body
        if prefix:
            return f"{lead}{linker}bis({text}){tail}"
        return f"{lead}{linker}di{unit_phrase(text, tail)}"
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
    units = [_unit_molecule(mol, atoms, root) for atoms, (_, root) in zip(arms, attachments)]
    if len({Chem.MolToSmiles(unit[0]) for unit in units}) != 1:
        return None
    try:
        _, _, parts = _select(*units[0])
    except UnsupportedStructure:
        return None
    if not all(_stereo_within(element, central) for element in stereo):
        raise UnsupportedStructure("stereodescriptors in a multiplicative name are not supported yet")
    prefix, body, tail, locant = parts[:4]
    lead = f"{locant},{locant}'-" if locant is not None else ""
    text = prefix + body
    if prefix:
        return f"{lead}{enclose(linker)}bis({text}){tail}"
    return f"{lead}{enclose(linker)}di{unit_phrase(text, tail)}"


def _stereo_within(element, atoms):
    kind, where, _ = element
    return where in atoms if kind == "atom" else all(a in atoms for a in where)


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
        units = [_unit_molecule(mol, atoms, root) for atoms, (_, root) in zip(arms, attachments)]
        if len({Chem.MolToSmiles(unit[0]) for unit in units}) != 1:
            continue
        unit, attach = units[0]
        try:
            _, _, parts = _select(unit, attach)
        except UnsupportedStructure:
            continue
        if spec is not None:
            linker = name_ring_component(mol, ring, attachments, [], None)
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
        if stereo:
            raise UnsupportedStructure("stereodescriptors in a multiplicative name are not supported yet")
        prefix, body, tail, locant = parts[:4]
        lead = f"{locant},{locant}'-" if locant is not None else ""
        text = prefix + body
        linker_text = f"({linker[0]})" if linker[1] else enclose(linker[0])
        if prefix:
            return f"{lead}{linker_text}bis({text}){tail}"
        return f"{lead}{linker_text}di{unit_phrase(text, tail)}"
    return None


def _stereo_arms(mol, stereo, parts_atoms, units):
    """(parts of the first unit, '(2R,2'S)-' text) for identical units whose stereo differs."""
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
            try:
                _, _, parts = _select(unit, attach, stereo=arm)
            except UnsupportedStructure as error:
                if str(error).startswith("a unit without a principal group"):
                    raise UnsupportedStructure("a unit without a principal group is not supported") from error
                raise
            entries, complete = _stereo_entries(arm, parts[4], ring_parent=parts[5])
            if not complete:
                raise UnsupportedStructure("stereodescriptors outside the unit parent are not supported by the chain engine yet")
            per_arm.append((parts, entries))
    finally:
        BRANCH_STEREO.reset(token)
    if len({tuple(parts[:4]) for parts, _ in per_arm}) != 1:
        raise UnsupportedStructure("the units name differently once their stereo is considered")
    if len(stereo) != sum(len(entries) for _, entries in per_arm):
        raise UnsupportedStructure("stereodescriptors outside the units are not supported by the chain engine yet")
    ordered = sorted(per_arm, key=lambda arm: tuple(0 if code in "RZr" else 1 for _, code in arm[1]))
    labels = [
        f"{locant}{chr(39) * i}{code}"
        for i, (_, entries) in enumerate(ordered)
        for locant, code in entries
    ]
    return per_arm[0][0], f"({','.join(labels)})-"


def _unit_molecule(mol, atoms, attach):
    editable = Chem.RWMol(mol)
    root = editable.GetAtomWithIdx(attach)
    root.SetAtomMapNum(1)
    cut = sum(1 for n in root.GetNeighbors() if n.GetIdx() not in atoms)
    if root.GetIsAromatic() and cut:
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


def _is_nitro_part(atom):
    """The charged atoms of a nitro, azido or isocyano group: N+ bonded to two oxygens (one O-)."""
    if _is_azide_part(atom) or is_halogen_oxo_part(atom.GetOwningMol(), atom):
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
    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        return len(oxygens) in (2, 3) and sum(o.GetFormalCharge() for o in oxygens) == -1 and atom.GetDegree() == 3
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
        if not oxo:
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
            and mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(others) != 1 or others[0].GetAtomicNum() == 7:
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


def _is_ester_like(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    if not _double_oxygens(mol, carbon):
        return False
    return any(
        n.GetAtomicNum() in (8, 7, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not _diacyl_chalcogen_chain(mol, carbon, n)
        and not (_acyloxy_amine_oxygen(mol, n, carbon) and any(c.GetAtomicNum() == 6 for c in atom.GetNeighbors()))
        and not _acyloxy_lambda_centre_oxygen(n, carbon)
        and not (n.GetAtomicNum() == 7 and (atom.GetTotalNumHs() == 1 or any(c.GetAtomicNum() == 6 for c in atom.GetNeighbors())) and _ring_nitrogen_acyl(mol, n, carbon))
        and not (n.GetAtomicNum() == 8 and _terminal_heteroatom(mol, n.GetIdx(), 1))
        and not (n.GetAtomicNum() == 7 and _terminal_heteroatom(mol, n.GetIdx(), 2))
        and not (n.GetAtomicNum() == 7 and _acyl_diazene_nitrogen(mol, n, carbon))
        for n in atom.GetNeighbors()
    )


def _terminal_chalcogen_hydride(atom):
    """The sulfur, selenium or tellurium of an -SH, -SeH or -TeH group on nitrogen (thiohydroxylamine analogues, P-68.3.1.1.1.6)."""
    return atom.GetAtomicNum() in (16, 34, 52) and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1 and not atom.GetFormalCharge()


def _substituted_amine_nitrogen(mol, atom):
    if any(n.GetAtomicNum() in MONONUCLEAR_HYDRIDES for n in atom.GetNeighbors()):
        return False
    if atom.GetAtomicNum() != 7 or (atom.GetFormalCharge() and not (AMINIUM.get() and atom.GetFormalCharge() == 1)) or atom.GetIsAromatic() or atom.IsInRing():
        return False
    if any(
        n.GetAtomicNum() not in (6, 8) and not is_oxo_nitrogen(mol, n) and not _terminal_chalcogen_hydride(n)
        for n in atom.GetNeighbors()
    ):
        return False
    if any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds()):
        return False
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(carbons) >= 2 or (len(carbons) == 1 and atom.GetDegree() >= 2)


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


def _amidrazone_group(mol, atom):
    """("hydrazonamide" | "imidohydrazide", owned nitrogens) of a -C(-N)=N-N or -C(=N)-N-N group on a carbon with at
    most one carbon neighbour (P-66.4.2.1); the owned atoms are the three nitrogens."""
    carbon = atom.GetIdx()
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    others = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 7]
    if len(nitrogens) != 2 or len(others) > 1 or any(n.GetAtomicNum() != 6 for n in others):
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


def _amidine_nitrogens(mol, atom):
    """(imino N, amino N) of a -C(=NR)-NR2 group (P-66.4.1.1) whose carbon has at most one carbon neighbour; the
    imino nitrogen carries hydrogen, a carbon group or an oxygen (an amidoxime, P-66.4.4)."""
    carbon = atom.GetIdx()
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    others = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 7]
    if len(nitrogens) != 2 or len(others) > 1 or any(n.GetAtomicNum() != 6 for n in others):
        return None
    imino = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    amino = [n for n in nitrogens if mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0]
    if len(imino) != 1 or len(amino) != 1:
        return None
    imino, amino = imino[0], amino[0]
    if any(n.GetFormalCharge() or n.IsInRing() or n.GetIsAromatic() for n in (imino, amino)):
        return None
    for nitrogen, allowed in ((imino, (6, 8, 9, 17, 35, 53)), (amino, (6, 9, 17, 35, 53))):
        for n in nitrogen.GetNeighbors():
            if n.GetIdx() == carbon:
                continue
            if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
                return None
            if n.GetAtomicNum() not in allowed or n.GetFormalCharge():
                return None
            if is_functional_carbon(mol, n.GetIdx()) and not (nitrogen is amino and _imidoyl_carbon(mol, n, nitrogen)):
                return None
    if imino.GetDegree() > 2 or amino.GetDegree() > 3:
        return None
    return imino.GetIdx(), amino.GetIdx()


def _acyclic_imine_nitrogen(mol, atom):
    """The =N- atom of a C=N group (P-62.3) whose carbon has only carbon neighbours and whose nitrogen carries
    hydrogen, a carbon group or an oxygen (an oxime or its ether, P-68.3.1.1.2); None otherwise."""
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
        n.GetAtomicNum() != 6 and not _nitro_or_nitroso(mol, n)
        for n in atom.GetNeighbors()
        if n.GetIdx() != nitrogen.GetIdx()
    ):
        return None
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == carbon:
            continue
        if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
        if n.GetAtomicNum() not in ((6,) if cationic else (6, 8)) or n.GetFormalCharge():
            return None
    return nitrogen.GetIdx()


def _nitro_or_nitroso(mol, atom):
    """The nitrogen of a nitro or nitroso group, which makes an oxime a nitrolic or nitrosolic acid (P-68.3.1.1.3)."""
    return atom.GetAtomicNum() == 7 and (
        is_nitro_nitrogen(mol, atom.GetIdx())
        or (atom.GetDegree() == 2 and any(n.GetAtomicNum() == 8 and n.GetDegree() == 1 for n in atom.GetNeighbors()))
    )


def _ring_imine_nitrogen(mol, nitrogen, ring_atom):
    """A neutral =N-H, =N-R, =N-OH or =N-OR on a saturated ring carbon (P-62.3, P-68.3.1.1)."""
    if nitrogen.GetFormalCharge() or nitrogen.GetDegree() > 2:
        return False
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == ring_atom:
            continue
        if n.GetAtomicNum() == 8:
            if n.GetFormalCharge() or n.GetDegree() > 2:
                return False
        elif n.GetAtomicNum() != 6:
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
        if not atom.IsInRing():
            continue
        r = atom.GetIdx()
        for n in atom.GetNeighbors():
            if n.IsInRing():
                continue
            z, i = n.GetAtomicNum(), n.GetIdx()
            order = mol.GetBondBetweenAtoms(r, i).GetBondTypeAsDouble()
            if z == 8 and order == 2.0 and n.GetDegree() == 1:
                found.append(("ketone", r, {i}))
            elif z in _CHALCOGEN_KETONE_CLASS and order == 2.0 and n.GetDegree() == 1 and not n.GetFormalCharge():
                found.append((_CHALCOGEN_KETONE_CLASS[z], r, {i}))
            elif z == 8 and _terminal_heteroatom(mol, i, 1):
                found.append(("alcohol", r, {i}))
            elif z == 8 and order == 1.0 and (hydroperoxy := _hydroperoxy_oxygen(mol, n, r)) is not None:
                found.append(("peroxol", r, {i, hydroperoxy}))
            elif z in (16, 34, 52) and _terminal_heteroatom(mol, i, 1):
                found.append(({16: "thiol", 34: "selenol", 52: "tellurol"}[z], r, {i}))
            elif z in (16, 34, 52) and order == 1.0:
                sulfonyl = _sulfonyl_group(mol, i, r)
                if sulfonyl is not None:
                    found.append((sulfonyl[0], r, sulfonyl[1]))
            elif z == 7 and _terminal_heteroatom(mol, i, 2):
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
            entry[1] = min(entry[1], name)
    primes = {}
    if group_count > 1:
        by_position = {}
        for group, (count, first, located) in groups.items():
            by_position.setdefault(_n_position(located, position_of), []).append((-count, first, group))
        for members in by_position.values():
            for prime, (_, _, group) in enumerate(sorted(members)):
                primes[group] = prime
    for name, compound, *locant in n_names:
        located = locant[0] if locant else "N"
        if isinstance(located, tuple):
            position = _n_position(located, position_of)
            if position is None:
                continue
            if group_count == 1:
                located = "N"
            else:
                mark = "'" * primes.get(located[2], 0) if len(located) > 2 else ""
                if len(located) > 5:
                    mark *= located[5]
                elif len(located) > 4:
                    mark *= 2
                located = f"{located[0]}{mark}{position}"
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
    """N-prefix names when the single imidic or hydrazonic acid group carries substituents on its =N atom."""
    carbon_centre = spec_from_key(cls).center == "C"
    centers = [
        carbon if carbon_centre else next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)), carbon)
        for carbon, owned in groups.get(cls, {}).items()
    ]
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            centers.append(
                next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in ((6,) if carbon_centre else (6, 16, 34, 52)))
            )
    found = [acid_group_at(mol, c) for c in dict.fromkeys(centers)]
    centers = list(dict.fromkeys(centers))
    substituted = [g for g in found if g is not None and g.n_substituents]
    if not substituted:
        return []
    if len(centers) != 1:
        raise UnsupportedStructure("several imidic acid groups with N-substitution are not handled by the chain engine")
    group = substituted[0]
    return [
        name_branch(graph, n, group.n_atom, halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for n in group.n_substituents
    ]


def _imidamide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls):
    """N, N' and N'' prefix entries of a sulfonimidamide-type group: N on the amino nitrogen, N' and N'' on the
    imido nitrogens (P-66.4.1.1); several groups cannot be told apart by these locants."""
    members = dict(groups.get(cls, {}))
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            members[next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)] = owned
    entries = []
    for carbon, owned in members.items():
        centre = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52))
        nitrogens = [a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
        amino = next(a for a in nitrogens if mol.GetBondBetweenAtoms(centre, a).GetBondTypeAsDouble() == 1.0)
        cited = {
            a: [
                name_branch(graph, n, a, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                for n in graph[a]
                if n != centre
            ]
            for a in nitrogens
        }
        imido = sorted(
            (a for a in nitrogens if a != amino),
            key=lambda a: (-len(cited[a]), sorted(alpha_sort_key(name) for name, _ in cited[a]), a),
        )
        for nitrogen, locant in [(amino, "N"), *zip(imido, ("N'", "N''"))]:
            entries.extend((name, compound, locant) for name, compound in cited[nitrogen])
    if entries and len(members) != 1:
        raise UnsupportedStructure("several imidamide groups with N-substitution are not handled by the chain engine")
    return entries


def _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls="hydrazide"):
    """N and N' prefix entries (name, compound, locant) of the hydrazide groups: N on the acylated nitrogen, N' on
    the terminal nitrogen; with several groups the locant carries the position of its group (P-66.3.3.2)."""
    members = {c: (owned, c) for c, owned in groups.get(cls, {}).items()}
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            key = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)
            members[key] = (owned, ring_atom)
    entries = []
    for carbon, (owned, anchor) in members.items():
        acyl = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)), carbon)
        alpha = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and acyl in graph[a])
        beta = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and a != alpha)
        rest = tuple(a for a in graph[carbon] if a not in owned)
        for nitrogen, partner, base in ((alpha, beta, "N"), (beta, alpha, "N'")):
            located = (base, anchor, carbon, rest, "paired") if len(members) > 1 else base
            for n in graph[nitrogen]:
                if n in (partner, acyl):
                    continue
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
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
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
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
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
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
        located = ("N", anchor, carbon, tuple(a for a in graph[carbon] if a not in owned))
        entries += [
            (*name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True), located)
            for n in subs
        ]
    return entries


def _amide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls="amide"):
    """N-prefix names when the single amide group is N-substituted, else []."""
    members = {c: (owned, c) for c, owned in groups.get(cls, {}).items()}
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            key = next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)
            members[key] = (owned, ring_atom)
    entries = []
    for carbon, (owned, anchor) in members.items():
        nitrogen = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
        center = next((a for a in graph[nitrogen] if a in owned), carbon)
        located = ("N", anchor, carbon, tuple(a for a in graph[carbon] if a not in owned)) if len(members) > 1 else "N"
        for n in graph[nitrogen]:
            if n != center:
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
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
    uncited_locants = length > 1 and not force and omits_all_locants(mol, chain_set | nitrogens, with_n, owned - nitrogens)
    if completely_substituted and not force and not n_names:
        prefix = single_site_prefixes(grouped)
    else:
        prefix = format_substituent_prefixes(with_n, omit_locants=length == 1 and not force and not n_names)
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
            prefix, body, tail = "", retained, ""
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
            prefix, body, tail = "", retained, ""
        else:
            body, tail = name_from_substituents(length, ene, yne, multiplied_word(count, "oic")), " acid"
            if attach is not None and length == 2 and count == 1 and not ene and not yne:
                body = "acetic"
    elif principal in ("peroxoic", "thioic", "imidic"):
        body, tail = name_from_substituents(length, ene, yne, multiplied_word(count, principal)), " acid"
    elif principal == "peroxol":
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, "peroxol"), suffix_locants, substituted=bool(grouped)
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
        word = {
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
        }[principal]
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, word), suffix_locants, force_own_locant=force,
            substituted=bool(grouped),
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
        _parent_configuration_rank(stereo, position_of),
        -total_count,
        locant_set,
        citation,
        _bonding_rank(mol, graph, chain_set, position_of),
        _n_group_positions(n_names, position_of),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (prefix, body, tail, reported_attach, position_of, False)


def _stereo_free(mol):
    return not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) and not any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
