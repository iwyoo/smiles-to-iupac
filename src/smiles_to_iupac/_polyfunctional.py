"""Acyclic chain parents carrying any mix of supported groups (P-41, P-44.1.1,
P-44.3, P-45): the most senior class present becomes the suffix and every
other group, ring or branch is cited as a substituent prefix through
`name_branch` (P-29.3.3, P-29.4). Rings may only be substituents; a principal
group on a ring is not handled here.
"""

import contextvars
from types import MappingProxyType
import re

from rdkit import Chem
from rdkit.Chem import CanonicalRankAtoms

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    ring_cycle,
    specified_stereo_elements,
    substituent_locant_set_and_citation,
)
from ._anion import ANION_PROP, anion_weight
from ._functional_prefixes import is_nitro_nitrogen
from ._hetero_prefixes import EXTENDED_PREFIXES, MONONUCLEAR_HYDRIDES, is_functional_carbon
from ._multiplicative import _bare_key
from ._multiplicative_text import PrimedLocant, enclose, unit_phrase
from ._multiplicative_ring import (
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
from ._acid_groups import acid_group_at
from ._acid_lexicon import carbo_suffix, chain_suffix, make_spec, rank_key, spec_from_key
from ._retained_acids import retained_chain_acid
from ._substituents import format_substituent_prefixes, name_branch

_SENIORITY = [
    "ide", "acid", "thioic", "peroxoic", "imidic", "sulfonic", "amide", "amidine", "sulfonamide", "hydrazide", "nitrile", "aldehyde", "ketone", "thione", "selone", "tellone", "alcohol", "peroxol",
    "thiol", "selenol", "tellurol", "amine", "imine",
]
_TERMINAL = {"acid", "thioic", "peroxoic", "imidic", "amide", "amidine", "hydrazide", "nitrile", "aldehyde"}
_PLAIN_ACID = make_spec("C", ["O"], ["O"])
_PLAIN_SULFONIC = make_spec("S", ["O", "O"], ["O"])


def _is_variant(name):
    return bool(name) and name.startswith("acid:")


def _acid_rank(name):
    spec = {"acid": _PLAIN_ACID, "sulfonic": _PLAIN_SULFONIC}.get(name) or spec_from_key(name)
    return rank_key(spec)


AMINIUM = contextvars.ContextVar("aminium", default=False)


def _principal_class(classes):
    if AMINIUM.get():
        return "amine" if "amine" in classes else None
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
_CHALCOGEN_KETONE_OK = {
    "acid", "thioic", "peroxoic", "imidic", "sulfonic", "amide", "amidine", "sulfonamide", "hydrazide", "nitrile", "aldehyde", "ketone",
    *_CHALCOGEN_KETONE_CLASS.values(),
}


def _chalcogen_ketone(mol, atom):
    """A carbon double-bonded to S, Se or Te (thione, selone, tellone)."""
    if atom.GetAtomicNum() != 6:
        return False
    return any(
        n.GetAtomicNum() in (16, 34, 52)
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )


def _plain_amide_nitrogen(mol, nitrogen, carbonyl):
    """An amide nitrogen carrying only carbon substituents or one hydroxy (a
    hydroxamic acid, P-65.1.3.4; no acyl group, so not an imide) -- named with
    'N-' prefixes on the amide parent."""
    if nitrogen.GetFormalCharge() or nitrogen.IsInRing() or nitrogen.GetIsAromatic():
        return False
    others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbonyl]
    hydroxy = [n for n in others if _terminal_heteroatom(mol, n.GetIdx(), 1) and n.GetAtomicNum() == 8]
    return bool(others) and len(hydroxy) <= 1 and all(
        n in hydroxy
        or (
            n.GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            and not _double_oxygens(mol, n.GetIdx())
            and not is_functional_carbon(mol, n.GetIdx())
        )
        for n in others
    )


def _ring_nitrogen_acyl(mol, nitrogen, carbonyl):
    """A neutral saturated ring nitrogen whose only acyl group is `carbonyl`:
    the carbonyl is a pseudoketone ('hidden amide', P-64.1.2.1(b), P-64.3.2)."""
    if nitrogen.GetFormalCharge() or nitrogen.GetIsAromatic() or not nitrogen.IsInRing() or nitrogen.GetDegree() != 3:
        return False
    return all(
        mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        and (n.GetIdx() == carbonyl or (n.GetAtomicNum() == 6 and not is_functional_carbon(mol, n.GetIdx())))
        for n in nitrogen.GetNeighbors()
    )


def _sulfonyl_group(mol, s_idx, attached):
    """("sulfonic" | "sulfonamide", owned atoms) for an S(=O)(=O)X group whose
    X is OH or an amine nitrogen and whose other neighbor is `attached`."""
    sulfur = mol.GetAtomWithIdx(s_idx)
    group = acid_group_at(mol, s_idx)
    if group is not None and any(n.GetIdx() == attached for n in sulfur.GetNeighbors()):
        if group.spec.kind[0] in ("S", "Se", "Te") and not (group.spec.kind == ("S", 2) and group.spec.plain):
            return group.spec.key, group.owned | {s_idx}
    if sulfur.GetAtomicNum() != 16 or sulfur.GetFormalCharge() or sulfur.GetDegree() != 4:
        return None
    oxygens = [
        n.GetIdx()
        for n in sulfur.GetNeighbors()
        if n.GetAtomicNum() == 8
        and n.GetDegree() == 1
        and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in sulfur.GetNeighbors() if n.GetIdx() not in oxygens and n.GetIdx() != attached]
    if len(oxygens) != 2 or len(rest) != 1:
        return None
    other = rest[0]
    if other.GetAtomicNum() == 8 and _terminal_heteroatom(mol, other.GetIdx(), 1):
        return "sulfonic", {s_idx, *oxygens, other.GetIdx()}
    if other.GetAtomicNum() == 7 and (
        _terminal_heteroatom(mol, other.GetIdx(), 2) or _plain_amide_nitrogen(mol, other, s_idx)
    ):
        return "sulfonamide", {s_idx, *oxygens, other.GetIdx()}
    return None


def _group_of(mol, carbon):
    """(class, atoms owned by the group) for a principal-capable group on
    `carbon`, else None. Raises on carbon-bound groups this engine cannot
    name (esters, acid halides, ...)."""
    atom = mol.GetAtomWithIdx(carbon)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
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
        others = [n for n in atom.GetNeighbors() if n.GetIdx() != nitrogens[0]]
        if len(others) != 1 or others[0].GetAtomicNum() != 6:
            raise UnsupportedStructure("a cyanide not bonded to carbon is not a nitrile")
        return "nitrile", {nitrogens[0]}
    for z, thione in _CHALCOGEN_KETONE_CLASS.items():
        found = [
            n.GetIdx()
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == z and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if found and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors() if n.GetIdx() != found[0]) and atom.GetDegree() == 3:
            return thione, {found[0]}
    oxygens = _double_oxygens(mol, carbon)
    if oxygens:
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        hetero = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 6 and n.GetIdx() != oxygens[0]]
        if any(
            mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() != 1.0
            for n in atom.GetNeighbors()
            if n.GetIdx() != oxygens[0]
        ):
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
            if other.GetAtomicNum() == 7 and len(carbon_neighbors) <= 1:
                beta = _hydrazide_beta_nitrogen(mol, other, carbon)
                if beta is not None:
                    return "hydrazide", {oxygens[0], other.GetIdx(), beta}
            if carbon_neighbors and other.GetAtomicNum() == 7 and _ring_nitrogen_acyl(mol, other, carbon):
                return "ketone", {oxygens[0]}
        return None
    amidine = _amidine_nitrogens(mol, atom)
    if amidine is not None:
        return "amidine", set(amidine)
    imine = _acyclic_imine_nitrogen(mol, atom)
    if imine is not None:
        return "imine", {imine}
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() == 16 and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0:
            sulfonyl = _sulfonyl_group(mol, n.GetIdx(), carbon)
            if sulfonyl is not None:
                return sulfonyl
    for z, hydrogens, name in ((8, 1, "alcohol"), (16, 1, "thiol"), (34, 1, "selenol"), (52, 1, "tellurol"), (7, 2, "amine")):
        for n in _single_neighbors(mol, carbon, z):
            if _terminal_heteroatom(mol, n, hydrogens):
                return name, {n}
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

    if has_linear_phane_shape(mol):
        return name_linear_phane(mol)
    split = split_isotopes(mol)
    if split is not None:
        return _name_isotopic(*split)
    cation = _aminium_base(mol)
    if cation is not None:
        return _name_aminium(cation)
    return _name_labelled(mol, {})


def _aminium_base(mol):
    """The mol with its ammonium nitrogen neutralised when the only charge is one N+ bonded to carbon and hydrogen
    (a cation outranks every acid, P-41), else None."""
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 1 or charged[0].GetFormalCharge() != 1 or charged[0].GetAtomicNum() != 7:
        return None
    nitrogen = charged[0]
    if (
        nitrogen.GetIsAromatic()
        or nitrogen.IsInRing()
        or nitrogen.GetDegree() + nitrogen.GetTotalNumHs() != 4
        or any(
            n.GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
            for n in nitrogen.GetNeighbors()
        )
    ):
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


def _name_aminium(base):
    token = AMINIUM.set(True)
    try:
        name = _name_labelled(base, {})
    finally:
        AMINIUM.reset(token)
    if not name.endswith(("amine", "aniline")):
        raise UnsupportedStructure("the ammonium cation is not named as an amine parent")
    return name[:-1] + "ium"


def _name_isotopic(clean, labels):
    if specified_stereo_elements(clean) is not None or any(
        a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms()
    ):
        raise UnsupportedStructure("stereodescriptors of an isotopically modified compound are not supported yet (P-82.4)")
    return _name_labelled(clean, labels)


def _name_labelled(mol, labels):
    from ._isotope_labels import descriptor
    from ._substituents import BRANCH_STEREO, ISOTOPE_LABELS

    stereo = _check_scope(mol) + [("isotope", atom, "") for atom in labels]
    context = {
        "atoms": {where: code for kind, where, code in stereo if kind == "atom"},
        "bonds": {where: code for kind, where, code in stereo if kind == "bond"},
        "used": set(),
    }
    token = BRANCH_STEREO.set(context if stereo else None)
    isotope_context = {"labels": labels, "consumed": set()}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if labels else None)
    try:
        multiplicative = _multiplicative_name(mol, stereo)
        if multiplicative is not None:
            if labels:
                raise UnsupportedStructure("isotopic modification of a multiplicative name is not supported yet")
            return multiplicative
        _, name, parts = _select(mol, stereo=stereo)
        if labels:
            in_parent = {a: e for a, e in labels.items() if a in parts[4]}
            if set(labels) - set(in_parent) - isotope_context["consumed"]:
                raise UnsupportedStructure("an isotopically modified atom outside the parent hydride is not supported yet")
            if in_parent:
                name = _with_descriptor(name, parts, descriptor(in_parent, parts[4], len(parts[4]) == 1))
        return _stereo_prefix(stereo, parts[4], ring_parent=parts[5], used=context["used"]) + name
    finally:
        BRANCH_STEREO.reset(token)
        ISOTOPE_LABELS.reset(isotope_token)


def _check_scope(mol):
    """The specified stereo elements as ("atom", idx, code) or
    ("bond", (a, b), code); [] when the molecule has none."""
    if mol.GetNumAtoms() > _MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("this molecule is out of scope for the polyfunctional chain engine")
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
                if ("bond", where) not in used:
                    complete = False
            else:
                entries.append((min(position_of[a], position_of[b]), code))
    return sorted(entries), complete


def _stereo_rank(stereo, position_of, ring_parent=False):
    entries, _ = _stereo_entries(stereo, position_of, ring_parent)
    isotopic = tuple(sorted(position_of[w] for k, w, _ in stereo or [] if k == "isotope" and w in position_of))
    return isotopic, tuple(0 if code in "RZr" else 1 for _, code in entries)


def _stereo_prefix(stereo, position_of, ring_parent=False, used=frozenset()):
    if not stereo:
        return ""
    entries, complete = _stereo_entries(stereo, position_of, ring_parent, used)
    if not complete:
        raise UnsupportedStructure("stereodescriptors outside the parent are not supported by the chain engine yet")
    if not entries:
        return ""
    return "(" + ",".join(f"{locant}{code}" for locant, code in entries) + ")-"


def _select(mol, attach=None, n_names=(), stereo=None):
    classes = {_group_of(mol, a.GetIdx())[0] for a in mol.GetAtoms() if _group_of(mol, a.GetIdx())} | {
        c for c, _, _ in _ring_occurrences(mol)
    }
    principal = _principal_class(classes)
    token = EXTENDED_PREFIXES.set(_is_acid_family(principal) if principal else False)
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
        return _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal, attach, n_names, stereo)
    finally:
        IDE_EXTRA.reset(token)


def _select_with_principal(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal, attach, n_names, stereo):
    if any(_chalcogen_ketone(mol, a) for a in mol.GetAtoms()) and principal not in _CHALCOGEN_KETONE_OK and not (principal and _is_variant(principal)):
        raise UnsupportedStructure("a thioketone-type group outranks the parents this engine can build here")
    for atom in mol.GetAtoms():
        if atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise UnsupportedStructure("isotopes and radicals are not supported by the polyfunctional chain engine")
        if atom.GetFormalCharge() and not _is_nitro_part(atom) and not _anionic_group_atom(mol, atom) and not (
            AMINIUM.get() and atom.GetAtomicNum() == 7
        ):
            raise UnsupportedStructure("charged atoms are not supported by the polyfunctional chain engine")
        if (
            atom.GetAtomicNum() == 6
            and not atom.IsInRing()
            and not _is_acid_family(principal)
            and not (AMINIUM.get() and principal in (None, "amine"))
            and principal not in ("peroxoic", "thioic", "imidic")
            and not (principal in ("amide", "sulfonamide", "hydrazide") and atom.GetIdx() in groups.get(principal, {}))
            and _is_ester_like(mol, atom.GetIdx())
        ):
            raise UnsupportedStructure("an ester outranks every parent this engine can build except an acid")
    if principal in (None, "amine") and any(_substituted_amine_nitrogen(mol, a) for a in mol.GetAtoms()):
        if attach is not None or n_names:
            raise UnsupportedStructure("N-substituted amines inside a unit are not handled by the chain engine")
        if stereo:
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

    if principal in ("amide", "sulfonamide"):
        amide_ns = _amide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, principal)
        if amide_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted amide inside a unit is not handled by the chain engine")
            n_names = amide_ns

    if principal == "hydrazide":
        hydrazide_ns = _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups)
        if hydrazide_ns:
            if attach is not None or n_names:
                raise UnsupportedStructure("an N-substituted hydrazide inside a unit is not handled by the chain engine")
            n_names = hydrazide_ns

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
    total_principal = len(anchors)
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
                or (_is_acid_family(principal) and _junior_end_group(mol, a.GetIdx()))
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
        carbo = _carbo_best(mol, graph, halogens, aromatic_atoms, principal, principal_atoms, stereo)
        if carbo is not None and -carbo[0][0] > chain_count:
            return _finish(carbo)
    if (
        _is_terminal(principal)
        and len(chain_best[2][4]) == 1
        and not (_is_variant(principal) and attach is None)
        and not chain_best[1].endswith("formic acid")
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
    spec = monocycle_spec(mol, ring)
    if spec is None:
        from ._diester_ring_diyl import evaluate_skeleton

        found = evaluate_skeleton(mol, graph, "ring", [ring], ring_set, [], set(), "")
        if found is None:
            raise UnsupportedStructure("this ring has no supported name")
        placed = found[2]
        name = _without_stereo(found[1])
        return (-len(roots), tuple(sorted(placed[r] for r, _ in roots)), name), ((0,), name, (None, None, None, 0, placed, True))
    if not roots:
        raise UnsupportedStructure("an unsubstituted ring is not a polyfunctional case")
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
    if len(entries) == 1 and spec.hetero is None:
        _, only_name, only_compound = entries[0]
        prefix_text = format_substituent_prefixes(
            {only_name: {"locants": [1], "compound": only_compound}}, omit_locants=True
        )
    else:
        prefix_text = _prefix_text(entries, best[1])
    name = _join(prefix_text, spec.parent)
    return (-len(roots), best[0][0], name), ((0,), name, (None, None, None, 0, best[1], True))


_GROUP_14 = (14, 32, 50, 82)


def _mononuclear_parent(mol, graph, halogens, aromatic_atoms, center):
    """A single Si, Ge, P, B, ... atom is the senior parent hydride when there is no principal group (P-44.1.2):
    'trimethyl(phenyl)silane', 'methoxy(trimethyl)silane'. On a Group 14 atom a hydroxy or amino group is the
    suffix (P-68.2): 'trimethylsilanol', '1,1,1-trimethylsilanamine'."""
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
    hydroxyls = [n for n in neighbors if _terminal_heteroatom(mol, n, 1) and mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
    amines = [
        n
        for n in neighbors
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 7
        and not mol.GetAtomWithIdx(n).GetFormalCharge()
        and not mol.GetAtomWithIdx(n).IsInRing()
        and all(mol.GetAtomWithIdx(m).GetAtomicNum() in (6, *MONONUCLEAR_HYDRIDES) for m in graph[n] if m != index)
    ]
    others = [n for n in neighbors if n not in hydroxyls and n not in amines]
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() not in (6, 8, *HALOGEN_PREFIXES) for n in others):
        return None
    if (hydroxyls or amines) and (z not in _GROUP_14 or (hydroxyls and amines) or len(amines) > 1):
        return None
    entries = [name_branch(graph, n, index, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in others]
    if hydroxyls:
        return format_mononuclear_prefixes(entries) + (
            stem[:-1] + "ol" if len(hydroxyls) == 1 else stem + multiplied_word(len(hydroxyls), "ol")
        )
    if amines:
        (nitrogen,) = amines
        grouped = group_substituents({1: entries} if entries else {})
        n_entries = [
            name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
            for n in graph[nitrogen]
            if n != index
        ]
        merged = _with_n_names(grouped, n_entries)
        return format_substituent_prefixes(merged) + stem[:-1] + "amine"
    return format_mononuclear_prefixes(entries) + stem


def _plain_parent(mol, graph, halogens, aromatic_atoms, stereo=None):
    """Parent without a principal group: the monocycle when there is one
    (P-44.1.2.2), else the longest chain, with every substituent a prefix."""
    from ._hydride_chain import name_hydride_chain

    chain_name = name_hydride_chain(mol, graph, halogens, aromatic_atoms)
    if chain_name is not None:
        return ((0,), chain_name, (None, None, None, 0, {}, False))
    centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not a.IsInRing()]
    if centers:
        named = [n for n in (_mononuclear_parent(mol, graph, halogens, aromatic_atoms, c) for c in centers) if n]
        if named:
            return ((0,), min(named), (None, None, None, 0, {}, False))
    if centers or any(
        b.GetBeginAtom().GetAtomicNum() == 7 and b.GetEndAtom().GetAtomicNum() == 7 and not b.IsInRing()
        for b in mol.GetBonds()
    ):
        raise UnsupportedStructure("a heteroatom hydride is the senior parent when there is no principal group (P-44.1.2.2)")
    ring_info = mol.GetRingInfo()
    rings = [r for r in ring_info.AtomRings()]
    if len(rings) >= 2:
        assembly = _assembly_parent(mol, graph, halogens, aromatic_atoms, None, [], stereo)
        if assembly is not None:
            return assembly[1]
    if rings and any(ring_info.NumAtomRings(a) != 1 for r in rings for a in r):
        return _fused_plain_parent(mol, graph, rings)
    if rings:
        if any(ring_info.NumAtomRings(a) != 1 for r in rings for a in r):
            raise UnsupportedStructure("several rings without a principal group are not named by the chain engine")
        if len(rings) >= 4:
            raise UnsupportedStructure("four or more rings need a phane or ring-assembly name")
        ranked = sorted(rings, key=lambda r: tuple(-x for x in _ring_rank(mol, r[0])))
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


def _system_rank(mol, system_rings, system_atoms):
    hetero = [mol.GetAtomWithIdx(a).GetSymbol() for a in system_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    ranks = [_HETERO_RANK.get(h, 99) for h in hetero]
    return (bool(hetero), "N" in hetero, -min(ranks, default=0), len(system_rings), len(system_atoms), len(hetero))


def _fused_plain_parent(mol, graph, rings):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    systems = []
    for ring in rings:
        system_rings, system_atoms = _system_of(mol, ring[0])
        if not any(set(system_atoms) == set(seen[1]) for seen in systems):
            systems.append((system_rings, system_atoms))
    ranked = sorted(systems, key=lambda s: _system_rank(mol, *s), reverse=True)
    if len(ranked) > 1 and _system_rank(mol, *ranked[0]) == _system_rank(mol, *ranked[1]):
        raise UnsupportedStructure("several equally senior ring systems need a multiplicative or assembly name")
    system_rings, system_atoms = ranked[0]
    if len(system_rings) > 1:
        _require_mancude_system(mol, system_atoms)
    found = evaluate_skeleton(mol, graph, "ring", system_rings, system_atoms, [], set(), "")
    if found is None:
        raise UnsupportedStructure("this fused ring system has no supported numbering")
    return ((0,), _without_stereo(found[1]), (None, None, None, 0, found[2], True))


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
    if not entries:
        raise UnsupportedStructure("an unsubstituted chain is not a polyfunctional case")
    grouped = group_substituents(entries)
    locant_set, total_count, citation = substituent_locant_set_and_citation(grouped)
    length = len(chain)
    single = total_count == 1
    prefix = format_substituent_prefixes(grouped, omit_locants=length == 1 or (length == 2 and (ene or yne) and single))
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
        -total_count,
        locant_set,
        citation,
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (None, None, None, 0, position_of, False)


_CARBO_WORDS = {
    "acid": "carboxylic acid",
    "amide": "carboxamide",
    "amidine": "carboximidamide",
    "hydrazide": "carbohydrazide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
}


def _carbo_best(mol, graph, halogens, aromatic_atoms, principal, principal_atoms, stereo):
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
                mol, graph, halogens, aromatic_atoms, chain, principal, group_carbons, owned, stereo
            )
            if candidate is not None and (best is None or candidate[0] < best[0]):
                best = candidate
    return best


def _evaluate_carbo(mol, graph, halogens, aromatic_atoms, chain, principal, group_carbons, owned, stereo):
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
    prefix = format_substituent_prefixes(grouped)
    body = name_from_substituents(length, ene, yne, multiplied_word(count, word), suffix_locants)
    name = prefix + body
    key = (
        -count,
        -length,
        tuple(suffix_locants),
        -(len(ene) + len(yne)),
        lowest_locant_set(ene + yne),
        -total_count,
        locant_set,
        citation,
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
    unsaturated = bool(specs) and specs[0].kind == "cycloalkene"
    for unprimed in (0, 1):
        for first in orientations(unprimed):
            for second in orientations(1 - unprimed):
                locants = {atom: (0, number) for atom, number in first.items()}
                locants.update({atom: (1, number) for atom, number in second.items()})
                ene = ()
                if unsaturated:
                    ene = (tuple(multiple_locants(specs[unprimed], first)[0]), tuple(multiple_locants(specs[1 - unprimed], second)[0]))
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
    if unsaturated and any(loc == len(rings[0]) for loc in best[2][0] + best[2][1]):
        raise UnsupportedStructure("a ring double bond closing the numbering (1(n) locant) is not supported in an assembly")
    if unsaturated and best[2][0] != best[2][1]:
        raise UnsupportedStructure("the rings of this assembly are not identical once numbered (P-28.7)")
    return best[1]


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
    if spec.kind == "cycloalkene":
        inside = {a: locants[a][1] for a in spec.cycle}
        parent = parent_text(spec, inside)
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
    if len(joins) != 1 or _bare_key(mol, set(own)) != _bare_key(mol, set(other)):
        return None
    specs = [spec_of(mol, r) for r in (own, other)]
    if any(sp is None or sp.kind == "pyrrole" for sp in specs) or specs[0].kind != specs[1].kind:
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
    base = _assembly_base(specs, locants, joins[0], elide=True, ylidene=ylidene)
    spot = locants[root]
    core = f"[{base}]-{spot[1]}{chr(39) * spot[0]}-yl"
    return (f"{prefix}-{core}" if prefix else core), True


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
    if len({_bare_key(mol, set(r)) for r in all_rings}) != 1:
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
            key = (-found[0], not ylidene, found[1][1])
            if best is None or key < best[0]:
                best = (key, (found[0], found[1]))
    return best[1] if best else None


def _pair_assembly(mol, graph, halogens, aromatic_atoms, principal, occurrences, rings):
    """The assembly name of two directly joined identical rings with every other part of `mol` cited as a
    substituent; None when the pair is not such an assembly or leaves a principal group outside it."""
    if set(rings[0]) & set(rings[1]):
        return None
    if occurrences and any(o[1] not in set(rings[0]) | set(rings[1]) for o in occurrences):
        return None
    specs = [spec_of(mol, r) for r in rings]
    if any(sp is None or sp.kind == "pyrrole" for sp in specs):
        return None
    if specs[0].kind != specs[1].kind or len(rings[0]) != len(rings[1]):
        return None
    if _bare_key(mol, set(rings[0])) != _bare_key(mol, set(rings[1])):
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
    if principal is None:
        core = _assembly_base(specs, locants, joins[0], elide=False, ylidene=ylidene)
    else:
        word = multiplied_word(count, _SUFFIX_WORDS[_RING_SUFFIX[principal]])
        base = _assembly_base(specs, locants, joins[0], elide=word[0] in "aeiouy", ylidene=ylidene)
        spots = ",".join(cite(locants[o[1]]) for o in sorted(occurrences, key=lambda o: _locant_order(locants[o[1]])))
        core = f"[{base}]-{spots}-{word}"
    name = f"{prefix}-{core}" if prefix else core
    return count, ((-count,), name, (None, None, None, 0, {a: PrimedLocant(*loc) for a, loc in locants.items()}, True)), joins[0]


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
    if len(member_rings) > 3 and any(len(r) not in (5, 6, 7) for r in member_rings):
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
    "amidine": "carboximidamide",
    "sulfonamide": "sulfonamide",
    "hydrazide": "carbohydrazide",
    "nitrile": "carbonitrile",
    "aldehyde": "carbaldehyde",
    "alcohol": "ol",
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
    if principal in ("ketone", *_CHALCOGEN_KETONE_CLASS.values()) and any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in attach):
        raise UnsupportedStructure("a ring-heteroatom oxide is not a ring ketone")
    found = evaluate_skeleton(mol, graph, "ring", rings, atoms, attach, blocked, _FUSED_SUFFIX[principal], n_names=n_names)
    if found is None:
        raise UnsupportedStructure("this fused ring system has no supported numbering")
    count = len(on_system)
    return count, ((-count,), _without_stereo(found[1]), (None, None, None, 0, found[2], True))


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
    "amidine": "amidine",
    "sulfonamide": "sulfonamide",
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
        return found[0], found[1], (0, (), (), found[1][1])
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
        and (spec.kind == "cycloalkane" or (spec.kind == "benzene" and suffix_name not in _RETAINED_BENZENE))
    ):
        word = _SUFFIX_WORDS[suffix_name]
        stem = spec.parent[:-1] if word[0] in "aeiouy" else spec.parent
        core = stem + word
    else:
        core, _ = _suffix_text(spec.parent, suffix_name, suffix_locants, spec)
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
    if len(amine_nitrogens) != 1 or primaries:
        raise UnsupportedStructure("several amine groups with N-substitution are not handled by the chain engine")
    nitrogen = amine_nitrogens[0]
    n_idx = nitrogen.GetIdx()
    if nitrogen.IsInRing() or nitrogen.GetTotalNumHs() > 1:
        raise UnsupportedStructure("a ring nitrogen is not an acyclic amine parent")
    neighbors = [n.GetIdx() for n in nitrogen.GetNeighbors()]
    if any(
        mol.GetAtomWithIdx(c).GetAtomicNum() != 6
        or _double_oxygens(mol, c)
        or is_functional_carbon(mol, c)
        or mol.GetBondBetweenAtoms(n_idx, c).GetBondTypeAsDouble() != 1.0
        for c in neighbors
    ):
        raise UnsupportedStructure("this nitrogen is not a plain amine nitrogen")
    arms = {c: _arm_atoms(graph, c, n_idx) for c in neighbors}
    if sum(len(a) for a in arms.values()) != len(set().union(*arms.values())) or n_idx in set().union(*arms.values()):
        raise UnsupportedStructure("a nitrogen closing a ring is not an acyclic amine parent")
    results = []
    for c in neighbors:
        others = [o for o in neighbors if o != c]
        n_names = [name_branch(graph, o, n_idx, halogens, aromatic_atoms, mol=mol, unsaturated=True) for o in others]
        parent, mapped = _amine_parent_molecule(mol, arms[c], c, n_idx)
        result = _select(parent, None, n_names)
        ring = parent.GetAtomWithIdx(mapped).IsInRing()
        results.append((ring, _ring_rank(parent, mapped), _chain_size(parent, mapped), result))
    results.sort(key=lambda r: (not r[0], tuple(-x for x in r[1]), -r[2], r[3][1]))
    return results[0][3]


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
            a.SetFormalCharge(0)
            a.SetNumExplicitHs(2)
            a.SetNoImplicit(True)
    Chem.SanitizeMol(parent)
    mapped = next(a.GetIdx() for a in parent.GetAtoms() if a.GetAtomMapNum() == 1)
    for a in parent.GetAtoms():
        a.SetAtomMapNum(0)
    return parent, mapped


def _ring_rank(mol, idx):
    ring = next((set(r) for r in mol.GetRingInfo().AtomRings() if idx in r), None)
    if ring is None:
        return (0, 0, 0)
    numbers = {mol.GetAtomWithIdx(a).GetAtomicNum() for a in ring}
    hetero_rank = 3 if 7 in numbers else 2 if 8 in numbers else 1 if numbers - {6} else 0
    ordered = [a for a in ring]
    unsaturation = sum(
        mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()
        for i, a in enumerate(ordered)
        for b in ordered[i + 1 :]
        if mol.GetBondBetweenAtoms(a, b) is not None
    )
    return (hetero_rank, len(ring), unsaturation)


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
                return True
    return False


_LINKER_WORDS = {8: "oxy", 16: "sulfanediyl", 34: "selanediyl"}
_DICHALCOGEN_WORDS = {8: "peroxy", 16: "disulfanediyl", 34: "diselanediyl"}


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
        _natural_product_linker_name(mol, graph, stereo)
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
    editable.GetAtomWithIdx(attach).SetAtomMapNum(1)
    for idx in sorted(set(range(mol.GetNumAtoms())) - atoms, reverse=True):
        editable.RemoveAtom(idx)
    unit = editable.GetMol()
    Chem.SanitizeMol(unit)
    attach_idx = next(a.GetIdx() for a in unit.GetAtoms() if a.GetAtomMapNum() == 1)
    canonical = Chem.Mol(unit)
    for a in canonical.GetAtoms():
        a.SetAtomMapNum(1 if a.GetIdx() == attach_idx else 0)
    return canonical, attach_idx


def _is_nitro_part(atom):
    """The charged atoms of a nitro group: N+ bonded to two oxygens (one O-)."""
    if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        return len(oxygens) == 2 and sum(o.GetFormalCharge() for o in oxygens) == -1 and atom.GetDegree() == 3
    if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (n,) = atom.GetNeighbors()
        return _is_nitro_part(n)
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


def _is_ester_like(mol, carbon):
    atom = mol.GetAtomWithIdx(carbon)
    if not _double_oxygens(mol, carbon):
        return False
    return any(
        n.GetAtomicNum() in (8, 7, 16, 9, 17, 35, 53) and mol.GetBondBetweenAtoms(carbon, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        and not (n.GetAtomicNum() == 7 and any(c.GetAtomicNum() == 6 for c in atom.GetNeighbors()) and _ring_nitrogen_acyl(mol, n, carbon))
        and not (n.GetAtomicNum() == 8 and _terminal_heteroatom(mol, n.GetIdx(), 1))
        and not (n.GetAtomicNum() == 7 and _terminal_heteroatom(mol, n.GetIdx(), 2))
        for n in atom.GetNeighbors()
    )


def _substituted_amine_nitrogen(mol, atom):
    if any(n.GetAtomicNum() in MONONUCLEAR_HYDRIDES for n in atom.GetNeighbors()):
        return False
    if atom.GetAtomicNum() != 7 or (atom.GetFormalCharge() and not (AMINIUM.get() and atom.GetFormalCharge() == 1)) or atom.GetIsAromatic() or atom.IsInRing():
        return False
    carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(carbons) >= 2


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
    for nitrogen, allowed in ((imino, (6, 8)), (amino, (6,))):
        for n in nitrogen.GetNeighbors():
            if n.GetIdx() == carbon:
                continue
            if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
                return None
            if n.GetAtomicNum() not in allowed or n.GetFormalCharge() or is_functional_carbon(mol, n.GetIdx()):
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
    if nitrogen.GetFormalCharge() or nitrogen.IsInRing() or nitrogen.GetIsAromatic() or nitrogen.GetDegree() > 2:
        return None
    if any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()):
        return None
    for n in nitrogen.GetNeighbors():
        if n.GetIdx() == carbon:
            continue
        if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            return None
        if n.GetAtomicNum() not in (6, 8) or n.GetFormalCharge():
            return None
    return nitrogen.GetIdx()


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
            elif z in (16, 34, 52) and _terminal_heteroatom(mol, i, 1):
                found.append(({16: "thiol", 34: "selenol", 52: "tellurol"}[z], r, {i}))
            elif z == 16 and order == 1.0:
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


def _with_n_names(grouped, n_names, position_of=None, group_count=2):
    """`grouped` plus the N-prefixes; a locant ("N", atom) of one of several groups reads N<position of atom>
    (P-66.1.1.4.2: 'N1,N5-dimethylpentanediamide')."""
    if not n_names:
        return grouped
    merged = {name: {"locants": list(info["locants"]), "compound": info["compound"]} for name, info in grouped.items()}
    for name, compound, *locant in n_names:
        located = locant[0] if locant else "N"
        if isinstance(located, tuple):
            if located[1] not in position_of:
                continue
            located = "N" if group_count == 1 else f"{located[0]}{position_of[located[1]]}"
        merged.setdefault(name, {"locants": [], "compound": compound})["locants"].append(located)
    return merged


def _n_group_positions(n_names, position_of):
    """Sorted numbers of the groups that carry N-prefixes, for choosing the numbering."""
    anchors = [entry[2][1] for entry in n_names if len(entry) > 2 and isinstance(entry[2], tuple)]
    return tuple(sorted(position_of[a] for a in anchors if a in position_of))


def _ring_prefix_text(entries, locants, n_names, group_count=2):
    grouped = {}
    for r, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locants[r])
    return format_substituent_prefixes(_with_n_names(grouped, n_names, locants, group_count)) if (grouped or n_names) else ""


def _imidic_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups, cls):
    """N-prefix names when the single imidic or hydrazonic acid group carries substituents on its =N atom."""
    centers = [
        next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (16, 34, 52)), carbon)
        for carbon, owned in groups.get(cls, {}).items()
    ]
    for group_cls, ring_atom, owned in ring_groups:
        if group_cls == cls:
            centers.append(next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() in (6, 16, 34, 52)))
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


def _hydrazide_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """N and N' prefix entries (name, compound, locant) of the hydrazide group: N on the acylated nitrogen, N' on
    the terminal nitrogen (P-66.3.3); several groups cannot be told apart by these locants."""
    members = dict(groups.get("hydrazide", {}))
    for cls, ring_atom, owned in ring_groups:
        if cls == "hydrazide":
            members[next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)] = owned
    entries = []
    for carbon, owned in members.items():
        alpha = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and carbon in graph[a])
        beta = next(a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7 and a != alpha)
        for nitrogen, partner, locant in ((alpha, beta, "N"), (beta, alpha, "N'")):
            for n in graph[nitrogen]:
                if n in (partner, carbon):
                    continue
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.append((name, compound, locant))
    if entries and len(members) != 1:
        raise UnsupportedStructure("several hydrazide groups with N-substitution are not handled by the chain engine")
    return entries


def _amidine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """N and N' prefix entries (name, compound, locant) of the amidine group: N on the amino nitrogen, N' on the
    imino nitrogen (P-66.4.1.1); several groups cannot be told apart by these locants."""
    members = dict(groups.get("amidine", {}))
    for cls, ring_atom, owned in ring_groups:
        if cls == "amidine":
            members[next((a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 6), ring_atom)] = owned
    entries = []
    for owned in members.values():
        for nitrogen in (a for a in owned if mol.GetAtomWithIdx(a).GetAtomicNum() == 7):
            atom = mol.GetAtomWithIdx(nitrogen)
            imino = any(
                b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 6 for b in atom.GetBonds()
            )
            for n in graph[nitrogen]:
                if mol.GetAtomWithIdx(n).GetAtomicNum() == 6 and _amidine_nitrogens(mol, mol.GetAtomWithIdx(n)) is not None:
                    continue
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.append((name, compound, "N'" if imino else "N"))
    if entries and len(members) != 1:
        raise UnsupportedStructure("several amidine groups with N-substitution are not handled by the chain engine")
    return entries


def _imine_n_names(mol, graph, halogens, aromatic_atoms, groups, ring_groups):
    """N-prefix names when the single ring imine group is substituted on its =N atom, else []."""
    members = [owned for cls, _, owned in ring_groups if cls == "imine"] + list(groups.get("imine", {}).values())
    substituted = []
    for owned in members:
        nitrogen = next(iter(owned))
        subs = [n for n in graph[nitrogen] if mol.GetBondBetweenAtoms(nitrogen, n).GetBondTypeAsDouble() == 1.0]
        if subs:
            substituted.append((nitrogen, subs))
    if not substituted:
        return []
    if len(members) != 1:
        raise UnsupportedStructure("several imine groups with N-substitution are not handled by the chain engine")
    nitrogen, subs = substituted[0]
    return [name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in subs]


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
        located = ("N", anchor) if len(members) > 1 else "N"
        for n in graph[nitrogen]:
            if n != center:
                name, compound = name_branch(graph, n, nitrogen, halogens, aromatic_atoms, mol=mol, unsaturated=True)
                entries.append((name, compound, located))
    anchors = [anchor for _, anchor in members.values()]
    if entries and len(set(anchors)) != len(anchors):
        raise UnsupportedStructure("amide groups on one atom need primed N locants, which are not supported yet")
    return entries


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
    force = (attach is not None and length != 1) or (principal == "ide" and bool(grouped) and length > 1)
    prefix = format_substituent_prefixes(
        _with_n_names(grouped, n_names, position_of, len(on_chain)), omit_locants=length == 1 and not force and not n_names
    )
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
    elif principal == "amide":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "amide"))
    elif principal == "amidine":
        body = name_from_substituents(length, ene, yne, multiplied_word(count, "imidamide"))
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
        }[principal]
        body = name_from_substituents(
            length, ene, yne, multiplied_word(count, word), suffix_locants, force_own_locant=force,
            substituted=bool(grouped),
        )
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
        -total_count,
        locant_set,
        citation,
        _n_group_positions(n_names, position_of),
        _stereo_rank(stereo, position_of),
        name,
    )
    return key, name, (prefix, body, tail, reported_attach, position_of, False)


def _stereo_free(mol):
    return not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) and not any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    )
