"""Substituent prefixes joined through a heteroatom or a functional carbon
(P-29.3.3, P-29.4, P-35, P-63.2, P-66): hydroxy, oxo, alkoxy/aryloxy, sulfanyl,
amino, nitro, cyano, formyl, carboxy, carbamoyl, alkoxycarbonyl and acyl.
`name_branch` calls `hetero_branch_name` before it walks a carbon chain.
"""

import contextvars
import re

from rdkit import Chem

from ._alkoxy import alkoxy_prefix
from ._common import nonstandard_bonding, UnsupportedStructure, alpha_sort_key, is_nitro_nitrogen, named_prefix
from ._free_valence import SUFFIX_OF_ORDER
from ._multiplicative_text import enclose
from ._numerals import alkane_name, multiplying_prefix
from ._pin import mark

_ALKOXY_STEMS = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy", "phenyl": "phenoxy"}
_SIMPLE_NAMES = {
    "methoxy", "ethoxy", "propoxy", "butoxy", "tert-butoxy", "phenoxy", "amino", "anilino", "hydroxy", "oxo",
    "nitro", "nitroso", "cyano", "sulfanyl", "formyl", "carboxy", "carbamoyl",
}
_MULTIPLE_TARGETS = {7, 8, 16, 34, 52}
# Acid-derived and chalcogen-chain prefixes are only valid under a principal acid group; elsewhere the
# groups they describe outrank the parent and the engine must decline rather than cite them as prefixes.
EXTENDED_PREFIXES = contextvars.ContextVar("extended_prefixes", default=False)
CHALCOGEN_PREFIXES = {16: "sulfanyl", 34: "selanyl", 52: "tellanyl"}


_SENIOR_TO_SELENOL = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=O)[OX2H1]",
        "[CX3](=O)[OX2][#6]",
        "[CX3](=O)[NX3]",
        "[CX3](=[S,Se])[NX3]",
        "[CX2]#[NX1]",
        "[CX3H1](=O)[#6]",
        "[#6][CX3](=O)[#6]",
        "[CX3](=[O,S,Se,Te])[NX2]=[NX2]",
        "[#6][CX3](=O)[O,S,Se,Te;X2][O,S,Se,Te;X2][O,S,Se,Te;X2]",
        "[#6][CX3](=O)[OX2][NX3;!R]",
        "[OX2H1][#6;!$([#6]=O)]",
        "[SX2H1][#6;!$([#6]=[O,S,Se,Te])]",
    )
]
AMINE_PATTERN = Chem.MolFromSmarts("[NX3;!$(N~[!#6;!#1;!#8;!#16]);!$(N-[#6]=[O,S,N])]-[CX4]")
_HYDRAZINE = Chem.MolFromSmarts("[NX3;!R]-[NX3;!R]")


def require_plain_chalcogen_kids(mol, z, kids):
    """Se/Te bonded to an acyl, carbamoyl, formyl or cyano carbon is a selenoate/selenocyanate-type group, not a
    plain selanyl/tellanyl prefix."""
    if z in (34, 52) and not EXTENDED_PREFIXES.get():
        for kid in kids:
            atom = mol.GetAtomWithIdx(kid)
            if atom.GetAtomicNum() in (14, 32, 50, 82) or (atom.GetAtomicNum() in (16, 34, 52) and atom.GetDegree() <= 2):
                continue
            if atom.GetAtomicNum() != 6 or any(
                b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52)
                for b in atom.GetBonds()
            ):
                raise UnsupportedStructure("an acyl, carbamoyl or cyano group on selenium/tellurium is not a selanyl prefix")


CATION_PARENT = contextvars.ContextVar("cation_parent", default=False)


def _has_senior_principal_group(mol):
    """A principal group senior to the hetero-hetero connection (hydroxylamine, hydrazine, peroxide classes) is
    present, so that connection is expressed as a prefix (P-41, P-29.4.1); a cationic parent outranks every group."""
    return (
        CATION_PARENT.get()
        or any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_SELENOL)
        or mol.HasSubstructMatch(AMINE_PATTERN)
        or mol.HasSubstructMatch(_HYDRAZINE)
    )


def _dichalcogenide_only(mol):
    """Every contiguous run of chalcogen atoms is a pair joining two carbon groups (a disulfide, diselenide, ditelluride or a
    mixed S-O, Se-S pair): two contiguous chalcogens are a prefix on a carbon parent, only three or more form a parent
    hydride of their own (P-68.4.1.1, P-63.3.2)."""
    hosts = (6, 5, 13, 14, 31, 32, 49, 50, 81, 82)
    if any(a.GetAtomicNum() not in (1, 8, 9, 17, 35, 53, 16, 34, 52, *hosts) for a in mol.GetAtoms()):
        return False
    pairs = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (8, 16, 34, 52):
            continue
        partners = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in (8, 16, 34, 52)]
        if not partners:
            continue
        carbons = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in hosts]
        if atom.GetDegree() != 2 or len(partners) != 1 or len(carbons) != 1 or partners[0].GetDegree() != 2:
            return False
        pairs += 1
    return pairs > 0


def _chain_prefix_allowed(mol):
    return EXTENDED_PREFIXES.get() or _has_senior_principal_group(mol) or _dichalcogenide_only(mol)


def require_senior_group(mol, z):
    """A free -SeH/-TeH/=Se/=Te is the principal group (selenol, tellurol, selone) unless a more senior class
    (acid, ester, amide, nitrile, aldehyde, ketone, alcohol, thiol) is present (P-41, P-63.1)."""
    if z in (34, 52) and not any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_SELENOL):
        raise UnsupportedStructure("a selenol, tellurol, selone or tellone is the principal group, not a prefix")


MONONUCLEAR_HYDRIDES = {
    5: ("borane", "boranyl", 3),
    13: ("alumane", "alumanyl", 3),
    14: ("silane", "silyl", 4),
    15: ("phosphane", "phosphanyl", 3),
    31: ("gallane", "gallanyl", 3),
    32: ("germane", "germyl", 4),
    33: ("arsane", "arsanyl", 3),
    49: ("indigane", "indiganyl", 3),
    50: ("stannane", "stannyl", 4),
    51: ("stibane", "stibanyl", 3),
    81: ("thallane", "thallanyl", 3),
    82: ("plumbane", "plumbyl", 4),
    83: ("bismuthane", "bismuthanyl", 3),
}
_SIMPLE_ACYLS = {"methanoyl", "ethanoyl", "propanoyl", "butanoyl", "benzoyl"}


def _carbonyl_oxygen(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return next(
        (
            b.GetOtherAtom(atom).GetIdx()
            for b in atom.GetBonds()
            if b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() == 8
        ),
        None,
    )


DIPOLAR_GROUPS = contextvars.ContextVar("dipolar_groups", default=False)


def is_functional_carbon(mol, idx):
    """A carbon that is a substituent group of its own (-C#N, -COOH, -COOR,
    -CONR2, ...) rather than a chain member. Ketone and aldehyde carbons stay
    in the chain and are cited as 'oxo' (P-64.3)."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    for b in atom.GetBonds():
        other = b.GetOtherAtom(atom)
        if b.GetBondTypeAsDouble() == 3.0 and other.GetAtomicNum() == 7:
            if DIPOLAR_GROUPS.get() and other.GetFormalCharge() == 1:
                continue
            return True
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() in _MULTIPLE_TARGETS:
            return any(
                n.GetIdx() != other.GetIdx() and n.GetAtomicNum() in _MULTIPLE_TARGETS and not _non_amino_nitrogen(n)
                for n in atom.GetNeighbors()
            )
    return False


def _onium_prefix(graph, root, order, others, halogens, aromatic_atoms, mol):
    """'trimethylphosphaniumyl', 'oxidodi(phenyl)phosphaniumyl': a cationic centre with organyl and oxido groups."""
    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    if not (
        z in _ONIUM_PREFIX_STEMS
        and atom.GetFormalCharge() == 1
        and order == 1.0
        and not atom.IsInRing()
        and atom.GetTotalValence() == _ONIUM_PREFIX_STEMS[z][1]
        and all(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 or _terminal_anion(mol, n) for n in others)
        and all(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 1.0 for n in others)
    ):
        return None
    from ._substituents import format_mononuclear_prefixes

    entries = _group_names(graph, mol, others, root, halogens, aromatic_atoms)
    return (format_mononuclear_prefixes(entries) if entries else "") + _ONIUM_PREFIX_STEMS[z][0] + "yl", bool(entries)


_YLIDYNIUM_STEMS = {7: "azaniumylidyne", 8: "oxidaniumylidyne", 16: "sulfaniumylidyne"}


def _terminal_anion(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetAtomicNum() in (8, 16) and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1


def _non_amino_nitrogen(atom):
    """A nitrogen of an N=N, nitro or nitroso group: the substituent of a formazan or nitrolic acid carbon is not an
    amino group."""
    return atom.GetAtomicNum() == 7 and any(
        b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8) for b in atom.GetBonds()
    )


def _enclose(name, compound):
    return enclose(name) if compound else name


_ACID_CENTRE_NUMBERS = frozenset({5, 15, 16, 33, 34, 51, 52})
POLYACID_SUBSTITUENT_REASON = (
    "the preferred prefix of a chain of acid centres is its skeletal replacement ('a') name (P-67.2.6)"
)


def phosphoryl_name(parts, group=("phosphono", "phosphoryl")):
    """'phosphono' for P(O)(OH)2, otherwise '(X)(Y)phosphoryl' with X, Y cited alphabetically (P-65.1.3.1)."""
    from ._substituents import multiplied_prefix

    if all(name == "hydroxy" for name, _ in parts):
        return group[0]
    if len(set(parts)) == 1:
        name, compound = parts[0]
        body = multiplied_prefix(2, name, compound)
    else:
        ordered = sorted(parts, key=lambda p: alpha_sort_key(p[0]))
        body = "".join(enclose(n) if c or (i and not ordered[i - 1][1]) else n for i, (n, c) in enumerate(ordered))
    return body + group[1]


def _phosphoryloxy(graph, phosphorus, oxygen, halogens, aromatic_atoms, mol):
    """'phosphonooxy' or '[(X)(Y)phosphoryl]oxy' for O-P(=O)(X)(Y)."""
    atom = mol.GetAtomWithIdx(phosphorus)
    if atom.GetFormalCharge() or atom.IsInRing():
        raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
    others = [n for n in graph[phosphorus] if n != oxygen]
    found = _pnictogen_oxo_group(graph, phosphorus, halogens, aromatic_atoms, mol, others)
    if found is None:
        raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
    group, compound = found
    if group == "phosphono":
        return "phosphonooxy", True
    if "phospho" in group.replace("phosphoryl", ""):
        mark(None, POLYACID_SUBSTITUENT_REASON)
    return _enclose(group, compound) + "oxy", True


def _alkoxy(rname, compound=False):
    return alkoxy_prefix(rname, compound)


def _compound(name):
    return name not in _SIMPLE_NAMES


def _acyl_name(mol, graph, carbon, from_atom, halogens=None, aromatic_atoms=None):
    """alkanoyl or benzoyl prefix for R-C(=O)-: an unbranched saturated chain
    whose carbons may carry substituents (2-aminoethanoyl), or phenyl."""
    from ._amino_acyl_group import amino_acyl_group
    from ._substituents import name_branch

    amino = amino_acyl_group(mol, graph, carbon, from_atom)
    if amino is not None:
        _mark_stereo_used(mol, graph, carbon, from_atom)
        return amino[0]
    halogens = halogens or {}
    aromatic_atoms = aromatic_atoms or frozenset()
    others = [n for n in graph[carbon] if n != from_atom and mol.GetBondBetweenAtoms(carbon, n).GetBondTypeAsDouble() == 1.0]
    (alkyl,) = others
    if mol.GetAtomWithIdx(alkyl).GetIsAromatic():
        name, _ = name_branch(graph, alkyl, carbon, mol=mol)
        if name == "phenyl":
            return "benzoyl"
        raise UnsupportedStructure("this aroyl group is not supported yet")
    chain = [carbon, alkyl]
    while True:
        node, parent = chain[-1], chain[-2]
        atom = mol.GetAtomWithIdx(node)
        if atom.GetAtomicNum() != 6 or atom.IsInRing() or is_functional_carbon(mol, node):
            raise UnsupportedStructure("this acyl group is not supported yet")
        if mol.GetBondBetweenAtoms(node, parent).GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure("an unsaturated acyl group is not supported yet")
        onward = [
            n
            for n in graph[node]
            if n != parent and mol.GetAtomWithIdx(n).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(n).IsInRing()
            and not is_functional_carbon(mol, n)
        ]
        if len(onward) > 1:
            raise UnsupportedStructure("a branched acyl group is not supported yet")
        if not onward:
            break
        chain.append(onward[0])
    entries = {}
    chain_set = set(chain)
    for position, atom_idx in enumerate(chain, start=1):
        for n in graph[atom_idx]:
            if n in chain_set or (atom_idx == carbon and n == from_atom):
                continue
            if atom_idx == carbon:
                continue
            name, compound = name_branch(graph, n, atom_idx, halogens, aromatic_atoms, mol=mol)
            entries.setdefault(position, []).append((name, compound))
    from ._common import group_substituents
    from ._retained_acids import retained_chain_acid
    from ._substituents import format_substituent_prefixes

    grouped = group_substituents(entries)
    from ._substituents import _branch_stereo_entries, _branch_stereo_prefix

    stereo = _branch_stereo_prefix(_branch_stereo_entries({a: i for i, a in enumerate(chain, start=1)}, record=True))
    retained = retained_chain_acid(grouped, len(chain), [], [], 1, "acyl")
    if retained is not None:
        return stereo + retained
    prefix = format_substituent_prefixes(grouped, omit_locants=len(chain) == 1) if entries else ""
    return stereo + prefix + alkane_name(len(chain))[:-1] + "oyl"


def _acyl_from_acid_name(mol, graph, carbon, from_atom):
    """'(9Z)-octadec-9-enoyl' from the name of the acid whose acyl group is rooted at `carbon` (P-65.1.7.1)."""
    from ._amino_acyl_group import amino_acyl_group
    from ._functional_prefixes import _acyl_prefix

    amino = amino_acyl_group(mol, graph, carbon, from_atom)
    if amino is not None:
        _mark_stereo_used(mol, graph, carbon, from_atom)
        return amino[0]
    atoms = {carbon}
    stack = [carbon]
    while stack:
        for n in graph[stack.pop()]:
            if n != from_atom and n not in atoms:
                atoms.add(n)
                stack.append(n)
    if from_atom in atoms or any(mol.GetAtomWithIdx(a).GetFormalCharge() for a in atoms):
        raise UnsupportedStructure("this acyl group is not supported yet")
    name = _acyl_prefix(mol, atoms, carbon, from_atom)
    _mark_stereo_used(mol, graph, carbon, from_atom)
    return name


def _mark_stereo_used(mol, graph, carbon, from_atom):
    """Record that the stereo elements inside the acyl group are cited by its name."""
    from ._substituents import BRANCH_STEREO

    atoms = {carbon}
    stack = [carbon]
    while stack:
        for n in graph[stack.pop()]:
            if n != from_atom and n not in atoms:
                atoms.add(n)
                stack.append(n)
    context = BRANCH_STEREO.get()
    if context:
        context["used"].update(("atom", a) for a in atoms)
        context["used"].update(
            ("bond", (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
        )


def _group_names(graph, mol, atoms, parent, halogens, aromatic_atoms):
    from ._substituents import name_branch

    out = []
    for a in atoms:
        if is_functional_carbon(mol, a) and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and any(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() == 2.0
            for n in graph[a]
        ):
            out.append(_functional_carbon(graph, a, parent, halogens, aromatic_atoms, mol))
        else:
            out.append(name_branch(graph, a, parent, halogens, aromatic_atoms, mol=mol))
    return out


def _ring_nitrogen_acyl(graph, x, root, halogens, aromatic_atoms, mol, infix):
    """Acyl prefix of a ring nitrogen bonded to the acyl centre: 'pyrrolidine-1-carbonyl' (P-65.1.3), '(azetidin-1-yl)sulfonyl' (P-65.3.2.3)."""
    from ._substituents import name_branch

    name, compound = name_branch(graph, x, root, halogens, aromatic_atoms, mol=mol)
    if infix == "sulfonyl":
        return (_enclose(name, compound) if compound else name) + "sulfonyl", True
    ending = re.search(r"-(\d+[a-z]?)-yl$", name)
    if ending is None:
        raise UnsupportedStructure("this ring nitrogen acyl group is not supported yet")
    return name[: ending.start()] + "e-" + ending.group(1) + "-" + infix, True


def _amino_stem(amino):
    return "phenyl" if amino == "anilino" else amino[: -len("amino")]


def _amino(names):
    if not names:
        return "amino"
    if names == [("phenyl", False)]:
        return "anilino"
    ordered = sorted(names, key=lambda item: alpha_sort_key(item[0]))
    if len(ordered) == 2 and ordered[0] == ordered[1]:
        name, compound = ordered[0]
        from ._substituents import multiplied_prefix

        return multiplied_prefix(2, name, compound) + "amino"
    parts = [_enclose(ordered[0][0], ordered[0][1])]
    parts += [f"({n})" for n, _ in ordered[1:]]
    return "".join(parts) + "amino"


_CHAIN_ELEMENTS = {8, 16, 34, 52}
_CHAIN_WORDS = {8: "oxy", 16: "sulfanyl", 34: "selanyl", 52: "tellanyl"}
_CHAIN_HYDRO = {8: "hydroxy", 16: "sulfanyl", 34: "selanyl", 52: "tellanyl"}


def _chalcogen_chain_group(graph, first, second, halogens, aromatic_atoms, mol):
    """-Z1-Z2...-H or -Z1-Z2...-R substituent groups of a run of chalcogen atoms (P-63.4.2): 'hydroperoxy',
    '(methylperoxy)', 'disulfanyl', '(methyltrisulfanyl)', '(sulfanyloxy)', '(hydroxysulfanyl)'."""
    from ._substituents import name_branch

    run = [first, second]
    while True:
        onward = [n for n in graph[run[-1]] if n != run[-2]]
        if len(onward) != 1 or mol.GetAtomWithIdx(onward[0]).GetAtomicNum() not in _CHAIN_ELEMENTS:
            break
        run.append(onward[0])
    tail = onward
    atoms = [mol.GetAtomWithIdx(i) for i in run]
    if (
        any(atom.GetFormalCharge() for atom in atoms)
        or len(tail) > 1
        or any(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble() != 1.0 for i, j in zip(run, run[1:]))
        or any(mol.GetBondBetweenAtoms(run[-1], n).GetBondTypeAsDouble() != 1.0 for n in tail)
    ):
        raise UnsupportedStructure("this chalcogen chain is not supported yet")
    organyl = None
    if tail:
        end = mol.GetAtomWithIdx(tail[0])
        hydride_end = end.GetAtomicNum() in MONONUCLEAR_HYDRIDES and not end.IsInRing()
        if (end.GetAtomicNum() != 6 and not hydride_end) or (
            end.GetAtomicNum() == 6 and is_functional_carbon(mol, tail[0]) and not _thioacyl(mol, tail[0])
        ):
            raise UnsupportedStructure("a functional group on a chalcogen chain is not supported yet")
        organyl = name_branch(graph, tail[0], run[-1], halogens, aromatic_atoms, mol=mol)
    elif not any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_SELENOL):
        raise UnsupportedStructure("a peroxol or its chalcogen analogue outranks an amine as the principal group")
    runs = [[atom.GetAtomicNum(), 0] for atom in atoms[:1]]
    for atom in atoms[1:]:
        if atom.GetAtomicNum() == runs[-1][0] and (atom.GetAtomicNum() == 8 or len(set(a.GetAtomicNum() for a in atoms)) == 1):
            runs[-1][1] += 1
        else:
            runs.append([atom.GetAtomicNum(), 1])
    runs[0][1] += 1
    inner = organyl
    for z, count in reversed(runs[1:]):
        if inner is None:
            inner = _run_hydro(z, count), False
        elif z == 8 and count == 1:
            inner = _alkoxy(*inner)
        else:
            inner = _enclose(*inner) + _run_word(z, count), True
    z, count = runs[0]
    if inner is None:
        return _run_hydro(z, count), False
    return _enclose(*inner) + _run_word(z, count), True


def _run_word(z, count):
    if z == 8:
        return {1: "oxy", 2: "peroxy"}.get(count) or multiplying_prefix(count) + "oxidanyl"
    return (multiplying_prefix(count) if count > 1 else "") + _CHAIN_WORDS[z]


def _run_hydro(z, count):
    if z == 8:
        return {1: "hydroxy", 2: "hydroperoxy"}.get(count) or multiplying_prefix(count) + "oxidanyl"
    return (multiplying_prefix(count) if count > 1 else "") + _CHAIN_WORDS[z]


def _thioacyl(mol, idx):
    """An acyl centre: a carbon with a terminal =O, =S, =Se or =Te, or a sulfonyl-type S, Se or Te with two terminal =O
    (P-65.1.7.2.3, P-66.1.1.4.3)."""
    atom = mol.GetAtomWithIdx(idx)
    terminal = [
        n
        for n in atom.GetNeighbors()
        if n.GetDegree() == 1 and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    if atom.GetAtomicNum() == 6:
        return any(n.GetAtomicNum() in (8, 16, 34, 52) for n in terminal)
    return atom.GetAtomicNum() in (16, 34, 52) and sum(n.GetAtomicNum() == 8 for n in terminal) == 2


def _chalcogen_formyl(mol, idx):
    """A -CH=S, -CH=Se or -CH=Te group, the analogue of formyl, which cannot join a chain from a substituent root."""
    atom = mol.GetAtomWithIdx(idx)
    return (
        atom.GetDegree() == 2
        and atom.GetTotalNumHs() == 1
        and any(
            n.GetAtomicNum() in (16, 34, 52)
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
    )


def _imidoyl_centre(mol, idx):
    """A carbon with a terminal =NH, or a sulfonyl-type S, Se or Te with two terminal =O or =NH of which at least one is
    =NH: the acyl group of an imidamide (P-66.4.1.3.5)."""
    atom = mol.GetAtomWithIdx(idx)
    terminal = [
        n
        for n in atom.GetNeighbors()
        if n.GetDegree() == 1 and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    imines = sum(n.GetAtomicNum() == 7 and n.GetTotalNumHs() == 1 and not n.GetFormalCharge() for n in terminal)
    if atom.GetAtomicNum() == 6:
        hydrazones = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 7
            and n.GetDegree() == 2
            and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
            and any(m.GetAtomicNum() == 7 and m.GetDegree() == 1 and m.GetTotalNumHs() == 2 for m in n.GetNeighbors())
        ]
        return (imines == 1 and len(terminal) == 1 and not hydrazones) or (not terminal and len(hydrazones) == 1)
    oxygens = sum(n.GetAtomicNum() == 8 for n in terminal)
    return atom.GetAtomicNum() in (16, 34, 52) and imines >= 1 and imines + oxygens == 2


PEROXY_PREFIXES = contextvars.ContextVar("peroxy_prefixes", default=False)


def hetero_branch_name(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """(name, is_compound) of a heteroatom- or functional-carbon-rooted
    substituent, or None when `root` is an ordinary carbon."""
    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    if z != 6 and _in_anionic_chain(mol, root, coming_from):
        from ._anion_chain import chain_prefix

        found = chain_prefix(mol, root, coming_from)
        if found is not None:
            return found
    if atom.HasProp("_anion"):
        if z == 6 and atom.IsInRing():
            from ._diester_ring_diyl import ring_substituent_name

            return ring_substituent_name(mol, graph, root, coming_from)
        if z == 6:
            return _carbon_anion_prefix(graph, root, coming_from, halogens, aromatic_atoms, mol)
        return _anionic_group(mol, root, coming_from)
    if atom.HasProp("_anion_word"):
        from ._anion_center import center_prefix

        return center_prefix(mol, root, coming_from)
    if z == 6 and not atom.HasProp("_anion") and _branch_has_anionic_carbon(graph, root, coming_from, mol):
        if atom.IsInRing():
            raise UnsupportedStructure("an anionic carbon beyond a ring substituent is not supported yet")
        return _anionic_chain_prefix(graph, root, coming_from, mol)
    if z == 6:
        from ._imidoyl_prefix import imidoyl_prefix, ketene_prefix

        named = ketene_prefix(mol, graph, root, coming_from) or imidoyl_prefix(mol, graph, root, coming_from)
        if named is not None:
            return named
        if is_functional_carbon(mol, root) or _carbonyl_oxygen(mol, root) is not None or (
            (EXTENDED_PREFIXES.get() or _chalcogen_formyl(mol, root)) and _thioacyl(mol, root)
        ):
            return _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol)
        return None
    if z in (15, 33) and any(mol.GetAtomWithIdx(n).HasProp("_anion_word") for n in graph[root]):
        return _oxoacid_anion_prefix(graph, root, coming_from, mol)
    if z == 5 and EXTENDED_PREFIXES.get() and not atom.GetFormalCharge():
        others = [n for n in graph[root] if n != coming_from]
        borono = _borono_prefix(mol, root, others)
        if borono is not None:
            return borono, borono != "borono"
    if z in _ONIUM_PREFIX_STEMS and atom.GetFormalCharge() == 1:
        onium = _onium_prefix(
            graph, root, mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble(),
            [n for n in graph[root] if n != coming_from], halogens, aromatic_atoms, mol,
        )
        if onium is not None:
            return onium
    if z in MONONUCLEAR_HYDRIDES and not atom.IsInRing():
        return _mononuclear_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if z in LAMBDA_CENTRE_STEMS and not atom.IsInRing():
        lambda_group = _lambda_centre_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
        if lambda_group is not None:
            return lambda_group
    if z in _HALOGEN_STEMS:
        named = halogen_oxo_prefix(mol, root, coming_from)
        if named is not None:
            return named, not named.startswith((_HALOGEN_STEMS[z], "per"))
    if atom.GetFormalCharge() == -1 and atom.GetDegree() == 1 and z in (8, 16):
        return {8: "oxido", 16: "sulfido"}[z], False
    if (
        z in _YLIDYNIUM_STEMS
        and atom.GetFormalCharge() == 1
        and mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() == 3.0
        and DIPOLAR_GROUPS.get()
        and all(_terminal_anion(mol, n) for n in graph[root] if n != coming_from)
    ):
        from ._substituents import format_mononuclear_prefixes

        entries = _group_names(graph, mol, [n for n in graph[root] if n != coming_from], root, halogens, aromatic_atoms)
        return (format_mononuclear_prefixes(entries) if entries else "") + _YLIDYNIUM_STEMS[z], True
    if atom.GetFormalCharge() and z != 7:
        raise UnsupportedStructure("a charged atom in a substituent group is not supported yet")
    if atom.IsInRing():
        return None
    others = [n for n in graph[root] if n != coming_from]
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    if z in _CYANATE_PREFIXES and order == 1.0 and len(others) == 1 and _is_cyanide_carbon(mol, others[0], root):
        return _CYANATE_PREFIXES[z], z != 8
    if z == 8:
        if order == 2.0:
            return "oxo", False
        if not others:
            return "hydroxy", False
        (other,) = others
        named = named_prefix(mol.GetAtomWithIdx(other))
        if named is not None:
            return ("ylooxidanyl" if named == "ylo" else named + "oxy"), True
        if mol.GetAtomWithIdx(other).GetAtomicNum() == 15 and any(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(other, n).GetBondTypeAsDouble() == 2.0
            for n in graph[other]
        ):
            return _phosphoryloxy(graph, other, root, halogens, aromatic_atoms, mol)
        if mol.GetAtomWithIdx(other).GetAtomicNum() in MONONUCLEAR_HYDRIDES:
            silyl, _ = _mononuclear_group(graph, other, root, halogens, aromatic_atoms, mol)
            return _enclose(silyl, True) + "oxy", True
        if mol.GetAtomWithIdx(other).GetAtomicNum() in _E_SYMBOL and _chalcogen_acyl_oxo(mol, other, root):
            acyl, acyl_compound = _sulfur_oxo_group(graph, other, root, halogens, aromatic_atoms, mol)
            return _enclose(acyl, acyl_compound) + "oxy", True
        if mol.GetAtomWithIdx(other).GetAtomicNum() in _CHAIN_ELEMENTS and _chain_prefix_allowed(mol):
            return _chalcogen_chain_group(graph, root, other, halogens, aromatic_atoms, mol)
        if (
            mol.GetAtomWithIdx(other).GetAtomicNum() == 8
            and order == 1.0
            and not mol.GetAtomWithIdx(other).GetFormalCharge()
            and mol.GetBondBetweenAtoms(root, other).GetBondTypeAsDouble() == 1.0
            and mol.GetAtomWithIdx(other).GetDegree() <= 2
            and PEROXY_PREFIXES.get()
        ):
            onward = [n for n in graph[other] if n != root]
            if not onward:
                return "hydroperoxy", False
            if mol.GetAtomWithIdx(onward[0]).GetAtomicNum() == 6:
                from ._substituents import name_branch

                rname, rcomp = name_branch(graph, onward[0], other, halogens, aromatic_atoms, mol=mol)
                return _enclose(rname, rcomp) + "peroxy", True
        if mol.GetAtomWithIdx(other).GetAtomicNum() != 6:
            if mol.GetAtomWithIdx(other).GetAtomicNum() == 8 or not _has_senior_principal_group(mol):
                raise UnsupportedStructure("this oxygen-linked group is not supported yet")
            from ._substituents import name_branch

            rname, rcomp = name_branch(graph, other, root, halogens, aromatic_atoms, mol=mol)
            return _enclose(rname, rcomp) + "oxy", True
        if is_functional_carbon(mol, other):
            (acyl,) = _group_names(graph, mol, [other], root, halogens, aromatic_atoms)
            return _enclose(acyl[0], acyl[1]) + "oxy", True
        from ._substituents import name_branch

        rname, rcomp = name_branch(graph, other, root, halogens, aromatic_atoms, mol=mol)
        return _alkoxy(rname, rcomp)
    if z in (16, 34, 52) and atom.GetDegree() > 2:
        return _sulfur_oxo_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if z in CHALCOGEN_PREFIXES:
        word = CHALCOGEN_PREFIXES[z]
        if not others:
            require_senior_group(mol, z)
        require_plain_chalcogen_kids(mol, z, others)
        if order == 2.0 and not others:
            return word[:-2] + "ylidene", False
        if order != 1.0 or len(others) > 1 or atom.GetDegree() > 2:
            raise UnsupportedStructure("this chalcogen-linked group is not supported yet")
        if not others:
            return word, False
        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() in MONONUCLEAR_HYDRIDES:
            silyl, silyl_compound = _mononuclear_group(graph, others[0], root, halogens, aromatic_atoms, mol)
            return _enclose(silyl, silyl_compound) + word, True
        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() in _E_SYMBOL and _chalcogen_acyl_oxo(mol, others[0], root):
            acyl, acyl_compound = _sulfur_oxo_group(graph, others[0], root, halogens, aromatic_atoms, mol)
            return _enclose(acyl, acyl_compound) + word, True
        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() in _CHAIN_ELEMENTS and _chain_prefix_allowed(mol):
            return _chalcogen_chain_group(graph, root, others[0], halogens, aromatic_atoms, mol)
        from ._substituents import name_branch

        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() != 6:
            if mol.GetAtomWithIdx(others[0]).GetAtomicNum() == z or not _has_senior_principal_group(mol):
                raise UnsupportedStructure("a chalcogen chain (disulfanyl, ...) is not supported yet")
            rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
            return _enclose(rname, rcomp) + word, True

        rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
        if re.fullmatch(r"[a-z]+animidoyl", rname):
            rcomp = False
        return _enclose(rname, rcomp) + word, True
    if z == 7 and atom.GetFormalCharge() == 1 and nitrogen_pseudohalide_prefix(mol, root, others, order) == "isocyano":
        return "isocyano", False
    if (
        z == 7
        and (
            (atom.GetFormalCharge() == 1 and atom.GetTotalValence() == 4)
            or (atom.HasProp("_cationic_amine") and not atom.GetFormalCharge())
        )
        and order == 1.0
        and all(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 or _terminal_anion(mol, n) for n in others)
    ):
        if any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() > 2.0 for n in others):
            raise UnsupportedStructure("this azaniumyl group is not supported yet")
        from ._substituents import format_mononuclear_prefixes

        entries = _group_names(graph, mol, others, root, halogens, aromatic_atoms)
        return (format_mononuclear_prefixes(entries) if entries else "") + "azaniumyl", bool(entries)
    onium = _onium_prefix(graph, root, order, others, halogens, aromatic_atoms, mol)
    if onium is not None:
        return onium
    if (
        z == 7
        and atom.GetFormalCharge() == 1
        and order == 2.0
        and atom.GetDegree() + atom.GetTotalNumHs() == 3
        and all(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 for n in others)
    ):
        if any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others):
            raise UnsupportedStructure("this azaniumylidene group is not supported yet")
        from ._substituents import format_mononuclear_prefixes

        entries = _group_names(graph, mol, others, root, halogens, aromatic_atoms)
        return (format_mononuclear_prefixes(entries) if entries else "") + "azaniumylidene", bool(entries)
    if z == 7:
        oxygens = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
        if len(oxygens) == len(others) and others:
            if len(oxygens) == 2 and is_nitro_nitrogen(mol, root):
                return "nitro", False
            if len(oxygens) == 1 and mol.GetBondBetweenAtoms(root, oxygens[0]).GetBondTypeAsDouble() == 2.0:
                return "nitroso", False
        pseudohalide = nitrogen_pseudohalide_prefix(mol, root, others, order)
        if pseudohalide is not None:
            return pseudohalide, False
        if order == 2.0 and atom.GetFormalCharge() == 1 and len(others) == 1 and mol.GetAtomWithIdx(coming_from).GetAtomicNum() == 6:
            far = mol.GetAtomWithIdx(others[0])
            if (
                far.GetAtomicNum() == 7
                and far.GetFormalCharge() == -1
                and far.GetDegree() == 1
                and mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() == 2.0
            ):
                return "diazo", False
        if EXTENDED_PREFIXES.get() and order == 2.0 and not atom.GetFormalCharge() and len(others) <= 1:
            if not others and atom.GetTotalNumHs() == 1:
                return "imino", False
            if others and not atom.GetTotalNumHs() and mol.GetAtomWithIdx(others[0]).GetAtomicNum() in (6, 8, 16):
                (named,) = _group_names(graph, mol, others, root, halogens, aromatic_atoms)
                return _enclose(*named) + "imino", True
        diazenyl = _diazenyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol, order, others)
        if diazenyl is not None:
            return diazenyl
        if order == 2.0 and len(others) == 1 and mol.GetAtomWithIdx(others[0]).GetAtomicNum() == 7:
            far = mol.GetAtomWithIdx(others[0])
            if (
                not far.GetFormalCharge()
                and not atom.GetFormalCharge()
                and mol.GetBondBetweenAtoms(root, far.GetIdx()).GetBondTypeAsDouble() == 1.0
                and far.GetTotalNumHs()
                + sum(mol.GetBondBetweenAtoms(far.GetIdx(), n).GetBondTypeAsDouble() for n in graph[far.GetIdx()] if n != root)
                == 2
            ):
                tail = [n for n in graph[far.GetIdx()] if n != root]
                if not tail:
                    return "hydrazinylidene", False
                from ._substituents import format_mononuclear_prefixes

                entries = _group_names(graph, mol, tail, far.GetIdx(), halogens, aromatic_atoms)
                prefix = f"({entries[0][0]})" if len(entries) == 1 and entries[0][1] else format_mononuclear_prefixes(entries)
                return prefix + "hydrazinylidene", True
        if order == 2.0 and not atom.GetFormalCharge() and len(others) <= 1 and _has_senior_principal_group(mol):
            if not others:
                return "imino", False
            if mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() == 1.0:
                from ._substituents import name_branch

                rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
                return _enclose(rname, rcomp) + "imino", True
        from ._chalcogenourea import is_oxo_nitrogen

        if order == 1.0 and not atom.GetFormalCharge():
            imidohydrazido = _hydrazido(graph, root, others, mol)
            if imidohydrazido is not None:
                return imidohydrazido
        if order == 1.0 and not atom.GetFormalCharge() and any(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 7 and not is_oxo_nitrogen(mol, mol.GetAtomWithIdx(n)) for n in others
        ):
            far = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 7]
            if not (len(others) == 1 and mol.GetAtomWithIdx(far[0]).GetDegree() == 1):
                return _chain_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
        if order == 1.0 and not atom.GetFormalCharge() and len(others) == 1 and mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() == 2.0:
            from ._substituents import name_branch

            rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
            return _enclose(rname, rcomp) + "amino", True
        aci = _aci_nitro_ligands(mol, root, others, order)
        if aci is not None:
            from ._substituents import format_mononuclear_prefixes, name_branch

            entries = [("oxo", False)] + [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in aci]
            return format_mononuclear_prefixes(entries) + "-λ5-azanylidene", True
        if order != 1.0 or atom.GetFormalCharge() or any(
            mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others
        ):
            raise UnsupportedStructure("this nitrogen-linked group is not supported yet")
        if len(others) == 1 and mol.GetAtomWithIdx(others[0]).GetAtomicNum() == 7:
            far = mol.GetAtomWithIdx(others[0])
            if far.GetDegree() == 1 and far.GetTotalNumHs() == 2 and not far.GetFormalCharge():
                return "hydrazinyl", False
        amido = _chalcogen_amido(graph, root, others, mol) or _hydrazido(graph, root, others, mol)
        if amido is not None:
            return amido
        if any(mol.GetAtomWithIdx(n).GetAtomicNum() not in (6,) + tuple(MONONUCLEAR_HYDRIDES) for n in others) and not (
            _has_senior_principal_group(mol)
        ):
            raise UnsupportedStructure("this nitrogen-linked group is not supported yet")
        if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in others):
            from ._substituents import name_branch

            names = []
            for n in others:
                z_n = mol.GetAtomWithIdx(n).GetAtomicNum()
                if z_n in MONONUCLEAR_HYDRIDES:
                    names.append(_mononuclear_group(graph, n, root, halogens, aromatic_atoms, mol))
                elif z_n == 6:
                    names.extend(_group_names(graph, mol, [n], root, halogens, aromatic_atoms))
                else:
                    names.append(name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol))
            name = _amino(names)
            return name, _compound(name)
        from ._anilino import anilino_prefix

        anilino = anilino_prefix(graph, mol, root, others, halogens, aromatic_atoms)
        if anilino is not None:
            return anilino[0], anilino[0] != "anilino"
        name = _amino(_group_names(graph, mol, others, root, halogens, aromatic_atoms))
        return name, _compound(name)
    raise UnsupportedStructure("this heteroatom-linked substituent is not supported yet")


def _side(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _acylated_nitrogens(fragment):
    """Carbon, sulfur, selenium or tellurium centres of `fragment` that carry a double-bonded atom and a nitrogen."""
    return sum(
        1
        for c in fragment.GetAtoms()
        if c.GetAtomicNum() in (6, 16, 34, 52)
        and any(n.GetAtomicNum() == 7 for n in c.GetNeighbors())
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(c).GetAtomicNum() != 6 for b in c.GetBonds())
    )


def _chalcogen_amido(graph, root, others, mol):
    """'acetamido', 'ethanethioamido', 'methanesulfonamido' for R-CO-NH-, R-CS-NH-, R-SO2-NH- and N-substituted
    forms: the final 'e' in the complete name of the amide becomes 'o' (P-66.1.1.4.3); None for any other nitrogen."""
    from ._polyfunctional import name_polyfunctional

    acyl = [n for n in others if _thioacyl(mol, n) or _imidoyl_centre(mol, n)]
    rest = [n for n in others if n not in acyl]
    if len(acyl) != 1 or len(rest) > 1 or any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in rest):
        return None
    if not rest:
        from ._amino_acyl_group import amino_acyl_group

        amino = amino_acyl_group(mol, graph, acyl[0], root)
        if amino is not None:
            return amino[0] + "amino", True
    if any(a.GetIsotope() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
        b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()
    ):
        return None
    keep = {root} | _side(graph, acyl[0], root)
    if rest:
        keep |= _side(graph, rest[0], root)
    fragment = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        fragment.RemoveAtom(index)
    fragment = fragment.GetMol()
    try:
        Chem.SanitizeMol(fragment)
        name = contextvars.Context().run(name_polyfunctional, fragment)
    except (UnsupportedStructure, ValueError):
        return None
    if not name.endswith("amide") or _acylated_nitrogens(fragment) > 1:
        return None
    prefix = name[:-1] + "o"
    return prefix, any(part in prefix for part in ("-", "ane", "benzene"))


def _hydrazido(graph, root, others, mol):
    """'ethanimidohydrazido', 'formohydrazido', 'hydrazinecarbohydrazido' for R-CO-NH-NH- and R-C(=NH)-NH-NH- joined
    through the terminal nitrogen (P-66.3.5.3, P-66.4.2.3.6): the final 'e' of the hydrazide name becomes 'o'; None for
    any other hydrazine nitrogen."""
    from ._carbonic_hydrazide import carbonic_hydrazide_name
    from ._polyfunctional import name_polyfunctional

    if len(others) != 1 or mol.GetAtomWithIdx(others[0]).GetAtomicNum() != 7 or mol.GetAtomWithIdx(root).IsInRing():
        return None
    alpha = others[0]
    onward = [n for n in graph[alpha] if n != root]
    if len(onward) != 1 or mol.GetAtomWithIdx(alpha).IsInRing() or mol.GetAtomWithIdx(onward[0]).GetAtomicNum() != 6:
        return None
    if not (_imidoyl_centre(mol, onward[0]) or _thioacyl(mol, onward[0])) or any(
        a.GetIsotope() or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()
    ):
        return None
    keep = {root, alpha} | _side(graph, onward[0], alpha)
    fragment = Chem.RWMol(mol)
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        fragment.RemoveAtom(index)
    fragment = fragment.GetMol()
    try:
        Chem.SanitizeMol(fragment)
        name = carbonic_hydrazide_name(fragment) or contextvars.Context().run(name_polyfunctional, fragment)
    except (UnsupportedStructure, ValueError):
        try:
            from .core import smiles_to_iupac

            name = contextvars.Context().run(smiles_to_iupac, Chem.MolToSmiles(fragment))
        except (UnsupportedStructure, ValueError):
            return None
    if not name.endswith("hydrazide") or name.endswith("dihydrazide") and name.startswith("dicarbonic"):
        return None
    return name[:-1] + "o", not name.startswith(("formo", "aceto", "benzo"))


_CYANATE_PREFIXES = {8: "cyanato", 16: "thiocyanato", 34: "selenocyanato", 52: "tellurocyanato"}


def _is_cyanide_carbon(mol, idx, chalcogen):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.GetDegree() != 2 or atom.GetFormalCharge():
        return False
    other = next(n for n in atom.GetNeighbors() if n.GetIdx() != chalcogen)
    return (
        other.GetAtomicNum() == 7
        and other.GetDegree() == 1
        and not other.GetFormalCharge()
        and mol.GetBondBetweenAtoms(idx, other.GetIdx()).GetBondTypeAsDouble() == 3.0
    )


_ISOCYANATE_PREFIXES = {8: "isocyanato", 16: "isothiocyanato", 34: "isoselenocyanato", 52: "isotellurocyanato"}


def nitrogen_pseudohalide_prefix(mol, root, others, order):
    """'azido' for -N=N(+)=N(-) and 'isocyanato' (-N=C=O, also S/Se/Te) for a neutral nitrogen singly bonded to its
    parent (P-61.7, P-61.11); None for any other nitrogen."""
    atom = mol.GetAtomWithIdx(root)
    if order != 1.0 or len(others) != 1:
        return None
    middle = mol.GetAtomWithIdx(others[0])
    if atom.GetFormalCharge() == 1:
        isocyanide = (
            mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() == 3.0
            and middle.GetAtomicNum() == 6
            and middle.GetDegree() == 1
            and middle.GetFormalCharge() == -1
        )
        return "isocyano" if isocyanide else None
    if atom.GetFormalCharge():
        return None
    if mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() != 2.0 or middle.GetDegree() != 2:
        return None
    ends = [n for n in middle.GetNeighbors() if n.GetIdx() != root]
    end = ends[0]
    if mol.GetBondBetweenAtoms(middle.GetIdx(), end.GetIdx()).GetBondTypeAsDouble() != 2.0 or end.GetDegree() != 1:
        return None
    if middle.GetAtomicNum() == 7 and middle.GetFormalCharge() == 1:
        return "azido" if end.GetAtomicNum() == 7 and end.GetFormalCharge() == -1 else None
    if middle.GetAtomicNum() == 6 and not middle.GetFormalCharge() and not end.GetFormalCharge():
        return _ISOCYANATE_PREFIXES.get(end.GetAtomicNum())
    return None


_HALOGEN_STEMS = {9: "fluor", 17: "chlor", 35: "brom", 53: "iod"}
_HALOGEN_CHALCOGENS = {8: "", 16: "thio", 34: "seleno", 52: "telluro"}


def _terminal_oxo_oxygens(atom, exclude):
    """The oxygens of `atom` that are only a =O or a -O(-): None when it also carries another kind of substituent."""
    oxygens = []
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == exclude:
            continue
        double = bond.GetBondTypeAsDouble() == 2.0 and not other.GetFormalCharge()
        anionic = bond.GetBondTypeAsDouble() == 1.0 and other.GetFormalCharge() == -1
        if other.GetAtomicNum() not in _HALOGEN_CHALCOGENS or other.GetDegree() != 1 or not (double or anionic):
            return None
        oxygens.append(other)
    return oxygens


def halogen_oxo_prefix(mol, root, coming_from):
    """'chlorosyl', 'chloryl' or 'perchloryl' (-XO, -XO2, -XO3; P-61.3.2.3, P-67.1.4.5) for a halogen bearing one to
    three terminal oxygens, else None."""
    atom = mol.GetAtomWithIdx(root)
    stem = _HALOGEN_STEMS.get(atom.GetAtomicNum())
    oxygens = _terminal_oxo_oxygens(atom, coming_from) if stem else None
    if not oxygens or len(oxygens) > 3:
        return None
    anions = sum(o.GetFormalCharge() for o in oxygens)
    if atom.GetFormalCharge() + anions != 0 or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        return None
    base = {1: f"{stem}osyl", 2: f"{stem}yl", 3: f"per{stem}yl"}[len(oxygens)]
    replaced = {}
    for o in oxygens:
        word = _HALOGEN_CHALCOGENS[o.GetAtomicNum()]
        if word:
            replaced[word] = replaced.get(word, 0) + 1
    from ._numerals import multiplying_prefix

    return "".join(
        (multiplying_prefix(count) if count > 1 else "") + word for word, count in sorted(replaced.items())
    ) + base


def is_halogen_oxo_part(mol, atom):
    """A charged atom of a halogen oxo group: the halogen cation or one of its oxide anions."""
    if atom.GetAtomicNum() in _HALOGEN_CHALCOGENS and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1:
        (halogen,) = atom.GetNeighbors()
        return halogen.GetAtomicNum() in _HALOGEN_STEMS and is_halogen_oxo_part(mol, halogen)
    if atom.GetAtomicNum() not in _HALOGEN_STEMS or atom.GetDegree() < 2:
        return False
    parents = [n for n in atom.GetNeighbors() if n.GetAtomicNum() not in _HALOGEN_CHALCOGENS or n.GetDegree() != 1]
    return len(parents) == 1 and halogen_oxo_prefix(mol, atom.GetIdx(), parents[0].GetIdx()) is not None


def _is_anionic_oxygen(atom):
    return atom.HasProp("_anion") or atom.GetFormalCharge() == -1


_ACYL_HALIDE_PREFIXES = {9: "carbonofluoridoyl", 17: "carbonochloridoyl", 35: "carbonobromidoyl", 53: "carbonoiodidoyl"}


def _is_amino_nitrogen(mol, atom_index, neighbor):
    atom = mol.GetAtomWithIdx(atom_index)
    return (
        atom.GetAtomicNum() == 7
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 2
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(atom_index, neighbor).GetBondTypeAsDouble() == 1.0
    )


def _is_plain_amidine(mol, graph, root, coming_from):
    """C(=NH)NH2 attached through `root` (P-66.4.1.3, carbamimidoyl)."""
    others = [n for n in graph[root] if n != coming_from]
    if len(others) != 2 or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        return False
    orders = sorted(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() for n in others)
    atoms = [mol.GetAtomWithIdx(n) for n in others]
    return (
        orders == [1.0, 2.0]
        and all(a.GetAtomicNum() == 7 and a.GetDegree() == 1 and not a.GetFormalCharge() for a in atoms)
        and sum(a.GetTotalNumHs() for a in atoms) == 3
    )


_NITRILE_OXIDE_YLIDENE = {8: "oxo", 16: "sulfanylidene", 34: "selanylidene", 52: "tellanylidene"}


def _nitrile_oxide_prefix(mol, nitrogen):
    """'(oxo-\u03bb5-azanylidyne)methyl' for the nitrogen of a nitrile oxide or its chalcogen analogue cited as a prefix
    (P-66.5.4.2), else None."""
    atom = mol.GetAtomWithIdx(nitrogen)
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() in _NITRILE_OXIDE_YLIDENE and n.GetDegree() == 1:
            return f"({_NITRILE_OXIDE_YLIDENE[n.GetAtomicNum()]}-\u03bb5-azanylidyne)methyl"
    return None


def _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol):
    from ._substituents import ISOTOPE_LABELS, _labelled_carboxy, name_branch

    if ISOTOPE_LABELS.get():
        carboxy = _labelled_carboxy(graph, root, coming_from, mol, ISOTOPE_LABELS.get())
        if carboxy is not None:
            return carboxy[0]

    atom = mol.GetAtomWithIdx(root)
    others = [n for n in graph[root] if n != coming_from]
    if mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() == 2.0 and len(others) in (1, 2) and not atom.IsInRing():
        amines = [mol.GetAtomWithIdx(n) for n in others]
        if all(
            a.GetAtomicNum() in (7, 8, 16) and not a.GetFormalCharge() and not a.IsInRing()
            and mol.GetBondBetweenAtoms(root, a.GetIdx()).GetBondTypeAsDouble() == 1.0
            for a in amines
        ):
            entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
            from ._substituents import format_mononuclear_prefixes

            return format_mononuclear_prefixes(entries) + "methylidene", True
    triple_n = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 3.0]
    if triple_n and len(others) == 1:
        nitrile_oxide = _nitrile_oxide_prefix(mol, triple_n[0])
        if nitrile_oxide is not None:
            return nitrile_oxide, True
        return "cyano", False
    from ._acid_prefixes import acid_group_prefix

    named = acid_group_prefix(mol, graph, root, coming_from, halogens, aromatic_atoms, name_branch)
    if named is not None:
        return named
    carbonyl = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0]
    if not carbonyl and _is_plain_amidine(mol, graph, root, coming_from):
        return "carbamimidoyl", False
    if len(carbonyl) != 1 or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")
    rest = [n for n in others if n != carbonyl[0]]
    if not rest:
        return "formyl", False
    (x,) = rest
    z = mol.GetAtomWithIdx(x).GetAtomicNum()
    if z in _ACYL_HALIDE_PREFIXES and mol.GetAtomWithIdx(x).GetDegree() == 1:
        return _ACYL_HALIDE_PREFIXES[z], False
    if z == 8:
        tail = [n for n in graph[x] if n != root]
        if not tail:
            return ("carboxylato" if _is_anionic_oxygen(mol.GetAtomWithIdx(x)) else "carboxy"), False
        if named_prefix(mol.GetAtomWithIdx(tail[0])) == "ylo":
            return "oxylcarbonyl", True
        rname, rcomp = name_branch(graph, tail[0], x, halogens, aromatic_atoms, mol=mol)
        return _enclose(*_alkoxy(rname, rcomp)) + "carbonyl", True
    if z == 7:
        if mol.GetAtomWithIdx(x).HasProp("_anion_word"):
            raise UnsupportedStructure("an anionic amide nitrogen is not named as a carbamoyl prefix")
        subs = [n for n in graph[x] if n != root]
        if not subs:
            return "carbamoyl", False
        if mol.GetAtomWithIdx(x).IsInRing():
            return _ring_nitrogen_acyl(graph, x, root, halogens, aromatic_atoms, mol, "carbonyl")
        if len(subs) == 1 and _is_amino_nitrogen(mol, subs[0], x):
            return "hydrazinecarbonyl", True
        if (
            len(subs) == 1
            and mol.GetAtomWithIdx(subs[0]).GetAtomicNum() == 7
            and mol.GetBondBetweenAtoms(x, subs[0]).GetBondTypeAsDouble() == 2.0
            and not mol.GetAtomWithIdx(subs[0]).GetFormalCharge()
        ):
            tail = [m for m in graph[subs[0]] if m != x]
            if not tail:
                return "diazenecarbonyl", True
            if len(tail) == 1 and mol.GetAtomWithIdx(tail[0]).GetAtomicNum() == 6:
                rname, rcomp = name_branch(graph, tail[0], subs[0], halogens, aromatic_atoms, mol=mol)
                return f"2-{_enclose(rname, rcomp)}diazene-1-carbonyl", True
        name = _amino(_group_names(graph, mol, subs, x, halogens, aromatic_atoms))
        return _amino_stem(name) + "carbamoyl", True
    if z == 6:
        try:
            name = _acyl_name(mol, graph, root, coming_from, halogens, aromatic_atoms)
        except UnsupportedStructure:
            name = _acyl_from_acid_name(mol, graph, root, coming_from)
        from ._retained_acids import is_compound_acyl

        return name, is_compound_acyl(name)
    raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")


_E_SYMBOL = {16: "S", 34: "Se", 52: "Te"}
_E_ACID_PREFIX = {
    ("S", 2): "sulfo",
    ("S", 1): "sulfino",
    ("Se", 2): "selenono",
    ("Se", 1): "selenino",
    ("Te", 2): "tellurono",
    ("Te", 1): "tellurino",
}
_REPLACEMENT_PREFIX = {"S": "thio", "Se": "seleno", "Te": "telluro"}
_REPLACEMENT_SYMBOL = {16: "S", 34: "Se", 52: "Te"}


ANIONIC_PREFIXES = {7: "azanidyl", 8: "oxido", 16: "sulfido", 34: "selenido", 52: "tellurido"}
_ONIUM_PREFIX_STEMS = {
    8: ("oxidanium", 3),
    16: ("sulfanium", 3),
    34: ("selanium", 3),
    52: ("telluranium", 3),
    15: ("phosphanium", 4),
    33: ("arsanium", 4),
    51: ("stibanium", 4),
}


def _in_anionic_chain(mol, root, coming_from):
    z = mol.GetAtomWithIdx(root).GetAtomicNum()
    seen, stack = {root}, [root]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            i = n.GetIdx()
            if i != coming_from and i not in seen and n.GetAtomicNum() == z:
                seen.add(i)
                stack.append(i)
    return len(seen) > 1 and any(mol.GetAtomWithIdx(i).HasProp("_anion_word") for i in seen)


def _branch_has_anionic_carbon(graph, root, coming_from, mol):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n == coming_from or n in seen:
                continue
            seen.add(n)
            stack.append(n)
    return any(i != root and mol.GetAtomWithIdx(i).HasProp("_anion") and mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in seen)


def _anionic_chain_prefix(graph, root, coming_from, mol):
    """'ethan-1-id-2-yl': the carbon chain's parent anion with the free valence cited after the ide center."""
    from ._polyfunctional import _arm_atoms, _select, _unit_molecule

    atoms = _arm_atoms(graph, root, coming_from)
    unit, attach = _unit_molecule(mol, atoms, root)
    _, _, parts = _select(unit, attach)
    prefix, body, tail, locant = parts[:4]
    if not body.endswith("ide"):
        raise UnsupportedStructure("this anionic substituent chain is not named as a parent anion yet")
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    ending = {1.0: "yl", 2.0: "ylidene", 3.0: "ylidyne"}[order]
    cited = f"-{locant}-" if locant is not None else ""
    return f"{prefix}{body[:-1]}{cited}{ending}", True


def _carbon_anion_prefix(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """methanidyl / methanediidyl (and -ylidene/-ylidyne) with its own substituents (P-72.6.3)."""
    from ._substituents import format_mononuclear_prefixes, name_branch

    atom = mol.GetAtomWithIdx(root)
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 and not is_functional_carbon(mol, n) for n in graph[root] if n != coming_from):
        return _anionic_chain_prefix(graph, root, coming_from, mol)
    charge = int(atom.GetProp("_anion"))
    entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in graph[root] if n != coming_from]
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    stem = "methane" if charge > 1 else "methan"
    word = {1: "id", 2: "diid"}[charge]
    ending = {1.0: "yl", 2.0: "ylidene", 3.0: "ylidyne"}[order]
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    return f"{prefixes}{stem}{word}{ending}", bool(entries)


def _oxoacid_anion_prefix(graph, root, coming_from, mol):
    """phosphonato / arsonato: every hydroxy of the -E(=O)(OH)2 group deprotonated (P-72.6.1)."""
    kids = [n for n in graph[root] if n != coming_from]
    marked = [n for n in kids if mol.GetAtomWithIdx(n).HasProp("_anion_word")]
    oxo = [
        n
        for n in kids
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 8
        and mol.GetAtomWithIdx(n).GetDegree() == 1
        and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0
    ]
    if len(kids) == 3 and len(marked) == 2 and len(oxo) == 1:
        return {15: "phosphonato", 33: "arsonato"}[mol.GetAtomWithIdx(root).GetAtomicNum()], False
    raise UnsupportedStructure("this partly deprotonated oxoacid substituent is not supported yet")


def _anionic_group(mol, root, coming_from):
    atom = mol.GetAtomWithIdx(root)
    z = atom.GetAtomicNum()
    if z == 7 and atom.GetDegree() == 2 and mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() == 1.0:
        from ._common import adjacency
        from ._substituents import format_mononuclear_prefixes, name_branch

        graph = adjacency(mol)
        entries = [
            name_branch(graph, n, root, {}, frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()), mol=mol, unsaturated=True)
            for n in graph[root]
            if n != coming_from
        ]
        return format_mononuclear_prefixes(entries) + "azanidyl", True
    if mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() == 1.0 and atom.GetDegree() == 1:
        if z in ANIONIC_PREFIXES:
            return ANIONIC_PREFIXES[z], False
    raise UnsupportedStructure("this anionic substituent group is not supported yet")


_CHALCOGEN_REPLACEMENT = {16: "thio", 34: "seleno", 52: "telluro"}


def _terminal_chalcogenol(mol, idx, center):
    atom = mol.GetAtomWithIdx(idx)
    return (
        atom.GetAtomicNum() in (8, 16, 34, 52)
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(center, idx).GetBondTypeAsDouble() == 1.0
    )


def _borono_prefix(mol, root, others):
    """'borono' for (HO)2B-; its chalcogen analogues use replacement prefixes: 'thioborono', 'diselenoborono' (P-68.1.4.2)."""
    if len(others) != 2 or not all(_terminal_chalcogenol(mol, n, root) for n in others):
        return None
    counts = {}
    for n in others:
        z = mol.GetAtomWithIdx(n).GetAtomicNum()
        if z != 8:
            counts[_CHALCOGEN_REPLACEMENT[z]] = counts.get(_CHALCOGEN_REPLACEMENT[z], 0) + 1
    words = "".join((multiplying_prefix(count) if count > 1 else "") + word for word, count in sorted(counts.items()))
    return words + "borono"


def _terminal_hydroxy(mol, idx, center):
    atom = mol.GetAtomWithIdx(idx)
    return (
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(center, idx).GetBondTypeAsDouble() == 1.0
    )


def _terminal_double_atom(mol, center, idx):
    """A terminal =O, =S, =Se, =Te or =NH on `center` (an oxo, thioxo or imido position of a sulfonic-type acyl group)."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetDegree() != 1 or atom.GetFormalCharge() or mol.GetBondBetweenAtoms(center, idx).GetBondTypeAsDouble() != 2.0:
        return False
    return atom.GetAtomicNum() in (8, 16, 34, 52) or (atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() == 1)


def _stereo_cited(name, atom):
    """`name` with the CIP descriptor of the stereogenic skeletal `atom` of the group, as in '(S)-methanesulfinyl' (P-93.3.4.1)."""
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    if not context or atom not in context["atoms"]:
        return name
    context["used"].add(("atom", atom))
    plain = re.sub(r"^\([RS]\)-", "", name)
    return f"({context['atoms'][atom]})-{plain}"


def _sulfur_oxo_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """Prefixes of the acyl groups of sulfonic, sulfinic, selenonic ... acids (P-65.3.2): 'sulfo', 'sulfamoyl',
    'benzenesulfonyl', 'methoxysulfonyl', 'chlorosulfinyl', 'trithiosulfo'."""
    from ._acid_lexicon import acyl_suffix, make_spec
    from ._substituents import name_branch

    atom = mol.GetAtomWithIdx(root)
    if mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0 or atom.GetFormalCharge():
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    center = _E_SYMBOL[atom.GetAtomicNum()]
    others = [n for n in graph[root] if n != coming_from]
    oxo = [
        (n, "NH" if mol.GetAtomWithIdx(n).GetAtomicNum() == 7 else mol.GetAtomWithIdx(n).GetSymbol())
        for n in others
        if _terminal_double_atom(mol, root, n)
    ]
    rest = [n for n in others if n not in {n for n, _ in oxo}]
    if len(rest) != 1 or len(oxo) not in (1, 2):
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    x = rest[0]
    zx = mol.GetAtomWithIdx(x).GetAtomicNum()
    acid_amide = zx == 7 and all(e == "O" for _, e in oxo)
    if not EXTENDED_PREFIXES.get() and not acid_amide and ((center != "S" and zx != 6) or any(e != "O" for _, e in oxo)):
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    symbols = [e for _, e in oxo]
    if zx == 8 and mol.GetAtomWithIdx(x).GetDegree() == 1 and mol.GetBondBetweenAtoms(root, x).GetBondTypeAsDouble() == 1.0:
        base = _E_ACID_PREFIX[(center, len(oxo))]
        replaced = [e for e in symbols if e != "O"]
        if not replaced:
            return ("sulfonato" if base == "sulfo" and _is_anionic_oxygen(mol.GetAtomWithIdx(x)) else base), False
        if len(set(replaced)) == 1 and replaced[0] in _REPLACEMENT_PREFIX:
            return {1: "", 2: "di", 3: "tri"}[len(replaced)] + _REPLACEMENT_PREFIX[replaced[0]] + base, True
        if not any(e in ("NH", "NNH2") for e in replaced):
            raise UnsupportedStructure("a mixed chalcogen acid group on a substituent is not supported yet")
    xatom = mol.GetAtomWithIdx(x)
    if (
        zx in _REPLACEMENT_SYMBOL
        and xatom.GetDegree() == 1
        and xatom.GetTotalNumHs() == 1
        and not xatom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(root, x).GetBondTypeAsDouble() == 1.0
    ):
        replaced = [*(e for e in symbols if e != "O"), _REPLACEMENT_SYMBOL[zx]]
        if len(set(replaced)) == 1 and replaced[0] in _REPLACEMENT_PREFIX:
            return {1: "", 2: "di", 3: "tri"}[len(replaced)] + _REPLACEMENT_PREFIX[replaced[0]] + _E_ACID_PREFIX[(center, len(oxo))], True
        if not any(e in ("NH", "NNH2") for e in replaced):
            raise UnsupportedStructure("a mixed chalcogen acid group on a substituent is not supported yet")
    if zx == 7 and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 7 for n in graph[x] if n != root):
        from ._hetero_carboxylic import hydrazine_acyl_prefix

        hydrazine = hydrazine_acyl_prefix(
            mol, graph, x, root, halogens, acyl_suffix(make_spec(center, symbols, ("O",)), chain=False, count=1)
        )
        if hydrazine is not None:
            return hydrazine, True
    pseudohalide = zx == 7 and any(b.GetBondTypeAsDouble() > 1.0 for b in mol.GetAtomWithIdx(x).GetBonds())
    if center == "S" and all(e == "O" for e in symbols) and len(oxo) == 2 and zx == 7 and not pseudohalide:
        subs = [n for n in graph[x] if n != root]
        if not subs:
            return "sulfamoyl", False
        if mol.GetAtomWithIdx(x).IsInRing():
            return _ring_nitrogen_acyl(graph, x, root, halogens, aromatic_atoms, mol, "sulfonyl")
        name = _amino(_group_names(graph, mol, subs, x, halogens, aromatic_atoms))
        return _amino_stem(name) + "sulfamoyl", True
    if zx != 6 and not EXTENDED_PREFIXES.get() and not acid_amide:
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    spec = make_spec(center, symbols, ("O",))
    acyl = acyl_suffix(spec, chain=False, count=1)
    if zx == 6:
        from ._functional_prefixes import _CARBOXYLIC_CLASS, _acyl_prefix

        try:
            return _stereo_cited(_acyl_prefix(mol, _subtree(graph, root, coming_from), root, coming_from), root), True
        except UnsupportedStructure:
            if not mol.HasSubstructMatch(_CARBOXYLIC_CLASS):
                raise
    if zx == 8 and any(mol.GetAtomWithIdx(n).GetAtomicNum() in _ACID_CENTRE_NUMBERS for n in graph[x] if n != root):
        from ._oxoacid_acyl import oxoacid_chain_group

        chain = oxoacid_chain_group(mol, graph, root, coming_from, halogens, aromatic_atoms)
        if chain is not None:
            return chain
        mark(None, POLYACID_SUBSTITUENT_REASON)
    z_name, z_compound = name_branch(graph, x, root, halogens, aromatic_atoms, mol=mol)
    located = "S-" if zx == 7 and "NH" in symbols and center == "S" else ""
    return located + (_enclose(z_name, z_compound) if z_compound else z_name) + acyl, True


def _subtree(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


_OXOACID_GROUP = {
    15: ("phosphono", "phosphoryl", "phosphinoyl"),
    33: ("arsono", "arsoryl", "arsinoyl"),
    51: ("stibono", "stiboryl", "stibinoyl"),
}


def _chalcogen_acyl_oxo(mol, center, attached):
    """Whether the S/Se/Te `center` bonded to `attached` carries a terminal =O/=S/=Se/=Te (a sulfonyl-type acyl group)."""
    return any(
        n.GetIdx() != attached and _terminal_double_atom(mol, center, n.GetIdx())
        for n in mol.GetAtomWithIdx(center).GetNeighbors()
    )


def _pnictogen_oxo_group(graph, root, halogens, aromatic_atoms, mol, others):
    """P(=O)(X)(Y)- as 'phosphono' (X = Y = hydroxy), '(X)(Y)phosphoryl', or '(R)(R')phosphinoyl' when both are carbon
    groups (P-67.1.4.1.1.3, P-67.1.4.1.3)."""
    from ._oxoacid_acyl import infix_acyl_name, oxoacid_chain_group
    from ._substituents import format_mononuclear_prefixes, name_branch

    attach = next((n for n in graph[root] if n not in others), None)
    if attach is not None:
        chain = oxoacid_chain_group(mol, graph, root, attach, halogens, aromatic_atoms)
        if chain is not None:
            return chain
        infix = infix_acyl_name(mol, graph, root, attach, halogens, aromatic_atoms)
        if infix is not None:
            return infix
    oxo = [
        n
        for n in others
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 8
        and mol.GetAtomWithIdx(n).GetDegree() == 1
        and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in others if n not in oxo]
    hydrogens = mol.GetAtomWithIdx(root).GetTotalNumHs()
    if hydrogens and len(oxo) == 1 and len(rest) + hydrogens == 2 and all(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 1.0 for n in rest):
        mono, acyl, carbon_acyl = _OXOACID_GROUP[mol.GetAtomWithIdx(root).GetAtomicNum()]
        if not rest:
            return carbon_acyl, False
        entry = name_branch(graph, rest[0], root, halogens, aromatic_atoms, mol=mol)
        stem = carbon_acyl if mol.GetAtomWithIdx(rest[0]).GetAtomicNum() == 6 else mono[:-1] + "oyl"
        return format_mononuclear_prefixes([entry]) + stem, True
    if len(oxo) != 1 or len(rest) != 2 or any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in rest):
        return None
    mono, acyl, carbon_acyl = _OXOACID_GROUP[mol.GetAtomWithIdx(root).GetAtomicNum()]
    parts = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in rest]
    both_carbon = all(mol.GetAtomWithIdx(n).GetAtomicNum() == 6 for n in rest)
    name = phosphoryl_name(parts, (mono, carbon_acyl if both_carbon else acyl))
    return name, name != mono


LAMBDA_CENTRE_STEMS = {
    16: ("sulfanyl", 2), 34: ("selanyl", 2), 52: ("tellanyl", 2), 17: ("chloranyl", 1), 35: ("bromanyl", 1), 53: ("iodanyl", 1),
}


def _lambda_centre_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """λ4-sulfanyl, λ3-iodanyl, ... a chalcogen or halogen with only single bonds whose bonding number exceeds the
    standard one: '(dihydroxy-λ3-iodanyl)', '[bis(acetyloxy)-λ3-iodanyl]' (P-68.4.3, P-68.5.1)."""
    from ._substituents import format_mononuclear_prefixes, name_branch

    atom = mol.GetAtomWithIdx(root)
    stem, standard = LAMBDA_CENTRE_STEMS[atom.GetAtomicNum()]
    bonding = atom.GetTotalValence()
    others = [n for n in graph[root] if n != coming_from]
    if (
        bonding <= standard
        or (bonding - standard) % 2
        or atom.GetFormalCharge()
        or atom.GetIsotope()
        or atom.GetNumRadicalElectrons()
        or any(b.GetBondTypeAsDouble() != 1.0 for b in atom.GetBonds())
    ):
        return None
    if bonding > standard and not _has_senior_principal_group(mol):
        raise UnsupportedStructure("a λ-bonded centre outranks the carbon parent unless a senior group is present")
    entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    return (prefix + "-" if prefix else "") + f"λ{bonding}-{stem}", True


def _mononuclear_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """silyl, germyl, phosphanyl, boranyl, ... with organyl substituents:
    '(trimethylsilyl)', '[dimethyl(phenyl)silyl]' (P-29.3.1)."""
    from ._substituents import format_mononuclear_prefixes, name_branch

    atom = mol.GetAtomWithIdx(root)
    _, base, valence = MONONUCLEAR_HYDRIDES[atom.GetAtomicNum()]
    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
        raise UnsupportedStructure("this mononuclear group is not supported yet")
    order = int(mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble())
    others = [n for n in graph[root] if n != coming_from]
    if order > 1:
        if any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 and not _oxo_atom(mol, root, n) for n in others):
            raise UnsupportedStructure("this mononuclear ylidene group is not supported yet")
        bonding = atom.GetTotalValence()
        if bonding > valence and (bonding - valence) % 2 or bonding > valence + 4:
            raise UnsupportedStructure("this mononuclear ylidene group is not supported yet")
        entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
        prefix = format_mononuclear_prefixes(entries) if entries else ""
        if bonding > valence:
            return (prefix + "-" if prefix else "") + f"λ{bonding}-" + base[:-2] + SUFFIX_OF_ORDER[order], True
        return prefix + base[:-2] + SUFFIX_OF_ORDER[order], bool(entries)
    if any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == atom.GetAtomicNum() and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 1.0
        for n in others
    ):
        return _chain_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if atom.GetAtomicNum() in _OXANE_ELEMENTS and any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 2 and not mol.GetAtomWithIdx(n).GetFormalCharge()
        and any(m != root and mol.GetAtomWithIdx(m).GetAtomicNum() == atom.GetAtomicNum() for m in graph[n])
        for n in others
    ):
        return _oxanyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if atom.GetAtomicNum() in _OXOACID_GROUP:
        found = _pnictogen_oxo_group(graph, root, halogens, aromatic_atoms, mol, others)
        if found is not None:
            return found
    bonding = atom.GetTotalValence()
    if bonding > valence and ((bonding - valence) % 2 or bonding > valence + 4):
        raise UnsupportedStructure("this mononuclear group carries a multiple bond")
    if atom.GetAtomicNum() == 5 and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in others) and not any(
        mol.HasSubstructMatch(query) for query in ACIDS_SENIOR_TO_BORON
    ):
        raise UnsupportedStructure("a boron group with a hydroxy or alkoxy substituent is a boronic or borinic acid (P-67.1.1)")
    if atom.GetAtomicNum() == 5:
        borono = _borono_prefix(mol, root, others)
        if borono is not None:
            return borono, borono != "borono"
    entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    if bonding > valence:
        return (prefix + "-" if prefix else "") + f"λ{bonding}-" + base, True
    return prefix + base, bool(entries)


ACIDS_SENIOR_TO_BORON = [
    Chem.MolFromSmarts(smarts) for smarts in ("[CX3](=O)[OX2H1]", "[#16,#34,#52;X3,X4](=O)[OX2H1]")
]


def _aci_nitro_ligands(mol, nitrogen, others, order):
    """The hydroxy or alkoxy oxygen of an aci-nitro group =N(O)OR (neutral or as the zwitterion C=N+(O-)OR), else None."""
    atom = mol.GetAtomWithIdx(nitrogen)
    if order != 2 or len(others) != 2 or atom.GetFormalCharge() not in (0, 1):
        return None
    neutral = atom.GetFormalCharge() == 0
    oxo = [n for n in others if (_oxo_atom(mol, nitrogen, n) if neutral else _oxide_atom(mol, n))]
    rest = [n for n in others if n not in oxo]
    if len(oxo) != 1 or len(rest) != 1:
        return None
    ether = mol.GetAtomWithIdx(rest[0])
    if ether.GetAtomicNum() != 8 or ether.GetFormalCharge() or mol.GetBondBetweenAtoms(nitrogen, rest[0]).GetBondTypeAsDouble() != 1.0:
        return None
    return rest


def _oxide_atom(mol, n):
    atom = mol.GetAtomWithIdx(n)
    return atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetFormalCharge() == -1


def _oxo_atom(mol, centre, n):
    atom = mol.GetAtomWithIdx(n)
    return (
        atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and not atom.GetFormalCharge()
        and mol.GetBondBetweenAtoms(centre, n).GetBondTypeAsDouble() == 2.0
    )


def _hydroxy_oxygen(mol, n):
    atom = mol.GetAtomWithIdx(n)
    return atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1 and not atom.GetFormalCharge()


def _chain_stem(z):
    return "azane" if z == 7 else MONONUCLEAR_HYDRIDES[z][0]


def _chain_paths(graph, chain_atoms, root):
    paths = []

    def extend(path):
        grew = False
        for n in graph[path[-1]]:
            if n in chain_atoms and n not in path:
                grew = True
                extend(path + [n])
        if not grew and root in path:
            paths.append(path)

    for start in chain_atoms:
        if sum(m in chain_atoms for m in graph[start]) <= 1:
            extend([start])
    return paths


DIAZENYL_PREFIX = contextvars.ContextVar("diazenyl_prefix", default=False)


def _diazenyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol, order, others):
    """R-N=N- attached through the nitrogen next to the parent (P-29.3.2.2, P-68.3.1.3): diazenyl with the far
    nitrogen's substituent cited as a prefix."""
    if (
        order != 1.0
        or len(others) != 1
        or mol.GetAtomWithIdx(root).GetFormalCharge()
        or not (DIAZENYL_PREFIX.get() or _has_senior_principal_group(mol))
    ):
        return None
    far = others[0]
    far_atom = mol.GetAtomWithIdx(far)
    if far_atom.GetAtomicNum() != 7 or far_atom.GetFormalCharge() or mol.GetBondBetweenAtoms(root, far).GetBondTypeAsDouble() != 2.0:
        return None
    tail = [n for n in graph[far] if n != root]
    if not tail:
        return "diazenyl", False
    if len(tail) != 1 or mol.GetBondBetweenAtoms(far, tail[0]).GetBondTypeAsDouble() != 1.0:
        return None
    from ._substituents import name_branch

    rname, rcomp = name_branch(graph, tail[0], far, halogens, aromatic_atoms, mol=mol)
    return _enclose(rname, rcomp) + "diazenyl", True


def _nitro_or_nitroso_nitrogen(mol, idx, attached_to):
    if is_nitro_nitrogen(mol, idx):
        return True
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetDegree() != 2 or atom.GetFormalCharge():
        return False
    oxygens = [n for n in atom.GetNeighbors() if n.GetIdx() != attached_to]
    return (
        len(oxygens) == 1
        and oxygens[0].GetAtomicNum() == 8
        and oxygens[0].GetDegree() == 1
        and mol.GetBondBetweenAtoms(idx, oxygens[0].GetIdx()).GetBondTypeAsDouble() == 2.0
    )


def _chain_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """disilanyl, triazan-1-yl, 3-silyltetrasilan-1-yl, 1-methyltetrasilan-1-yl: a homogeneous heteroatom chain
    attached through one of its atoms; the longest chain through the free valence is the parent (P-29.4.1, P-44.3)."""
    from ._common import group_substituents, substituent_locant_set_and_citation
    from ._numerals import multiplying_prefix
    from ._substituents import format_substituent_prefixes, name_branch

    z = mol.GetAtomWithIdx(root).GetAtomicNum()
    stem = _chain_stem(z)
    chain_atoms = {root}
    stack = [root]
    while stack:
        current = stack.pop()
        for n in graph[current]:
            if (
                n != coming_from
                and n not in chain_atoms
                and mol.GetAtomWithIdx(n).GetAtomicNum() == z
                and mol.GetBondBetweenAtoms(current, n).GetBondTypeAsDouble() in (1.0, 2.0)
                and not (z == 7 and _nitro_or_nitroso_nitrogen(mol, n, current))
            ):
                chain_atoms.add(n)
                stack.append(n)
    if any(mol.GetAtomWithIdx(a).IsInRing() or mol.GetAtomWithIdx(a).GetFormalCharge() for a in chain_atoms):
        raise UnsupportedStructure("this heteroatom chain substituent is not supported yet")
    paths = _chain_paths(graph, chain_atoms, root)
    longest = max(len(p) for p in paths)
    best = None
    for walk in (p for p in paths if len(p) == longest):
        subs = {}
        for i, atom in enumerate(walk):
            for n in graph[atom]:
                if n in walk or n == coming_from:
                    continue
                subs.setdefault(i + 1, []).append(name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol))
        grouped = group_substituents(subs)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        attach = walk.index(root) + 1
        ene = tuple(i + 1 for i in range(len(walk) - 1) if mol.GetBondBetweenAtoms(walk[i], walk[i + 1]).GetBondTypeAsDouble() == 2.0)
        lam = {i + 1: n for i, a in enumerate(walk) if (n := nonstandard_bonding(mol.GetAtomWithIdx(a)))}
        key = ((attach,), ene, tuple(sorted(lam)), tuple(-lam[p] for p in sorted(lam)), -sum(len(v) for v in subs.values()), locant_set, citation)
        if best is None or key < best[0]:
            best = (key, grouped, attach, ene, lam)
    _, grouped, attach, ene, lam = best
    if len(ene) > 1:
        raise UnsupportedStructure("a heteroatom chain with several double bonds is not supported yet")
    lam_text = ",".join(f"{p}λ{lam[p]}" for p in sorted(lam))
    if lam_text and not ene:
        word = multiplying_prefix(longest) + stem[:-1]
        base = f"{lam_text}-{word}-{attach}-yl"
        prefix = format_substituent_prefixes(grouped) if grouped else ""
        return prefix + ("-" if prefix else "") + base, True
    if ene:
        base = f"{multiplying_prefix(longest)}{stem[:-3]}-{ene[0]}-en-{attach}-yl"
        prefix = format_substituent_prefixes(grouped) if grouped else ""
        return prefix + base, True
    word = multiplying_prefix(longest) + stem[:-1]
    base = word + "yl" if longest == 2 and attach == 1 else f"{word}-{attach}-yl"
    if z in _GROUP_13_CHAIN and not ene:
        base = f"{word}({longest + 2})-{attach}-yl"
    ylidene_only = bool(grouped) and all(name.endswith("ylidene") for name in grouped)
    if z == 7 and longest == 2 and attach == 1:
        # P-29.3.1: a substituted hydrazinyl cites its free valence ('2-phenylhydrazin-1-yl'); an ylidene group can only
        # sit on the second nitrogen, so no locant is cited (P-66.3.6)
        base = "hydrazinyl" if not grouped or ylidene_only else "hydrazin-1-yl"
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if base == "hydrazinyl" and ylidene_only:
        prefix = format_substituent_prefixes(grouped, omit_all=True)
    return prefix + base, bool(prefix) or "-" in base


_OXANE_ELEMENTS = {5, 13, 14, 31, 32, 49, 50, 81, 82}
_GROUP_13_CHAIN = {5, 13, 31, 49, 81}


def _oxanyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """disiloxanyl, diboroxanyl, trisiloxan-1-yl: an unbranched E-O-E... chain attached through a terminal atom E (P-68.1.2)."""
    from ._common import group_substituents
    from ._hydride_chain import _alternating_parent
    from ._substituents import format_substituent_prefixes, name_branch

    element = mol.GetAtomWithIdx(root).GetAtomicNum()
    walk, previous = [root], coming_from
    while True:
        bridges = [
            n
            for n in graph[walk[-1]]
            if n != previous
            and mol.GetAtomWithIdx(n).GetAtomicNum() == 8
            and any(m != walk[-1] and mol.GetAtomWithIdx(m).GetAtomicNum() == element for m in graph[n])
        ]
        if not bridges:
            break
        if len(bridges) > 1:
            raise UnsupportedStructure("a branched oxane substituent is not supported yet")
        (centre,) = [m for m in graph[bridges[0]] if m != walk[-1]]
        if mol.GetAtomWithIdx(bridges[0]).GetDegree() != 2:
            raise UnsupportedStructure("this oxane substituent is not supported yet")
        walk += [bridges[0], centre]
        previous = bridges[0]
    chain = set(walk)
    subs = {}
    for i, atom in enumerate(walk):
        for n in graph[atom]:
            if n in chain or n == coming_from:
                continue
            if mol.GetAtomWithIdx(n).GetAtomicNum() not in (6,) and n not in halogens:
                raise UnsupportedStructure("this oxane substituent carries something other than organyl groups")
            subs.setdefault(i + 1, []).append(name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol))
    grouped = group_substituents(subs)
    count = (len(walk) + 1) // 2
    parent = _alternating_parent(element, 8, count)
    base = parent[:-1] + "yl" if count == 2 else parent[:-1] + "-1-yl"
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return prefix + base, bool(prefix) or "-" in base
