"""Substituent prefixes joined through a heteroatom or a functional carbon
(P-29.3.3, P-29.4, P-35, P-63.2, P-66): hydroxy, oxo, alkoxy/aryloxy, sulfanyl,
amino, nitro, cyano, formyl, carboxy, carbamoyl, alkoxycarbonyl and acyl.
`name_branch` calls `hetero_branch_name` before it walks a carbon chain.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, alpha_sort_key
from ._free_valence import SUFFIX_OF_ORDER
from ._multiplicative_text import enclose
from ._numerals import alkane_name

_ALKOXY_STEMS = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy", "phenyl": "phenoxy"}
_SIMPLE_NAMES = {
    "methoxy", "ethoxy", "propoxy", "butoxy", "tert-butoxy", "phenoxy", "amino", "anilino", "hydroxy", "oxo",
    "nitro", "nitroso", "cyano", "sulfanyl", "formyl", "carboxy", "carbamoyl",
}
_MULTIPLE_TARGETS = {7, 8, 16}
CHALCOGEN_PREFIXES = {16: "sulfanyl", 34: "selanyl", 52: "tellanyl"}


_SENIOR_TO_SELENOL = [
    Chem.MolFromSmarts(smarts)
    for smarts in (
        "[CX3](=O)[OX2H1]",
        "[CX3](=O)[OX2][#6]",
        "[CX3](=O)[NX3]",
        "[CX2]#[NX1]",
        "[CX3H1](=O)[#6]",
        "[#6][CX3](=O)[#6]",
        "[OX2H1][#6;!$([#6]=O)]",
        "[SX2H1][#6;!$([#6]=[O,S,Se,Te])]",
    )
]


def require_plain_chalcogen_kids(mol, z, kids):
    """Se/Te bonded to an acyl, carbamoyl, formyl or cyano carbon is a selenoate/selenocyanate-type group, not a
    plain selanyl/tellanyl prefix."""
    if z in (34, 52):
        for kid in kids:
            atom = mol.GetAtomWithIdx(kid)
            if atom.GetAtomicNum() != 6 or any(
                b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(atom).GetAtomicNum() in (7, 8, 16, 34, 52)
                for b in atom.GetBonds()
            ):
                raise UnsupportedStructure("an acyl, carbamoyl or cyano group on selenium/tellurium is not a selanyl prefix")


def _has_senior_principal_group(mol):
    """A principal group senior to the hetero-hetero connection (hydroxylamine, hydrazine, peroxide classes) is
    present, so that connection is expressed as a prefix (P-41, P-29.4.1)."""
    return any(mol.HasSubstructMatch(query) for query in _SENIOR_TO_SELENOL)


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
            return True
        if b.GetBondTypeAsDouble() == 2.0 and other.GetAtomicNum() in _MULTIPLE_TARGETS:
            return any(
                n.GetIdx() != other.GetIdx() and n.GetAtomicNum() in _MULTIPLE_TARGETS for n in atom.GetNeighbors()
            )
    return False


def _enclose(name, compound):
    return enclose(name) if compound else name


def phosphoryl_name(parts):
    """'phosphono' for P(O)(OH)2, otherwise '(X)(Y)phosphoryl' with X, Y cited alphabetically (P-65.1.3.1)."""
    from ._numerals import multiplying_prefix

    if all(name == "hydroxy" for name, _ in parts):
        return "phosphono"
    if len(set(parts)) == 1:
        name, compound = parts[0]
        body = multiplying_prefix(2, compound=compound) + (enclose(name) if compound else name)
    else:
        ordered = sorted(parts, key=lambda p: alpha_sort_key(p[0]))
        body = "".join(enclose(n) if c else n for n, c in ordered)
    return body + "phosphoryl"


def _phosphoryloxy(graph, phosphorus, oxygen, halogens, aromatic_atoms, mol):
    """'phosphonooxy' or '[(X)(Y)phosphoryl]oxy' for O-P(=O)(OX)(OY)."""
    from ._substituents import name_branch

    atom = mol.GetAtomWithIdx(phosphorus)
    others = [n for n in graph[phosphorus] if n != oxygen]
    terminal = [
        n for n in others
        if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 1
        and mol.GetBondBetweenAtoms(phosphorus, n).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in others if n not in terminal]
    if (
        atom.GetFormalCharge() or atom.IsInRing() or len(terminal) != 1 or len(rest) != 2
        or any(mol.GetAtomWithIdx(n).GetAtomicNum() != 8 or mol.GetBondBetweenAtoms(phosphorus, n).GetBondTypeAsDouble() != 1.0 for n in rest)
    ):
        raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
    parts = [name_branch(graph, n, phosphorus, halogens, aromatic_atoms, mol=mol) for n in rest]
    if any("phospho" in name for name, _ in parts):
        raise UnsupportedStructure("a polyphosphate chain substituent is not supported yet")
    group = phosphoryl_name(parts)
    return ("phosphonooxy" if group == "phosphono" else enclose(group) + "oxy"), True


def _alkoxy(rname):
    if rname.startswith("("):
        return enclose(rname) + "oxy"
    if rname == "tert-butyl":
        return "tert-butoxy"
    for stem, short in _ALKOXY_STEMS.items():
        if rname.endswith(stem) and not rname.endswith("cyclo" + stem):
            return rname[: -len(stem)] + short
    if rname[0].isdigit():
        return enclose(rname) + "oxy"
    return rname + "oxy"


def _compound(name):
    return name not in _SIMPLE_NAMES


def _acyl_name(mol, graph, carbon, from_atom, halogens=None, aromatic_atoms=None):
    """alkanoyl or benzoyl prefix for R-C(=O)-: an unbranched saturated chain
    whose carbons may carry substituents (2-aminoethanoyl), or phenyl."""
    from ._substituents import name_branch

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
    from ._substituents import format_substituent_prefixes

    prefix = format_substituent_prefixes(group_substituents(entries), omit_locants=len(chain) == 1) if entries else ""
    return prefix + alkane_name(len(chain))[:-1] + "oyl"


def _acyl_from_acid_name(mol, graph, carbon, from_atom):
    """'(9Z)-octadec-9-enoyl' from the name of the acid whose acyl group is rooted at `carbon` (P-65.1.7.1)."""
    from ._functional_prefixes import _acyl_prefix
    from ._substituents import BRANCH_STEREO

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
    context = BRANCH_STEREO.get()
    if context:
        context["used"].update(("atom", a) for a in atoms)
        context["used"].update(
            ("bond", (b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
            for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
        )
    return name


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


def _amino(names):
    if not names:
        return "amino"
    if names == [("phenyl", False)]:
        return "anilino"
    ordered = sorted(names, key=lambda item: alpha_sort_key(item[0]))
    if len(ordered) == 2 and ordered[0] == ordered[1]:
        name, compound = ordered[0]
        return (f"bis({name})" if compound else "di" + name) + "amino"
    parts = [_enclose(ordered[0][0], ordered[0][1])]
    parts += [f"({n})" for n, _ in ordered[1:]]
    return "".join(parts) + "amino"


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
        if is_functional_carbon(mol, root) or _carbonyl_oxygen(mol, root) is not None:
            return _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol)
        return None
    if z in (15, 33) and any(mol.GetAtomWithIdx(n).HasProp("_anion_word") for n in graph[root]):
        return _oxoacid_anion_prefix(graph, root, coming_from, mol)
    if z in MONONUCLEAR_HYDRIDES:
        return _mononuclear_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if atom.GetFormalCharge() and z != 7:
        raise UnsupportedStructure("a charged atom in a substituent group is not supported yet")
    if atom.IsInRing():
        return None
    others = [n for n in graph[root] if n != coming_from]
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    if z == 8:
        if order == 2.0:
            return "oxo", False
        if not others:
            return "hydroxy", False
        (other,) = others
        if mol.GetAtomWithIdx(other).GetAtomicNum() == 15 and any(
            mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(other, n).GetBondTypeAsDouble() == 2.0
            for n in graph[other]
        ):
            return _phosphoryloxy(graph, other, root, halogens, aromatic_atoms, mol)
        if mol.GetAtomWithIdx(other).GetAtomicNum() in MONONUCLEAR_HYDRIDES:
            silyl, _ = _mononuclear_group(graph, other, root, halogens, aromatic_atoms, mol)
            return _enclose(silyl, True) + "oxy", True
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

        rname, _ = name_branch(graph, other, root, halogens, aromatic_atoms, mol=mol)
        name = _alkoxy(rname)
        return name, _compound(name)
    if z == 16 and atom.GetDegree() > 2:
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
            silyl, _ = _mononuclear_group(graph, others[0], root, halogens, aromatic_atoms, mol)
            return _enclose(silyl, True) + word, True
        from ._substituents import name_branch

        if mol.GetAtomWithIdx(others[0]).GetAtomicNum() != 6:
            if mol.GetAtomWithIdx(others[0]).GetAtomicNum() == z or not _has_senior_principal_group(mol):
                raise UnsupportedStructure("a chalcogen chain (disulfanyl, ...) is not supported yet")
            rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
            return _enclose(rname, rcomp) + word, True

        rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
        return _enclose(rname, rcomp) + word, True
    if z == 7 and atom.GetFormalCharge() == 1 and order == 1.0 and atom.GetDegree() + atom.GetTotalNumHs() == 4:
        if any(
            mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0
            for n in others
        ):
            raise UnsupportedStructure("this azaniumyl group is not supported yet")
        from ._substituents import format_mononuclear_prefixes

        entries = _group_names(graph, mol, others, root, halogens, aromatic_atoms)
        return (format_mononuclear_prefixes(entries) if entries else "") + "azaniumyl", bool(entries)
    if z == 7:
        oxygens = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
        if len(oxygens) == len(others) and others:
            if len(oxygens) == 2 and atom.GetTotalDegree() == 3:
                return "nitro", False
            if len(oxygens) == 1 and mol.GetBondBetweenAtoms(root, oxygens[0]).GetBondTypeAsDouble() == 2.0:
                return "nitroso", False
        diazenyl = _diazenyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol, order, others)
        if diazenyl is not None:
            return diazenyl
        if order == 2.0 and len(others) == 1 and mol.GetAtomWithIdx(others[0]).GetAtomicNum() == 7:
            far = mol.GetAtomWithIdx(others[0])
            if far.GetDegree() == 1 and far.GetTotalNumHs() == 2 and not far.GetFormalCharge() and not atom.GetFormalCharge():
                return "hydrazinylidene", False
        if order == 2.0 and not atom.GetFormalCharge() and len(others) <= 1 and _has_senior_principal_group(mol):
            if not others:
                return "imino", False
            if mol.GetBondBetweenAtoms(root, others[0]).GetBondTypeAsDouble() == 1.0:
                from ._substituents import name_branch

                rname, rcomp = name_branch(graph, others[0], root, halogens, aromatic_atoms, mol=mol)
                return _enclose(rname, rcomp) + "imino", True
        if order == 1.0 and not atom.GetFormalCharge() and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 7 for n in others):
            far = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 7]
            if not (len(others) == 1 and mol.GetAtomWithIdx(far[0]).GetDegree() == 1):
                return _chain_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
        if order != 1.0 or atom.GetFormalCharge() or any(
            mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others
        ):
            raise UnsupportedStructure("this nitrogen-linked group is not supported yet")
        if len(others) == 1 and mol.GetAtomWithIdx(others[0]).GetAtomicNum() == 7:
            far = mol.GetAtomWithIdx(others[0])
            if far.GetDegree() == 1 and far.GetTotalNumHs() == 2 and not far.GetFormalCharge():
                return "hydrazinyl", False
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
        name = _amino(_group_names(graph, mol, others, root, halogens, aromatic_atoms))
        return name, _compound(name)
    raise UnsupportedStructure("this heteroatom-linked substituent is not supported yet")


def _functional_carbon(graph, root, coming_from, halogens, aromatic_atoms, mol):
    from ._substituents import name_branch

    atom = mol.GetAtomWithIdx(root)
    others = [n for n in graph[root] if n != coming_from]
    triple_n = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 3.0]
    if triple_n and len(others) == 1:
        return "cyano", False
    carbonyl = [n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0]
    if len(carbonyl) != 1 or mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")
    rest = [n for n in others if n != carbonyl[0]]
    if not rest:
        return "formyl", False
    (x,) = rest
    z = mol.GetAtomWithIdx(x).GetAtomicNum()
    if z == 8:
        tail = [n for n in graph[x] if n != root]
        if not tail:
            return ("carboxylato" if mol.GetAtomWithIdx(x).HasProp("_anion") else "carboxy"), False
        rname, _ = name_branch(graph, tail[0], x, halogens, aromatic_atoms, mol=mol)
        return _alkoxy(rname) + "carbonyl", True
    if z == 7:
        if mol.GetAtomWithIdx(x).HasProp("_anion_word"):
            raise UnsupportedStructure("an anionic amide nitrogen is not named as a carbamoyl prefix")
        subs = [n for n in graph[x] if n != root]
        if not subs:
            return "carbamoyl", False
        name = _amino(_group_names(graph, mol, subs, x, halogens, aromatic_atoms))
        return name[: -len("amino")] + "carbamoyl", True
    if z == 6:
        try:
            name = _acyl_name(mol, graph, root, coming_from, halogens, aromatic_atoms)
        except UnsupportedStructure:
            name = _acyl_from_acid_name(mol, graph, root, coming_from)
        return name, any(ch.isdigit() for ch in name) or "(" in name or "[" in name
    raise UnsupportedStructure("this carbonyl-derived substituent is not supported yet")


ANIONIC_PREFIXES = {7: "azanidyl", 8: "oxido", 16: "sulfido", 34: "selenido", 52: "tellurido"}


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


def _sulfur_oxo_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """-SO3H (sulfo), -SO2NR2 (sulfamoyl), -SO2R (R-sulfonyl), -SOR (R-sulfinyl)."""
    from ._substituents import name_branch

    atom = mol.GetAtomWithIdx(root)
    if mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble() != 1.0 or atom.GetFormalCharge():
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    others = [n for n in graph[root] if n != coming_from]
    oxygens = [
        n for n in others if mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() == 2.0
    ]
    rest = [n for n in others if n not in oxygens]
    if len(rest) != 1 or len(oxygens) not in (1, 2):
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    x = rest[0]
    zx = mol.GetAtomWithIdx(x).GetAtomicNum()
    if len(oxygens) == 2 and zx == 8 and mol.GetAtomWithIdx(x).GetDegree() == 1:
        return ("sulfonato" if mol.GetAtomWithIdx(x).HasProp("_anion") else "sulfo"), False
    if len(oxygens) == 2 and zx == 7:
        subs = [n for n in graph[x] if n != root]
        if not subs:
            return "sulfamoyl", False
        name = _amino(_group_names(graph, mol, subs, x, halogens, aromatic_atoms))
        return name[: -len("amino")] + "sulfamoyl", True
    if zx != 6:
        raise UnsupportedStructure("this sulfur-linked group is not supported yet")
    rname, _ = name_branch(graph, x, root, halogens, aromatic_atoms, mol=mol)
    if rname == "phenyl":
        stem = "benzene"
    elif rname.endswith("yl") and rname[:-2] in ("meth", "eth", "prop", "but", "pent", "hex", "hept", "oct"):
        stem = rname[:-2] + "ane"
    else:
        raise UnsupportedStructure("this sulfonyl group is not supported yet")
    return stem + ("sulfonyl" if len(oxygens) == 2 else "sulfinyl"), True


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
        if any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others) or len(others) > valence - order:
            raise UnsupportedStructure("this mononuclear ylidene group is not supported yet")
        entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
        prefix = format_mononuclear_prefixes(entries) if entries else ""
        return prefix + base[:-2] + SUFFIX_OF_ORDER[order], bool(entries)
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() == atom.GetAtomicNum() for n in others):
        return _chain_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if atom.GetAtomicNum() == 14 and any(
        mol.GetAtomWithIdx(n).GetAtomicNum() == 8 and mol.GetAtomWithIdx(n).GetDegree() == 2 and not mol.GetAtomWithIdx(n).GetFormalCharge()
        and any(m != root and mol.GetAtomWithIdx(m).GetAtomicNum() == 14 for m in graph[n])
        for n in others
    ):
        return _siloxanyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol)
    if len(others) > valence - 1 or any(mol.GetBondBetweenAtoms(root, n).GetBondTypeAsDouble() != 1.0 for n in others):
        raise UnsupportedStructure("this mononuclear group carries a multiple bond")
    if atom.GetAtomicNum() == 5 and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in others):
        raise UnsupportedStructure("a boron group with a hydroxy or alkoxy substituent is a boronic or borinic acid (P-67.1.1)")
    entries = [name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol) for n in others]
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    return prefix + base, bool(entries)


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


def _diazenyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol, order, others):
    """R-N=N- attached through the nitrogen next to the parent (P-29.3.2.2, P-68.3.1.3): diazenyl with the far
    nitrogen's substituent cited as a prefix."""
    if order != 1.0 or len(others) != 1 or mol.GetAtomWithIdx(root).GetFormalCharge() or not _has_senior_principal_group(mol):
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
        for n in graph[stack.pop()]:
            if n != coming_from and n not in chain_atoms and mol.GetAtomWithIdx(n).GetAtomicNum() == z:
                chain_atoms.add(n)
                stack.append(n)
    if any(
        mol.GetAtomWithIdx(a).IsInRing() or mol.GetAtomWithIdx(a).GetFormalCharge() or any(
            mol.GetBondBetweenAtoms(a, m).GetBondTypeAsDouble() != 1.0 for m in graph[a] if m in chain_atoms or m == coming_from
        )
        for a in chain_atoms
    ):
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
                if mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure("a multiple bond on a heteroatom chain substituent is not supported yet")
                subs.setdefault(i + 1, []).append(name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol))
        grouped = group_substituents(subs)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        attach = walk.index(root) + 1
        key = ((attach,), -sum(len(v) for v in subs.values()), locant_set, citation)
        if best is None or key < best[0]:
            best = (key, grouped, attach)
    _, grouped, attach = best
    word = multiplying_prefix(longest) + stem[:-1]
    base = word + "yl" if longest == 2 and attach == 1 else f"{word}-{attach}-yl"
    if z == 7 and longest == 2 and attach == 1:
        base = "hydrazinyl"
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return prefix + base, bool(prefix) or "-" in base


def _siloxanyl_group(graph, root, coming_from, halogens, aromatic_atoms, mol):
    """disiloxanyl, trisiloxan-1-yl: an unbranched Si-O-Si... chain attached through a terminal silicon."""
    from ._common import group_substituents
    from ._numerals import multiplying_prefix
    from ._substituents import format_substituent_prefixes, name_branch

    walk, previous = [root], coming_from
    while True:
        bridges = [
            n
            for n in graph[walk[-1]]
            if n != previous
            and mol.GetAtomWithIdx(n).GetAtomicNum() == 8
            and any(m != walk[-1] and mol.GetAtomWithIdx(m).GetAtomicNum() == 14 for m in graph[n])
        ]
        if not bridges:
            break
        if len(bridges) > 1:
            raise UnsupportedStructure("a branched siloxane substituent is not supported yet")
        (silicon,) = [m for m in graph[bridges[0]] if m != walk[-1]]
        if mol.GetAtomWithIdx(bridges[0]).GetDegree() != 2:
            raise UnsupportedStructure("this siloxane substituent is not supported yet")
        walk += [bridges[0], silicon]
        previous = bridges[0]
    chain = set(walk)
    subs = {}
    for i, atom in enumerate(walk):
        for n in graph[atom]:
            if n in chain or n == coming_from:
                continue
            if mol.GetAtomWithIdx(n).GetAtomicNum() not in (6,) and n not in halogens:
                raise UnsupportedStructure("this siloxane substituent carries something other than organyl groups")
            subs.setdefault(i + 1, []).append(name_branch(graph, n, atom, halogens, aromatic_atoms, mol=mol))
    grouped = group_substituents(subs)
    silicon = (len(walk) + 1) // 2
    base = multiplying_prefix(silicon) + "siloxan" + ("yl" if silicon == 2 else "-1-yl")
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    return prefix + base, bool(prefix) or "-" in base
