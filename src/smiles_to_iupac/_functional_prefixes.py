"""Prefix names for heteroatom groups and acyl groups hanging off a skeleton (P-63, P-65.1.7), computed
bottom-up so `name_branch`/`_name_acyl_part` can cite them as plain leaves: hydroxy, oxo, amino,
(alkyl)amino, nitro, cyano, sulfanyl, alkoxy, alkylsulfanyl and acyl/acyloxy (acid name '-ic acid' to '-yl').
"""

import re

from rdkit import Chem

from ._alkoxy import alkoxy_prefix
from ._anilino import anilino_prefix
from ._amino_acid import SYSTEMATIC_ACID_PROBE
from ._common import HALOGEN_PREFIXES, UnsupportedStructure, alpha_sort_key, is_nitro_nitrogen
from ._hetero_prefixes import (
    DIPOLAR_GROUPS,
    POLYACID_SUBSTITUENT_REASON,
    ANIONIC_PREFIXES,
    CHALCOGEN_PREFIXES,
    MONONUCLEAR_HYDRIDES,
    azoxy_atoms,
    _functional_carbon,
    IMINE_PARENT,
    _has_senior_principal_group,
    _imidoyl_centre,
    _thioacyl,
    diazo_prefix,
    nitrogen_pseudohalide_prefix,
    phosphoryl_name,
    require_plain_chalcogen_kids,
)
from ._imidoyl_prefix import ketene_prefix
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix
from ._pin import mark
from ._retained_acids import is_compound_acyl
from ._substituents import FORCED_BRANCH_NAMES, _terminal_amide, name_branch

_NATIVE_ROOTS = frozenset({6, 7, 8, 9, 16, 17, 34, 35, 52, 53})


def _anionic_branch(mol, root, branch):
    """A substituent whose anionic centre is named by the anion prefixes (P-72.6.3), not by the functional-group
    prefixes: a nitrogen or a carbon without substituents at the root, or a centre beyond the root."""
    root_atom = mol.GetAtomWithIdx(root)
    if root_atom.HasProp("_anion"):
        return root_atom.GetAtomicNum() == 7 or root_atom.GetDegree() == 1
    return any(mol.GetAtomWithIdx(a).HasProp("_anion") for a in branch)


def _acylamino_nitrogen(mol, root, parent):
    """A neutral acyclic nitrogen joined singly to `parent` that carries an acyl group (acetamido, sulfonamido,
    ethanimidamido; P-66.1.1.4.3) and at most one carbon group, or two acyl groups (N-acetylacetamido; P-66.1.2)."""
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() or atom.IsInRing() or atom.GetDegree() > 3:
        return False
    if mol.GetBondBetweenAtoms(root, parent).GetBondTypeAsDouble() != 1.0:
        return False
    kids = [n for n in atom.GetNeighbors() if n.GetIdx() != parent]
    acyl = [n for n in kids if _thioacyl(mol, n.GetIdx()) or _imidoyl_centre(mol, n.GetIdx())]
    return len(acyl) in (1, 2) and len(kids) <= 2 and all(n in acyl or n.GetAtomicNum() == 6 for n in kids)


def _nitrogen_on_ring_nitrogen(mol, root, parent):
    """A neutral acyclic nitrogen joined singly to a ring nitrogen: its group (diazenyl, hydrazinyl, ...) is one prefix of
    the ring (P-68.3.1.2, P-68.3.1.3)."""
    atom = mol.GetAtomWithIdx(root)
    return (
        atom.GetAtomicNum() == 7
        and not atom.GetFormalCharge()
        and not atom.IsInRing()
        and mol.GetAtomWithIdx(parent).GetAtomicNum() == 7
        and mol.GetAtomWithIdx(parent).IsInRing()
        and mol.GetBondBetweenAtoms(root, parent).GetBondTypeAsDouble() == 1.0
    )


def _acyl_hydrazine_sulfur(mol, root):
    """A sulfonyl-type S, Se or Te group whose singly bonded nitrogen carries a second nitrogen (hydrazinesulfinyl)."""
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() not in (16, 34, 52):
        return False
    has_oxo = any(
        n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(root, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )
    return has_oxo and any(
        n.GetAtomicNum() == 7 and any(m.GetAtomicNum() == 7 for m in n.GetNeighbors() if m.GetIdx() != root)
        for n in atom.GetNeighbors()
    )


_enclose = enclose


def _is_compound(name):
    return any(ch.isdigit() for ch in name) or "(" in name


def _bond_order(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def nitro_atoms(mol):
    found = set()
    for atom in mol.GetAtoms():
        if is_nitro_nitrogen(mol, atom.GetIdx()):
            found.add(atom.GetIdx())
            found.update(n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8)
    return found


def nitroso_atoms(mol):
    """The nitrogen and oxygen of every -N=O group bonded through the nitrogen to one other atom."""
    found = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetFormalCharge():
            continue
        oxygens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(oxygens) == 1:
            found.update((atom.GetIdx(), oxygens[0].GetIdx()))
    return found


_FREE_CARBOXY = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3;!R](=[O,S])[O,N,S,F,Cl,Br,I]),$([CX3H1;!R]=O)]")
_MULTIPLIED_ACID = re.compile(r"(?:di|tri|tetra)carboxylic acid$|(?:di|tri|tetra)oic acid$")
_CARBOXYLIC_CLASS = Chem.MolFromSmarts("[$([CX3](=[O,S,N])[N,O,F,Cl,Br,I]),$([CX2]#N)]")


def _restore_cation(editable, index):
    """The cationic amine that a surrogate name stands for: the diazonium group N#N+ or the ammonium group."""
    atom = editable.GetAtomWithIdx(index)
    atom.SetFormalCharge(1)
    atom.SetNoImplicit(True)
    if atom.HasProp("_diazonium_surrogate"):
        atom.SetNumExplicitHs(0)
        editable.AddBond(index, editable.AddAtom(Chem.Atom(7)), Chem.BondType.TRIPLE)
    else:
        atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)


def _acyl_prefix(mol, subtree, root, parent):
    """'-yl' name of the acyl group rooted at `root`, from its parent acid's name."""
    from .core import smiles_to_iupac

    from ._isotope_labels import restore_labels
    from ._substituents import ISOTOPE_LABELS

    rw = Chem.RWMol(mol)
    context = ISOTOPE_LABELS.get()
    centre = {n.GetIdx() for n in mol.GetAtomWithIdx(root).GetNeighbors() if n.GetDegree() == 1}
    labelled = {a: e for a, e in context["labels"].items() if a in subtree and a not in centre} if context else {}
    restore_labels(rw, mol, labelled)
    new_o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(root, new_o, Chem.BondType.SINGLE)
    cationic_amines = [a for a in subtree if mol.GetAtomWithIdx(a).HasProp("_cationic_amine") and not mol.GetAtomWithIdx(a).GetFormalCharge()]
    for index in cationic_amines:
        _restore_cation(rw, index)
    # an ammonium group in the acid outranks it as a parent, so the acid is named as its anion beside the cation
    zwitterion = bool(cationic_amines) or any(mol.GetAtomWithIdx(a).GetFormalCharge() > 0 for a in subtree)
    if zwitterion:
        rw.GetAtomWithIdx(new_o).SetFormalCharge(-1)
    for idx in sorted(set(range(mol.GetNumAtoms())) - subtree, reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    try:
        Chem.SanitizeMol(sub)
    except Chem.MolSanitizeException as error:
        raise UnsupportedStructure("the acyl group cannot be cut from its ring system as a valid molecule") from error
    root_atom = mol.GetAtomWithIdx(root)
    own_acid = None
    if root_atom.GetAtomicNum() != 6 and sub.HasSubstructMatch(_CARBOXYLIC_CLASS):
        own_acid = _imido_acid_class(sub, sorted(subtree).index(root))
        if own_acid is None:
            raise UnsupportedStructure("a carboxylic acid outranks the sulfur or phosphorus acid, which is then a prefix, not an acyl group")
    from ._acid_derivatives import acyl_name

    from ._polyfunctional import RING_CENTER

    token = SYSTEMATIC_ACID_PROBE.set(True)
    labels_token = ISOTOPE_LABELS.set(None)
    centre_token = RING_CENTER.set(False)
    acid = Chem.Mol(sub)
    acid.GetAtomWithIdx(sorted(subtree).index(root)).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    adjusted = _match_whole_labels(mol, acid, sorted(subtree))
    others = [match[0] for match in acid.GetSubstructMatches(_FREE_CARBOXY) if match[0] != sorted(subtree).index(root)]
    try:
        if own_acid is not None:
            probe = _name_with_principal(acid, own_acid)
        elif others:
            probe = _name_with_free_carboxy(acid, others) or smiles_to_iupac(Chem.MolToSmiles(acid))
            if probe != "oxalic acid" and _MULTIPLIED_ACID.search(probe):
                raise UnsupportedStructure("a further carboxylic acid group would become a second acyl group of the prefix")
        else:
            probe = _name_without_cation_context(acid) if not zwitterion else smiles_to_iupac(Chem.MolToSmiles(acid))
    finally:
        RING_CENTER.reset(centre_token)
        ISOTOPE_LABELS.reset(labels_token)
        SYSTEMATIC_ACID_PROBE.reset(token)
    if zwitterion and probe.endswith("ate"):
        probe = probe[: -len("ate")] + "ic acid"
    name = acyl_name(probe)
    _cite_stereo_of_acyl(mol, acid, sorted(subtree), adjusted, name)
    if labelled:
        context["consumed"].update(labelled)
    # HOOC-CO- keeps one acid group, so the retained 'oxalyl' (-CO-CO-) becomes 'oxalo' (P-65.1.2.2.3)
    return "oxalo" if name == "oxalyl" else name


def _name_without_cation_context(acid):
    """The name of an acid that holds no cation, outside the cationic-parent context of the molecule it was cut from."""
    from ._polyfunctional import AMINIUM, FORCED_PRINCIPAL
    from .core import smiles_to_iupac

    aminium_token, forced_token = AMINIUM.set(False), FORCED_PRINCIPAL.set(None)
    try:
        return smiles_to_iupac(Chem.MolToSmiles(acid))
    finally:
        FORCED_PRINCIPAL.reset(forced_token)
        AMINIUM.reset(aminium_token)


def _name_with_free_carboxy(acid, others):
    """The acid name in which the carboxylic acid groups at `others` stay 'carboxy' prefixes, so that the acyl group
    made from it has one acyl position only (P-65.1.7); None when the structure has no such name."""
    from ._polyfunctional import AMINIUM, FORCED_PRINCIPAL, JUNIOR_GROUP, name_polyfunctional

    for index in others:
        acid.GetAtomWithIdx(index).SetBoolProp(JUNIOR_GROUP, True)
    aminium_token, forced_token = AMINIUM.set(False), FORCED_PRINCIPAL.set(None)
    try:
        return name_polyfunctional(acid)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_PRINCIPAL.reset(forced_token)
        AMINIUM.reset(aminium_token)


def _imido_acid_class(sub, root):
    """Class key of the sulfonimidic-type acid on `root` of the cut-out group, else None: the =NR position is lost when
    such a group is cited as a prefix by the sulfonyl-type fallback (P-65.3.2.2)."""
    from ._acid_groups import acid_group_at

    group = acid_group_at(sub, root)
    if group is None or not any(x in ("NH", "NNH2") for x in group.spec.oxo):
        return None
    return group.spec.key


def _name_with_principal(sub, principal):
    from ._polyfunctional import FORCED_PRINCIPAL, name_polyfunctional

    token = FORCED_PRINCIPAL.set(principal)
    try:
        name = name_polyfunctional(sub)
    finally:
        FORCED_PRINCIPAL.reset(token)
    if name is None:
        raise UnsupportedStructure("the acid group is not a principal group of its fragment")
    return name


def _match_whole_labels(mol, acid, kept):
    """The cut-out acid gets the CIP labels its centres and double bonds have in the whole molecule, so that its name
    cites those: a label can change or be lost, because the arms of a centre become identical, once the rest of the
    molecule is cut away (P-92.1.4, P-93.5.1). The elements that had no label of their own are returned."""
    from ._common import BOND_STEREO_OVERRIDE, STEREO_OVERRIDE, assign_cip_labels, match_cip_labels
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    if context is None:
        return set()
    in_acyl = {index: new for new, index in enumerate(kept)}
    match_cip_labels(mol, acid, {new: index for index, new in in_acyl.items()})
    assign_cip_labels(acid)
    adjusted = set()
    for atom, code in context["atoms"].items():
        if atom in in_acyl and code in ("R", "S") and not acid.GetAtomWithIdx(in_acyl[atom]).HasProp("_CIPCode"):
            acid.GetAtomWithIdx(in_acyl[atom]).SetProp(STEREO_OVERRIDE, code)
            adjusted.add(("atom", atom))
    for (a, b), code in context["bonds"].items():
        bond = acid.GetBondBetweenAtoms(in_acyl[a], in_acyl[b]) if a in in_acyl and b in in_acyl else None
        if bond is not None and code in ("E", "Z") and not bond.HasProp("_CIPCode"):
            bond.SetProp(BOND_STEREO_OVERRIDE, code)
            adjusted.add(("bond", (a, b)))
    return adjusted


def _cite_stereo_of_acyl(mol, acid, kept, adjusted, name):
    """The acid name of an acyl group carries the descriptors of the whole molecule, which are cited with the acyl
    prefix."""
    from ._substituents import BRANCH_STEREO

    context = BRANCH_STEREO.get()
    if context is None:
        return
    in_acyl = {index: new for new, index in enumerate(kept)}
    atoms = [a for a in context["atoms"] if a in in_acyl and mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
    bonds = [(a, b) for a, b in context["bonds"] if a in in_acyl and b in in_acyl]
    if len(adjusted) > len(re.findall(r"(?<=[(,])\d+[\u2032']*[RSEZ](?=[,)])", name)):
        raise UnsupportedStructure("a stereo element of the acyl group is not cited by the name of its acid")
    for atom in atoms:
        if ("atom", atom) in adjusted or acid.GetAtomWithIdx(in_acyl[atom]).HasProp("_CIPCode"):
            context["used"].add(("atom", atom))
    for a, b in bonds:
        twin = acid.GetBondBetweenAtoms(in_acyl[a], in_acyl[b])
        if ("bond", (a, b)) in adjusted or (twin is not None and twin.HasProp("_CIPCode")):
            context["used"].add(("bond", (a, b)))


def _multiplied_amino(children, tail):
    parts = sorted(children, key=lambda c: alpha_sort_key(c[0]))
    names = [(_enclose(n) if c else n) for n, c in parts]
    if len(set(parts)) == 1 and len(parts) > 1:
        n, c = parts[0]
        return multiplying_prefix(len(parts), compound=c) + (_enclose(n) if c else n) + tail
    first, rest = names[0], names[1:]
    return first + "".join(_enclose(r) if not r.startswith(("(", "[")) else r for r in rest) + tail


def _phosphoryl_parts(mol, node, parent, kids, named, bond_order):
    """Substituent names of a P(=O)(X)(Y) group bonded through oxygen, from the already-named oxygen children."""
    atom = mol.GetAtomWithIdx(node)
    if atom.GetFormalCharge() != 0 or atom.GetDegree() != 4:
        raise UnsupportedStructure("a charged or non-tetracoordinate phosphorus substituent is not supported yet")
    oxo = [
        k
        for k in kids
        if bond_order(mol, node, k) == 2.0 and mol.GetAtomWithIdx(k).GetDegree() == 1 and mol.GetAtomWithIdx(k).GetAtomicNum() == 8
    ]
    rest = [k for k in kids if k not in oxo]
    if len(oxo) != 1 or len(rest) != 2 or any(bond_order(mol, node, k) != 1.0 for k in rest):
        raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
    parts = []
    for k in rest:
        if k not in named or mol.GetAtomWithIdx(k).GetAtomicNum() != 8:
            raise UnsupportedStructure("this phosphorus-bearing substituent is not supported yet")
        if "phospho" in named[k][0]:
            mark(None, POLYACID_SUBSTITUENT_REASON)
        parts.append(named[k])
    return parts


def functional_names(mol, graph, seeds, blocked, halogens, aromatic_atoms=frozenset(), chain_seeds=False):
    """({atom: (name, is_compound)}, {atom: display string}, covered atoms) for every functional group rooted at an
    atom reachable outward from `seeds` ((parent, root) pairs) without entering `blocked`."""
    parent_of = {}
    order = []
    for parent, root in seeds:
        if root in parent_of or root in blocked:
            continue
        parent_of[root] = parent
        queue = [root]
        while queue:
            node = queue.pop(0)
            order.append(node)
            for n in graph[node]:
                if n not in parent_of and n not in blocked and n != parent_of[node]:
                    parent_of[n] = node
                    queue.append(n)

    named = {}
    shown = dict(halogens)
    children_of = {}
    for node in order:
        children_of.setdefault(parent_of[node], []).append(node)

    def subtree(node):
        found = {node}
        stack = [node]
        while stack:
            for child in children_of.get(stack.pop(), []):
                found.add(child)
                stack.append(child)
        return found

    nitro = nitro_atoms(mol)
    nitroso = nitroso_atoms(mol)
    root_set = set() if chain_seeds else {root for _, root in seeds}

    from ._diester_ring_diyl import _system_of

    ring_info = mol.GetRingInfo()
    ring_entries = set()
    for node in order:
        if ring_info.NumAtomRings(node) and node not in nitro:
            if parent_of[node] not in _system_of(mol, node)[1]:
                ring_entries.add(node)
    skip = set()
    for node in ring_entries:
        skip |= subtree(node) - {node}
    delegated = {}
    for _, root in seeds:
        if (
            root in parent_of
            and root not in skip
            and (
                mol.GetAtomWithIdx(root).GetAtomicNum() not in _NATIVE_ROOTS
                or _anionic_branch(mol, root, subtree(root))
                or _acyl_hydrazine_sulfur(mol, root)
                or _acylamino_nitrogen(mol, root, parent_of[root])
                or _nitrogen_on_ring_nitrogen(mol, root, parent_of[root])
            )
        ):
            try:
                delegated[root] = name_branch(graph, root, parent_of[root], shown, aromatic_atoms, mol=mol)
            except UnsupportedStructure:
                continue
            skip |= subtree(root) - {root}

    def child_name(child, via):
        if child in named:
            return named[child]
        if mol.GetAtomWithIdx(child).GetAtomicNum() == 6:
            return name_branch(graph, child, via, shown, aromatic_atoms, mol=mol)
        raise UnsupportedStructure("this substituent group is not supported yet")

    def record(node, name, compound):
        named[node] = (name, compound)
        shown[node] = _enclose(name) if compound else name

    pseudohalides = {}
    for node in order:
        if node in skip or node not in parent_of or mol.GetAtomWithIdx(node).GetAtomicNum() != 7:
            continue
        kids = [n for n in graph[node] if n != parent_of[node]]
        order_to_parent = _bond_order(mol, node, parent_of[node])
        prefix = nitrogen_pseudohalide_prefix(mol, node, kids, order_to_parent)
        if prefix is None and mol.GetAtomWithIdx(parent_of[node]).GetAtomicNum() == 6 and diazo_prefix(mol, node, kids, order_to_parent):
            prefix = "diazo"
        if prefix is not None:
            pseudohalides[node] = prefix
            skip |= subtree(node) - {node}

    azoxy_groups = {}
    for node in order:
        if node in skip or node not in parent_of or _bond_order(mol, node, parent_of[node]) != 1.0:
            continue
        atom = mol.GetAtomWithIdx(node)
        if atom.GetAtomicNum() == 7 and (
            azoxy_atoms(mol, atom) is not None
            or any(azoxy_atoms(mol, mol.GetAtomWithIdx(n)) for n in graph[node] if n != parent_of[node])
        ):
            azoxy_groups[node] = name_branch(graph, node, parent_of[node], shown, aromatic_atoms, mol=mol)
            skip |= subtree(node) - {node}

    phosphoryl_nodes = set()
    for node in reversed(order):
        if node in skip:
            continue
        if node in azoxy_groups:
            record(node, *azoxy_groups[node])
            continue
        if node in pseudohalides:
            record(node, pseudohalides[node], False)
            continue
        if node in ring_entries:
            name, compound = name_branch(graph, node, parent_of[node], shown, aromatic_atoms, mol=mol)
            record(node, name, compound)
            continue
        forced = FORCED_BRANCH_NAMES.get()
        if forced and mol.GetNumAtoms() == forced[0] and node in forced[1]:
            record(node, *forced[1][node])
            continue
        if node in delegated:
            record(node, *delegated[node])
            continue
        atom = mol.GetAtomWithIdx(node)
        z = atom.GetAtomicNum()
        parent = parent_of[node]
        kids = [n for n in graph[node] if n != parent]
        if z in HALOGEN_PREFIXES:
            continue
        if z == 15 and mol.GetAtomWithIdx(parent).GetAtomicNum() == 8:
            from ._oxoacid_acyl import infix_acyl_name

            replaced = infix_acyl_name(mol, graph, node, parent, shown, aromatic_atoms)
            if replaced is not None:
                phosphoryl_nodes.add(node)
                record(node, *replaced)
                continue
            if any(mol.GetAtomWithIdx(k).GetAtomicNum() == 6 for k in kids):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol, unsaturated=True))
                phosphoryl_nodes.add(node)
                continue
            group = phosphoryl_name(_phosphoryl_parts(mol, node, parent, kids, named, _bond_order))
            phosphoryl_nodes.add(node)
            record(node, group, group != "phosphono")
            continue
        if node in nitro:
            if is_nitro_nitrogen(mol, node):
                record(node, "nitro", False)
            continue
        if node in nitroso:
            if z == 7:
                record(node, "nitroso", False)
            continue
        if (
            atom.HasProp("_anion") or (DIPOLAR_GROUPS.get() and atom.GetFormalCharge() == -1 and z in (8, 16))
        ) and z in ANIONIC_PREFIXES and not kids and _bond_order(mol, node, parent) == 1.0:
            record(node, ANIONIC_PREFIXES[z], False)
            continue
        if z == 8 and node not in named:
            if _bond_order(mol, node, parent) == 2.0 and not kids:
                record(node, "oxo", False)
            elif not kids:
                record(node, "hydroxy", False)
            elif len(kids) == 1 and kids[0] in phosphoryl_nodes:
                rname, rcompound = named[kids[0]]
                record(node, "phosphonooxy" if rname == "phosphono" else _enclose(rname) + "oxy", True)
            elif len(kids) == 1 and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() in (6, 9, 17, 35, 53):
                rname, rcompound = child_name(kids[0], node)
                record(node, *alkoxy_prefix(rname, rcompound))
            elif (
                len(kids) == 1
                and not atom.IsInRing()
                and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() in (7, 8, 16, 34, 52, *MONONUCLEAR_HYDRIDES)
                and _bond_order(mol, node, kids[0]) == 1.0
                and _bond_order(mol, node, parent) == 1.0
            ):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol, unsaturated=True))
            else:
                raise UnsupportedStructure("this oxygen-bearing substituent is not supported yet")
        elif z in CHALCOGEN_PREFIXES and atom.GetTotalValence() > 2 and not atom.IsInRing():
            record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol, unsaturated=True))
        elif z in CHALCOGEN_PREFIXES:
            word = CHALCOGEN_PREFIXES[z]
            require_plain_chalcogen_kids(mol, z, kids)
            if not kids and _bond_order(mol, node, parent) == 2.0:
                record(node, word[:-2] + "ylidene", False)
            elif not kids:
                record(node, word, False)
            elif len(kids) == 1 and atom.GetDegree() == 2 and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() == 6:
                rname, rcompound = child_name(kids[0], node)
                name = (_enclose(rname) if rcompound else rname) + word
                record(node, name, True)
            else:
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
        elif z == 7:
            bond = _bond_order(mol, node, parent)
            cationic = atom.GetFormalCharge() == 1 and atom.GetTotalValence() == 4
            if (cationic or (atom.HasProp("_cationic_amine") and not atom.GetFormalCharge())) and not atom.IsInRing():
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
                continue
            if bond == 3.0 and not kids:
                continue
            if bond == 2.0 and len(kids) == 1 and named.get(kids[0]) == ("amino", False) and _bond_order(mol, node, kids[0]) == 1.0:
                record(node, "hydrazinylidene", False)
                continue
            if (
                bond == 2.0
                and len(kids) == 1
                and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() == 7
                and _bond_order(mol, node, kids[0]) == 1.0
                and not mol.GetAtomWithIdx(kids[0]).IsInRing()
            ):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
                continue
            if (
                bond == 2.0
                and len(kids) <= 1
                and not atom.GetFormalCharge()
                and all(_bond_order(mol, node, k) == 1.0 for k in kids)
                and (_has_senior_principal_group(mol) or IMINE_PARENT.get())
            ):
                if not kids:
                    record(node, "imino", False)
                else:
                    rname, rcompound = child_name(kids[0], node)
                    record(node, (_enclose(rname) if rcompound else rname) + "imino", True)
                continue
            if (
                bond == 1.0
                and mol.GetAtomWithIdx(parent).GetAtomicNum() == 7
                and len(kids) == 1
                and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() == 6
                and _bond_order(mol, node, kids[0]) == 2.0
            ):
                continue
            if (
                bond == 1.0
                and not atom.IsInRing()
                and len(kids) == 1
                and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() in (6, 7)
                and _bond_order(mol, node, kids[0]) == 2.0
            ):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
                continue
            if (
                bond == 2.0
                and mol.GetAtomWithIdx(parent).GetAtomicNum() == 7
                and not atom.IsInRing()
                and len(kids) == 1
                and mol.GetAtomWithIdx(kids[0]).GetAtomicNum() == 6
                and _bond_order(mol, node, kids[0]) == 1.0
            ):
                continue
            if (
                bond == 2.0
                and z == 7
                and not atom.IsInRing()
                and mol.GetAtomWithIdx(parent).GetAtomicNum() == 6
                and all(_bond_order(mol, node, k) == 1.0 for k in kids)
            ):
                continue
            if bond != 1.0 or any(_bond_order(mol, node, k) != 1.0 for k in kids):
                raise UnsupportedStructure("an imine/azo/nitroso-type substituent is not supported yet")
            if not atom.IsInRing() and any(
                mol.GetAtomWithIdx(k).GetAtomicNum() == 7 and not mol.GetAtomWithIdx(k).IsInRing() for k in kids
            ):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
                continue
            if not atom.IsInRing() and any(mol.GetAtomWithIdx(k).GetAtomicNum() in MONONUCLEAR_HYDRIDES for k in kids):
                record(node, *name_branch(graph, node, parent, shown, aromatic_atoms, mol=mol))
                continue
            if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 and k not in named for k in kids):
                raise UnsupportedStructure("a hetero-substituted nitrogen substituent is not supported yet")
            anilino = anilino_prefix(graph, mol, node, kids, shown, aromatic_atoms) if kids else None
            if anilino is not None:
                record(node, anilino[0], anilino[0] != "anilino")
            elif not kids:
                record(node, "amino", False)
            else:
                children = [child_name(k, node) for k in kids]
                name = _multiplied_amino(children, "amino")
                record(node, name, True)
        elif z == 6 and node not in named:
            ketene = ketene_prefix(mol, graph, node, parent)
            if ketene is not None:
                record(node, *ketene)
                continue
            carbonyl = [k for k in kids if mol.GetAtomWithIdx(k).GetAtomicNum() == 8 and not mol.GetAtomWithIdx(k).GetFormalCharge() and _bond_order(mol, node, k) == 2.0]
            triple_n = [k for k in kids if mol.GetAtomWithIdx(k).GetAtomicNum() == 7 and _bond_order(mol, node, k) == 3.0]
            if triple_n and len(kids) == 1 and mol.GetAtomWithIdx(triple_n[0]).GetFormalCharge() != 1:
                record(node, "cyano", False)
            elif carbonyl and (node in root_set or mol.GetAtomWithIdx(parent).GetAtomicNum() != 6):
                others = [k for k in kids if k not in carbonyl]
                if len(carbonyl) != 1 or len(others) > 1:
                    raise UnsupportedStructure("this carbonyl substituent is not supported yet")
                if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 for k in others):
                    record(node, *_functional_carbon(graph, node, parent, shown, aromatic_atoms, mol))
                    continue
                name = "formyl" if not others else _acyl_prefix(mol, subtree(node), node, parent)
                record(node, name, is_compound_acyl(name))
            elif carbonyl and not _terminal_amide(graph, mol, node, parent):
                others = [k for k in kids if k not in carbonyl]
                if any(mol.GetAtomWithIdx(k).GetAtomicNum() != 6 for k in others):
                    record(node, *_functional_carbon(graph, node, parent, shown, aromatic_atoms, mol))
    covered = set()
    for node in named:
        covered |= subtree(node)
    return named, shown, covered
