"""Secondary and tertiary amides (R-CO)2NH, (R-SO2)2NH, (R-CO)3N (P-66.1.2): N-acyl derivatives of the senior primary
amide, e.g. N-acetylbenzamide, N,N-di(cyclohexanecarbonyl)cyclohexanecarboxamide. The other acyl groups are stood in
for by methyl groups whose names are replaced by the acyl prefixes, so the amide is named with its usual N-prefixes."""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._diester_ring_diyl import _system_of
from ._dipolar import _smiles_with_order
from ._isotope_labels import restore_labels, split_isotopes
from ._ring_system_seniority import ring_seniority_key
from ._substituents import FORCED_BRANCH_NAMES, ISOTOPE_LABELS, name_branch


def _acyl_roots(mol, nitrogen):
    """Neighbours of `nitrogen` that are the carbon of a carboxylic acyl group or the sulfur of a sulfonyl group."""
    roots = []
    for n in nitrogen.GetNeighbors():
        if mol.GetBondBetweenAtoms(nitrogen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        others = [m for m in n.GetNeighbors() if m.GetIdx() != nitrogen.GetIdx()]
        doubly = [m for m in others if m.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(n.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() == 2.0]
        rest = [m for m in others if m not in doubly]
        if n.GetAtomicNum() == 6 and len(doubly) == 1 and not n.IsInRing():
            if (not rest and n.GetTotalNumHs() == 1) or (len(rest) == 1 and rest[0].GetAtomicNum() == 6):
                roots.append(n.GetIdx())
        elif n.GetAtomicNum() == 16 and len(doubly) == 2 and len(rest) == 1 and rest[0].GetAtomicNum() == 6 and not n.IsInRing():
            roots.append(n.GetIdx())
    return roots


def _side(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _rank(mol, graph, root, nitrogen):
    """Smaller is senior: carboxamides before sulfonamides, ring parents before chains, heterocycles before
    carbocycles, then the larger ring or the longer chain."""
    atom = mol.GetAtomWithIdx(root)
    sulfur = atom.GetAtomicNum() == 16
    carbons = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    if not carbons:
        return (int(sulfur), 2, (), 0)
    start = carbons[0]
    side = _side(graph, start, root)
    if mol.GetAtomWithIdx(start).IsInRing():
        return (int(sulfur), 0, ring_seniority_key(mol, _system_of(mol, start)[1]), 0)
    length = len([a for a in side if mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(a).IsInRing()])
    return (int(sulfur), 1, (), -length)


def donor_names(mol, graph, donors, halogens, aromatic, labels=None):
    """{root: name} of the acyl donor groups, each cited with the descriptors of its own stereo elements and the nuclides
    of `labels` on its atoms; raises UnsupportedStructure when one of them is left uncited."""
    from ._common import assign_cip_labels

    from ._substituents import BRANCH_STEREO

    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    names = {}
    inside = set().union(*(_side(graph, root, nitrogen) for root, nitrogen in donors))
    isotope_context = {"labels": {a: e for a, e in (labels or {}).items() if a in inside}, "consumed": set(), "mol": mol}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if isotope_context["labels"] else None)
    try:
        for root, nitrogen in donors:
            side = _side(graph, root, nitrogen)
            context = {
                "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in side and a.HasProp("_CIPCode")},
                "bonds": {
                    (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
                    for b in probe.GetBonds()
                    if b.GetBeginAtomIdx() in side and b.GetEndAtomIdx() in side and b.HasProp("_CIPCode")
                },
                "used": set(),
            }
            token = BRANCH_STEREO.set(context)
            try:
                names[root] = name_branch(graph, root, nitrogen, halogens, aromatic, mol=mol)
            finally:
                BRANCH_STEREO.reset(token)
            if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
                ("bond", b) not in context["used"] for b in context["bonds"]
            ):
                raise UnsupportedStructure("a stereo element of an acyl group is not cited by its name")
    finally:
        ISOTOPE_LABELS.reset(isotope_token)
    if set(isotope_context["labels"]) - isotope_context["consumed"]:
        raise UnsupportedStructure("a nuclide of an acyl group is not cited by its name")
    return names


def forced_for_unlabelled(stand_in, position, names):
    """(atom count, {atom index: name}) of the stand-in without its isotopic hydrogen atoms, in the SMILES order."""
    heavy = sorted(position[a.GetIdx()] for a in stand_in.GetAtoms() if a.GetAtomicNum() != 1)
    return len(heavy), {heavy.index(position[atom]): name for atom, name in names.items()}


def modified_atoms_beyond(graph, labels, root, nitrogen):
    """Negative count of the nuclides in the acyl group at `root`: of equally senior parents the more modified is senior
    (P-44.4.1.11.1)."""
    return -sum(
        (1 if entry["skeleton"] else 0) + sum(entry["H"].values())
        for atom, entry in labels.items()
        if atom in _side(graph, root, nitrogen)
    )


def unlabelled_with_labels(mol):
    """(mol without its nuclide labels, labels), or (None, {}) when a nuclide sits on a heteroatom or the structure has
    stereo elements."""
    try:
        split = split_isotopes(mol)
    except UnsupportedStructure:
        return None, {}
    if split is None:
        return mol, {}
    clean, labels, _ = split
    stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms()) or any(
        b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) for b in clean.GetBonds()
    )
    if stereo or any(clean.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in labels):
        return None, {}
    return clean, labels


def methylate_donor(editable, donor, nitrogen, removed):
    """Turns the acyl atom `donor` into a methyl-type carbon bonded to `nitrogen`, dropping its bonds to the atoms `removed`;
    it gets its hydrogens back even when a nuclide label made it a bracket atom."""
    atom = editable.GetAtomWithIdx(donor)
    atom.SetAtomicNum(6)
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)
    for bond in list(atom.GetBonds()):
        if bond.GetOtherAtomIdx(donor) in removed:
            editable.RemoveBond(donor, bond.GetOtherAtomIdx(donor))
    editable.AddBond(nitrogen, donor, Chem.BondType.SINGLE)


def diacylamine_name(mol):
    """The name of a molecule with an acyclic nitrogen bearing two or three acyl groups, else None."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    mol, labels = unlabelled_with_labels(mol)
    if mol is None:
        return None
    nitrogens = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 7 and not a.IsInRing() and not a.GetIsAromatic() and len(_acyl_roots(mol, a)) >= 2
    ]
    if len(nitrogens) != 1:
        return None
    nitrogen = nitrogens[0]
    graph = adjacency(mol)
    roots = _acyl_roots(mol, nitrogen)
    parent = min(
        roots,
        key=lambda r: (_rank(mol, graph, r, nitrogen.GetIdx()), modified_atoms_beyond(graph, labels, r, nitrogen.GetIdx()), r),
    )
    donors = [r for r in roots if r != parent]
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    try:
        names = donor_names(mol, graph, [(r, nitrogen.GetIdx()) for r in donors], halogens, aromatic, labels)
    except UnsupportedStructure:
        return None
    editable = Chem.RWMol(mol)
    removed = set()
    for r in donors:
        side = _side(graph, r, nitrogen.GetIdx())
        removed |= side - {r}
        editable.RemoveBond(nitrogen.GetIdx(), r)
    keep = sorted(removed, reverse=True)
    for donor in donors:
        methylate_donor(editable, donor, nitrogen.GetIdx(), removed)
    restore_labels(editable, mol, labels, skip=removed | set(donors))
    for idx in keep:
        editable.RemoveAtom(idx)
    shift = 0
    old_to_new = {}
    for i in range(mol.GetNumAtoms()):
        if i in removed:
            shift += 1
        else:
            old_to_new[i] = i - shift
    stand_in = editable.GetMol()
    try:
        Chem.SanitizeMol(stand_in)
    except Chem.rdchem.MolSanitizeException:
        return None
    smiles, position = _smiles_with_order(stand_in)
    forced = forced_for_unlabelled(stand_in, position, {old_to_new[r]: names[r] for r in donors})
    from .core import smiles_to_iupac

    token = FORCED_BRANCH_NAMES.set(forced)
    try:
        name = smiles_to_iupac(smiles)
    except UnsupportedStructure:
        return None
    finally:
        FORCED_BRANCH_NAMES.reset(token)
    if not name.endswith("amide") or any(f"{donor_name}" not in name for donor_name, _ in names.values()):
        return None
    return name
