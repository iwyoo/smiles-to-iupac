"""Ester names built from their two parts (P-65.6.3.2.1): the alkyl/aryl group
cited as a substituent group, then the anion of the acid part. The molecule is
split at the O-alkyl bond; the acid fragment is named as a free acid and its
'-ic acid' ending becomes '-ate', so every acid shape the engine can name also
works as an ester acyl part.
"""

import re

from rdkit import Chem

from ._carbonic_amide import _carbonic_carbon
from ._common import assign_cip_labels
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, halogen_substituents
from ._isotope_labels import restore_labels
from ._numerals import multiplied_ester_word, multiplying_prefix
from ._polycarbonic import chain_numbering
from ._substituents import BRANCH_STEREO, name_branch

_ESTER = Chem.MolFromSmarts("[CX3;!R](=O)[OX2;!R][#6]")
_PEROXY_ESTER = Chem.MolFromSmarts("[CX3;!R](=O)[OX2;!R][OX2;!R][#6]")
_CARBONIC_LINK = Chem.MolFromSmarts("[C;!R]O[C;!R]")
_FREE_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
_HYDROGEN_WORDS = {1: "hydrogen", 2: "dihydrogen", 3: "trihydrogen"}


def _carbonyl_with_heteroatom(mol, idx):
    """A carbonyl carbon bonded through a single bond to oxygen, sulfur, selenium or tellurium: an acid, ester or anhydride,
    which ranks above an ester; an amide or an acid halide (classes 11 and 10, P-41) is a prefix of the alkyl part."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    has_double_o = any(
        n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )
    return has_double_o and any(
        n.GetAtomicNum() in (8, 16, 34, 52) and mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in atom.GetNeighbors()
    )


def _insert_before_ending(anion, text):
    for ending in ("oate", "sulfonate", "sulfinate", "ate"):
        if anion.endswith(ending):
            return anion[: -len(ending)] + text + ending
    raise UnsupportedStructure("the ester ending is not delimited")


def _anion_name(acid_name):
    if not acid_name.endswith(" acid") and " " not in acid_name:
        from ._amino_acid_derivative import _anion_stem

        stem = _anion_stem(acid_name)
        if stem is not None:
            return stem
    if not acid_name.endswith(" acid") or " " in acid_name[: -len(" acid")].strip():
        raise UnsupportedStructure("the acid part of this ester is not named as a plain acid")
    stem = acid_name[: -len(" acid")]
    if not stem.endswith("ic"):
        raise UnsupportedStructure("the acid part of this ester has an unexpected acid name ending")
    return stem[:-2] + "ate"


def _branch_atoms(graph, root, blocked):
    seen, stack = {root}, [root]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def name_ester_by_parts(mol) -> str:
    from ._isotope_labels import split_isotopes

    split = split_isotopes(mol)
    if split is None:
        return _name_ester_parts(mol, {})
    clean, labels, index_map = split
    chiral = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms())
    return _name_ester_parts(clean, labels, (mol, index_map) if chiral else None)


def _acid_with_labels(original, index_map, removed, matches, group_labelled):
    """The free acid cut from the labelled molecule itself, which keeps the configuration of a centre that the nuclide
    makes stereogenic (P-82.4); `removed` are the atoms of the alkyl parts in the numbering of the unlabelled molecule and
    `group_labelled` the ester-group atoms whose nuclides are cited in the ester word instead (P-82.6.4)."""
    inverse = {new: old for old, new in index_map.items()}
    dropped = set()
    for atom in removed:
        dropped.add(inverse[atom])
        dropped |= {n.GetIdx() for n in original.GetAtomWithIdx(inverse[atom]).GetNeighbors() if n.GetAtomicNum() == 1}
    for atom in group_labelled:
        dropped |= {
            n.GetIdx() for n in original.GetAtomWithIdx(inverse[atom]).GetNeighbors() if n.GetAtomicNum() == 1 and n.GetIsotope()
        }
    editable = Chem.RWMol(original)
    for atom in group_labelled:
        editable.GetAtomWithIdx(inverse[atom]).SetIsotope(0)
    for _, _, ester_oxygen, _ in matches:
        oxygen = editable.GetAtomWithIdx(inverse[ester_oxygen])
        oxygen.SetNumExplicitHs(1)
        oxygen.SetNoImplicit(True)
    for idx in sorted(dropped, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    Chem.SanitizeMol(acid)
    return acid


def _ester_locants(mol, acid, removed, matches, named_by_arm):
    """{alkyl name: sorted locants of its ester groups} in the substitutive name of the acid, or {} when the ester carbons
    are equivalent in the acid and the bare name is unambiguous (P-65.6.3.2.2, P-14.3.4.4)."""
    from ._polyfunctional import _select

    kept = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    new_index = {old: new for new, old in enumerate(kept)}
    carbons = [new_index[acyl] for acyl, _, _, _ in matches]
    plain = Chem.Mol(acid)
    ranks = Chem.CanonicalRankAtoms(plain, breakTies=False)
    numbering = chain_numbering(acid)
    if numbering:
        ends = [c for c, locant in numbering.items() if locant in (1, 2 * len(numbering) - 1)]
        bare = Chem.Mol(acid)
        for atom in bare.GetAtoms():
            atom.SetIsotope(0)
        bare_ranks = Chem.CanonicalRankAtoms(bare, breakTies=False)
        if len({bare_ranks[c] for c in ends}) > 1 or len({ranks[c] for c in ends}) == 1:
            return {}
        if set(ends) <= set(carbons) and len(set(named_by_arm.values())) == 1:
            return {}
        found = {}
        for arm, carbon in enumerate(carbons):
            found.setdefault(named_by_arm[arm], []).append(str(numbering[carbon]))
        return {name: sorted(locs, key=lambda t: (len(t), t)) for name, locs in found.items()}
    if len({ranks[c] for c in carbons}) == 1:
        return {}
    for atom in acid.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    try:
        selected = _select(acid)
    except UnsupportedStructure:
        selected = None
    if selected is None:
        raise UnsupportedStructure("the locants of the ester groups are not available")
    parts = selected[2]
    position = {acid.GetAtomWithIdx(i).GetAtomMapNum() - 1: loc for i, loc in parts[4].items()}
    found = {}
    for arm, carbon in enumerate(carbons):
        if carbon in position:
            locant = position[carbon]
        else:
            attach = [n.GetIdx() for n in acid.GetAtomWithIdx(carbon).GetNeighbors() if n.GetIdx() in position]
            if len(attach) != 1:
                raise UnsupportedStructure("the locants of the ester groups are not available")
            locant = position[attach[0]]
        found.setdefault(named_by_arm[arm], []).append(str(locant))
    return {name: sorted(locs, key=lambda t: (len(t), t)) for name, locs in found.items()}


def _arms_in_parent(mol, matches):
    """The ester groups whose acid carbon is a suffix of the parent acid; the others are cited as prefixes of a
    monoester (P-65.6.3.3.2.2.2). None when the parent cannot be determined."""
    from ._polyfunctional import _select

    graph = adjacency(mol)
    removed = set()
    for _, _, ester_oxygen, alkyl_carbon in matches:
        removed |= _branch_atoms(graph, alkyl_carbon, ester_oxygen)
    if any(acyl in removed for acyl, _, _, _ in matches):
        return None
    editable = Chem.RWMol(mol)
    for _, _, ester_oxygen, _ in matches:
        oxygen = editable.GetAtomWithIdx(ester_oxygen)
        oxygen.SetNumExplicitHs(1)
        oxygen.SetNoImplicit(True)
    kept = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    acid = editable.GetMol()
    try:
        Chem.SanitizeMol(acid)
        for atom in acid.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        selected = _select(acid)
    except (UnsupportedStructure, ValueError):
        return None
    if selected is None:
        return None
    parts = selected[2]
    if not (parts[1] or "").endswith(("oic", "carboxylic")):
        return None
    new_index = {old: new for new, old in enumerate(kept)}
    position = {acid.GetAtomWithIdx(i).GetAtomMapNum() - 1 for i in parts[4]}
    carbons = [new_index[acyl] for acyl, _, _, _ in matches]
    if any(c in position for c in carbons):
        return [arm for arm, c in enumerate(carbons) if c in position]
    attached = [[n.GetIdx() for n in acid.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 6] for c in carbons]
    return [arm for arm, nbrs in enumerate(attached) if any(n in position for n in nbrs)]


def _require_ester_with_its_free_acid(mol, matches):
    """A free carboxylic acid outranks the ester (P-41), so the monoester 'ethyl hydrogen ...ate' of one acid parent names
    an ester only when its carbonyl carbon is joined to a free acid carbon through carbon atoms alone; an ester on another
    part of the molecule is an alkoxy-oxo prefix of the acid (P-65.6.3.2.3)."""
    free = {m[0] for m in mol.GetSubstructMatches(_FREE_ACID)}
    if not free:
        return
    carbons = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
    for acyl_carbon, *_ in matches:
        reached, stack = {acyl_carbon}, [acyl_carbon]
        while stack:
            current = stack.pop()
            for n in mol.GetAtomWithIdx(current).GetNeighbors():
                onward = [n.GetIdx()] if n.GetIdx() in carbons else []
                if n.GetAtomicNum() == 8 and n.GetDegree() == 2:
                    onward = [m.GetIdx() for m in n.GetNeighbors() if m.GetIdx() != current and _carbonic_link(mol, current, m.GetIdx())]
                for idx in onward:
                    if idx not in reached:
                        reached.add(idx)
                        stack.append(idx)
        if not reached & free:
            raise UnsupportedStructure("an ester that is not on the carbon skeleton of the free acid is an alkoxy-oxo prefix")


def _carbonic_link(mol, first, second):
    """Whether two carbons are the units of a polycarbonic acid joined through one oxygen."""
    return _carbonic_carbon(mol.GetAtomWithIdx(first)) and _carbonic_carbon(mol.GetAtomWithIdx(second))


def _polycarbonic(mol):
    return any(_carbonic_link(mol, m[0], m[2]) for m in mol.GetSubstructMatches(_CARBONIC_LINK))


def _ester_matches(mol):
    """(acyl carbon, carbonyl oxygen, ester oxygen, alkyl carbon) of each ester group; the ester oxygen of a peroxy ester
    is the one bonded to the alkyl carbon, so the peroxoic acid is what the split leaves."""
    peroxy = [(m[0], m[1], m[3], m[4]) for m in mol.GetSubstructMatches(_PEROXY_ESTER)]
    esters = [m for m in mol.GetSubstructMatches(_ESTER) if not _carbonic_link(mol, m[0], m[3])]
    return (*esters, *peroxy)


def _name_ester_parts(mol, labels, original=None) -> str:
    from ._substituents import ISOTOPE_LABELS
    from .core import smiles_to_iupac

    matches = _ester_matches(mol)
    if not matches:
        raise UnsupportedStructure("an acyclic ester group is required for part-wise ester naming")
    if len(matches) > 1 and not labels:
        in_parent = _arms_in_parent(mol, matches)
        if in_parent:
            matches = tuple(matches[arm] for arm in in_parent)
    _require_ester_with_its_free_acid(mol, matches)
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    arms = []
    for acyl_carbon, _, ester_oxygen, alkyl_carbon in matches:
        atoms = _branch_atoms(graph, alkyl_carbon, ester_oxygen)
        if acyl_carbon in atoms:
            raise UnsupportedStructure("a ring-closing ester (lactone) is not named part-wise")
        if any(_carbonyl_with_heteroatom(mol, a) for a in atoms):
            raise UnsupportedStructure("an acid or ester group inside the alkyl part outranks this ester")
        if any(mol.GetAtomWithIdx(a).GetIsotope() or mol.GetAtomWithIdx(a).GetNumRadicalElectrons() for a in atoms):
            raise UnsupportedStructure("isotopes and radicals in the alkyl part are not supported")
        arms.append(atoms)
    removed = set().union(*arms)
    if sum(len(a) for a in arms) != len(removed):
        raise UnsupportedStructure("the alkyl parts share atoms, so these esters are of a polyol, not a polyacid")

    ester_labels = {}
    free_hydrogens = []
    free_hydrogen_atoms = []
    labels = dict(labels)
    if labels:
        bridging = {ester_oxygen: i for i, (_, _, ester_oxygen, _) in enumerate(matches)}
        polycarbonic = _polycarbonic(mol)
        carbonyl = {
            m.GetIdx()
            for acyl_carbon, _, _, _ in matches
            for m in mol.GetAtomWithIdx(acyl_carbon).GetNeighbors()
            if m.GetAtomicNum() == 8 and m.GetIdx() not in bridging and not polycarbonic
        }
        for idx in [i for i in labels if i in bridging or i in carbonyl]:
            entry = labels[idx]
            if idx in carbonyl and not entry["skeleton"] and sum(entry["H"].values()) == 1 and mol.GetAtomWithIdx(idx).GetTotalNumHs() == 1:
                free_hydrogens.extend(entry["H"])
                free_hydrogen_atoms.append(idx)
                del labels[idx]
                continue
            if entry["H"] or not entry["skeleton"]:
                continue
            ester_labels[idx] = (entry["skeleton"], bridging.get(idx))
            del labels[idx]

    if original is not None:
        acid = _acid_with_labels(*original, removed, matches, [*ester_labels, *free_hydrogen_atoms])
    else:
        editable = Chem.RWMol(mol)
        for _, _, ester_oxygen, _ in matches:
            oxygen = editable.GetAtomWithIdx(ester_oxygen)
            oxygen.SetNumExplicitHs(1)
            oxygen.SetNoImplicit(True)
        restore_labels(editable, mol, labels, removed)
        for idx in sorted(removed, reverse=True):
            editable.RemoveAtom(idx)
        acid = editable.GetMol()
        Chem.SanitizeMol(acid)
    if len(Chem.GetMolFrags(acid)) != 1:
        raise UnsupportedStructure("the acid parts of these esters are separate groups")
    anion = _anion_name(smiles_to_iupac(Chem.MolToSmiles(acid)))

    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    context = {
        "atoms": {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in removed and a.HasProp("_CIPCode")},
        "bonds": {
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
            for b in probe.GetBonds()
            if b.GetBeginAtomIdx() in removed and b.GetEndAtomIdx() in removed and b.HasProp("_CIPCode")
        },
        "used": set(),
    }
    named = {}
    named_by_arm = {}
    arm_label = {i: nuclide for nuclide, i in ester_labels.values() if i is not None}
    token = BRANCH_STEREO.set(context)
    isotope_context = {"labels": {a: e for a, e in labels.items() if a in removed}, "consumed": set(), "mol": mol}
    isotope_token = ISOTOPE_LABELS.set(isotope_context if isotope_context["labels"] else None)
    try:
        for arm, (acyl_carbon, _, ester_oxygen, alkyl_carbon) in enumerate(matches):
            name, compound = name_branch(graph, alkyl_carbon, ester_oxygen, halogens, mol=mol)
            locant = (f"{arm_label[arm][:-1]}O" if arm in arm_label else "O") if ester_labels else ""
            entry = named.setdefault((name, locant), [0, compound, []])
            entry[0] += 1
            named_by_arm[arm] = name
    finally:
        BRANCH_STEREO.reset(token)
        ISOTOPE_LABELS.reset(isotope_token)
    if set(isotope_context["labels"]) - isotope_context["consumed"]:
        raise UnsupportedStructure("an isotopically modified atom of the alkyl part is not cited by any supported name")
    if any(("atom", a) not in context["used"] for a in context["atoms"]) or any(
        ("bond", b) not in context["used"] for b in context["bonds"]
    ):
        raise UnsupportedStructure("a stereo element of the alkyl part is not cited by any supported name")
    if ester_labels:
        counts = {}
        for nuclide, _ in ester_labels.values():
            counts[nuclide] = counts.get(nuclide, 0) + 1
        text = "(" + ",".join(f"{n}{c}" for n, c in sorted(counts.items())) + ")"
        stem = re.search(r"(?:di|tri|tetra)?carbonate$", anion)
        anion = anion[: stem.start()] + text + anion[stem.start():] if stem else _insert_before_ending(anion, text)
    grouped = {}
    for (name, locant), (count, compound, _) in named.items():
        grouped.setdefault(name, []).append((locant, count, compound))
    located = {}
    if (len(grouped) > 1 or _polycarbonic(mol)) and not ester_labels:
        located = _ester_locants(mol, acid, removed, matches, named_by_arm)
    parts = []
    for name in sorted(grouped, key=alpha_sort_key):
        if name in located:
            locants = located[name]
            compound = grouped[name][0][2]
            word = multiplying_prefix(len(locants), compound=compound) if len(locants) > 1 else ""
            text = f"{word}({name})" if compound and len(locants) > 1 else word + name
            parts.append(f"{','.join(locants)}-{text}")
            continue
        for locant, count, compound in sorted(grouped[name]):
            text = name if count == 1 else (
                f"{multiplied_ester_word(count, name, compound=compound)}({name})" if compound else f"{multiplied_ester_word(count, name, compound=compound)}{name}"
            )
            if locant:
                text = f"{','.join([locant] * count)}-{text}"
            parts.append(text)
    free = len(acid.GetSubstructMatches(_FREE_ACID)) - len(matches)
    if free > 0:
        if free not in _HYDROGEN_WORDS:
            raise UnsupportedStructure("too many free acid groups beside the esters")
        word = _HYDROGEN_WORDS[free]
        if free_hydrogens:
            if free != 1 or len(free_hydrogens) != 1:
                raise UnsupportedStructure("isotopic modification of several acid hydrogens is not supported yet")
            word = f"({free_hydrogens[0]})hydrogen"
        parts.append(word)
    elif free_hydrogens:
        raise UnsupportedStructure("an isotopically modified hydroxy hydrogen of this ester is not supported yet")
    return " ".join(parts + [anion])
