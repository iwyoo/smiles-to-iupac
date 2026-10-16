"""Naming of urea (H2N-C(=O)-NH2) and its N-substituted derivatives, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'urea' is a
  retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- the same way 'carbamic acid' (`_carbamate.py`)
  is retained rather than derived. Confirmed via PubChem structure match:
  `NC(=O)N` -> "urea". Unlike every other module in this project, urea has
  no parent-hydride chain to number and no chain-length logic of its own.
- Each nitrogen may carry 0, 1, or 2 plain, unbranched, saturated alkyl
  substituents, cited as 'N-'-prefixed substituents directly ahead of
  'urea', mirroring `_amide.py`'s/`_carbamate.py`'s own N-substitution
  citation -- confirmed directly from the Blue Book's own text
  (the Blue Book): 'N-methyl-N-nitrosourea
  (PIN)' cites two different substituents on the SAME nitrogen both as
  'N-' (no prime needed -- there's only one substituted nitrogen to name),
  and '(i) The symbols N,N' are used for the 'unprimed' parent
  structure...' establishes that a SECOND, distinct nitrogen's
  substituents are cited with a primed 'N'-' instead (the same convention
  `_common.py`'s multi-nitrogen relatives use, e.g.
  'N,N'-methylenediethanamine (PIN)'). PubChem's own generated names for
  these use numeric locants instead (e.g. "1,3-dimethylurea"), a
  structure-only match, not a naming-convention one -- this project
  follows the Blue Book's letter-locant convention, per the source text
  above.

Each nitrogen carries any number of substituents named through `name_branch` (acyclic, ring, halogenated or
unsaturated); the nitrogen with more substituents takes the unprimed locant, otherwise the one holding the
alphanumerically first substituent (P-14.3.5, P-14.5.2), so 'N,N'-dimethyl-N-phenylurea' and
'N-phenyl-N'-(pyridin-2-yl)urea' come from one citation rule shared with `_chalcogenourea.py`. A ring fused to the
urea core (hydantoin) and groups senior to urea (acids, esters, amides) are out of scope; hydroxy, alkoxy, amino, nitrile and similar groups are prefixes.

- Semicarbazide (H2N-NH-C(=O)-NH2, P-68.3.1.4): PubChem structure match
  confirms `NC(=O)NN` -> "aminourea" -- the unsubstituted parent is named
  as urea carrying a plain 'amino' substituent, not with its own retained
  name. Scope here is deliberately limited to that single unsubstituted
  case; any carbon substituent alongside the amino nitrogen is out of
  scope (raise `UnsupportedStructure`) -- PubChem's own examples for that
  combination (`CNC(=O)NN` -> "1-amino-3-methylurea") use a numeric-locant
  style this project doesn't otherwise follow for urea (see above), so
  mixing amino with alkyl substitution needs its own follow-up scoping
  pass rather than being folded in here.
"""

from rdkit import Chem

from ._chalcogen_imide import ylidene_group_atoms
from ._chalcogenourea import is_urea_substituent_root, n_prefix, n_substituent_names, unprimed_first, unprimed_first_labelled
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, group_substituents
from ._isotope_labels import nitrogen_hydrogen_descriptor, split_isotopes
from ._parent_choice import senior_parent_name
from ._substituents import format_substituent_prefixes


_CHALCOGENS = (8, 16, 34, 52)
_UREA_WORD = {8: "urea", 16: "thiourea", 34: "selenourea", 52: "tellurourea"}
_AMIDE_ENDING = {8: "carboxamide", 16: "carbothioamide", 34: "carboselenoamide", 52: "carbotelluroamide"}


def _urea_cores(mol, atomic_num=8):
    """Every (carbon_idx, (nitrogen1_idx, nitrogen2_idx)) for a urea carbonyl carbon and its two nitrogens (a carbon
    with exactly one double-bonded, terminal chalcogen and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic() or atom.IsInRing():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == atomic_num]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(oxygens) != 1 or len(nitrogens) != 2:
            continue
        (oxygen,) = oxygens
        if oxygen.GetDegree() != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
            continue
        if all(n.IsInRing() for n in nitrogens):
            continue
        if any(
            b.GetBondTypeAsDouble() == 2.0
            and b.GetOtherAtom(n).GetAtomicNum() != 6
            and not ylidene_group_atoms(mol, b.GetOtherAtom(n).GetIdx(), n.GetIdx())
            for n in nitrogens
            for b in n.GetBonds()
        ):
            continue
        if any(
            b.GetBondTypeAsDouble() != 1.0 and not (b.GetBeginAtomIdx() == atom.GetIdx() or b.GetEndAtomIdx() == atom.GetIdx())
            for n in nitrogens
            for b in n.GetBonds()
            if b.GetOtherAtom(n).GetAtomicNum() == 7 and b.GetOtherAtom(n).GetIdx() != atom.GetIdx()
        ):
            continue
        if any(
            not is_urea_substituent_root(mol, nn) and not _is_hydrazine_tail(mol, nn, n.GetIdx())
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        tails = [
            [nn for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx() and _is_hydrazine_tail(mol, nn, n.GetIdx())]
            for n in nitrogens
        ]
        if max(len(t) for t in tails) > 1:
            continue
        yield atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())


def _urea_core(mol, atomic_num=8):
    return next(_urea_cores(mol, atomic_num), None)


def _is_sulfonyl_sulfur(atom):
    return atom.GetAtomicNum() == 16 and atom.GetDegree() == 4 and sum(
        1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1
    ) == 2


def _is_hydrazine_tail(mol, atom, exclude_idx, sulfonyl=False):
    """True if `atom` is the second nitrogen of a semicarbazide: bonded to the nitrogen at `exclude_idx` and to
    carbon groups (or a sulfonyl group with `sulfonyl`) only, or to one ylidene carbon (a semicarbazone)."""
    if atom.GetAtomicNum() != 7 or atom.GetIsAromatic() or atom.IsInRing():
        return False
    if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return False
    if mol.GetBondBetweenAtoms(atom.GetIdx(), exclude_idx).GetBondTypeAsDouble() != 1.0:
        return False
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != exclude_idx]
    orders = [mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in others]
    return all(n.GetAtomicNum() == 6 or (sulfonyl and _is_sulfonyl_sulfur(n)) for n in others) and (
        orders == [2.0] or 2.0 not in orders and 3.0 not in orders
    )


def has_urea_shape(mol) -> bool:
    return _urea_core(mol) is not None


def senior_chalcogen(mol):
    """The chalcogen of the most senior urea group present: urea above thiourea, selenourea and tellurourea (P-41)."""
    return next((z for z in _CHALCOGENS if _urea_core(mol, z) is not None), None)


def _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx):
    """idx of the extra terminal amino nitrogen attached to one of urea's
    two core nitrogens (the semicarbazide shape), or None if neither core
    nitrogen carries one."""
    for n_idx in (n1_idx, n2_idx):
        nitrogen = mol.GetAtomWithIdx(n_idx)
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() in (carbon_idx,):
                continue
            if _is_hydrazine_tail(mol, neighbor, n_idx):
                return neighbor.GetIdx()
    return None


def _labelled_urea_name(mol, labels):
    """N-substituted urea with hydrogen nuclides on its nitrogens: the labelled nitrogen takes the unprimed locant when
    the two are otherwise equivalent (P-82.2.5, P-14.3.5)."""
    core = _urea_core(mol)
    if core is None or len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("only a plain urea core can carry nuclide labels on its nitrogens here")
    carbon_idx, (n1_idx, n2_idx) = core
    if set(labels) - {n1_idx, n2_idx} or _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx) is not None:
        raise UnsupportedStructure("only the two nitrogens of a plain urea can carry nuclide labels here")
    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 8)
    n1_names, n2_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, oxygen_idx, n1_idx, n2_idx}, (n1_idx, n2_idx), carbon_idx, junior_groups=True
    )
    first = unprimed_first_labelled(n1_idx, n2_idx, n1_names, n2_names, labels)
    (unprimed, primed), (unprimed_names, primed_names) = (
        ((n1_idx, n2_idx), (n1_names, n2_names)) if first else ((n2_idx, n1_idx), (n2_names, n1_names))
    )
    locant = {unprimed: "N", primed: "N'"}
    nuclides = nitrogen_hydrogen_descriptor(mol, {n: (locant[n], entry) for n, entry in sorted(labels.items(), key=lambda i: locant[i[0]])})
    positions = {"N": unprimed_names, "N'": primed_names}
    prefixes = format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))
    return f"{prefixes}{nuclides}urea"


def name_urea(mol, atomic_num=8) -> str:
    split = split_isotopes(mol)
    if split is not None:
        if atomic_num != 8:
            raise UnsupportedStructure("nuclide labels on a chalcogenourea are not supported yet")
        return _labelled_urea_name(*split[:2])
    return name_cores(mol, atomic_num)


def name_cores(mol, atomic_num):
    """Name with the senior urea group of the chalcogen `atomic_num` as the parent (P-45.2.1, P-45.5)."""
    cores = list(_urea_cores(mol, atomic_num))
    if not cores:
        raise UnsupportedStructure(
            "no urea (H2N-C(=O)-NH2 or an N-substituted derivative) shape "
            "found; this module only handles urea and simple N-substituted "
            "ureas"
        )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    if len(cores) == 1:
        return _core_name(mol, atomic_num, cores[0])[1]
    by_carbon = {core[0]: core for core in cores}
    return senior_parent_name(mol, list(by_carbon), lambda carbon: _core_name(mol, atomic_num, by_carbon[carbon]))


def _junior_cores(mol, carbon_idx, nitrogens):
    """Carbons of the other urea groups that are not acyl groups on the nitrogens of the parent urea, so they rank
    equal to or below it (P-41) and are cited as N-substituents."""
    return {
        c
        for z in _CHALCOGENS
        for c, _ in _urea_cores(mol, z)
        if c != carbon_idx and not {n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors()} & nitrogens
    }


def _core_name(mol, atomic_num, core):
    """(number of N-substituents, name) with `core` as the parent; every other urea group is an N-substituent."""
    carbon_idx, (n1_idx, n2_idx) = core
    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == atomic_num)
    amino_nitrogen_idx = _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx)

    core_atoms = {carbon_idx, oxygen_idx, n1_idx, n2_idx}
    graph = adjacency(mol)
    tails = [next((t for t in graph[n] if t != carbon_idx and _is_hydrazine_tail(mol, mol.GetAtomWithIdx(t), n)), None) for n in (n1_idx, n2_idx)]
    if all(t is not None for t in tails):
        core_atoms.update(tails)
        units = ((n1_idx, tails[0]), (n2_idx, tails[1]))
        names = n_substituent_names(mol, graph, core_atoms, (*units[0], *units[1]), carbon_idx, junior_groups=True)
        return sum(map(len, names)), _hydrazinecarbohydrazide_name(
            (names[0], names[1]), (names[2], names[3]), _HYDRAZIDE_ENDING[atomic_num]
        )
    if amino_nitrogen_idx is not None:
        core_atoms.add(amino_nitrogen_idx)
        alpha, amide = (
            (n1_idx, n2_idx) if amino_nitrogen_idx in graph[n1_idx] else (n2_idx, n1_idx)
        )
        names = n_substituent_names(
            mol, graph, core_atoms, (amide, alpha, amino_nitrogen_idx), carbon_idx, junior_groups=True
        )
        return sum(map(len, names)), _hydrazinecarboxamide_name(*names, _AMIDE_ENDING[atomic_num])
    n1_names, n2_names = n_substituent_names(
        mol, graph, core_atoms, (n1_idx, n2_idx), carbon_idx, junior_groups=True,
        junior_cores=_junior_cores(mol, carbon_idx, {n1_idx, n2_idx}),
    )
    return len(n1_names) + len(n2_names), f"{n_prefix(n1_names, n2_names)}{_UREA_WORD[atomic_num]}"


def nitrogen_locants(mol):
    """({nitrogen: 'N' | "N'"}, the name with every N locant cited) of an N-substituted urea, None for any other shape."""
    core = _urea_core(mol)
    if core is None or len(Chem.GetMolFrags(mol)) > 1:
        return None
    carbon_idx, (n1_idx, n2_idx) = core
    if _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx) is not None:
        return None
    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 8)
    n1_names, n2_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, oxygen_idx, n1_idx, n2_idx}, (n1_idx, n2_idx), carbon_idx, junior_groups=True
    )
    first = unprimed_first(n1_names, n2_names)
    (unprimed, primed), (unprimed_names, primed_names) = (
        ((n1_idx, n2_idx), (n1_names, n2_names)) if first else ((n2_idx, n1_idx), (n2_names, n1_names))
    )
    positions = {"N": unprimed_names, "N'": primed_names}
    prefixes = format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))
    return {unprimed: "N", primed: "N'"}, f"{prefixes}urea"


def _hydrazinecarboxamide_name(amide_names, alpha_names, beta_names, ending="carboxamide"):
    """Substituted hydrazinecarboxamide (P-66.1.1.1.1.3, P-66.3.5): N on the amide nitrogen, 1 and 2 on the
    hydrazine nitrogens."""
    positions = {"N": amide_names, 1: alpha_names, 2: beta_names}
    grouped = group_substituents({k: v for k, v in positions.items() if v})
    if not grouped:
        return f"hydrazine{ending}"
    numbered = any(isinstance(loc, int) for info in grouped.values() for loc in info["locants"])
    return f"{format_substituent_prefixes(grouped)}hydrazine{'-1-' if numbered else ''}{ending}"


_HYDRAZIDE_ENDING = {8: "carbohydrazide", 16: "carbothiohydrazide", 34: "carboselenohydrazide", 52: "carbotellurohydrazide"}


def _locant_rank(locant):
    return (0, locant) if isinstance(locant, int) else (1, len(locant))


def _hydrazinecarbohydrazide_name(unit_a, unit_b, ending):
    """Hydrazine-1-carbohydrazide (P-68.3.1.2.6): the hydrazine unit with more substituents is the parent chain (1 on
    the nitrogen bonded to the carbonyl carbon, 2 on the other), the other unit is the hydrazide (N, N'); then the
    lowest locants, then the alphanumerically first prefix at the lowest locant."""
    options = []
    for parent, hydrazide in ((unit_a, unit_b), (unit_b, unit_a)):
        positions = {1: parent[0], 2: parent[1], "N": hydrazide[0], "N'": hydrazide[1]}
        grouped = group_substituents({k: v for k, v in positions.items() if v})
        by_name = sorted(grouped, key=alpha_sort_key)
        key = (
            -(len(parent[0]) + len(parent[1])),
            sorted(_locant_rank(loc) for info in grouped.values() for loc in info["locants"]),
            [sorted(_locant_rank(loc) for loc in grouped[name]["locants"]) for name in by_name],
        )
        options.append((key, grouped))
    grouped = min(options, key=lambda option: option[0])[1]
    if not grouped:
        return f"hydrazine{ending}"
    numbered = any(isinstance(loc, int) for info in grouped.values() for loc in info["locants"])
    return f"{format_substituent_prefixes(grouped)}hydrazine{'-1-' if numbered else ''}{ending}"
