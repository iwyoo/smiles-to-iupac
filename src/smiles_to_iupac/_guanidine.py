"""Naming of guanidine (HN=C(NH2)2) and its N-substituted derivatives on
all three nitrogens, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'guanidine'
  is a retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- the same way 'urea' (`_urea.py`) and
  'thiourea' (`_thiourea.py`) are retained rather than derived. Confirmed
  via PubChem structure match: `NC(=N)N` -> "guanidine". Like urea,
  guanidine has no parent-hydride chain to number and no chain-length
  logic of its own.
- P-66.4.1.2.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  guanidine's three nitrogens (one imino =NH, two amino -NH2) use the
  locants N, N', and N'' in preferred IUPAC names -- confirmed worked
  examples: "N,N,N′,N′-tetramethyl-N′′-phenylguanidine (PIN)" (both amino
  nitrogens fully substituted as N/N', the imino nitrogen substituted as
  N'') and "N,N′-dimethylguanidine (PIN) (not N,N′′-dimethylguanidine)"
  (one substituent on each of the two amino nitrogens, imino nitrogen
  left unsubstituted) -- the latter is structurally identical to
  `_urea.py`'s own N,N'- case, just with guanidine's extra (here,
  unsubstituted) imino nitrogen. P-66.4.1.2.1.2 also notes a
  minimum-number-of-primes tie-break for "when the position of the double
  bond is unknown" -- that ambiguity doesn't arise here at all, since an
  RDKit molecule's C=N double bond position is always explicit, so which
  nitrogen is 'imino' vs 'amino' is read directly off the graph, never
  guessed.
- Confirmed via PubChem structure match (PubChem itself uses the
  deprecated numeral locants, not the N/N'/N'' PIN style, so these are
  structure-only checks): `CNC(=N)N` -> "2-methylguanidine",
  `CN(C)C(=N)N` -> "1,1-dimethylguanidine", `CNC(=N)NC` ->
  "1,2-dimethylguanidine" (this last one is the same structure as the
  Blue Book's own "N,N′-dimethylguanidine (PIN)" worked example above),
  `CN=C(N)N` -> "2-methylguanidine" (imino nitrogen alone, same numeral
  locant PubChem uses for the amino case -- structure-only confirms the
  shape, not which letter this project assigns), `CCN=C(N)N` ->
  "2-ethylguanidine", `CN=C(NC)N` -> "1,2-dimethylguanidine" (imino +
  one amino nitrogen substituted).

Every nitrogen carries any substituent `name_branch` can name (acyclic, ring, halogenated, unsaturated): the amino
nitrogen with more substituents takes the unprimed locant, then the one holding the alphanumerically first
substituent, and the imino nitrogen takes N''. A ring fused to the guanidine core is out of scope.
"""

from rdkit import Chem

from ._chalcogenourea import (
    c_n_configuration,
    is_core_substituent_root,
    is_urea_substituent_root,
    n_substituent_names,
    unprimed_first,
    unprimed_first_labelled,
)
from ._common import UnsupportedStructure, adjacency, alpha_sort_key, group_substituents, multiplied_word
from ._isotope_labels import nitrogen_hydrogen_descriptor, split_isotopes
from ._substituents import format_substituent_prefixes
from ._urea import _is_hydrazine_tail


def _is_guanidine_substituent_root(mol, atom):
    return is_core_substituent_root(mol, atom) or (atom.GetAtomicNum() != 7 and is_urea_substituent_root(mol, atom))


def _carbonyl_carbon(mol, atom):
    """A carbon double-bonded to a chalcogen: an acyl group on a guanidine nitrogen makes the amide the parent
    (P-41, P-66.4.1.2.1.3)."""
    return atom.GetAtomicNum() == 6 and any(
        n.GetAtomicNum() in (8, 16, 34, 52) and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in atom.GetNeighbors()
    )


def _guanidine_core(mol, hydrazine=False):
    """(carbon_idx, imino_idx, (amino1_idx, amino2_idx)) for the guanidine
    carbon and its three nitrogens, or None if the molecule isn't shaped
    like a guanidine core at all (a carbon with exactly one double-bonded
    nitrogen and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides). With `hydrazine`, one amino nitrogen is the first nitrogen of a hydrazine
    instead (P-68.3.1.2.4)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic() or atom.IsInRing():
            continue
        nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
        if len(nitrogens) != 3:
            continue
        imino = [
            n
            for n in nitrogens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        amino = [
            n
            for n in nitrogens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(imino) != 1 or len(amino) != 2:
            continue
        (imino_n,) = imino
        if imino_n.GetFormalCharge() != 0 or imino_n.GetIsotope() != 0:
            continue
        if imino_n.GetDegree() == 1:
            if imino_n.GetTotalNumHs() != 1:
                continue
        elif imino_n.GetDegree() == 2:
            other = [n for n in imino_n.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
            if imino_n.GetTotalNumHs() != 0 or not _is_guanidine_substituent_root(mol, other[0]):
                continue
        else:
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in amino):
            continue
        tails = [
            nn.GetIdx()
            for n in amino
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx() and _is_hydrazine_tail(mol, nn, n.GetIdx(), sulfonyl=True)
        ]
        if len(tails) != hydrazine:
            continue
        if any(
            not _is_guanidine_substituent_root(mol, nn) and nn.GetIdx() not in tails
            for n in amino
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        if any(
            _carbonyl_carbon(mol, nn)
            for n in (*amino, imino_n)
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        return atom.GetIdx(), imino_n.GetIdx(), (amino[0].GetIdx(), amino[1].GetIdx())
    return None


def has_guanidine_shape(mol) -> bool:
    return _guanidine_core(mol) is not None or _guanidine_core(mol, hydrazine=True) is not None


def has_carbonic_amide_parent(mol) -> bool:
    """A urea, thiourea or guanidine group: a class that outranks the element hydride whose atom is only a substituent
    on its nitrogen (P-41, P-44.1.1)."""
    from ._urea import has_urea_shape

    return has_urea_shape(mol) or has_guanidine_shape(mol)


_AMIDINE_CARBON = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")


def _amidine_groups(mol, hydrazine):
    """(hydrazine nitrogen, carbon, imino nitrogen, amino nitrogen) of every amidine group on a hydrazine nitrogen."""
    groups = []
    for h in hydrazine:
        for carbon in mol.GetAtomWithIdx(h).GetNeighbors():
            if (
                carbon.GetAtomicNum() != 6
                or carbon.GetDegree() != 3
                or carbon.IsInRing()
                or carbon.GetIsAromatic()
                or carbon.GetFormalCharge()
                or mol.GetBondBetweenAtoms(h, carbon.GetIdx()).GetBondTypeAsDouble() != 1.0
            ):
                continue
            nitrogens = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != h]
            orders = {n.GetIdx(): mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in nitrogens}
            imino = [n for n in nitrogens if orders[n.GetIdx()] == 2.0]
            amino = [n for n in nitrogens if orders[n.GetIdx()] == 1.0]
            substituent_roots = [
                nn for n in (*imino, *amino) for nn in n.GetNeighbors() if nn.GetIdx() != carbon.GetIdx()
            ]
            if any(not _is_guanidine_substituent_root(mol, nn) or _carbonyl_carbon(mol, nn) for nn in substituent_roots):
                continue
            if len(imino) == 1 and len(amino) == 1:
                groups.append((h, carbon.GetIdx(), imino[0].GetIdx(), amino[0].GetIdx()))
    return groups


def _hydrazine_polycarboximidamide_name(mol, graph, hydrazine, groups):
    """Hydrazine-1,2-dicarboximidamide and its geminal forms: the amino nitrogen of a group is N and its imino nitrogen
    N', both with the position of the group as superscript; the second of two geminal groups continues with N'' and
    N''', the more substituted group taking the lower locants (P-66.4.1.4.1, P-66.4.1.4.2)."""
    core = set(hydrazine)
    for _, carbon, imino, amino in groups:
        core |= {carbon, imino, amino}
    nitrogens = (*hydrazine, *(n for _, _, imino, amino in groups for n in (amino, imino)))
    names = n_substituent_names(mol, graph, core, nitrogens, groups[0][1], junior_groups=True, hydrazide_junior=True)
    on_hydrazine = dict(zip(hydrazine, names))
    on_group = [(names[2 + 2 * i], names[3 + 2 * i]) for i in range(len(groups))]
    rank = {}
    for h in hydrazine:
        here = sorted(
            (i for i, group in enumerate(groups) if group[0] == h),
            key=lambda i: (-sum(map(len, on_group[i])), min((alpha_sort_key(name) for name, _ in (*on_group[i][0], *on_group[i][1])), default="")),
        )
        rank.update({i: r for r, i in enumerate(here)})
    best = None
    for numbering in ((1, 2), (2, 1)):
        position = dict(zip(hydrazine, numbering))
        located = [(position[h], entry) for h in hydrazine for entry in on_hydrazine[h]]
        n_located = [
            ((position[group[0]], primes), entry)
            for i, group in enumerate(groups)
            for nitrogen_names, primes in zip(on_group[i], (2 * rank[i], 2 * rank[i] + 1))
            for entry in nitrogen_names
        ]
        key = (
            sorted(position[group[0]] for group in groups),
            sorted(p for p, _ in located),
            sorted(locant for locant, _ in n_located),
            sorted((entry[0], str(locant)) for locant, entry in (*located, *n_located)),
        )
        if best is None or key < best[0]:
            best = (key, located, n_located)
    _, located, n_located = best
    by_locant = {}
    for locant, entry in located:
        by_locant.setdefault(locant, []).append(entry)
    for (p, primes), entry in n_located:
        mark = "'" * primes
        by_locant.setdefault(f"N{mark}{p}", []).append(entry)
    prefixes = format_substituent_prefixes(group_substituents(by_locant)) if by_locant else ""
    suffix_locants = ",".join(map(str, best[0][0]))
    return f"{prefixes}hydrazine-{suffix_locants}-{multiplied_word(len(groups), 'carboximidamide')}"


def _hydrazinecarboximidamide_name(mol, core):
    """Hydrazine-1-carboximidamide with N on the amino nitrogen of the amidine, N' on its imino nitrogen and 1 and 2 on
    the hydrazine nitrogens (P-66.4.1.1, P-68.3.1.2.4)."""
    carbon_idx, imino_idx, amino = core
    graph = adjacency(mol)
    alpha = next(
        n for n in amino if any(_is_hydrazine_tail(mol, mol.GetAtomWithIdx(t), n, sulfonyl=True) for t in graph[n] if t != carbon_idx)
    )
    amide = next(n for n in amino if n != alpha)
    tail = next(
        t for t in graph[alpha] if t != carbon_idx and _is_hydrazine_tail(mol, mol.GetAtomWithIdx(t), alpha, sulfonyl=True)
    )
    groups = _amidine_groups(mol, (alpha, tail))
    if any(match[0] not in {group[1] for group in groups} for match in mol.GetSubstructMatches(_AMIDINE_CARBON)):
        raise UnsupportedStructure("an amidine group away from the hydrazine nitrogens is not named with it (P-66.4.1.1)")
    if len(groups) > 1:
        return _hydrazine_polycarboximidamide_name(mol, graph, (alpha, tail), groups)
    names = n_substituent_names(
        mol,
        graph,
        {carbon_idx, imino_idx, alpha, amide, tail},
        (amide, imino_idx, alpha, tail),
        carbon_idx,
        junior_groups=True,
        hydrazide_junior=True,
    )
    grouped = group_substituents({k: v for k, v in zip(("N", "N'", 1, 2), names) if v})
    numbered = any(isinstance(loc, int) for info in grouped.values() for loc in info["locants"])
    return f"{_configuration_prefix(mol, carbon_idx)}{format_substituent_prefixes(grouped) if grouped else ''}hydrazine{'-1-' if numbered else ''}carboximidamide"


def _labelled_guanidine_name(mol, labels):
    """Guanidine with hydrogen nuclides on its nitrogens: the descriptor with their N, N', N'' locants stands in front of
    the retained name (P-82.2.5, P-66.4.1.2.1.1); of two equivalent amino nitrogens the labelled one takes the unprimed
    locant (P-14.3.5)."""
    core = _guanidine_core(mol)
    if core is None or len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("only a plain guanidine core can carry nuclide labels on its nitrogens here")
    carbon_idx, imino_idx, (n1_idx, n2_idx) = core
    if set(labels) - {n1_idx, n2_idx, imino_idx}:
        raise UnsupportedStructure("only the three nitrogens of a guanidine can carry nuclide labels here")
    n1_names, n2_names, imino_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, imino_idx, n1_idx, n2_idx}, (n1_idx, n2_idx, imino_idx), carbon_idx, junior_groups=True
    )
    first = unprimed_first_labelled(n1_idx, n2_idx, n1_names, n2_names, labels)
    (unprimed, primed), (unprimed_names, primed_names) = (
        ((n1_idx, n2_idx), (n1_names, n2_names)) if first else ((n2_idx, n1_idx), (n2_names, n1_names))
    )
    locant = {unprimed: "N", primed: "N'", imino_idx: "N''"}
    nuclides = nitrogen_hydrogen_descriptor(mol, {n: (locant[n], entry) for n, entry in sorted(labels.items(), key=lambda i: locant[i[0]])})
    prefixes = format_substituent_prefixes(group_substituents({k: v for k, v in {"N": unprimed_names, "N'": primed_names, "N''": imino_names}.items() if v}))
    return f"{prefixes}{nuclides}guanidine"


def _configuration_prefix(mol, carbon_idx):
    """'(Z)-' for the specified C=N bond: the name holds no locant for it (P-93.4.2.1.3)."""
    found = c_n_configuration(mol, [carbon_idx])
    return f"({found[0][1]})-" if found else ""


def name_guanidine(mol) -> str:
    split = split_isotopes(mol)
    if split is not None:
        return _labelled_guanidine_name(*split[:2])
    core = _guanidine_core(mol)
    if core is None and len(Chem.GetMolFrags(mol)) == 1:
        hydrazine_core = _guanidine_core(mol, hydrazine=True)
        if hydrazine_core is not None:
            return _hydrazinecarboximidamide_name(mol, hydrazine_core)
    if core is None:
        raise UnsupportedStructure(
            "no guanidine (HN=C(NH2)2 or an N-/N'-/N''-substituted "
            "derivative) shape found; this module only handles guanidine "
            "and simple nitrogen-substituted guanidines"
        )
    carbon_idx, imino_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(
        _is_guanidine_carbon(mol, n.GetIdx())
        for idx in (n1_idx, n2_idx, imino_idx)
        for n in mol.GetAtomWithIdx(idx).GetNeighbors()
        if n.GetIdx() != carbon_idx
    ):
        raise UnsupportedStructure("a condensed guanidine (biguanide, imidodicarbonimidic diamide) is named on its own parent")

    n1_names, n2_names, imino_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, imino_idx, n1_idx, n2_idx}, (n1_idx, n2_idx, imino_idx), carbon_idx, junior_groups=True
    )
    unprimed, primed = _amino_assignment(n1_names, n2_names)
    positions = {"N": unprimed, "N'": primed, "N''": imino_names}
    grouped = group_substituents({k: v for k, v in positions.items() if v})
    return f"{_configuration_prefix(mol, carbon_idx)}{format_substituent_prefixes(grouped)}guanidine"


def _is_guanidine_carbon(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
    return (
        atom.GetAtomicNum() == 6
        and not atom.IsInRing()
        and atom.GetDegree() == 3
        and len(nitrogens) == 3
        and sorted(mol.GetBondBetweenAtoms(idx, n.GetIdx()).GetBondTypeAsDouble() for n in nitrogens) == [1.0, 1.0, 2.0]
    )


def nitrogen_locants(mol):
    """({nitrogen: 'N' | "N'" | "N''"}, the name with every N locant cited) of a guanidine, None for any other shape."""
    core = _guanidine_core(mol)
    if core is None or len(Chem.GetMolFrags(mol)) > 1:
        return None
    carbon_idx, imino_idx, (n1_idx, n2_idx) = core
    n1_names, n2_names, imino_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, imino_idx, n1_idx, n2_idx}, (n1_idx, n2_idx, imino_idx), carbon_idx, junior_groups=True
    )
    first = unprimed_first(n1_names, n2_names)
    (unprimed, primed), (unprimed_names, primed_names) = (
        ((n1_idx, n2_idx), (n1_names, n2_names)) if first else ((n2_idx, n1_idx), (n2_names, n1_names))
    )
    positions = {"N": unprimed_names, "N'": primed_names, "N''": imino_names}
    prefixes = format_substituent_prefixes(group_substituents({k: v for k, v in positions.items() if v}))
    return {unprimed: "N", primed: "N'", imino_idx: "N''"}, f"{prefixes}guanidine"


def _amino_assignment(n1_names, n2_names):
    """The amino nitrogen with more substituents takes the unprimed locant, then the one whose substituent comes
    first alphanumerically (P-66.4.1.2.1.2: the minimum number of primes)."""
    return (n1_names, n2_names) if unprimed_first(n1_names, n2_names) else (n2_names, n1_names)
