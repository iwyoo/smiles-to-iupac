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

Scope, deliberately narrow, mirroring `_urea.py`'s/`_thiourea.py`'s own:
substituents landing on the two amino nitrogens (one nitrogen with one or
two, using the same 'N-'/'N,N-di' citation established there; an identical
single substituent on each of the two different amino nitrogens using
'N,N'-di...'; or one DIFFERENT substituent on each, the alphanumerical
order (P-14.5.2, ignoring italicized prefixes like 'tert-') deciding
which becomes 'N-' and which 'N''-', same rule as `_urea.py`/
`_thiourea.py` -- PubChem structure match: `CCNC(=N)NC` ->
'1-ethyl-2-methylguanidine', CID 17814701), the imino nitrogen (at most
one substituent -- it only has one open valence beyond its C=N double
bond), or both at once, combined and alphabetized per the
tetramethyl-phenyl worked example above. Each N-substituent's own name is
built with `name_branch` (P-29 PIN style, fixed project-wide by PR #237;
mirrors `_urea.py`'s/`_thiourea.py`'s identical fix, PR #332/#333) -- a
branched N-substituent is supported (e.g. 'N-tert-butylguanidine', CID
12830400, a retained non-compound name), and a *compound* one is always
parenthesized -- 'N-(propan-2-yl)guanidine', not PubChem's own raw
'N-propan-2-ylguanidine' (CID 11491919), same correction as `_urea.py`
(see that module's docstring for the Blue Book citations), applied per
substituent (single or grouped) in `_combine_prefixes` (e.g.
'N,N'-di(propan-2-yl)guanidine', CID 198192; 'N,N'-ditert-butylguanidine',
CID 23103888). Explicitly out of scope (raise
`UnsupportedStructure`): a different substituent *count* on each amino
nitrogen (no confirmed worked example settles that locant tie-break), an
unsaturated N-substituent, a ring-bearing N-substituent other than a
single plain (unsubstituted) benzene ring, and a ring-fused guanidine.

- A plain benzene ring bonded directly to a nitrogen is cited as
  'phenyl', mirroring `_urea.py`'s/`_thiourea.py`'s identical fix
  (PR #382/#383) and the tetramethyl-phenyl worked example already cited
  above -- PubChem structure match: `c1ccccc1NC(=N)N` -> "2-phenylguanidine"
  (N-phenylguanidine), `c1ccccc1NC(=N)Nc1ccccc1` -> "1,2-diphenylguanidine"
  (N,N'-diphenylguanidine), `CNC(=N)Nc1ccccc1` -> "2-methyl-1-phenylguanidine"
  (N-methyl-N'-phenylguanidine). A *substituted* phenyl ring or a second
  substituent sharing that same nitrogen is out of scope.
"""

from collections import defaultdict

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._numerals import multiplying_prefix
from ._substituents import alpha_sort_key, name_branch


def _guanidine_core(mol):
    """(carbon_idx, imino_idx, (amino1_idx, amino2_idx)) for the guanidine
    carbon and its three nitrogens, or None if the molecule isn't shaped
    like a guanidine core at all (a carbon with exactly one double-bonded
    nitrogen and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
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
            if imino_n.GetTotalNumHs() != 0 or other[0].GetAtomicNum() != 6:
                continue
        else:
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in amino):
            continue
        if any(nn.GetAtomicNum() != 6 for n in amino for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx()):
            continue
        return atom.GetIdx(), imino_n.GetIdx(), (amino[0].GetIdx(), amino[1].GetIdx())
    return None


def has_guanidine_shape(mol) -> bool:
    return _guanidine_core(mol) is not None


def _n_substituent_carbons(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(
        n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx
    )


def _substituent_names(full_graph, nitrogen_idx, substituent_carbons, aromatic_atoms=frozenset(), mol=None):
    return [
        name_branch(full_graph, c, nitrogen_idx, {}, aromatic_atoms, mol=mol)
        for c in substituent_carbons
    ]


def _plain_phenyl_substituent_atoms(mol, graph, roots):
    """Union of ring atoms for every plain, unsubstituted benzene ring in
    `mol` that hangs directly off one of `roots` (a guanidine nitrogen's
    substituent-carbon neighbors) with no other exocyclic attachment --
    i.e. a lone 'phenyl' N-substituent, as opposed to a fused or
    otherwise-substituted ring. Mirrors `_urea.py`'s identical helper."""
    atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = set(ring)
        if not is_plain_benzene_ring(mol, ring_atoms):
            continue
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            continue
        ring_atom, _ = attachment
        if ring_atom in roots:
            atoms |= ring_atoms
    return atoms


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        reached, _ = bfs(carbon_graph, root)
        atoms.update(reached)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def _prime_rank(letter):
    return letter.count("'")


def _di_name(name, is_compound):
    return f"({name})" if is_compound else name


def _combine_prefixes(letters_and_entries):
    """Combine (letter, substituent name, is_compound) triples into a
    single citation, grouping identical substituent names under a shared
    multiplying prefix and alphanumerically ordering groups by name
    (P-14.5.2, via `alpha_sort_key`), per the Blue Book's own
    "N,N,N′,N′-tetramethyl-N′′-phenylguanidine" worked example. A
    multiplied compound name (has its own locant) is parenthesized to
    avoid ambiguity (e.g. 'N,N'-di(propan-2-yl)'); a single occurrence or
    a retained (non-compound) name is not."""
    groups = defaultdict(list)
    compound_of = {}
    for letter, name, is_compound in letters_and_entries:
        groups[name].append(letter)
        compound_of[name] = is_compound
    parts = []
    for name in sorted(groups, key=alpha_sort_key):
        letters = sorted(groups[name], key=lambda l: (_prime_rank(l), l))
        if len(letters) == 1:
            prefix_name = _di_name(name, compound_of[name])
        else:
            prefix_name = multiplying_prefix(len(letters)) + _di_name(name, compound_of[name])
        parts.append(f"{','.join(letters)}-{prefix_name}")
    return "-".join(parts)


def name_guanidine(mol) -> str:
    core = _guanidine_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "no guanidine (HN=C(NH2)2 or an N-/N'-/N''-substituted "
            "derivative) shape found; this module only handles guanidine "
            "and simple nitrogen-substituted guanidines"
        )
    carbon_idx, imino_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    n1_carbons = _n_substituent_carbons(mol, n1_idx, carbon_idx)
    n2_carbons = _n_substituent_carbons(mol, n2_idx, carbon_idx)
    imino_carbons = _n_substituent_carbons(mol, imino_idx, carbon_idx)

    full_graph = adjacency(mol)
    phenyl_atoms = _plain_phenyl_substituent_atoms(
        mol, full_graph, n1_carbons + n2_carbons + imino_carbons
    )
    if phenyl_atoms and (
        (any(c in phenyl_atoms for c in n1_carbons) and len(n1_carbons) > 1)
        or (any(c in phenyl_atoms for c in n2_carbons) and len(n2_carbons) > 1)
    ):
        raise UnsupportedStructure(
            "a phenyl N-substituent alongside another substituent on the "
            "same nitrogen is not supported yet"
        )

    if mol.GetRingInfo().NumRings() > 0:
        all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
        if all_ring_atoms - phenyl_atoms:
            raise UnsupportedStructure(
                "a ring-fused guanidine or a ring N-substituent other than "
                "a plain, unsubstituted benzene ring is out of scope for "
                "this module"
            )

    carbon_graph = carbon_adjacency(mol)
    n1_chain_atoms = _substituent_chain_atoms(carbon_graph, n1_carbons)
    n2_chain_atoms = _substituent_chain_atoms(carbon_graph, n2_carbons)
    imino_chain_atoms = _substituent_chain_atoms(carbon_graph, imino_carbons)
    known_atoms = (
        {carbon_idx, imino_idx, n1_idx, n2_idx}
        | n1_chain_atoms
        | n2_chain_atoms
        | imino_chain_atoms
    )
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in known_atoms:
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the "
                "guanidine core and its plain N-alkyl substituents is not "
                "supported yet"
            )

    _reject_unsaturated_substituents(mol, n1_chain_atoms - phenyl_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms - phenyl_atoms)
    _reject_unsaturated_substituents(mol, imino_chain_atoms - phenyl_atoms)

    aromatic_atoms = frozenset(phenyl_atoms)
    n1_names = _substituent_names(full_graph, n1_idx, n1_carbons, aromatic_atoms, mol=mol)
    n2_names = _substituent_names(full_graph, n2_idx, n2_carbons, aromatic_atoms, mol=mol)
    imino_names = _substituent_names(full_graph, imino_idx, imino_carbons, aromatic_atoms, mol=mol)

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1:
            raise UnsupportedStructure(
                "a different substituent count on each of guanidine's two "
                "amino nitrogens is not supported yet (no confirmed "
                "worked example settles the locant tie-break for that "
                "case)"
            )
        first, second = sorted((n1_names[0], n2_names[0]), key=lambda e: alpha_sort_key(e[0]))
        amino_entries = [("N", *first), ("N'", *second)]
    else:
        letter = "N" if n1_names else "N'"
        names = n1_names or n2_names
        amino_entries = [(letter, name, is_compound) for name, is_compound in names]

    imino_entries = [("N''", name, is_compound) for name, is_compound in imino_names]

    entries = amino_entries + imino_entries
    if not entries:
        return "guanidine"
    return f"{_combine_prefixes(entries)}guanidine"
