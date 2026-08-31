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
two, using the same 'N-'/'N,N-di' citation established there, or an
identical single substituent on each of the two different amino nitrogens
using 'N,N'-di...'), the imino nitrogen (at most one substituent -- it
only has one open valence beyond its C=N double bond), or both at once,
combined and alphabetized per the tetramethyl-phenyl worked example
above. Explicitly out of scope (raise `UnsupportedStructure`): two
DIFFERENT substituents split across the two amino nitrogens, a different
substituent count on each amino nitrogen (no confirmed worked example
settles the locant tie-break for either case), a
branched/unsaturated/ring-bearing N-substituent, and a ring-fused
guanidine.
"""

from collections import defaultdict

from rdkit import Chem

from ._common import UnsupportedStructure, carbon_adjacency, linear_branch, non_single_bonds
from ._numerals import alkyl_name, multiplying_prefix


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


def _substituent_names(carbon_graph, substituent_carbons):
    names = []
    for c in substituent_carbons:
        length = linear_branch(carbon_graph, c, None)
        if length is None:
            raise UnsupportedStructure("a branched N-substituent is not supported yet")
        names.append(alkyl_name(length))
    return names


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        previous, current = None, root
        while current is not None:
            atoms.add(current)
            neighbors = [n for n in carbon_graph[current] if n != previous]
            previous, current = current, (neighbors[0] if neighbors else None)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def _prime_rank(letter):
    return letter.count("'")


def _combine_prefixes(letters_and_names):
    """Combine (letter, substituent name) pairs into a single citation,
    grouping identical substituent names under a shared multiplying
    prefix and alphabetizing groups by name (P-14.5.2), per the Blue
    Book's own "N,N,N′,N′-tetramethyl-N′′-phenylguanidine" worked
    example."""
    groups = defaultdict(list)
    for letter, name in letters_and_names:
        groups[name].append(letter)
    parts = []
    for name in sorted(groups):
        letters = sorted(groups[name], key=lambda l: (_prime_rank(l), l))
        prefix_name = name if len(letters) == 1 else multiplying_prefix(len(letters)) + name
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
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a ring-fused guanidine is out of scope for this module"
        )

    n1_carbons = _n_substituent_carbons(mol, n1_idx, carbon_idx)
    n2_carbons = _n_substituent_carbons(mol, n2_idx, carbon_idx)
    imino_carbons = _n_substituent_carbons(mol, imino_idx, carbon_idx)

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

    _reject_unsaturated_substituents(mol, n1_chain_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms)
    _reject_unsaturated_substituents(mol, imino_chain_atoms)

    n1_names = _substituent_names(carbon_graph, n1_carbons)
    n2_names = _substituent_names(carbon_graph, n2_carbons)
    imino_names = _substituent_names(carbon_graph, imino_carbons)

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1 or n1_names[0] != n2_names[0]:
            raise UnsupportedStructure(
                "different substituents, or a different substituent count, "
                "split across guanidine's two amino nitrogens is not "
                "supported yet (no confirmed worked example settles the "
                "locant tie-break in that case)"
            )
        amino_entries = [("N", n1_names[0]), ("N'", n2_names[0])]
    else:
        letter = "N" if n1_names else "N'"
        names = n1_names or n2_names
        amino_entries = [(letter, name) for name in names]

    imino_entries = [("N''", name) for name in imino_names]

    entries = amino_entries + imino_entries
    if not entries:
        return "guanidine"
    return f"{_combine_prefixes(entries)}guanidine"
