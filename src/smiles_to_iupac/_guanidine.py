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

from ._chalcogenourea import is_core_substituent_root, n_substituent_names
from ._common import UnsupportedStructure, adjacency, group_substituents
from ._substituents import alpha_sort_key, format_substituent_prefixes


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
            if imino_n.GetTotalNumHs() != 0 or not is_core_substituent_root(mol, other[0]):
                continue
        else:
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in amino):
            continue
        if any(
            not is_core_substituent_root(mol, nn) for n in amino for nn in n.GetNeighbors() if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        return atom.GetIdx(), imino_n.GetIdx(), (amino[0].GetIdx(), amino[1].GetIdx())
    return None


def has_guanidine_shape(mol) -> bool:
    return _guanidine_core(mol) is not None


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

    n1_names, n2_names, imino_names = n_substituent_names(
        mol, adjacency(mol), {carbon_idx, imino_idx, n1_idx, n2_idx}, (n1_idx, n2_idx, imino_idx), carbon_idx
    )
    unprimed, primed = _amino_assignment(n1_names, n2_names)
    positions = {"N": unprimed, "N'": primed, "N''": imino_names}
    grouped = group_substituents({k: v for k, v in positions.items() if v})
    return f"{format_substituent_prefixes(grouped)}guanidine"


def _amino_assignment(n1_names, n2_names):
    """The amino nitrogen with more substituents takes the unprimed locant, then the one whose substituent comes
    first alphanumerically (P-66.4.1.2.1.2: the minimum number of primes)."""
    if len(n1_names) != len(n2_names):
        return (n1_names, n2_names) if len(n1_names) > len(n2_names) else (n2_names, n1_names)
    key = lambda names: min((alpha_sort_key(name) for name, _ in names), default="")
    first, second = sorted((n1_names, n2_names), key=key)
    return first, second
