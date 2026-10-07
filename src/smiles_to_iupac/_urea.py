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
urea core (hydantoin) and characteristic groups outside it are out of scope.

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

from ._chalcogenourea import n_prefix, n_substituent_names
from ._common import UnsupportedStructure, adjacency


def _urea_core(mol):
    """(carbon_idx, (nitrogen1_idx, nitrogen2_idx)) for the urea carbonyl
    carbon and its two nitrogens, or None if the molecule isn't shaped like
    a urea core at all (a carbon with exactly one double-bonded, terminal
    oxygen and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
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
        if any(
            nn.GetAtomicNum() != 6 and not _is_terminal_amino_nitrogen(nn, n.GetIdx())
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        amino_neighbors = [
            nn
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx() and _is_terminal_amino_nitrogen(nn, n.GetIdx())
        ]
        if len(amino_neighbors) > 1:
            continue
        return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    return None


def _is_terminal_amino_nitrogen(atom, exclude_idx):
    """True if `atom` is a plain terminal -NH2 nitrogen (semicarbazide's
    extra nitrogen) bonded only to the nitrogen at `exclude_idx`."""
    if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
        return False
    if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return False
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    return neighbors == [exclude_idx]


def has_urea_shape(mol) -> bool:
    return _urea_core(mol) is not None


def _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx):
    """idx of the extra terminal amino nitrogen attached to one of urea's
    two core nitrogens (the semicarbazide shape), or None if neither core
    nitrogen carries one."""
    for n_idx in (n1_idx, n2_idx):
        nitrogen = mol.GetAtomWithIdx(n_idx)
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() in (carbon_idx,):
                continue
            if _is_terminal_amino_nitrogen(neighbor, n_idx):
                return neighbor.GetIdx()
    return None


def name_urea(mol) -> str:
    core = _urea_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "no urea (H2N-C(=O)-NH2 or an N-substituted derivative) shape "
            "found; this module only handles urea and simple N-substituted "
            "ureas"
        )
    carbon_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 8)
    amino_nitrogen_idx = _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx)

    core_atoms = {carbon_idx, oxygen_idx, n1_idx, n2_idx}
    if amino_nitrogen_idx is not None:
        core_atoms.add(amino_nitrogen_idx)
    n1_names, n2_names = n_substituent_names(mol, adjacency(mol), core_atoms, n1_idx, n2_idx, carbon_idx)

    if amino_nitrogen_idx is not None:
        if n1_names or n2_names:
            raise UnsupportedStructure(
                "a semicarbazide (amino-substituted urea nitrogen) "
                "combined with a plain N-alkyl substituent is not "
                "supported yet -- PubChem's own examples for that "
                "combination use a numeric-locant style this module "
                "doesn't otherwise follow for urea"
            )
        return "aminourea"

    return f"{n_prefix(n1_names, n2_names)}urea"
