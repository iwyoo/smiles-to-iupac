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

from ._chalcogenourea import is_urea_substituent_root, n_prefix, n_substituent_names, unprimed_first
from ._common import UnsupportedStructure, adjacency, group_substituents
from ._substituents import format_substituent_prefixes


_AMIDE_ENDING = {8: "carboxamide", 16: "carbothioamide", 34: "carboselenoamide", 52: "carbotelluroamide"}


def _urea_core(mol, atomic_num=8):
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
        amino_neighbors = [
            nn
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx() and _is_hydrazine_tail(mol, nn, n.GetIdx())
        ]
        if len(amino_neighbors) > 1:
            continue
        return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    return None


def _is_hydrazine_tail(mol, atom, exclude_idx):
    """True if `atom` is the second nitrogen of a semicarbazide: bonded to the nitrogen at `exclude_idx` and to
    carbon groups only, or to one ylidene carbon (a semicarbazone)."""
    if atom.GetAtomicNum() != 7 or atom.GetIsAromatic() or atom.IsInRing():
        return False
    if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return False
    if mol.GetBondBetweenAtoms(atom.GetIdx(), exclude_idx).GetBondTypeAsDouble() != 1.0:
        return False
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != exclude_idx]
    orders = [mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in others]
    return all(n.GetAtomicNum() == 6 for n in others) and (orders == [2.0] or 2.0 not in orders and 3.0 not in orders)


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
            if _is_hydrazine_tail(mol, neighbor, n_idx):
                return neighbor.GetIdx()
    return None


def name_urea(mol, atomic_num=8) -> str:
    core = _urea_core(mol, atomic_num)
    if core is None:
        raise UnsupportedStructure(
            "no urea (H2N-C(=O)-NH2 or an N-substituted derivative) shape "
            "found; this module only handles urea and simple N-substituted "
            "ureas"
        )
    carbon_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == atomic_num)
    amino_nitrogen_idx = _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx)

    core_atoms = {carbon_idx, oxygen_idx, n1_idx, n2_idx}
    graph = adjacency(mol)
    if amino_nitrogen_idx is not None:
        core_atoms.add(amino_nitrogen_idx)
        alpha, amide = (
            (n1_idx, n2_idx) if amino_nitrogen_idx in graph[n1_idx] else (n2_idx, n1_idx)
        )
        amide_names, alpha_names, beta_names = n_substituent_names(
            mol, graph, core_atoms, (amide, alpha, amino_nitrogen_idx), carbon_idx, junior_groups=True
        )
        return _hydrazinecarboxamide_name(amide_names, alpha_names, beta_names, _AMIDE_ENDING[atomic_num])
    n1_names, n2_names = n_substituent_names(mol, graph, core_atoms, (n1_idx, n2_idx), carbon_idx, junior_groups=True)
    return f"{n_prefix(n1_names, n2_names)}urea"


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
