"""Naming of a 'hidden amide' pseudoketone -- an acyl group (R-CO-) bonded
directly to the nitrogen of an otherwise-plain, saturated monocyclic ring
(e.g. 1-acetylpiperidine), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-64.1.2.1(b) (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  a carbonyl group bonded to a heteroatom of a ring or ring system is a
  'pseudoketone'; when that ring heteroatom is nitrogen, "the compound has
  been called an unexpressed or 'hidden amide'" -- confirmed worked
  example: '1-(piperidin-1-yl)ethan-1-one (PIN)' / '1-acetylpiperidine
  (a hidden amide)'.
- P-64.3.2: such compounds are named substitutively with the 'one' suffix
  on the acyl chain, citing the ring as a substituent-group prefix (e.g.
  'piperidin-1-yl') -- "This method is preferred to that using acyl
  groups" (i.e. not 'N-acetylpiperidine'-style). Second confirmed worked
  example: '1-(piperidin-1-yl)propan-1-one (PIN)' (= 1-propanoylpiperidine).
- The ring substituent name reuses `_hetero_monocyclic.saturated_ring_name`
  (already used by `_ketone.py`'s hetero-ring-ketone path) -- the ring's
  sole nitrogen is always locant 1 for a single-heteroatom ring (P-22.2.1),
  so the substituent name is always '<stem>-1-yl' (e.g. 'piperidine' ->
  'piperidin-1-yl'), cited as a compound (parenthesized) prefix at
  position 1 of the acyl chain -- the same position the '-one' suffix's
  own locant lands on, confirmed by every worked example above.
- The acyl chain otherwise reuses this project's ordinary chain-naming
  machinery (`ordered_chain`/`name_branch`/`format_substituent_prefixes`),
  restricted to a plain unbranched, saturated chain with halogen
  substituents allowed (mirrors `_carboxylic_acid.py`'s/`_ester.py`'s
  acyl-side first-pass scope) -- e.g. 'ClCC(=O)N1CCCCC1' ->
  '2-chloro-1-(piperidin-1-yl)ethan-1-one' (PubChem CID 222312's own
  '2-chloro-1-piperidin-1-ylethanone', with this project's usual
  explicit-locant/parenthesization convention applied, same divergence
  already documented in `_carbamate.py`/`_ester.py`).

`core.py` must route to this module *before* both
`_ketone.py`'s `has_hetero_ring_ketone_shape` (which claims any
single-heteroatom saturated ring alongside any oxygen in the molecule,
without checking the ring itself actually contains a ketone, and would
otherwise misname this shape's acyl branch as an "N-alkyl substituent" via
`_linear_alkyl_substituent`) and `_amide.py`'s `has_amide_shape` (whose
`name_amide` unconditionally rejects any non-benzene ring in the
molecule) -- neither module supports this shape today.

Explicitly out of scope (raise `UnsupportedStructure`):
- R = H (a formyl group on the ring nitrogen) -- a distinct suffix
  construction, 'ring-1-carbaldehyde' (PubChem CID 17429's own
  'piperidine-1-carbaldehyde' for 'O=CN1CCCCC1'), not a pseudoketone at
  all; out of scope here.
- A branched or unsaturated acyl chain.
- A ring with more than one heteroatom (e.g. morpholine/piperazine --
  attaching via the nitrogen there is locant 4, not 1, and needs separate
  verification), a heteroatom other than nitrogen, a ring size other than
  5/6/7 (`saturated_ring_name`'s own P-22.2.1 scope, mirroring
  `_ketone.py`'s identical `_HETERO_RING_SIZES` restriction), any
  substituent on the ring itself, or any ring unsaturation.
- Any other heteroatom, charged/isotopically modified atom, or
  multi-fragment structure.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    non_single_bonds,
    ordered_chain,
)
from ._hetero_monocyclic import saturated_ring_name
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_RING_SIZES = (5, 6, 7)


def _is_carbonyl_carbon(mol, carbon_atom):
    return any(
        o.GetAtomicNum() == 8
        and o.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        for o in carbon_atom.GetNeighbors()
    )


def _hidden_amide_shape(mol):
    """(ring_nitrogen_idx, ring_size, acyl_carbon_idx) if `mol` is a
    single plain, unsubstituted, saturated monocyclic ring (5/6/7-membered)
    with one nitrogen heteroatom whose sole exocyclic bond is a single
    bond to an acyl (carbonyl) carbon -- else None."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = set(ring_info.AtomRings()[0])
    if len(ring_atoms) not in _RING_SIZES:
        return None
    if any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring_atoms):
        return None
    heteroatoms = [a for a in ring_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    (n_idx,) = heteroatoms
    n_atom = mol.GetAtomWithIdx(n_idx)
    if n_atom.GetAtomicNum() != 7 or n_atom.GetFormalCharge() != 0 or n_atom.GetIsotope() != 0:
        return None
    exo = [n for n in n_atom.GetNeighbors() if n.GetIdx() not in ring_atoms]
    if len(exo) != 1:
        return None
    (acyl_atom,) = exo
    if acyl_atom.GetAtomicNum() != 6:
        return None
    if mol.GetBondBetweenAtoms(n_idx, acyl_atom.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    if not _is_carbonyl_carbon(mol, acyl_atom):
        return None
    for a in ring_atoms:
        if a == n_idx:
            continue
        atom = mol.GetAtomWithIdx(a)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if any(n.GetIdx() not in ring_atoms for n in atom.GetNeighbors()):
            return None
        if atom.GetTotalNumHs() != 2:
            return None
    if any(a in ring_atoms and b in ring_atoms for a, b, _ in non_single_bonds(mol)):
        return None
    return n_idx, len(ring_atoms), acyl_atom.GetIdx()


def has_hidden_amide_shape(mol) -> bool:
    return _hidden_amide_shape(mol) is not None


def name_hidden_amide_ketone(mol) -> str:
    shape = _hidden_amide_shape(mol)
    if shape is None:
        raise UnsupportedStructure(
            "no hidden-amide pseudoketone (acyl group on a plain "
            "saturated monocyclic ring nitrogen) shape found"
        )
    n_idx, ring_size, acyl_atom_idx = shape

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    acyl_atom = mol.GetAtomWithIdx(acyl_atom_idx)
    (carbonyl_oxygen,) = [n.GetIdx() for n in acyl_atom.GetNeighbors() if n.GetAtomicNum() == 8]

    graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    chain = ordered_chain(carbon_graph, acyl_atom_idx, n_idx, set())
    if chain is None:
        raise UnsupportedStructure("a branched acyl chain is not supported yet")
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a formyl group (R = H) attached to a ring nitrogen uses a "
            "separate 'ring-1-carbaldehyde' suffix construction, out of "
            "scope for this module (see P-64.1.2.1)"
        )

    chain_set = set(chain)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])
    halogens = halogen_substituents(mol)
    allowed_atoms = ring_atoms | chain_set | {carbonyl_oxygen}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in allowed_atoms:
            continue
        if idx in halogens:
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetIdx() in chain_set:
                continue
        raise UnsupportedStructure(
            "heteroatoms/substituents other than the acyl chain's halogens "
            "are not supported yet for this hidden-amide pseudoketone path"
        )

    other_non_single = [
        b
        for b in non_single_bonds(mol)
        if not ({b[0], b[1]} == {acyl_atom_idx, carbonyl_oxygen})
    ]
    if other_non_single:
        raise UnsupportedStructure("unsaturation in the acyl chain is not supported yet")

    element = "N"
    stem = saturated_ring_name(element, ring_size)
    if stem is None:
        raise UnsupportedStructure(
            f"no retained/Hantzsch-Widman name for a {ring_size}-membered "
            f"{element}-heteroatom saturated ring (P-22.2.1)"
        )
    ring_yl = f"{stem[:-1]}-1-yl"

    substituents = {1: [(ring_yl, True)]}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [x for x in graph[atom] if x not in chain_set and x not in (n_idx, carbonyl_oxygen)]
        for root in branch_roots:
            entry = name_branch(graph, root, atom, halogens, mol=mol)
            substituents.setdefault(position, []).append(entry)

    grouped = group_substituents(substituents)
    chain_length = len(chain)
    stem_name = alkane_name(chain_length)
    base = stem_name[:-1]
    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{base}-1-one"
