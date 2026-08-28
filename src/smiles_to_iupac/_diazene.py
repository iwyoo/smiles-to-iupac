"""Naming of simple diazenes/azo compounds (HN=NH bearing 0-2 unbranched,
saturated alkyl substituents, one per nitrogen), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-68.3.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  azo compounds R-N=N-R' are preferably named substitutively as
  substituted derivatives of the parent hydride 'diazene' (HN=NH) --
  mirroring `_phosphane.py`'s own "central atom is the parent hydride,
  alkyl groups are substituent prefixes" structure, e.g. 'dimethyldiazene'
  (not an amine-style name). Confirmed via PubChem PUG REST: CID 123195
  (`N=N`) -> "diazene", CID 123421 (`CN=N`) -> "methyldiazene", CID 10421
  (`CN=NC`) -> "dimethyldiazene", CID 526060 (`CN=NCC`) ->
  "ethyl(methyl)diazene", CID 13183 (`CCN=NCC`) -> "diethyldiazene".
- Unlike phosphane's single trivalent phosphorus (0-3 substituents),
  diazene's two nitrogens each have only one bond free after their N=N
  double bond, so each carries at most one substituent -- but the same
  P-16.5.1.3.1 parenthesization rule and ordinary multiplying-prefix rule
  `_phosphane.py` already implements for its own mononuclear-parent
  substituent list apply unchanged here (confirmed for the two-different-
  substituents case above, 'ethyl(methyl)diazene' -- alphabetically-first
  unparenthesized, the other parenthesized).
- No locants are ever needed (not just P-14.3.4.2(a)'s "always 1"
  omission, but a genuine absence of ambiguity): diazene's two nitrogens
  are interchangeable by the molecule's own symmetry, so a single
  substituent or a differing pair has only one possible structure
  regardless of which nitrogen is nominally "first" -- confirmed by every
  PubChem example above citing no locants at all.
- P-35.2.1: a halogen substituent on a carbon branch (not directly on a
  diazene nitrogen) coexists freely, reusing `_hydrazine.py`'s own
  `name_branch`/halogen-aware branch-forking pre-check. Confirmed via
  PubChem PUG REST: CID 76139658 (`ClCCN=N`) -> "2-chloroethyldiazene"
  (the sole-substituent case concatenates directly, no parens, even
  though the substituent itself is a compound name -- matching PubChem's
  own auto-generated name exactly), CID 56629798 (`ClCCN=NCCCl`) ->
  "bis(2-chloroethyl)diazene" (two identical compound substituents use
  the compound 'bis' multiplying prefix, P-14.2.2, rather than the plain
  'di' `format_mononuclear_prefixes` uses for simple substituents).
  **Deliberately conservative scope**: the case of two *different*
  substituents where at least one is a compound (halogen-bearing) name is
  not verified against any worked example and is rejected explicitly
  rather than guessed at -- only (a) zero or one substituent, and (b) two
  identical substituents, are supported when a compound name is involved;
  two different plain (non-compound) substituents keep working exactly as
  before via the existing `format_mononuclear_prefixes` path.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any atom other than the two diazene nitrogens, carbon, halogen, and
  hydrogen.
- Two different substituents where at least one is a compound
  (halogen-bearing) name -- unverified, see above.
- A branched (real carbon fork) or unsaturated substituent, an aromatic
  substituent, or any ring anywhere in the molecule.
- More than one N=N unit (bis(azo) compounds), azoxy compounds (an N-oxide
  of this group, P-68.3.1.3.3), or hydrazine-type N-N single-bonded
  compounds (P-68.3.1.2, a structurally unrelated parent).
- Charged or isotopically modified atoms.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, non_single_bonds
from ._numerals import multiplying_prefix
from ._substituents import format_mononuclear_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _is_unbranched_ignoring_halogens(graph, root, coming_from, halogens):
    """True iff the branch's carbon skeleton (halogen leaves set aside) is
    a straight, non-forking chain (see `_hydrazine.py`'s identical
    helper)."""
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous and n not in halogens]
        if len(neighbors) > 1:
            return False
        if not neighbors:
            return True
        previous, current = current, neighbors[0]


def _diazene_nitrogens(mol):
    """The two nitrogens of a plain diazene skeleton: N=N (double bond),
    each nitrogen degree <= 2 (the double bond plus at most one other
    single-bonded neighbor), formal charge 0 -- or None if `mol` isn't
    shaped this way."""
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(nitrogens) != 2:
        return None
    n1, n2 = nitrogens
    bond = mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 2.0:
        return None
    if n1.GetDegree() > 2 or n2.GetDegree() > 2:
        return None
    if n1.GetFormalCharge() != 0 or n2.GetFormalCharge() != 0:
        return None
    return n1, n2


def has_diazene_shape(mol) -> bool:
    return _diazene_nitrogens(mol) is not None


def name_diazene(mol) -> str:
    nitrogens = _diazene_nitrogens(mol)
    if nitrogens is None:
        raise UnsupportedStructure("no plain diazene (N=N) skeleton found; this module only handles diazenes")
    n1, n2 = nitrogens

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the diazene's own two nitrogens "
                "(P-68.3.1.3.2) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic substituent is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    n1_idx, n2_idx = n1.GetIdx(), n2.GetIdx()
    if any(a not in (n1_idx, n2_idx) and b not in (n1_idx, n2_idx) for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturation in a substituent is out of scope for this module"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    substituents = []
    for n_idx in (n1_idx, n2_idx):
        (root,) = [n for n in graph[n_idx] if n not in (n1_idx, n2_idx)] or (None,)
        if root is None:
            continue
        if root in halogens:
            raise UnsupportedStructure(
                "a halogen bonded directly to a diazene nitrogen is out "
                "of scope for this module (see module docstring)"
            )
        if not _is_unbranched_ignoring_halogens(graph, root, n_idx, halogens):
            raise UnsupportedStructure("a branched substituent is out of scope for this module")
        substituents.append(name_branch(graph, root, n_idx, halogens))

    if not substituents:
        return "diazene"
    if len(substituents) == 1:
        (name, _), = substituents
        return name + "diazene"

    (name_a, compound_a), (name_b, compound_b) = substituents
    if name_a == name_b:
        if compound_a:
            return f"{multiplying_prefix(2, compound=True)}({name_a})diazene"
        return multiplying_prefix(2) + name_a + "diazene"
    if compound_a or compound_b:
        raise UnsupportedStructure(
            "two different substituents where at least one is a compound "
            "(halogen-bearing) name is unverified and out of scope for "
            "this module (see module docstring)"
        )
    return format_mononuclear_prefixes([name_a, name_b]) + "diazene"
