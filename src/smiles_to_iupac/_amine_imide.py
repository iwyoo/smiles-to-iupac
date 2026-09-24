"""Naming of amine imides (R3N+-N(-)-R, P-74.2.1.3), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-74.2.1.3 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf,
  ~2938-2951): an amine imide is named as a zwitterion based on hydrazine
  "in order not to break the nitrogen chain" -- the -1-charged nitrogen is
  always cited as the '-ide' suffix at locant 1, the +1-charged nitrogen
  as the '-ium' suffix at locant 2 (fixed by which nitrogen carries which
  charge, not a free numbering choice the way plain `_hydrazine.py`'s own
  two interchangeable nitrogens are). Confirmed worked example
  `1,2,2-trimethylhydrazin-2-ium-1-ide (PIN)`: RDKit-sanitizable as
  `C[N-][N+](C)(C)C` -- N1 (-1 charge, degree 2: the N-N bond plus one
  methyl) / N2 (+1 charge, degree 4: the N-N bond plus three methyls).
- Reuses `_hydrazine.py`'s own substituent-collection helpers
  (`_substituent_names`/`_group`) directly rather than duplicating them --
  the substituent shape (plain alkyl/halogenated-alkyl/phenyl branches
  off either nitrogen) is identical to the neutral parent, only the
  charge pattern and the fixed (not tie-broken) locant assignment differ.
  `_hydrazine.py`'s own `_candidate_key`/direction tie-break isn't needed
  here since the anion/cation charge already fixes which nitrogen is
  locant 1 vs 2.
- Suffix assembly: `"hydrazin"` (hydrazine's own stem, trailing 'e'
  elided since `-ium` starts with a vowel, mirroring `_ammonium.py`'s
  identical `amine`->`aminium` elision) + `f"-2-ium-1-ide"`, substituent
  locants 1/2 assigned to whichever nitrogen is the anion/cation
  respectively.

Scope: exactly one N-N single bond, one nitrogen -1 charged (degree <= 2:
the N-N bond plus at most one substituent), the other +1 charged (degree
<= 4: the N-N bond plus up to three substituents), net molecule charge 0,
substituents restricted to `_hydrazine.py`'s own existing scope (plain
alkyl/halogenated-alkyl chains, a plain unsubstituted phenyl ring).

Explicitly out of scope (raise `UnsupportedStructure`): nitrile imide
(P-74.2.2.2.1.1 -- the cationic nitrogen's own substituent is a
triple-bonded 'ylidyne' carbon fragment there, not a plain alkyl group;
a separate, unresearched mechanism), phosphine imide, more than one
+1/-1 nitrogen pair, any ring/unsaturated substituent beyond
`_hydrazine.py`'s own scope, any heteroatom beyond the two nitrogens, a
halogen bonded directly to either nitrogen, or any other charged or
isotopically modified atom.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents, non_single_bonds
from ._hydrazine import _group, _substituent_names
from ._substituents import format_substituent_prefixes

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _amine_imide_nitrogens(mol):
    """(anion_n, cation_n) if `mol` has P-74.2.1.3's amine-imide charge
    pattern -- an N-N single bond, one nitrogen -1 charged (degree <= 2),
    the other +1 charged (degree <= 4) -- else None."""
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(nitrogens) != 2:
        return None
    n1, n2 = nitrogens
    bond = mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    charges = {n1.GetFormalCharge(), n2.GetFormalCharge()}
    if charges != {-1, 1}:
        return None
    anion, cation = (n1, n2) if n1.GetFormalCharge() == -1 else (n2, n1)
    if anion.GetDegree() > 2 or cation.GetDegree() > 4:
        return None
    return anion, cation


def has_amine_imide_shape(mol) -> bool:
    return _amine_imide_nitrogens(mol) is not None


def name_amine_imide(mol) -> str:
    nitrogens = _amine_imide_nitrogens(mol)
    if nitrogens is None:
        raise UnsupportedStructure(
            "no amine-imide (N-N single bond, one -1/one +1 charged "
            "nitrogen) shape found; this module only handles amine imides"
        )
    anion, cation = nitrogens

    aromatic_atoms = frozenset(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the amine imide's own two "
                "nitrogens and halogen substituents are not supported yet"
            )
        if atom.GetIdx() not in (anion.GetIdx(), cation.GetIdx()) and atom.GetFormalCharge() != 0:
            raise UnsupportedStructure("more than one charged nitrogen pair is not supported yet")
        if atom.GetIsotope() != 0:
            raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    if any(
        not (a in aromatic_atoms and b in aromatic_atoms)
        for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure("an unsaturated substituent is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    names_anion = _substituent_names(graph, anion.GetIdx(), cation.GetIdx(), halogens, aromatic_atoms, mol=mol)
    names_cation = _substituent_names(graph, cation.GetIdx(), anion.GetIdx(), halogens, aromatic_atoms, mol=mol)

    entries = [(1, name, compound) for name, compound in names_anion] + [
        (2, name, compound) for name, compound in names_cation
    ]
    grouped = _group(entries)
    return format_substituent_prefixes(grouped) + "hydrazin-2-ium-1-ide"
