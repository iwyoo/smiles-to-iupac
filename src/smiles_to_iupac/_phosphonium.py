"""Naming of simple phosphonium cations (the '-phosphanium' suffix,
P-73.1.1.2), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium' --
  the same rule `_ammonium.py` already applies to amines, just for
  phosphorus's own parent hydride 'phosphane' (PH3) instead of nitrogen's
  'azane'. Since `_phosphane.py`'s own substitutive naming already covers
  0-3 plain alkyl substituents on phosphorus (unlike `_amine.py`, which
  only covers a single primary substituent), this module reuses that
  scope wholesale: a phosphonium cation is neutralized to the equivalent
  phosphane (formal charge 0, one fewer hydrogen), named via
  `name_simple_phosphane`, and the resulting name's terminal 'e' is
  replaced with 'ium'.
- Confirmed via PubChem structure match: `[PH4+]` -> "phosphanium",
  `C[PH3+]` -> "methylphosphanium", `C[PH2+]C` -> "dimethylphosphanium",
  `C[PH+](C)C` -> "trimethylphosphanium". A fourth substituent
  (`C[P+](C)(C)C` -> "tetramethylphosphanium", a genuine quaternary
  phosphonium salt) is NOT reachable by the neutralize-then-rename
  approach above -- a neutral phosphorus atom cannot carry four
  substituents at all, so there is no phosphane name to derive it from.
  Instead, mirroring `_sulfonium.py`'s own direct-substituent-construction
  approach (used there for the same "no neutral counterpart" reason), a
  quaternary phosphonium's name is built directly from its four carbon
  substituents via `format_mononuclear_prefixes` + the '-phosphanium'
  suffix. Confirmed via PubChem structure match for
  `C[P+](C)(C)C` -> "tetramethylphosphanium"; `CC[P+](C)(C)C`'s structure
  is named "ethyltri(methyl)phosphanium" here (the multiplying prefix
  sits outside the parentheses, per P-16.5.1.3.1's own text and the Blue
  Book's "chlorodi(methyl)borane (PIN)"/"ethyldi(methyl)phosphane (PIN)"
  worked examples, `tmp/bluebook/P6.txt`/`P1.html` -- see
  `_phosphane.py`'s docstring for the full derivation), not PubChem's own
  raw "ethyl(trimethyl)phosphanium".
- A plain, unsubstituted benzene ring bonded directly to phosphorus is
  cited as a 'phenyl' substituent here too. The degree 0-3 case already
  gets this for free via `_phosphane.py`'s own identical extension (PR
  #397) through the neutralize-then-rename path above; the quaternary
  degree-4 case needs its own copy of that same detection (mirroring
  `_phosphane.py`'s `plain_phenyl_substituent_atoms`) since it's built
  directly rather than reusing `name_simple_phosphane`. Confirmed via
  PubChem PUG REST: `c1ccccc1[P+](c1ccccc1)(c1ccccc1)c1ccccc1` ->
  "tetraphenylphosphanium" (CID 164912).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any phosphonium phosphorus not shaped like PH4+, a phosphorus bonded to
  1-3 carbons (with the remaining valence as hydrogens), or a phosphorus
  bonded to exactly 4 carbons -- e.g. formal charge other than +1, more
  than one charged atom, isotopic modification, a halogen or other
  heteroatom substituent, a branched/unsaturated substituent, or an
  aromatic substituent other than a plain, unsubstituted phenyl, or a
  ring other than a plain phenyl substituent (inherited unchanged from
  `_phosphane.py`'s own scope for the degree 0-3 case via neutralization;
  independently verified for the quaternary degree-4 case above).
- P-92 stereocenters:
  unlike `_ammonium.py`'s nitrogen (which inverts too fast to be a real
  stereocenter), a phosphonium phosphorus with three or four distinct
  substituents is itself a genuine, configurationally stable stereocenter
  -- confirmed via RDKit `FindPotentialStereo` on `C[PH+](CC)CCC` (degree
  3) as well as the quaternary degree-4 case (both `@`/`@@` currently
  collapse to the same, silently wrong name). Since substituents here are
  always plain unbranched alkyl (never a stereocenter on their own), the
  phosphorus is the only possible stereocenter, and this project has no
  established heteroatom-centered stereodescriptor convention (same
  policy as `_sulfonium.py`/`_sulfinic_acid.py`) -- so a specified
  stereocenter is explicitly rejected rather than silently dropped.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    linear_branch,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
    specified_stereocenters,
)
from ._numerals import alkyl_name
from ._phosphane import name_simple_phosphane
from ._substituents import format_mononuclear_prefixes


def has_phosphonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged phosphorus
    shaped like a genuine phosphonium (PH4+, a phosphorus singly bonded to
    1-3 carbons with the rest hydrogens, or a phosphorus singly bonded to
    4 carbons). Used by `core.py` to route here before `_phosphane.py`,
    which rejects any charged atom outright."""
    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() == 1
    ]
    if len(charged_phosphorus) != 1:
        return False
    phosphorus = charged_phosphorus[0]
    if phosphorus.GetIsotope() != 0:
        return False
    degree = phosphorus.GetDegree()
    if degree > 4 or phosphorus.GetTotalNumHs() + degree != 4:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in phosphorus.GetNeighbors()
    )


def name_phosphonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_phosphorus = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() != 0
    ]
    (phosphorus,) = charged_phosphorus
    if phosphorus.GetFormalCharge() != 1 or phosphorus.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "phosphonium phosphorus is supported (P-73.1.1.2)"
        )
    if specified_stereocenters(mol) is not None:
        # The phosphonium phosphorus is itself a genuine, configurationally
        # stable stereocenter in virtually every real R3HP+/R4P+ molecule
        # (module docstring), and this project has no established way to
        # cite a heteroatom-centered stereodescriptor -- explicitly reject
        # rather than silently drop the marker (P-92).
        raise UnsupportedStructure(
            "a specified stereocenter (the phosphonium phosphorus itself) "
            "is not supported yet (see P-92, module docstring)"
        )
    degree = phosphorus.GetDegree()
    if any(n.GetAtomicNum() != 6 for n in phosphorus.GetNeighbors()):
        raise UnsupportedStructure(
            "a phosphonium substituent other than carbon is out of scope "
            "for this module"
        )
    if any(
        mol.GetBondBetweenAtoms(phosphorus.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
        for n in phosphorus.GetNeighbors()
    ):
        raise UnsupportedStructure("the phosphonium phosphorus must be singly bonded to each substituent")

    if degree == 4:
        return _name_quaternary_phosphonium(mol, phosphorus)

    neutral_rw = Chem.RWMol(mol)
    neutral_phosphorus = neutral_rw.GetAtomWithIdx(phosphorus.GetIdx())
    neutral_phosphorus.SetFormalCharge(0)
    neutral_phosphorus.SetNoImplicit(True)
    neutral_phosphorus.SetNumExplicitHs(3 - degree)
    neutral_mol = neutral_rw.GetMol()
    Chem.SanitizeMol(neutral_mol)

    phosphane_name = name_simple_phosphane(neutral_mol)
    return phosphane_name[:-1] + "ium"


def _name_quaternary_phosphonium(mol, phosphorus) -> str:
    graph = adjacency(mol)
    roots = set(graph[phosphorus.GetIdx()])
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots)

    other_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() != phosphorus.GetIdx()]
    for atom in other_atoms:
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms other than the phosphonium phosphorus itself "
                "are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetIsAromatic() and atom.GetIdx() not in phenyl_atoms:
            raise UnsupportedStructure(
                "an aromatic substituent other than a plain, unsubstituted "
                "phenyl group is out of scope for this module"
            )
    all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if all_ring_atoms - phenyl_atoms:
        raise UnsupportedStructure(
            "a ring other than a plain phenyl substituent directly on "
            "phosphorus is out of scope for this module"
        )
    non_ring_unsaturation = [
        b for b in non_single_bonds(mol) if b[0] not in phenyl_atoms and b[1] not in phenyl_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure("an unsaturated substituent is out of scope for this module")

    substituent_names = []
    for root in roots:
        if root in phenyl_atoms:
            substituent_names.append(("phenyl", False))
            continue
        length = linear_branch(graph, root, phosphorus.GetIdx())
        if length is None:
            raise UnsupportedStructure("a branched substituent is out of scope for this module")
        substituent_names.append((alkyl_name(length), False))

    return format_mononuclear_prefixes(substituent_names) + "phosphanium"
