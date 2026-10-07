"""Naming of simple sulfonium cations (the '-sulfanium' suffix,
P-73.1.1.2, bearing 0-3 substituent groups named by `name_branch`),
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium'.
  For sulfur this is 'sulfane' (H2S) -- unlike `_ammonium.py`'s own
  '-aminium' derivation, this module can't reuse an existing neutral
  module's output the same way: a neutralized R-SH2+ is just R-SH, an
  ordinary thiol already named by `_thiol.py`'s own suffix-style
  '-thiol' convention ('methanethiol'), which has no textual relationship
  to 'methylsulfanium' at all. Sulfane's own retained-name substitutive
  style ('methylsulfane' for a hypothetical neutral CH3-SH2, mirroring
  `_phosphane.py`'s identical 'methylphosphane' pattern for PH3) is a
  different, otherwise-unused naming branch that only ever surfaces
  through this cation -- so this module builds the '-sulfanium' name
  directly from the charged sulfur's own substituents, without ever
  constructing (or exposing to `core.py`'s general dispatch) a neutral
  'sulfane' molecule or name.
- Confirmed via PubChem structure match: `[SH3+]` -> "sulfanium",
  `C[SH2+]` -> "methylsulfanium", `C[SH+]C` -> "dimethylsulfanium",
  `C[S+](C)C` -> "trimethylsulfanium". A charged sulfur's cation valence
  is 3, the same as phosphane's own neutral valence -- so this module's
  substituent-count range (0-3) and prefix formatting
  (`format_mononuclear_prefixes`, P-16.5.1.3.1's parenthesization rule)
  are lifted directly from `_phosphane.py`.
- A plain, unsubstituted benzene ring bonded directly to the sulfonium
  sulfur is cited as a 'phenyl' substituent, mixed freely with alkyl
  substituents -- confirmed via PubChem PUG REST: `c1ccccc1[SH2+]` ->
  "phenylsulfanium" (CID 12099100), `c1ccccc1[S+](c1ccccc1)c1ccccc1` ->
  "triphenylsulfanium" (CID 61344). Mirrors `_phosphonium.py`'s
  identical quaternary-phosphonium extension (PR #398); the sulfonium
  cation's own valence of 3 never triggers a lambda-convention label,
  same as that module's degree 0-3 phosphane-reuse path. Detection
  reuses the same `plain_phenyl_substituent_atoms` pattern
  (`is_plain_benzene_ring` + `ring_chain_attachment`).

Explicitly out of scope (raise `UnsupportedStructure`):
- A substituent that carries a principal characteristic group (oxo,
  hydroxy, carboxy, amino, cyano, ...), a halogen, or a heteroatom other
  than an ether-type chalcogen link; substituents are named by
  `name_branch`.
- P-92 stereocenters:
  unlike `_oxonium.py`'s/`_ammonium.py`'s nitrogen/oxygen (which invert too
  fast to be a real stereocenter), a sulfonium sulfur with three distinct
  substituents is itself a genuine, configurationally stable stereocenter,
  confirmed via RDKit `FindPotentialStereo` on `C[S+](CC)CCC` (both `@`/`@@`
  currently collapse to the same, silently wrong name). Since substituents
  here are always plain unbranched alkyl (never a stereocenter on their
  own), the sulfur is the only possible stereocenter, and this project has
  no established heteroatom-centered stereodescriptor convention (same
  policy as `_sulfinic_acid.py`) -- so a specified stereocenter is
  explicitly rejected rather than silently dropped.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, specified_stereocenters
from ._onium_prefixes import onium_name


def has_sulfonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged sulfur shaped
    like a genuine sulfonium (SH3+, or a sulfur singly bonded to 1-3
    carbons with the rest hydrogens). Used by `core.py` to route here
    before any other branch, none of which recognize a charged atom."""
    charged_sulfurs = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16 and atom.GetFormalCharge() == 1
    ]
    if len(charged_sulfurs) != 1:
        return False
    sulfur = charged_sulfurs[0]
    if sulfur.GetIsotope() != 0:
        return False
    degree = sulfur.GetDegree()
    if degree > 3 or sulfur.GetTotalNumHs() + degree != 3:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in sulfur.GetNeighbors()
    )


def name_sulfonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_sulfurs = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16 and atom.GetFormalCharge() != 0
    ]
    (sulfur,) = charged_sulfurs
    if sulfur.GetFormalCharge() != 1 or sulfur.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "sulfonium sulfur is supported (P-73.1.1.2)"
        )
    if specified_stereocenters(mol) is not None:
        # The sulfonium sulfur is itself a genuine, configurationally
        # stable stereocenter in virtually every real R3S+/R2HS+ molecule
        # (module docstring), and this project has no established way to
        # cite a heteroatom-centered stereodescriptor -- explicitly reject
        # rather than silently drop the marker (P-92).
        raise UnsupportedStructure(
            "a specified stereocenter (the sulfonium sulfur itself) is "
            "not supported yet (see P-92, module docstring)"
        )
    if sulfur.GetDegree() > 3:
        raise UnsupportedStructure(
            "a sulfur atom with more than three substituents is not a "
            "sulfonium"
        )
    if any(n.GetAtomicNum() != 6 for n in sulfur.GetNeighbors()):
        raise UnsupportedStructure(
            "a sulfonium substituent other than carbon is out of scope "
            "for this module"
        )
    if any(
        mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
        for n in sulfur.GetNeighbors()
    ):
        raise UnsupportedStructure("the sulfonium sulfur must be singly bonded to each substituent")

    return onium_name(mol, sulfur, "sulfanium")
