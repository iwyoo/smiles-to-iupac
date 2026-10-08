"""Naming of simple oxonium cations (the '-oxidanium' suffix,
P-73.1.1.2, bearing 0-3 unbranched, saturated alkyl and/or plain phenyl
substituents), per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-73.1.1.2 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/PDF/P7.pdf):
  a cation formed by adding a hydron to a parent hydride is named by
  changing the parent hydride name's terminal 'e' to the suffix 'ium'.
  For oxygen this is 'oxidane' (H2O) -- the same situation `_sulfonium.py`
  already documented for sulfur's own 'sulfane': this module can't reuse
  an existing neutral module's output, since a neutralized R-OH2+ is just
  ordinary R-OH/R-O-R', already named by `_alcohol.py`'s/`_ether.py`'s own
  unrelated conventions ('methanol'/'methoxymethane'), which have no
  textual relationship to 'methyloxidanium' at all. Oxidane's own
  retained-name substitutive style ('methyloxidane' for a hypothetical
  neutral CH3-OH2, mirroring `_phosphane.py`'s 'methylphosphane' pattern
  for PH3) is a different, otherwise-unused naming branch that only ever
  surfaces through this cation -- so this module builds the '-oxidanium'
  name directly from the charged oxygen's own substituents, exactly the
  way `_sulfonium.py` does for sulfur.
- Confirmed via PubChem structure match: `[OH3+]` -> "oxidanium",
  `C[OH2+]` -> "methyloxidanium", `CC[OH2+]` -> "ethyloxidanium",
  `C[OH+]C` -> "dimethyloxidanium", `C[O+](C)C` -> "trimethyloxidanium".
  A charged oxygen's cation valence is 3, the same as phosphane's own
  neutral valence -- so this module's substituent-count range (0-3) and
  prefix formatting (`format_mononuclear_prefixes`, P-16.5.1.3.1's
  parenthesization rule) are lifted directly from `_phosphane.py`/
  `_sulfonium.py`.
- A plain, unsubstituted benzene ring bonded directly to the oxonium
  oxygen is cited as a 'phenyl' substituent, mixed freely with alkyl
  substituents -- confirmed via PubChem PUG REST: `c1ccccc1[OH2+]` ->
  "phenyloxidanium" (CID 5152889),
  `c1ccccc1[O+](c1ccccc1)c1ccccc1` -> "triphenyloxidanium" (CID
  3474029). Mirrors `_sulfonium.py`'s identical extension (PR #399);
  the oxonium cation's own valence of 3 never triggers a
  lambda-convention label. Detection reuses the same
  `plain_phenyl_substituent_atoms` pattern (`is_plain_benzene_ring` +
  `ring_chain_attachment`).

Explicitly out of scope (raise `UnsupportedStructure`):
- A branched or unsaturated substituent, or an aromatic substituent other
  than a plain, unsubstituted phenyl (a substituted or heteroaromatic
  ring, e.g.), or a ring other than a plain phenyl substituent directly
  on oxygen.
- A halogen substituent directly on oxygen (unverified for this cation,
  unlike `_phosphane.py`'s own confirmed halophosphane case).
- Any oxonium oxygen not shaped like OH3+ or an oxygen bonded to 1-3
  carbons (with the remaining valence as hydrogens) -- e.g. formal charge
  other than +1, more than one charged atom, isotopic modification.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._onium_prefixes import onium_name


def has_oxonium_shape(mol) -> bool:
    """True if the molecule contains exactly one +1-charged oxygen shaped
    like a genuine oxonium (OH3+, or an oxygen singly bonded to 1-3
    carbons with the rest hydrogens). Used by `core.py` to route here
    before any other branch, none of which recognize a charged atom."""
    charged_oxygens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 1
    ]
    if len(charged_oxygens) != 1 or sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge()) != 1:
        return False
    oxygen = charged_oxygens[0]
    if oxygen.GetIsotope() != 0:
        return False
    degree = oxygen.GetDegree()
    if degree > 3 or oxygen.GetTotalNumHs() + degree != 3:
        return False
    return all(
        n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(oxygen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        for n in oxygen.GetNeighbors()
    )


def name_oxonium(mol) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    charged_oxygens = [
        atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() != 0
    ]
    (oxygen,) = charged_oxygens
    if oxygen.GetFormalCharge() != 1 or oxygen.GetIsotope() != 0:
        raise UnsupportedStructure(
            "only a single, singly-charged, non-isotopically-modified "
            "oxonium oxygen is supported (P-73.1.1.2)"
        )
    if oxygen.GetDegree() > 3:
        raise UnsupportedStructure(
            "an oxygen atom with more than three substituents is not an "
            "oxonium"
        )
    if any(n.GetAtomicNum() != 6 for n in oxygen.GetNeighbors()):
        raise UnsupportedStructure(
            "an oxonium substituent other than carbon is out of scope "
            "for this module"
        )
    if any(
        mol.GetBondBetweenAtoms(oxygen.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
        for n in oxygen.GetNeighbors()
    ):
        raise UnsupportedStructure("the oxonium oxygen must be singly bonded to each substituent")

    return onium_name(mol, oxygen, "oxidanium")
