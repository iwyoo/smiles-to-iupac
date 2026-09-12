"""Naming of the seven 1989 IUPAC steroid parent ring hydrides (Rule 2.1
and Rule 3S-2.2/2.3/2.4 Table 1, https://iupac.qmul.ac.uk/steroid/3S02a.html,
the document Blue Book P-6 defers to for steroid parent names), each by
hardcoded retained-name recognition of its exact unsubstituted structure.

Each skeleton in this lineage is the previous one plus one defined
extension, built on the plain perhydrocyclopenta[a]phenanthrene ring system
(three fused six-membered rings plus one five-membered ring, angularly
ortho-fused):

- gonane (Rule 2.1): the bare ring system, no angular methyls, no C17 side
  chain. PubChem CID 6857523.
- androstane (Rule 3S-2.3): gonane + angular methyls at C10/C13. CID 6857536.
- estrane (Rule 3S-2.2): gonane + only the C13 angular methyl (no C10
  methyl). CID 5460658.
- pregnane (Rule 3S-2.4, Table 1): androstane + a plain two-carbon ethyl
  side chain at C17 (C20-C21). CID 439513.
- cholane: pregnane's side chain extended to five carbons (C20-C24) with a
  methyl branch at C20. CID 6857459.
- cholestane: cholane's side chain extended by a symmetric gem-dimethyl
  terminus (C25-C27). CID 6857534.
- ergostane: cholestane's side chain extended by one more methyl branch at
  C24 (C28). CID 6857535.

Each entry below was independently confirmed by removing the next
skeleton's side-chain/methyl carbons from its canonical SMILES and
re-canonicalizing, reproducing the previous skeleton's own canonical key
exactly (e.g. removing androstane's two angular methyls reproduces gonane's
key; removing ergostane's nine-carbon C17 side chain reproduces
androstane's key).

Since each module's only job is recognizing one exact unsubstituted parent,
a whole-molecule canonical-SMILES match is sufficient: any substituent,
extra methyl, longer side chain, or unsaturation changes the canonical
SMILES and so is correctly left unmatched, falling through to
`_polycyclic.py`'s general von Baeyer engine (which already names the bare
gonane skeleton as "tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane" if nothing
here intercepts it first -- confirmed by direct testing).

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the input
mol, so any stereo-specified ring-fusion or side-chain input naturally
fails to match these stereo-free canonical keys and falls through the same
way -- none of Rule 2.1/3S-2.2/3S-2.3 specify ring-fusion stereochemistry,
and Table 1's cholane/cholestane/ergostane side-chain stereocenters (C20,
and C24 for ergostane) are likewise left to the general engine once
specified.
"""

from rdkit import Chem

_PARENT_HYDRIDES = {
    "C1CCCC2CCC3C(C12)CCC4C3CCC4": "gonane",
    "CC12CCCC1C3CCC4CCCCC4(C3CC2)C": "androstane",
    "CC12CCCC1C1CCC3CCCCC3C1CC2": "estrane",
    "CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "pregnane",
    "CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholane",
    "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholestane",
    "CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "ergostane",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for smiles, name in _PARENT_HYDRIDES.items()}


def has_steroid_parent_hydride_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_steroid_parent_hydride(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
