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
- campestane: ergostane's C24 epimer -- identical constitution (a
  5,6-dimethylheptan-2-yl side chain), opposite configuration at that
  carbon. CID 6857532.
- poriferastane: a 5-ethyl-6-methylheptan-2-yl side chain (ergostane's
  extra C24 methyl replaced by an ethyl group). CID 6857528.
- stigmastane: poriferastane's C24 epimer -- identical constitution,
  opposite configuration at that carbon. CID 6857438.

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
mol, so a stereo-specified ring-fusion or side-chain input only matches if
it's the exact natural configuration listed below (Rule 3S-2.2/2.3/2.4's
own worked structures, PubChem's own isomeric SMILES for each CID above)
-- any other stereoisomer (a ring-fusion epimer, a partially-specified
input, or a side-chain epimer at C20/C24 not matching the natural series)
still falls through to the general engine unmatched, since it produces a
different canonical SMILES string. None of Rule 2.1/3S-2.2/3S-2.3 assign
ring-fusion stereochemistry beyond this one natural configuration per
skeleton -- a differently configured stereoisomer (e.g. 5-beta) needs its
own descriptor and lookup row, out of scope here (tracked separately).

Each skeleton's own PubChem CID above sometimes leaves the C5 ring-fusion
stereocenter unspecified in its own isomeric SMILES record (a PubChem
data-depiction quirk, not a chemistry fact -- C5 is a genuine stereocenter
in every skeleton here) -- confirmed for gonane/androstane/estrane/
cholestane/ergostane/campestane/poriferastane/stigmastane (only pregnane
and cholane's own CID records happen to already be fully specified). Each
of those 8 skeletons therefore also has a second natural-configuration
entry below, from that skeleton's own dedicated "5alpha-<name>" PubChem
name/CID, which does fully specify C5 -- additive alongside the
already-registered partially-specified one, so a real input matching
either still resolves correctly.
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

_PARENT_HYDRIDES = {
    "C1CCCC2CCC3C(C12)CCC4C3CCC4": "gonane",
    "CC12CCCC1C3CCC4CCCCC4(C3CC2)C": "androstane",
    "CC12CCCC1C1CCC3CCCCC3C1CC2": "estrane",
    "CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "pregnane",
    "CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholane",
    "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholestane",
    "CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "ergostane",
    # Natural-configuration isomeric SMILES, one per skeleton above (same
    # PubChem CIDs), additive alongside the stereo-free entries -- these
    # are the shape essentially every real-world instance of these
    # compounds actually has.
    "C1CCC2CC[C@H]3[C@@H]4CCC[C@H]4CC[C@@H]3[C@H]2C1": "gonane",
    "C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@]4([C@H]3CC2)C": "androstane",
    "C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@H]4[C@H]3CC2": "estrane",
    "CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C": "pregnane",
    "CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C": "cholane",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "cholestane",
    "C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "ergostane",
    "C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "campestane",
    "CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C": "poriferastane",
    "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C": "stigmastane",
    # Fully C5-specified natural-configuration entries (see module
    # docstring), from each skeleton's own "5alpha-<name>" PubChem CID.
    "C1CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H]2CC[C@H]4[C@H]3CCC4": "gonane",
    "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C": "androstane",
    "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@H]4[C@H]3CC2": "estrane",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "cholestane",
    "C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "ergostane",
    "C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "campestane",
    "CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C": "poriferastane",
    "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C": "stigmastane",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for smiles, name in _PARENT_HYDRIDES.items()}
assert len(_CANONICAL_TO_NAME) == len(_PARENT_HYDRIDES), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row (see module docstring)"
)



# Androstane's own bare-skeleton constitution (no stereo), in the same atom
# order as every stereo-specified androstane entry above (confirmed by
# substructure match: parsing this string always yields an identity
# mapping onto each specified entry's atom indices, since adding/removing
# only stereo markers never reorders RDKit's own SMILES atom parse order).
# `_ANDROSTANE_FUSION_ATOMS` are its six ring-fusion stereocenters (the two
# quaternary angular-methyl carbons, C10/C13, plus the four ring-fusion CH
# carbons, C5/C8/C9/C14); `_ANDROSTANE_C5_ATOM` (index 9) is C5 specifically
# -- the one ring-fusion stereocenter whose configuration distinguishes the
# 5-alpha (natural) and 5-beta series, confirmed via `rdCIPLabeler`: it's
# 'R' in the natural-configuration entry above and 'S' in real 5-beta-
# androstane (etiocholane, PubChem CID 6857462), while every other fusion
# stereocenter's CIP label is unchanged between the two.
_ANDROSTANE_PLAIN = Chem.MolFromSmiles("CC12CCCC1C3CCC4CCCCC4(C3CC2)C")
_ANDROSTANE_FUSION_ATOMS = (1, 5, 6, 9, 14, 15)
_ANDROSTANE_C5_ATOM = 9
_ANDROSTANE_NATURAL_CIP = {1: "S", 5: "S", 6: "S", 9: "R", 14: "S", 15: "S"}


def _androstane_5beta_name(mol):
    """If `mol`'s constitution (ignoring stereo) is exactly the androstane
    skeleton, and every one of its six ring-fusion stereocenters is
    specified and matches the natural configuration except C5 (which is
    the opposite, non-natural configuration), return "5-beta-androstane".
    Otherwise return None -- constitution mismatch, a stereocenter left
    unspecified, or any fusion stereocenter other than C5 not matching
    the natural configuration (a different diastereomer entirely, not
    simply "the 5-beta epimer") are all out of scope here, same as the
    module's existing exact-match-only policy for every other case."""
    stripped = Chem.Mol(mol)
    Chem.RemoveStereochemistry(stripped)
    if Chem.MolToSmiles(stripped) != Chem.MolToSmiles(_ANDROSTANE_PLAIN):
        return None
    match = mol.GetSubstructMatch(_ANDROSTANE_PLAIN)
    if not match:
        return None
    rdCIPLabeler.AssignCIPLabels(mol)
    codes = {}
    for plain_idx in _ANDROSTANE_FUSION_ATOMS:
        atom = mol.GetAtomWithIdx(match[plain_idx])
        if not atom.HasProp("_CIPCode"):
            return None
        codes[plain_idx] = atom.GetProp("_CIPCode")
    mismatched = [idx for idx in _ANDROSTANE_FUSION_ATOMS if codes[idx] != _ANDROSTANE_NATURAL_CIP[idx]]
    if mismatched == [_ANDROSTANE_C5_ATOM]:
        return "5-beta-androstane"
    return None


def has_steroid_parent_hydride_name(mol) -> bool:
    if Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME:
        return True
    return _androstane_5beta_name(mol) is not None


def name_steroid_parent_hydride(mol) -> str:
    key = Chem.MolToSmiles(mol)
    if key in _CANONICAL_TO_NAME:
        return _CANONICAL_TO_NAME[key]
    return _androstane_5beta_name(mol)
