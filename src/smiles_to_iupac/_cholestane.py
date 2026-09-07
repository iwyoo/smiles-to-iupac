"""Naming of cholestane, the steroid parent hydride one step up from
cholane, by hardcoded retained-name recognition of its exact unsubstituted
structure:

Cholestane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.4,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document
`_gonane.py`/`_androstane.py`/`_estrane.py`/`_pregnane.py`/`_cholane.py`
defer to, Table 1) as the hydrocarbon with methyl groups at both C-10 and
C-13 (like the rest of this family) plus an eight-carbon C17 side chain,
C20-C27, with a methyl branch at C20 and a gem-dimethyl (isopropyl-style)
terminus at C25 -- i.e. androstane (`_androstane.py`) plus cholane's own
five-carbon C20-branched side chain (`_cholane.py`) extended by three more
carbons ending in two symmetric terminal methyls (C26/C27).

As with the rest of this family, since this module's only job is
recognizing one exact unsubstituted parent, a whole-molecule
canonical-SMILES match is sufficient: any further substituent or
unsaturation changes the canonical SMILES and so is correctly left
unmatched, falling through to `_polycyclic.py`'s general von Baeyer engine.

Like cholane's C20 (which bears a methyl branch and so is a genuine
stereocenter per Table 1), cholestane's C20 is also a stereocenter, but
C25 is not: its two methyl termini (C26/C27) are symmetric, so C25 has no
stereocenter the way `_cholane.py`'s module docstring explains for C20 --
confirmed by PubChem CID 6857534's own IUPACName citing only a single
side-chain stereocenter, "(2R)-6-methylheptan-2-yl" (the "2R" is C20; "6-
methyl" cites C25's substituent but no locant stereodescriptor for C25
itself). `Chem.MolToSmiles` only emits a stereo marker when the input mol
has a stereocenter specified, so a plain, stereo-unspecified input still
canonicalizes to this module's stereo-free key and matches, the same way
it does for the rest of this family; any stereo-specified input (this C20,
or the ring fusions) naturally fails to match and falls through unchanged.

Structure cross-checked: PubChem CID 6857534 ("cholestane")'s
ConnectivitySMILES is 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C'
(MolecularFormula C27H48). Removing that CID's eight-carbon side chain
(the 'CC(C)CCCC(C)' branch, C20-C27) and re-canonicalizing reproduces
`_androstane.py`'s own `_ANDROSTANE_CANONICAL` exactly, confirming the
side chain sits at C17 and the ring system plus both angular methyls are
unchanged from the rest of this family.

Explicitly out of scope: any substituent beyond the plain C20-C27 side
chain and the two angular methyls, any longer/branched side chain, any
unsaturated derivative, and stereo-specified ring-fusion or C20 input --
all fall through to `core.py`'s other dispatch branches unchanged."""

from rdkit import Chem

_CHOLESTANE_SMILES = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C"
_CHOLESTANE_CANONICAL = Chem.CanonSmiles(_CHOLESTANE_SMILES)
_CHOLESTANE_NAME = "cholestane"


def has_cholestane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _CHOLESTANE_CANONICAL


def name_cholestane(mol) -> str:
    return _CHOLESTANE_NAME
