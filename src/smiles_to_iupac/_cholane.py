"""Naming of cholane, the steroid parent hydride one step up from pregnane,
by hardcoded retained-name recognition of its exact unsubstituted structure:

Cholane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.4,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document
`_gonane.py`/`_androstane.py`/`_estrane.py`/`_pregnane.py` defer to, Table 1)
as the hydrocarbon with methyl groups at both C-10 and C-13 (like androstane/
pregnane) plus a five-carbon C17 side chain, C20-C24, with a methyl branch at
C20 -- i.e. androstane (`_androstane.py`) plus a plain 2-methylpentan... side
chain (the same isohexyl-style branch also found at C17 of cholesterol and
related sterols), unlike pregnane's plain two-carbon ethyl side chain.

As with `_gonane.py`/`_androstane.py`/`_estrane.py`/`_pregnane.py`, since
this module's only job is recognizing one exact unsubstituted parent, a
whole-molecule canonical-SMILES match is sufficient: any further substituent
or unsaturation changes the canonical SMILES and so is correctly left
unmatched, falling through to `_polycyclic.py`'s general von Baeyer engine.

Unlike pregnane's C20 (which bears two hydrogens and so is not itself a
stereocenter), cholane's C20 bears a methyl branch and so, per Table 1, is a
genuine stereocenter (cholane's own definition leaves it unspecified,
covering both epimers per Rule 3S-2.4's own convention for the parent
hydride name). `Chem.MolToSmiles` only emits a stereo marker (@/@@) when the
input mol has that stereocenter specified, so a plain, stereo-unspecified
input still canonicalizes to this module's stereo-free key and matches, the
same way ring-fusion-unspecified input does for the other modules in this
family; any stereo-specified input (this C20, or the ring fusions) naturally
fails to match and falls through unchanged.

Structure cross-checked: PubChem CID 6857459 ("cholane")'s
ConnectivitySMILES is 'CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C' (MolecularFormula
C24H42, IUPACName '...-17-[(2R)-pentan-2-yl]-...-cyclopenta[a]phenanthrene' --
not 'cholane' itself, like the other modules in this family, only usable to
verify the structure). Removing that CID's five-carbon side chain (the
'CCCC(C)' branch, C20-C24 plus the C20 methyl) and re-canonicalizing
reproduces `_androstane.py`'s own `_ANDROSTANE_CANONICAL` exactly, confirming
the side chain sits at C17 and the ring system plus both angular methyls are
unchanged from androstane/pregnane.

Explicitly out of scope: any substituent beyond the plain C20-C24 side chain
and the two angular methyls, any longer/branched side chain (cholestane and
up), any unsaturated derivative, and stereo-specified ring-fusion or C20
input -- all fall through to `core.py`'s other dispatch branches unchanged."""

from rdkit import Chem

_CHOLANE_SMILES = "CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C"
_CHOLANE_CANONICAL = Chem.CanonSmiles(_CHOLANE_SMILES)
_CHOLANE_NAME = "cholane"


def has_cholane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _CHOLANE_CANONICAL


def name_cholane(mol) -> str:
    return _CHOLANE_NAME
