"""Naming of ergostane, the steroid parent hydride one step up from
cholestane, by hardcoded retained-name recognition of its exact
unsubstituted structure:

Ergostane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.4,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document
`_gonane.py`/`_androstane.py`/`_estrane.py`/`_pregnane.py`/`_cholane.py`/
`_cholestane.py` defer to, Table 1) as the hydrocarbon with methyl groups
at both C-10 and C-13 (like the rest of this family) plus a nine-carbon
C17 side chain, C20-C28 -- cholestane's own eight-carbon side chain
(`_cholestane.py`) with one extra methyl branch at C24.

As with the rest of this family, since this module's only job is
recognizing one exact unsubstituted parent, a whole-molecule
canonical-SMILES match is sufficient: any further substituent or
unsaturation changes the canonical SMILES and so is correctly left
unmatched, falling through to `_polycyclic.py`'s general von Baeyer engine.

Structure cross-checked: PubChem CID 6857535 ("ergostane")'s
ConnectivitySMILES is
'CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C' (MolecularFormula
C28H50). Removing that CID's nine-carbon side chain (the
'CC(C)C(C)CCC(C)' branch, C20-C28) and re-canonicalizing reproduces
`_androstane.py`'s own `_ANDROSTANE_CANONICAL` exactly, confirming the
side chain sits at C17 and the ring system plus both angular methyls are
unchanged from the rest of this family.

Explicitly out of scope: any substituent beyond the plain C20-C28 side
chain and the two angular methyls, any longer/branched side chain, any
unsaturated derivative, and stereo-specified ring-fusion or side-chain
input (C20/C24 are both genuine stereocenters per Table 1, but -- same as
the rest of this family -- a plain stereo-unspecified input still matches
this module's stereo-free canonical key; a stereo-specified one falls
through unchanged) -- all fall through to `core.py`'s other dispatch
branches unchanged."""

from rdkit import Chem

_ERGOSTANE_SMILES = "CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C"
_ERGOSTANE_CANONICAL = Chem.CanonSmiles(_ERGOSTANE_SMILES)
_ERGOSTANE_NAME = "ergostane"


def has_ergostane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _ERGOSTANE_CANONICAL


def name_ergostane(mol) -> str:
    return _ERGOSTANE_NAME
