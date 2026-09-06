"""Naming of pregnane, the steroid parent hydride one step up from
androstane, by hardcoded retained-name recognition of its exact
unsubstituted structure:

Pregnane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.4,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document
`_gonane.py`/`_androstane.py`/`_estrane.py` defer to, Table 1) as the
hydrocarbon with methyl groups at both C-10 and C-13 (like androstane) plus
a plain ethyl side chain at C-17 -- i.e. androstane (`_androstane.py`) plus
a two-carbon C17 side chain (C20/C21), unlike androstane which has no C17
side chain at all.

As with `_gonane.py`/`_androstane.py`/`_estrane.py`, since this module's
only job is recognizing one exact unsubstituted parent, a whole-molecule
canonical-SMILES match is sufficient: any further substituent or
unsaturation changes the canonical SMILES and so is correctly left
unmatched, falling through to `_polycyclic.py`'s general von Baeyer engine.

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the
input mol, so any stereo-specified input (ring-fusion or the C17 side
chain's own C20, which per Table 1 is only a real stereocenter for longer
side chains than pregnane's plain ethyl -- pregnane's C20 bears two
hydrogens, so it is not itself a stereocenter) naturally fails to match
this module's stereo-free canonical key and falls through the same way.

Structure cross-checked: PubChem CID 439513 ("pregnane")'s
ConnectivitySMILES is 'CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C'
(MolecularFormula C21H36, IUPACName
'...-17-ethyl-10,13-dimethyl-...-cyclopenta[a]phenanthrene' -- not
'pregnane' itself, like the other modules in this family, only usable to
verify the structure). Removing that CID's terminal two-carbon (ethyl)
side chain and re-canonicalizing reproduces `_androstane.py`'s own
`_ANDROSTANE_CANONICAL` exactly, confirming the side chain sits at C17 and
the ring system plus both angular methyls are unchanged from androstane.

Explicitly out of scope: any substituent beyond the plain ethyl side
chain and the two angular methyls, any longer/branched/stereo-specified
side chain (cholane and up, which do have a real C20 stereocenter -- Table
1), any unsaturated derivative, and stereo-specified ring-fusion input --
all fall through to `core.py`'s other dispatch branches unchanged."""

from rdkit import Chem

_PREGNANE_SMILES = "CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C"
_PREGNANE_CANONICAL = Chem.CanonSmiles(_PREGNANE_SMILES)
_PREGNANE_NAME = "pregnane"


def has_pregnane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _PREGNANE_CANONICAL


def name_pregnane(mol) -> str:
    return _PREGNANE_NAME
