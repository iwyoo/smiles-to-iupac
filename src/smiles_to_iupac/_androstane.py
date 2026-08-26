"""Naming of androstane, the steroid parent hydride one step up from gonane,
by hardcoded retained-name recognition of its exact unsubstituted structure:

Androstane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.3,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document `_gonane.py`
defers to) as "the hydrocarbon with methyl groups at C-10 and C-13 but
without a side chain at C-17" -- i.e. gonane (`_gonane.py`) plus the two
angular methyls (C18 on C13, C19 on C10) that most real steroids (androgens
etc.) actually have, still with no C17 side chain (present from pregnane
up).

As with `_gonane.py`, since this module's only job is recognizing one exact
unsubstituted parent, a whole-molecule canonical-SMILES match is sufficient:
any further substituent, side chain, or unsaturation changes the canonical
SMILES and so is correctly left unmatched, falling through to
`_polycyclic.py`'s general von Baeyer engine.

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the input
mol, so any stereo-specified ring-fusion/angular-methyl input naturally
fails to match this module's stereo-free canonical key and falls through
the same way -- exactly like `_gonane.py`, ring-fusion stereochemistry
(3S-2.2) is a separate rule, out of scope here.

Structure cross-checked: the test SMILES below, canonicalized, is identical
to PubChem CID 6857536's ConnectivitySMILES
('CC12CCCC1C3CCC4CCCCC4(C3CC2)C') canonicalized the same way (PubChem's own
computed IUPACName for that CID is a von Baeyer/fused-ring name,
'...-10,13-dimethyl-...-tetradecahydro-1H-cyclopenta[a]phenanthrene', not
'androstane' -- like `_gonane.py`, not usable for verifying the name
itself, only the structure). Removing CID 6857536's two degree-1 (methyl)
carbons and re-canonicalizing the remaining skeleton reproduces
`_gonane.py`'s own `_GONANE_CANONICAL` exactly, confirming the two methyls
sit at the correct angular positions (C10/C13) and not some other pair of
ring-fusion carbons.

Explicitly out of scope: any substituent beyond the two angular methyls, any
C17 side chain (pregnane and up), estrane (only one angular methyl, C19,
with an aromatic A ring -- a structurally distinct parent, not a substituent
variant of androstane), any unsaturated derivative, and stereo-specified
input -- all fall through to `core.py`'s other dispatch branches unchanged.
"""

from rdkit import Chem

_ANDROSTANE_SMILES = "CC12CCCC1C3CCC4CCCCC4(C3CC2)C"
_ANDROSTANE_CANONICAL = Chem.CanonSmiles(_ANDROSTANE_SMILES)
_ANDROSTANE_NAME = "androstane"


def has_androstane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _ANDROSTANE_CANONICAL


def name_androstane(mol) -> str:
    return _ANDROSTANE_NAME
