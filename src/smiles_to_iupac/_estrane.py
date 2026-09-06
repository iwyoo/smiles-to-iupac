"""Naming of estrane, the steroid parent hydride with only the C-13 angular
methyl, by hardcoded retained-name recognition of its exact unsubstituted
structure:

Estrane is defined by the 1989 IUPAC steroid nomenclature (Rule 3S-2.2,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the same document `_gonane.py`/
`_androstane.py` defer to) as "the hydrocarbon with a methyl group at C-13
but without a methyl group at C-10 and without a side chain at C-17" --
i.e. gonane (`_gonane.py`) plus only the C13 angular methyl (C18), unlike
androstane (`_androstane.py`) which has both the C13 and C10 (C19) methyls.

As with `_gonane.py`/`_androstane.py`, since this module's only job is
recognizing one exact unsubstituted parent, a whole-molecule canonical-
SMILES match is sufficient: any further substituent, side chain, or
unsaturation changes the canonical SMILES and so is correctly left
unmatched, falling through to `_polycyclic.py`'s general von Baeyer engine.

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the input
mol, so any stereo-specified ring-fusion/angular-methyl input naturally
fails to match this module's stereo-free canonical key and falls through
the same way -- exactly like `_gonane.py`/`_androstane.py`, ring-fusion
stereochemistry (3S-2.2's own diagram) is a separate rule, out of scope
here.

Structure cross-checked: the test SMILES below, canonicalized, is identical
to PubChem CID 5460658's ("estrane") IsomericSMILES with stereo markers
stripped, canonicalized the same way. Independently re-derived by removing
`_androstane.py`'s C10 methyl (the angular methyl attached to the
six-six-ring-fusion atom, as opposed to C13's five-six-ring-fusion atom)
from its canonical structure and re-canonicalizing -- both routes produce
the identical canonical SMILES, confirming the remaining methyl sits at the
real C13 position and not some other ring-fusion atom.

Explicitly out of scope: any substituent beyond the one angular methyl, any
C17 side chain (pregnane and up), androstane (has the extra C10 methyl --
a structurally distinct parent, not a substituent variant of this one), any
unsaturated derivative, and stereo-specified input -- all fall through to
`core.py`'s other dispatch branches unchanged.
"""

from rdkit import Chem

_ESTRANE_SMILES = "CC12CCCC1C1CCC3CCCCC3C1CC2"
_ESTRANE_CANONICAL = Chem.CanonSmiles(_ESTRANE_SMILES)
_ESTRANE_NAME = "estrane"


def has_estrane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _ESTRANE_CANONICAL


def name_estrane(mol) -> str:
    return _ESTRANE_NAME
