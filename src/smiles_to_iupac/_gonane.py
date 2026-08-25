"""Naming of gonane, the fundamental steroid parent ring hydride, by
hardcoded retained-name recognition of its exact unsubstituted structure:

Gonane is defined by the 1989 IUPAC steroid nomenclature (Rule 2.1,
https://iupac.qmul.ac.uk/steroid/3S02a.html, the document Blue Book P-6
defers to for steroid parent names) as "the parent tetracyclic hydrocarbon
without methyl groups at C-10 and C-13 and without a side chain at C-17" --
i.e. the plain perhydrocyclopenta[a]phenanthrene skeleton (three fused
six-membered rings plus one five-membered ring, angularly ortho-fused,
C17H28) with none of the angular methyls (C18/C19, present from androstane
up) or the C17 side chain (present from pregnane up).

As with `_cyclophane.py`/`_fullerene.py`/`_peri_fused_aromatic.py`, since
this module's only job is recognizing one exact unsubstituted parent, a
whole-molecule canonical-SMILES match is sufficient: any substituent,
methyl, side chain, or unsaturation changes the canonical SMILES and so is
correctly left unmatched, falling through to `_polycyclic.py`'s general
von Baeyer engine (which already names this exact skeleton as
"tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane" if this module doesn't
intercept it first -- confirmed by direct testing).

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the
input mol, so any stereo-specified ring-fusion input naturally fails to
match this module's stereo-free canonical key and falls through the same
way -- gonane's own definition (Rule 2.1) does not specify ring-fusion
stereochemistry (that's a separate rule, 3S-2.2), so stereo-specified
steroid naming is out of scope here.

Structure cross-checked: the test SMILES below, canonicalized, is
identical to PubChem CID 6857523's ConnectivitySMILES
(`C1CCC2C(C1)CCC3C2CCC4C3CCC4`) canonicalized the same way (PubChem's own
computed IUPACName for that CID is a von Baeyer/fused-ring name,
"...hexadecahydro-1H-cyclopenta[a]phenanthrene", not "gonane" -- like the
phane modules, it isn't usable for verifying the name itself, only the
structure).

Explicitly out of scope: any substituent, any angular methyl (androstane
and up), any C17 side chain (pregnane and up), any unsaturated derivative
(e.g. estrane's aromatic A ring), and stereo-specified input -- all fall
through to `core.py`'s other dispatch branches unchanged.
"""

from rdkit import Chem

_GONANE_SMILES = "C1CCCC2CCC3C(C12)CCC4C3CCC4"
_GONANE_CANONICAL = Chem.CanonSmiles(_GONANE_SMILES)
_GONANE_NAME = "gonane"


def has_gonane_name(mol) -> bool:
    return Chem.MolToSmiles(mol) == _GONANE_CANONICAL


def name_gonane(mol) -> str:
    return _GONANE_NAME
