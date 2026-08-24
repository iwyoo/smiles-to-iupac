"""Naming of saturated monocyclic rings (3- to 7-membered) containing
exactly one heteroatom (O, S, or N) and no substituents, using
Hantzsch-Widman-system names, per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-22.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf,
  Table 22.1): the Hantzsch-Widman stem for a saturated ring is
  '-irane'/'-irene' (3), '-etane' (4), '-olane' (5), '-ane' (6), '-epane'
  (7) for O and S ('oxa'/'thia' + stem), giving oxirane, oxetane,
  oxolane, oxane, oxepane and thiirane, thietane, thiolane, thiane,
  thiepane.
- For N, the saturated-ring stem is '-iridine' (3), '-etidine' (4),
  '-olidine' (5), '-inane' (6), '-epane' (7) ('aza' + stem) -- but four of
  these five (all but 7-membered azepane) are retained names that are
  themselves the preferred IUPAC name in place of the literal 'aza' +
  stem form: aziridine, azetidine, pyrrolidine, and piperidine (not
  azirane/azetane/azolidine/azinane) -- confirmed via IUPAC's P-22.2.1
  Table 2.3 listing piperidine/pyrrolidine (among others) as retained
  names that are PINs. Aziridine/azetidine happen to coincide with the
  literal stem-based construction; pyrrolidine/piperidine don't.
- Since this module's only job is recognizing the exact unsubstituted
  parent for a fixed, small (element, ring size) table -- no locants to
  assign, no substituent numbering -- an exact whole-molecule canonical-
  SMILES match against each name's structure is both sufficient and
  simplest, mirroring `_peri_fused_aromatic.py`'s approach.

Formulas cross-checked (all well-known compounds): oxirane C2H4O,
piperidine C5H11N, thiane C5H10S, etc. -- see `_RETAINED_NAME_SMILES`.

Explicitly out of scope: any substituent, any unsaturation (P-22's
indicated-hydrogen/partial-saturation forms are a separate future task),
two or more heteroatoms (dioxane, morpholine, etc.), heteroatoms other
than O/S/N, and ring sizes outside 3-7. `has_hetero_monocyclic_name`
returns False for all of these, so `core.py`'s existing dispatch (which
already rejects heteroatoms outside a few specific recognized shapes)
continues to raise `UnsupportedStructure` for them, unchanged.
"""

from rdkit import Chem

_RETAINED_NAME_SMILES = {
    ("O", 3): ("oxirane", "C1CO1"),
    ("O", 4): ("oxetane", "C1CCO1"),
    ("O", 5): ("oxolane", "C1CCCO1"),
    ("O", 6): ("oxane", "C1CCCCO1"),
    ("O", 7): ("oxepane", "C1CCCCCO1"),
    ("S", 3): ("thiirane", "C1CS1"),
    ("S", 4): ("thietane", "C1CCS1"),
    ("S", 5): ("thiolane", "C1CCCS1"),
    ("S", 6): ("thiane", "C1CCCCS1"),
    ("S", 7): ("thiepane", "C1CCCCCS1"),
    ("N", 3): ("aziridine", "C1CN1"),
    ("N", 4): ("azetidine", "C1CCN1"),
    ("N", 5): ("pyrrolidine", "C1CCCN1"),
    ("N", 6): ("piperidine", "C1CCCCN1"),
    ("N", 7): ("azepane", "C1CCCCCN1"),
}
_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name for name, smiles in _RETAINED_NAME_SMILES.values()
}


def has_hetero_monocyclic_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_hetero_monocyclic(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
