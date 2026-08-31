"""Naming of saturated and mancude (maximally unsaturated) monocyclic
rings containing exactly one heteroatom (O, S, Se, Te, or N) and no
substituents, using Hantzsch-Widman-system and retained names, per the
IUPAC 2013 Recommendations ("the Blue Book"):

Saturated rings (3- to 7-membered):
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

Mancude (aromatic) rings, single heteroatom, 5- and 6-membered only
(P-22.2.1 Table 2.2): furan/thiophene/selenophene/tellurophene (5-membered
O/S/Se/Te) and pyridine (6-membered N) are retained names that are PINs
outright; pyrrole (5-membered N) needs an indicated-hydrogen prefix --
its PIN is '1H-pyrrole', not bare 'pyrrole' (confirmed via Table 2.2 and
PubChem). Two-heteroatom mancude rings (imidazole, oxazole, pyrazine,
etc., also listed in Table 2.2) and the 6-membered O/S/Se/Te rings
(pyran/thiopyran/selenopyran/telluropyran, which need an indicated-
hydrogen prefix themselves since they aren't fully mancude with a single
chalcogen) are out of scope -- separate future tasks.

Since this module's only job is recognizing the exact unsubstituted
parent for a fixed, small (element, ring size, saturation) table -- no
locants to assign, no substituent numbering -- an exact whole-molecule
canonical-SMILES match against each name's structure is both sufficient
and simplest, mirroring `_peri_fused_aromatic.py`'s approach.

Formulas cross-checked (all well-known compounds): oxirane C2H4O,
piperidine C5H11N, thiane C5H10S, furan C4H4O, pyridine C5H5N, etc. --
see `_RETAINED_NAME_SMILES` and `_MANCUDE_NAME_SMILES`.

Explicitly out of scope: any substituent, partially-saturated indicated-
hydrogen forms other than 1H-pyrrole, two or more heteroatoms (dioxane,
morpholine, imidazole, etc.), heteroatoms other than O/S/Se/Te/N, and
ring sizes outside the tables above. `has_hetero_monocyclic_name` returns
False for all of these, so `core.py`'s existing dispatch (which already
rejects heteroatoms outside a few specific recognized shapes) continues
to raise `UnsupportedStructure` for them, unchanged.
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
_MANCUDE_NAME_SMILES = {
    ("O", 5): ("furan", "c1ccoc1"),
    ("S", 5): ("thiophene", "c1ccsc1"),
    ("Se", 5): ("selenophene", "c1cc[se]c1"),
    ("Te", 5): ("tellurophene", "c1cc[te]c1"),
    ("N", 5): ("1H-pyrrole", "c1cc[nH]c1"),
    ("N", 6): ("pyridine", "c1ccncc1"),
}
_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name
    for name, smiles in (*_RETAINED_NAME_SMILES.values(), *_MANCUDE_NAME_SMILES.values())
}


def has_hetero_monocyclic_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_hetero_monocyclic(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
