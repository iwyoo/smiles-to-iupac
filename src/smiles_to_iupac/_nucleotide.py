"""Naming of the 7 retained nucleotide names (P-106.1, Chapter P-10,
https://iupac.qmul.ac.uk/BlueBook/P10.html, `tmp/bluebook/P10.txt` lines
4504-4535): the phosphoric-acid monoester of each of the 7 already-
recognized nucleosides (`_nucleoside.py`), at the ring-primed position
P-106.1 itself prescribes -- 5' for six of them, 3' for xanthosine alone
-- names as "5'-adenylic acid", "5'-guanylic acid", "5'-inosinic acid",
"3'-xanthylic acid", "5'-cytidylic acid", "5'-thymidylic acid", and
"5'-uridylic acid".

Matched by a whole-molecule canonical-SMILES lookup, mirroring
`_nucleoside.py`'s own identical retained-name-by-exact-structure pattern:
each entry here is that module's own nucleoside SMILES with a plain,
unsubstituted -OPO(OH)2 group esterifying the one hydroxyl P-106.1 names
(the 5'-CH2-OH for six of them, the 3'-OH for xanthosine), confirmed
against AMP's real structure (PubChem CID 6083's isomeric SMILES
canonicalizes identically to the 'adenosine' entry below). A different
phosphorylation position, a salt, or any other deviation changes the
canonical SMILES and so is correctly left unmatched, falling through to
whatever generic module (most commonly `_phosphate.py`) would otherwise
apply.

P-106.2 (`tmp/bluebook/P10.txt` lines 4536-4553) names the same 7
nucleosides' unbranched, fully-protonated di-/triphosphate esters as a
phrase instead of a retained '-ylic acid' name -- e.g. "adenosine
5'-(trihydrogen diphosphate)" (ADP, PubChem CID 6022), "uridine
5'-(tetrahydrogen triphosphate)" (worked example, line 4541) -- built
here by chaining one/two more '-O-P(=O)(OH)-' units onto each
monophosphate entry's own terminal phosphate group (the hydrogen count
is just the phosphate-chain length plus one, per the Blue Book's own
worked examples). Any branched, substituted, or non-fully-protonated
phosphate chain, or a chain longer than triphosphate, is out of scope --
still falls through unmatched, same as an unmatched monophosphate above.
"""

from rdkit import Chem

_MONOPHOSPHATE_TAIL = "P(=O)(O)O"
_DIPHOSPHATE_TAIL = "P(=O)(O)OP(=O)(O)O"
_TRIPHOSPHATE_TAIL = "P(=O)(O)OP(=O)(O)OP(=O)(O)O"

_NUCLEOTIDE_SMILES = {
    "5'-adenylic acid": "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N",
    "5'-guanylic acid": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N",
    "5'-inosinic acid": "C1=NC2=C(C(=O)N1)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O",
    "3'-xanthylic acid": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)NC(=O)NC2=O",
    "5'-cytidylic acid": "C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O",
    "5'-thymidylic acid": "CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)COP(=O)(O)O)O",
    "5'-uridylic acid": "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O",
}

# {monophosphate name -> (nucleoside display name, ring-primed position)}
# for P-106.2's own phrase-name construction below.
_DI_TRIPHOSPHATE_BASE = {
    "5'-adenylic acid": ("adenosine", "5'"),
    "5'-guanylic acid": ("guanosine", "5'"),
    "5'-inosinic acid": ("inosine", "5'"),
    "3'-xanthylic acid": ("xanthosine", "3'"),
    "5'-cytidylic acid": ("cytidine", "5'"),
    "5'-thymidylic acid": ("thymidine", "5'"),
    "5'-uridylic acid": ("uridine", "5'"),
}

for _mono_name, _mono_smiles in list(_NUCLEOTIDE_SMILES.items()):
    assert _mono_smiles.count(_MONOPHOSPHATE_TAIL) == 1
    _base_name, _position = _DI_TRIPHOSPHATE_BASE[_mono_name]
    _NUCLEOTIDE_SMILES[f"{_base_name} {_position}-(trihydrogen diphosphate)"] = (
        _mono_smiles.replace(_MONOPHOSPHATE_TAIL, _DIPHOSPHATE_TAIL)
    )
    _NUCLEOTIDE_SMILES[f"{_base_name} {_position}-(tetrahydrogen triphosphate)"] = (
        _mono_smiles.replace(_MONOPHOSPHATE_TAIL, _TRIPHOSPHATE_TAIL)
    )

_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name for name, smiles in _NUCLEOTIDE_SMILES.items()
}
assert len(_CANONICAL_TO_NAME) == len(_NUCLEOTIDE_SMILES), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row"
)


def has_nucleotide_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_nucleotide(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
