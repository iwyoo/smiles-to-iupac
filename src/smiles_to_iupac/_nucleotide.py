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
phosphorylation position, a diphosphate/triphosphate (P-106.2), a salt, or
any other deviation changes the canonical SMILES and so is correctly left
unmatched, falling through to whatever generic module (most commonly
`_phosphate.py`) would otherwise apply.
"""

from rdkit import Chem

_NUCLEOTIDE_SMILES = {
    "5'-adenylic acid": "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N",
    "5'-guanylic acid": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N=C(NC2=O)N",
    "5'-inosinic acid": "C1=NC2=C(C(=O)N1)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O",
    "3'-xanthylic acid": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)NC(=O)NC2=O",
    "5'-cytidylic acid": "C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O",
    "5'-thymidylic acid": "CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)COP(=O)(O)O)O",
    "5'-uridylic acid": "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O",
}
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
