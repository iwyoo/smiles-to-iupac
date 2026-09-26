"""Naming of the 7 retained nucleoside names (P-105.1, Chapter P-10,
https://iupac.qmul.ac.uk/BlueBook/P10.html, `tmp/bluebook/P10.txt` lines
4406-4436): "The following names are retained: adenosine, guanosine,
inosine, xanthosine, cytidine, thymidine, uridine."

Matched by a whole-molecule canonical-SMILES lookup, mirroring
`_steroid_parent_hydrides.py`'s/`_fullerene.py`'s/`_inositol.py`'s
identical retained-name-by-exact-structure pattern: a differently
configured stereoisomer, a base attached via O instead of N, or any other
deviation changes the canonical SMILES and so is correctly left
unmatched, falling through to whatever generic module would otherwise
apply. P-105.2 (substituted nucleosides) and P-106 (nucleotides) are
out of scope here.
"""

from rdkit import Chem

_NUCLEOSIDE_SMILES = {
    "adenosine": "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N",
    "guanosine": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=C(NC2=O)N",
    "inosine": "C1=NC2=C(C(=O)N1)N=CN2[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O",
    "xanthosine": "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)NC(=O)NC2=O",
    "cytidine": "C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
    "thymidine": "CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CO)O",
    "uridine": "C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O",
}
_CANONICAL_TO_NAME = {
    Chem.CanonSmiles(smiles): name for name, smiles in _NUCLEOSIDE_SMILES.items()
}
assert len(_CANONICAL_TO_NAME) == len(_NUCLEOSIDE_SMILES), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row"
)


def has_nucleoside_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_nucleoside(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
