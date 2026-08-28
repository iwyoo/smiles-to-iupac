"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
the carbocyclic (non-heteroatom) ring of indole or 1-benzofuran -- a base
component that is itself already a fused bicyclic retained-name system,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-25.3.1.3 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf,
  full text confirmed accessible 2026-08-29 via
  https://iupac.qmul.ac.uk/BlueBook/P2.html after several earlier sessions'
  access failures -- see tasks/general-fusion-naming.md's own log):
  "Isomers are distinguished by lettering, continuously, each peripheral
  side of the parent component ... using the italic letters a, b, c, etc.,
  beginning with a for the side numbered '1,2', b for '2,3' etc." This
  module actually walks indole's/1-benzofuran's own fixed retained-name
  numbering (N1/O1, C2, C3, C3a, C4, C5, C6, C7, C7a -- all periphery
  atoms, no interior atom, since it's a plain two-ring ortho-fusion) to
  compute the letter for each of the three periphery bonds this module
  covers (C4-C5 = 'e', C5-C6 = 'f', C6-C7 = 'g'), rather than a hardcoded
  per-compound lookup table like `_heteroaromatic_fused.py`'s -- this is
  the general algorithm this project's own roadmap had been blocked on
  for lack of primary-source access to this exact section.
- Confirmed via PubChem (structure, not just autoname -- see caveat
  below): CID 98617 (`C1=CC=C2C(=C1)C=CC3=C2NC=C3`) -> '1H-benzo[g]indole
  (PIN)' (fusion bond C6-C7, letter 'g'; the indole base needs its own
  indicated hydrogen '1H-' regardless of where a benzo ring fuses onto
  it, the same P-25.7.1.3 requirement `_heteroaromatic_fused.py` already
  documents for bare indole), CID 9192
  (`C1=CC=C2C(=C1)C=CC3=C2C=CO3`) -> 'benzo[e][1]benzofuran' (fusion bond
  C4-C5, letter 'e'), CID 67475 (`C1=CC=C2C(=C1)C=CC3=C2OC=C3`) ->
  'benzo[g][1]benzofuran' (fusion bond C6-C7, letter 'g' -- same letter
  as the indole case above, a useful independent cross-check that the
  algorithm's periphery walk is base-hetero-atom-agnostic). Both letters
  found this way (via `GetSubstructMatches` locating which two base-role
  atoms are adjacent to the four new ring atoms) match the confirmed
  names exactly. CAUTION: an earlier, different PubChem CID (170319311,
  bare "benzo[g]indole", no '1H-') was considered in an earlier session
  pass and is a malformed database entry (molecular formula C12H7N, two
  hydrogens short of the real C12H9N structure, with a broken non-
  aromatic N=C=C pyrrole ring) -- do not reuse that CID; CID 98617 is the
  correct, properly aromatic reference.
- When 1-benzofuran (not indole) is the base component within a larger
  fusion name, its retained name must be qualified as '[1]benzofuran'
  (not bare 'benzofuran') to distinguish it from the '2-benzofuran'
  ("isobenzofuran") isomer -- confirmed directly by both benzofuran-based
  PubChem names above, which both carry the '[1]' qualifier that bare,
  unfused '1-benzofuran' (`_heteroaromatic_fused.py`) does not need on
  its own.
- The middle bond of the same periphery stretch, C5-C6 ('f'), is included
  by direct analogy (same mechanism, no atom-identity dependence) but has
  no PubChem-confirmed compound of its own -- a reviewed, not
  independently verified, extension.

Scope, deliberately narrow (see tasks/general-fusion-naming.md,
tasks/polycyclic-component-fusion-naming.md): exactly one plain
unsubstituted benzo ring ortho-fused onto indole's or 1-benzofuran's own
6-membered (carbocyclic) ring, at one of the three periphery bonds not
touching a ring-fusion atom or the base's own heteroatom (C4-C5, C5-C6,
C6-C7 -- letters 'e', 'f', 'g'). Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any other periphery bond, in particular one touching a
  fusion atom (C3a, C7a -- an 'ortho- and peri-fused' system needing a
  different citation mechanism per P-25.3.1.3) or the base's own
  heteroatom (a different, unverified shape).
- Any substituent, any heteroatom in the new ring, any non-benzo attached
  component, or more than one extra ring.
- Thiophene/selenophene/tellurophene analogues of benzofuran (a separate,
  unverified extension).
"""

from rdkit import Chem

from ._common import UnsupportedStructure

# Atom-role -> index in each reference SMILES (both share the same
# topology, indole's [nH]/benzofuran's o sit at role '<het>1'; verified
# directly via RDKit bond inspection, see module docstring).
_INDOLE_REF = Chem.MolFromSmiles("c1ccc2[nH]ccc2c1")
_BENZOFURAN_REF = Chem.MolFromSmiles("c1ccc2occc2c1")
_ROLE_TO_IDX = {"C5": 0, "C6": 1, "C7": 2, "C7a": 3, "1": 4, "C2": 5, "C3": 6, "C3a": 7, "C4": 8}

# P-25.3.1.3: periphery bonds lettered continuously from the '1,2' side --
# only the three bonds this module supports (none touching a fusion atom
# or the heteroatom) are listed.
_SUPPORTED_FUSION_BONDS = {
    frozenset(("C4", "C5")): "e",
    frozenset(("C5", "C6")): "f",
    frozenset(("C6", "C7")): "g",
}

_BASES = {
    "indole": (_INDOLE_REF, "1H-benzo[{letter}]indole"),
    "benzofuran": (_BENZOFURAN_REF, "benzo[{letter}][1]benzofuran"),
}


def _find_fusion_letter(mol, reference):
    matches = mol.GetSubstructMatches(reference, useChirality=False)
    if len(matches) != 1:
        return None
    (match,) = matches
    idx_to_role = {v: k for k, v in _ROLE_TO_IDX.items()}
    target_to_role = {match[i]: idx_to_role[i] for i in range(len(match))}
    core_atoms = set(match)
    extra_atoms = set(range(mol.GetNumAtoms())) - core_atoms
    if len(extra_atoms) != 4:
        return None
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(a).GetIsAromatic() for a in extra_atoms):
        return None

    fusion_roles = set()
    for a in extra_atoms:
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            if n.GetIdx() in target_to_role:
                fusion_roles.add(target_to_role[n.GetIdx()])
            elif n.GetIdx() not in extra_atoms:
                return None
    if len(fusion_roles) != 2:
        return None
    return _SUPPORTED_FUSION_BONDS.get(frozenset(fusion_roles))


def _find_core(mol):
    if mol.GetNumAtoms() != 13:
        return None
    if mol.GetRingInfo().NumRings() != 3:
        return None
    if any(not atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return None
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if atom.GetAtomicNum() not in (6, 7, 8):
            return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    for base_name, (reference, template) in _BASES.items():
        heteroatom_num = 7 if base_name == "indole" else 8
        heteroatoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == heteroatom_num]
        others = [a for a in mol.GetAtoms() if a.GetAtomicNum() not in (6, heteroatom_num)]
        if len(heteroatoms) != 1 or others:
            continue
        letter = _find_fusion_letter(mol, reference)
        if letter is not None:
            return template.format(letter=letter)
    return None


def has_polycyclic_component_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_polycyclic_component_fusion(mol) -> str:
    name = _find_core(mol)
    if name is None:
        raise UnsupportedStructure(
            "this tricyclic system is not a supported benzo-fused "
            "indole/1-benzofuran shape (see P-25.3.1.3)"
        )
    return name
