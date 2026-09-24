"""Naming of a hydrate adduct -- one connected organic-compound fragment
plus one or more separate water molecules in the same SMILES -- per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-14.8 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf,
  ~1161-1192, "Adducts"): components are cited in order of P-41
  seniority, joined by an em dash ('—'), with the proportion cited
  afterward as '(m/n)' in parentheses -- confirmed directly from the
  primary source's own worked example, `ethanol—pyridine (1/1) (PIN)`.
  P-14.8's own text explicitly extends this to hydrates: "Solvates,
  including hydrates..., are treated as adducts and preferred IUPAC names
  must use the notation for denoting the proportion of components
  described above." Water is always the lowest-seniority component in a
  hydrate, so the format simplifies to a fixed pattern: `<compound
  name>—water (1/n)` for n waters -- n=1 is always cited, never
  omitted, per the worked example's own `(1/1)` for a 1:1 adduct.
- PubChem's own auto-generated name for a hydrate (e.g. oxalic acid
  dihydrate, CID 61373 -> "oxalic acid;dihydrate") uses a different,
  non-PIN semicolon-joined convention, not this module's em-dash format
  -- not used for verification here (mirrors a pattern already seen
  elsewhere in this project, e.g. `[60]fullerene`'s own PIN vs PubChem's
  non-PIN algorithmic name); verification instead relies on the primary
  source's own literal worked examples.
- The non-water fragment's own name is obtained by recursing into this
  project's own top-level `smiles_to_iupac` (passed in as `namer` to
  avoid a circular import with `core.py`), mirroring
  `_hydrohalide_salt.py`'s identical mechanism.

Scope: exactly one connected non-water fragment (any organic molecule
this project can already name standalone) plus one or more separate `O`
(water) fragments, each a single, uncharged, isotopically-unmodified
oxygen atom with 2 implicit hydrogens and no bonds.

Explicitly out of scope: more than one non-water fragment, a non-water
solvate (P-14.8.1's general adduct case, e.g. `ethanol—pyridine`), any
mixed organic-inorganic adduct (P-14.8.2), a charged fragment (already
`_salt.py`'s territory).
"""

from rdkit import Chem


def _is_bare_water(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return (
        atom.GetAtomicNum() == 8
        and atom.GetFormalCharge() == 0
        and atom.GetIsotope() == 0
        and atom.GetTotalNumHs() == 2
    )


def _split_compound_and_waters(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    water_frags = [frag for frag in frags if _is_bare_water(frag)]
    if not water_frags:
        return None
    non_water_frags = [frag for frag in frags if not _is_bare_water(frag)]
    if len(non_water_frags) != 1:
        return None
    (compound_frag,) = non_water_frags
    if any(atom.GetFormalCharge() != 0 for atom in compound_frag.GetAtoms()):
        return None
    return compound_frag, len(water_frags)


def has_hydrate_adduct_shape(mol) -> bool:
    return _split_compound_and_waters(mol) is not None


def name_hydrate_adduct(mol, namer) -> str:
    compound_frag, water_count = _split_compound_and_waters(mol)
    compound_name = namer(Chem.MolToSmiles(compound_frag))
    return f"{compound_name}—water (1/{water_count})"
