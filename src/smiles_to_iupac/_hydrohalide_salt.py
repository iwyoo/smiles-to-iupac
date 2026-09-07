"""Naming of a two-fragment SMILES where one fragment is a bare hydrogen
halide (HF/HCl/HBr/HI, written as a lone halogen atom) and the other is
any neutral fragment this project can otherwise name, per PubChem's own
naming convention (Blue Book P-77.1.3(3) covers the traditional
'hydrochloride'-style general-nomenclature form, but not this exact
semicolon-joined presentation):

- PubChem PUG REST confirms the pattern is mechanical, not restricted to
  organic bases: 'Nc1ccccc1.Cl' -> 'aniline;hydrochloride',
  'CN.Cl'/'CN.Br'/'CN.I'/'CN.F' -> 'methanamine;hydrochloride'/
  'hydrobromide'/'hydroiodide'/'hydrofluoride', 'CNC.Cl' ->
  'N-methylmethanamine;hydrochloride', and even a non-basic fragment like
  'CCO.Cl' -> 'ethanol;hydrochloride'. The other fragment's own name is
  used unchanged (no '-ium' cation conversion, unlike the Blue Book PIN
  'anilinium chloride') and joined with a literal semicolon.
- The other fragment's name is obtained by recursing into this project's
  own top-level `smiles_to_iupac` (passed in as `namer` to avoid a
  circular import with `core.py`) -- any fragment already supported by
  this library works here for free, and an unsupported one surfaces the
  same `UnsupportedStructure` it always would standalone. Termination is
  guaranteed since the recursed-into fragment always has fewer atoms than
  the original two-fragment molecule.

Explicitly out of scope: more than one halide fragment (e.g. a
dihydrochloride), an ionized SMILES ('[NH3+]'/[Cl-]' rather than the
neutral 'N'/'Cl' atoms), a charged base fragment, and salts of any acid
other than a hydrohalic one (sulfates etc. -- a separate, larger axis)."""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES


def _halide_word(atomic_num):
    return "hydro" + HALOGEN_PREFIXES[atomic_num][:-1] + "ide"


def _is_bare_hydrogen_halide(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return (
        atom.GetAtomicNum() in HALOGEN_PREFIXES
        and atom.GetFormalCharge() == 0
        and atom.GetIsotope() == 0
        and atom.GetTotalNumHs() == 1
    )


def _split_base_and_halide(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 2:
        return None
    for i, halide_frag in enumerate(frags):
        if not _is_bare_hydrogen_halide(halide_frag):
            continue
        (base_frag,) = frags[:i] + frags[i + 1 :]
        if any(atom.GetFormalCharge() != 0 for atom in base_frag.GetAtoms()):
            return None
        return base_frag, halide_frag.GetAtomWithIdx(0).GetAtomicNum()
    return None


def has_hydrohalide_salt_shape(mol) -> bool:
    return _split_base_and_halide(mol) is not None


def name_hydrohalide_salt(mol, namer) -> str:
    base_frag, atomic_num = _split_base_and_halide(mol)
    base_name = namer(Chem.MolToSmiles(base_frag))
    return f"{base_name};{_halide_word(atomic_num)}"
