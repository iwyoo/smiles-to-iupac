"""Naming of simple hydrazines (H2N-NH2 bearing 0-2 unbranched, saturated
alkyl substituents per nitrogen), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-68.3.1.2.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  "Hydrazine is a retained name describing the structure H2N-NH2; it is a
  preselected name... As a parent hydride, hydrazine is numbered by
  numerical locants, 1 and 2... Hydrazine is substituted by hydrocarbyl
  groups and characteristic groups expressed by suffixes and prefixes."
  Unlike `_silane_chain.py`'s variable-length Si-Si chain, hydrazine is a
  single, fixed two-nitrogen parent hydride -- structurally the same
  "skeleton itself is the parent, alkyl groups are substituent prefixes"
  shape as `_phosphane.py`/`_diazene.py`, just with a single N-N bond (so
  each nitrogen can carry up to 2 substituents, unlike diazene's N=N,
  which leaves only 1 free valence per nitrogen). Confirmed worked
  example: '1,1-dimethylhydrazine (PIN)'. Further confirmed via PubChem
  PUG REST, which agrees with the PIN in every case checked (an unusually
  complete match for this project): CID 713 (`NN`) -> "hydrazine", CID
  ... (`CNN`) -> "methylhydrazine", (`CN(C)N`) -> "1,1-dimethylhydrazine",
  (`CNNC`) -> "1,2-dimethylhydrazine", (`CN(C)NC`) ->
  "1,1,2-trimethylhydrazine", (`CN(C)N(C)C`) -> "1,1,2,2-tetramethylhydrazine",
  (`CCN(CC)N`) -> "1,1-diethylhydrazine".
- P-35.2.1: a halogen substituent on a carbon branch (not directly on a
  hydrazine nitrogen) coexists freely, reusing `_amide.py`/`_hydrazide.py`'s
  own `name_branch`/`halogen_substituents` pattern -- a straight
  (unbranched, once halogens are set aside) carbon chain bearing one or
  more halogens is named as a compound substituent the same way, e.g.
  '2-chloroethyl'. Confirmed via PubChem PUG REST: CID 19348410
  (`ClCCNN`) -> "2-chloroethylhydrazine", CID 165881
  (`ClCCN(N)CCCl`) -> "1,1-bis(2-chloroethyl)hydrazine" (two identical
  compound substituents combine with the ordinary 'bis' multiplying
  prefix, same mechanism `format_substituent_prefixes` already uses
  elsewhere for a compound substituent), CID 140409285 (`ClCNN`) ->
  "chloromethylhydrazine" (the sole-substituent case omits its own
  hydrazine locant regardless of the substituent itself being a compound
  name, same as the plain-alkyl sole-substituent case below). A real
  carbon branch (e.g. isopropyl) is supported too, built with
  `name_branch` (P-29 PIN style, fixed project-wide by PR #237; mirrors
  `_diazene.py`'s identical fix, PR #338) -- confirmed via PubChem:
  CID 52789 (`CC(C)NN`) -> "propan-2-ylhydrazine", CID 18459757
  (`CC(C)NNC`) -> "1-methyl-2-propan-2-ylhydrazine". **Not** independently
  re-verified: two identical branched substituents on the same nitrogen
  (e.g. `CC(C)N(N)C(C)C`) -- PubChem gives "1,1-di(propan-2-yl)hydrazine"
  (CID 70201) for that shape, but `format_substituent_prefixes`'s
  existing `compound=True` path (already confirmed correct for the
  chloroethyl case, 'bis(2-chloroethyl)') would instead produce
  '1,1-bis(propan-2-yl)hydrazine' -- the same di-vs-bis PubChem
  inconsistency `_diazene.py` already ran into (see that module's own
  docstring) and deliberately left unresolved; not touched here either.
- A single substituent (on either nitrogen) needs no locant -- hydrazine's
  two nitrogens are interchangeable by symmetry when only one substituent
  is present, the same "genuine absence of ambiguity" `_diazene.py`
  documents for its own single-substituent case (e.g. 'methylhydrazine',
  not '1-methylhydrazine').
- Two or more substituents do need locants (unlike diazene, which never
  needs them at all): with up to 2 substituents per nitrogen, '1,1-' and
  '1,2-' are genuinely different structures. The locant-assignment logic
  here is new (there are only two candidate numbering directions -- N1 as
  locant 1 or N2 as locant 1 -- so the existing lowest-locant-set/citation-
  order tie-break pattern used elsewhere in this project is reused, just
  simplified to two candidates instead of a general search); substituent
  grouping/formatting reuses `_substituents.py`'s ordinary
  `format_substituent_prefixes` (ordinary locant-carrying prefixes, not
  `_phosphane.py`'s special no-locant P-16.5.1.3.1 parenthesization rule,
  which only applies when locants are omitted entirely).
- A plain, unsubstituted benzene ring bonded directly to a hydrazine
  nitrogen is cited as a 'phenyl' substituent -- `name_branch`'s own
  `aromatic_atoms` parameter already recognizes this shape (mirrors
  `_diazene.py`'s identical fix, PR #401). Confirmed via PubChem PUG
  REST: `c1ccccc1NN` -> "phenylhydrazine" (CID 7516),
  `c1ccccc1NNc1ccccc1` -> "1,2-diphenylhydrazine" (CID 31222),
  `c1ccccc1NNC` -> "1-methyl-2-phenylhydrazine" (CID 10197813).

Explicitly out of scope (raise `UnsupportedStructure`):
- Any atom other than the two hydrazine nitrogens, carbon, hydrogen, and
  a halogen substituent on a carbon branch.
- A halogen bonded directly to a hydrazine nitrogen (a different, more
  complex shape than a plain P-35.2.1 carbon-branch substituent; not
  covered here).
- An unsaturated substituent, an aromatic substituent other than a plain,
  unsubstituted phenyl (a substituted or heteroaromatic ring, e.g.), or a
  ring other than a plain phenyl substituent.
- Hydrazine derivatives with their own suffix/prefix mechanism: hydrazone
  (P-68.3.1.2.2), azine (P-68.3.1.2.3), semicarbazide (P-68.3.1.2.4),
  hydrazide (R-CO-NH-NH2, P-66.3).
- Charged or isotopically modified atoms.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
)
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _hydrazine_nitrogens(mol):
    """The two nitrogens of a plain hydrazine skeleton: N-N (single bond),
    each nitrogen degree <= 3 (the N-N bond plus at most two other
    single-bonded neighbors), formal charge 0 -- or None if `mol` isn't
    shaped this way."""
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(nitrogens) != 2:
        return None
    n1, n2 = nitrogens
    bond = mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 1.0:
        return None
    if n1.GetDegree() > 3 or n2.GetDegree() > 3:
        return None
    if n1.GetFormalCharge() != 0 or n2.GetFormalCharge() != 0:
        return None
    return n1, n2


def has_hydrazine_shape(mol) -> bool:
    return _hydrazine_nitrogens(mol) is not None


def _substituent_names(graph, n_idx, other_n_idx, halogens, aromatic_atoms, mol=None):
    names = []
    for root in graph[n_idx]:
        if root == other_n_idx:
            continue
        if root in halogens:
            raise UnsupportedStructure(
                "a halogen bonded directly to a hydrazine nitrogen is out "
                "of scope for this module (see module docstring)"
            )
        names.append(name_branch(graph, root, n_idx, halogens, aromatic_atoms, mol=mol))
    return names


def _group(entries):
    grouped = {}
    for locant, name, compound in entries:
        info = grouped.setdefault(name, {"locants": [], "compound": compound})
        info["locants"].append(locant)
    return grouped


def _candidate_key(grouped):
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return locant_set, citation_locants


def name_hydrazine(mol) -> str:
    nitrogens = _hydrazine_nitrogens(mol)
    if nitrogens is None:
        raise UnsupportedStructure(
            "no plain hydrazine (N-N single bond) skeleton found; this "
            "module only handles hydrazines"
        )
    n1, n2 = nitrogens

    aromatic_atoms = frozenset(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the hydrazine's own two nitrogens "
                "(P-68.3.1.2.1) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    n1_idx, n2_idx = n1.GetIdx(), n2.GetIdx()
    if any(
        not (a in aromatic_atoms and b in aromatic_atoms)
        for a, b, _ in non_single_bonds(mol)
    ):
        raise UnsupportedStructure("an unsaturated substituent is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    names_n1 = _substituent_names(graph, n1_idx, n2_idx, halogens, aromatic_atoms, mol=mol)
    names_n2 = _substituent_names(graph, n2_idx, n1_idx, halogens, aromatic_atoms, mol=mol)

    total = len(names_n1) + len(names_n2)
    if total == 0:
        return "hydrazine"
    if total == 1:
        ((name, _),) = names_n1 + names_n2
        return name + "hydrazine"

    candidates = []
    for first, second in ((names_n1, names_n2), (names_n2, names_n1)):
        entries = [(1, name, compound) for name, compound in first] + [
            (2, name, compound) for name, compound in second
        ]
        grouped = _group(entries)
        candidates.append((_candidate_key(grouped), grouped))
    _, best_grouped = min(candidates, key=lambda candidate: candidate[0])
    return format_substituent_prefixes(best_grouped) + "hydrazine"
