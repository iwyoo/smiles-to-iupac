"""Naming of urea (H2N-C(=O)-NH2) and its N-substituted derivatives, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'urea' is a
  retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- the same way 'carbamic acid' (`_carbamate.py`)
  is retained rather than derived. Confirmed via PubChem structure match:
  `NC(=O)N` -> "urea". Unlike every other module in this project, urea has
  no parent-hydride chain to number and no chain-length logic of its own.
- Each nitrogen may carry 0, 1, or 2 plain, unbranched, saturated alkyl
  substituents, cited as 'N-'-prefixed substituents directly ahead of
  'urea', mirroring `_amide.py`'s/`_carbamate.py`'s own N-substitution
  citation -- confirmed directly from the Blue Book's own text
  (`tmp/bluebook/P6.txt` lines 630, 1373-1378): 'N-methyl-N-nitrosourea
  (PIN)' cites two different substituents on the SAME nitrogen both as
  'N-' (no prime needed -- there's only one substituted nitrogen to name),
  and '(i) The symbols N,N' are used for the 'unprimed' parent
  structure...' establishes that a SECOND, distinct nitrogen's
  substituents are cited with a primed 'N'-' instead (the same convention
  `_common.py`'s multi-nitrogen relatives use, e.g.
  'N,N'-methylenediethanamine (PIN)'). PubChem's own generated names for
  these use numeric locants instead (e.g. "1,3-dimethylurea"), a
  structure-only match, not a naming-convention one -- this project
  follows the Blue Book's letter-locant convention, per the source text
  above.

Scope, deliberately narrow: substituents landing on a single nitrogen (one
or two, using the same 'N-'/'N,N-di' citation as `_amide.py`), an
identical single substituent on each of the two different nitrogens
(symmetric 'N,N'-di...' citation), or one DIFFERENT substituent on each of
the two nitrogens -- confirmed via the Blue Book's own PIN worked example
'N-[1-cyano-3-(methylsulfanyl)propyl]-N'-methylurea': the alphanumerical
order (P-14.5.2, locants and italicized prefixes like 'tert-' ignored --
via `alpha_sort_key`) decides which substituent becomes 'N-' and which
becomes 'N''-', e.g. 'N-ethyl-N'-methylurea' (PubChem structure match,
numeric-locant style: `CCNC(=O)NC` -> '1-ethyl-3-methylurea', CID 206567).
Each N-substituent's own name is built with `name_branch` (P-29 PIN
style, fixed project-wide by PR #237; mirrors `_carbamate.py`'s identical
fix, PR #328/#331) -- a branched N-substituent is supported (e.g.
'N-tert-butylurea', CID 14233, a retained non-compound name), and a
*compound* one (has its own locant, e.g. 'propan-2-yl') is always
parenthesized -- 'N-(propan-2-yl)urea', not PubChem's own raw
'N-propan-2-ylurea' (CID 12725), per P-16.5.1.5's/P-66.1.6.1.3.1's own
worked examples ('N-(2-chloroethyl)propan-1-amine (PIN)',
'N-(butan-2-yl)selenourea (PIN)', `tmp/bluebook/P1.html`/`P6a.txt`) --
mirrors `_carbamate.py`'s identical correction. An identical-pair
'N,N-di'/'N,N'-di' name is parenthesized only when the substituent name
is compound (e.g. 'N,N'-di(propan-2-yl)urea', CID 20084) and not when
it's a retained name (e.g. 'N,N'-ditert-butylurea', CID 21420).
Explicitly out of scope (raise `UnsupportedStructure`): a
different substituent *count* on each nitrogen (e.g. one with two
substituents, the other with one -- no confirmed worked example settles
that locant tie-break), an unsaturated N-substituent, a ring-bearing
N-substituent other than a single plain (unsubstituted) benzene ring, a
ring-fused urea (e.g. hydantoin), and thiourea (the sulfur analogue).

- A plain benzene ring bonded directly to one nitrogen is cited as
  'phenyl', mirroring every other alkyl case above ('N-phenylurea',
  PubChem structure match: `NC(=O)Nc1ccccc1` -> "phenylurea", CID 6145 --
  this project's letter-locant convention over PubChem's numeric one, same
  deviation as the plain-alkyl case). Combines with the existing
  symmetric/asymmetric machinery unchanged (`N,N'-diphenylurea`, CID 7595;
  `N-methyl-N'-phenylurea`, CID 13880). A *substituted* phenyl ring (e.g.
  `(4-methylphenyl)urea`, CID 12148) or a second substituent sharing that
  same nitrogen is out of scope.

- Semicarbazide (H2N-NH-C(=O)-NH2, P-68.3.1.4): PubChem structure match
  confirms `NC(=O)NN` -> "aminourea" -- the unsubstituted parent is named
  as urea carrying a plain 'amino' substituent, not with its own retained
  name. Scope here is deliberately limited to that single unsubstituted
  case; any carbon substituent alongside the amino nitrogen is out of
  scope (raise `UnsupportedStructure`) -- PubChem's own examples for that
  combination (`CNC(=O)NN` -> "1-amino-3-methylurea") use a numeric-locant
  style this project doesn't otherwise follow for urea (see above), so
  mixing amino with alkyl substitution needs its own follow-up scoping
  pass rather than being folded in here.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    is_plain_benzene_ring,
    non_single_bonds,
    ring_chain_attachment,
)
from ._substituents import alpha_sort_key, name_branch


def _urea_core(mol):
    """(carbon_idx, (nitrogen1_idx, nitrogen2_idx)) for the urea carbonyl
    carbon and its two nitrogens, or None if the molecule isn't shaped like
    a urea core at all (a carbon with exactly one double-bonded, terminal
    oxygen and two singly-bonded nitrogens, each nitrogen bonded only to
    that carbon and 0-2 carbons besides)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(oxygens) != 1 or len(nitrogens) != 2:
            continue
        (oxygen,) = oxygens
        if oxygen.GetDegree() != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if any(mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in nitrogens):
            continue
        if any(n.GetFormalCharge() != 0 or n.GetIsotope() != 0 for n in nitrogens):
            continue
        if any(
            nn.GetAtomicNum() != 6 and not _is_terminal_amino_nitrogen(nn, n.GetIdx())
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx()
        ):
            continue
        amino_neighbors = [
            nn
            for n in nitrogens
            for nn in n.GetNeighbors()
            if nn.GetIdx() != atom.GetIdx() and _is_terminal_amino_nitrogen(nn, n.GetIdx())
        ]
        if len(amino_neighbors) > 1:
            continue
        return atom.GetIdx(), (nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    return None


def _is_terminal_amino_nitrogen(atom, exclude_idx):
    """True if `atom` is a plain terminal -NH2 nitrogen (semicarbazide's
    extra nitrogen) bonded only to the nitrogen at `exclude_idx`."""
    if atom.GetAtomicNum() != 7 or atom.GetIsAromatic():
        return False
    if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
        return False
    neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
    return neighbors == [exclude_idx]


def has_urea_shape(mol) -> bool:
    return _urea_core(mol) is not None


def _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx):
    """idx of the extra terminal amino nitrogen attached to one of urea's
    two core nitrogens (the semicarbazide shape), or None if neither core
    nitrogen carries one."""
    for n_idx in (n1_idx, n2_idx):
        nitrogen = mol.GetAtomWithIdx(n_idx)
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() in (carbon_idx,):
                continue
            if _is_terminal_amino_nitrogen(neighbor, n_idx):
                return neighbor.GetIdx()
    return None


def _n_substituent_carbons(mol, nitrogen_idx, carbon_idx):
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    return tuple(
        n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx
    )


def _substituent_names(full_graph, nitrogen_idx, substituent_carbons, aromatic_atoms=frozenset(), mol=None):
    return [name_branch(full_graph, c, nitrogen_idx, {}, aromatic_atoms, mol=mol) for c in substituent_carbons]


def _plain_phenyl_substituent_atoms(mol, graph, roots):
    """Union of ring atoms for every plain, unsubstituted benzene ring in
    `mol` that hangs directly off one of `roots` (a urea nitrogen's
    substituent-carbon neighbors) with no other exocyclic attachment --
    i.e. a lone 'phenyl' N-substituent, as opposed to a fused or
    otherwise-substituted ring."""
    atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = set(ring)
        if not is_plain_benzene_ring(mol, ring_atoms):
            continue
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            continue
        ring_atom, _ = attachment
        if ring_atom in roots:
            atoms |= ring_atoms
    return atoms


def _substituent_chain_atoms(carbon_graph, substituent_carbons):
    atoms = set()
    for root in substituent_carbons:
        reached, _ = bfs(carbon_graph, root)
        atoms.update(reached)
    return atoms


def _reject_unsaturated_substituents(mol, atoms):
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def _di_name(name, is_compound):
    return f"({name})" if is_compound else name


def _n_letter_entry(letter, name, is_compound):
    return f"{letter}-({name})" if is_compound else f"{letter}-{name}"


def _n_prefix(letter, entries):
    if not entries:
        return ""
    if len(entries) == 1:
        (name, is_compound), = entries
        return _n_letter_entry(letter, name, is_compound)
    (name_a, compound_a), (name_b, compound_b) = entries
    if name_a == name_b:
        return f"{letter},{letter}-di{_di_name(name_a, compound_a)}"
    (a, ca), (b, cb) = sorted(entries, key=lambda e: alpha_sort_key(e[0]))
    return f"{_n_letter_entry(letter, a, ca)}-{_n_letter_entry(letter, b, cb)}"


def name_urea(mol) -> str:
    core = _urea_core(mol)
    if core is None:
        raise UnsupportedStructure(
            "no urea (H2N-C(=O)-NH2 or an N-substituted derivative) shape "
            "found; this module only handles urea and simple N-substituted "
            "ureas"
        )
    carbon_idx, (n1_idx, n2_idx) = core

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (oxygen_idx,) = (n.GetIdx() for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 8)
    amino_nitrogen_idx = _semicarbazide_amino_nitrogen(mol, n1_idx, n2_idx, carbon_idx)
    n1_carbons = tuple(c for c in _n_substituent_carbons(mol, n1_idx, carbon_idx) if c != amino_nitrogen_idx)
    n2_carbons = tuple(c for c in _n_substituent_carbons(mol, n2_idx, carbon_idx) if c != amino_nitrogen_idx)

    full_graph = adjacency(mol)
    phenyl_atoms = _plain_phenyl_substituent_atoms(mol, full_graph, n1_carbons + n2_carbons)
    if phenyl_atoms and (
        (any(c in phenyl_atoms for c in n1_carbons) and len(n1_carbons) > 1)
        or (any(c in phenyl_atoms for c in n2_carbons) and len(n2_carbons) > 1)
    ):
        raise UnsupportedStructure(
            "a phenyl N-substituent alongside another substituent on the "
            "same nitrogen is not supported yet"
        )

    if mol.GetRingInfo().NumRings() > 0:
        all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
        if all_ring_atoms - phenyl_atoms:
            raise UnsupportedStructure(
                "a ring-fused urea (e.g. hydantoin) or a ring N-substituent "
                "other than a plain, unsubstituted benzene ring is out of "
                "scope for this module"
            )

    carbon_graph = carbon_adjacency(mol)
    n1_chain_atoms = _substituent_chain_atoms(carbon_graph, n1_carbons)
    n2_chain_atoms = _substituent_chain_atoms(carbon_graph, n2_carbons)
    known_atoms = {carbon_idx, oxygen_idx, n1_idx, n2_idx} | n1_chain_atoms | n2_chain_atoms
    if amino_nitrogen_idx is not None:
        known_atoms.add(amino_nitrogen_idx)
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in known_atoms:
            raise UnsupportedStructure(
                "a heteroatom or other characteristic group outside the "
                "urea core and its plain N-alkyl substituents is not "
                "supported yet"
            )

    _reject_unsaturated_substituents(mol, n1_chain_atoms - phenyl_atoms)
    _reject_unsaturated_substituents(mol, n2_chain_atoms - phenyl_atoms)

    n1_names = _substituent_names(full_graph, n1_idx, n1_carbons, frozenset(phenyl_atoms), mol=mol)
    n2_names = _substituent_names(full_graph, n2_idx, n2_carbons, frozenset(phenyl_atoms), mol=mol)

    if amino_nitrogen_idx is not None:
        if n1_names or n2_names:
            raise UnsupportedStructure(
                "a semicarbazide (amino-substituted urea nitrogen) "
                "combined with a plain N-alkyl substituent is not "
                "supported yet -- PubChem's own examples for that "
                "combination use a numeric-locant style this module "
                "doesn't otherwise follow for urea"
            )
        return "aminourea"

    if not n1_names and not n2_names:
        return "urea"

    if n1_names and n2_names:
        if len(n1_names) != 1 or len(n2_names) != 1:
            raise UnsupportedStructure(
                "a different substituent count on each of urea's two "
                "nitrogens is not supported yet (no confirmed worked "
                "example settles the locant tie-break for that case)"
            )
        (name_a, compound_a), (name_b, compound_b) = n1_names[0], n2_names[0]
        if name_a == name_b:
            return f"N,N'-di{_di_name(name_a, compound_a)}urea"
        (first, first_compound), (second, second_compound) = sorted(
            (n1_names[0], n2_names[0]), key=lambda e: alpha_sort_key(e[0])
        )
        first_entry = _n_letter_entry("N", first, first_compound)
        second_entry = _n_letter_entry("N'", second, second_compound)
        return f"{first_entry}-{second_entry}urea"

    names = n1_names or n2_names
    return f"{_n_prefix('N', names)}urea"
