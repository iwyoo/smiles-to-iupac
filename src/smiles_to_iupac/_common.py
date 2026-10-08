"""Shared graph utilities and validation for saturated hydrocarbon parent
hydrides, used by both the acyclic (`_acyclic.py`) and monocyclic (`_cyclic.py`)
naming modules.

- P-35.2.1 (Chapter P-3, https://iupac.qmul.ac.uk/BlueBook/PDF/P3.pdf): 'fluoro',
  'chloro', 'bromo', and 'iodo' are the preselected substituent prefixes for
  -F, -Cl, -Br, and -I respectively. These halogen atoms are always
  monovalent, non-skeletal substituents (P-44.3's "skeletal atoms" always
  means carbon in this repository) — never part of a parent hydride's
  counted chain/ring — so `carbon_adjacency` below lets chain/ring-skeleton
  search ignore them while substituent detection still finds them.
- P-92 (Chapter P-9, https://iupac.qmul.ac.uk/BlueBook/P9.html):
  `specified_stereocenters` below delegates all CIP priority-rule
  computation (atomic number, duplicate-atom treatment, mass number,
  pseudoasymmetry, ...) to RDKit's
  `rdCIPLabeler` rather than reimplementing P-92's rules directly -- this
  project's own contribution is only formatting the resulting label(s)
  into a name, not computing them. P-91.3: when two or more are cited
  together, they're joined in ascending locant order, comma-separated,
  in one parenthesized group, e.g. '(2S,3S)-3-chloro-2-hydroxybutanoic
  acid' (a Blue Book worked example) -- confirmed directly from the
  primary source text, not derived from this project's own reasoning.
"""

import re
from contextvars import ContextVar

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler

from ._numerals import alkane_name, numerical_term

_NUCLIDE = r"(?:[\d,]+-)?\d+[A-Z][a-z]?\d*"
_LEADING_LOCANTS_RE = re.compile(rf"^\x01?(?:[\d,\-]+(?:\((?!{_NUCLIDE}(?:,{_NUCLIDE})*\)[A-Za-z]))?)?")
_ITALIC_PREFIX_RE = re.compile(r"^(tert|sec|iso)-")
_LEADING_ISOTOPE_RE = re.compile(rf"^\x01?\({_NUCLIDE}(?:,{_NUCLIDE})*\)(?=[A-Za-z])")
_LEADING_STEREO_RE = re.compile(r"^\([\dRSEZrsez,' ]+\)-")
_LEADING_ANOMER_RE = re.compile(r"^[αβ]-[DL]-")


class UnsupportedStructure(NotImplementedError):
    pass


HALOGEN_PREFIXES = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALIDE_WORDS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_ALLOWED_ATOMIC_NUMS = {6, *HALOGEN_PREFIXES}

ENE_BOND_ORDER = 2.0
YNE_BOND_ORDER = 3.0


def validate_atoms_and_bonds(mol):
    """Structure-independent checks shared by every parent hydride kind: a
    single-fragment, all-carbon skeleton optionally bearing monovalent
    halogen substituents (P-35.2.1). Bond order (all single, or exactly one
    double/triple bond) and ring shape (none, one simple ring, or more) are
    checked separately by each naming module, since what's allowed there
    differs (see `_acyclic.py`, `_cyclic.py`, `_unsaturated.py`)."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than halogen substituents (F, Cl, Br, I; "
                "see P-35.2.1) are not supported yet (see P-21.2.3, skeletal "
                "replacement nomenclature)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
        elif atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute (see P-44.3)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )


def validate_allowed_atoms(mol, heteroatom_message, group_checks, aromatic_ring_atoms=frozenset()):
    """The same "every atom's atomic number is in this module's own
    allowlist, and no charge/isotope" validation loop shared by every
    two-characteristic-group-coexistence module (e.g. `_aldehyde_amine.py`,
    `_thiol_amine.py`): carbon (tracked for a final "has any carbon" check,
    and rejected if aromatic outside `aromatic_ring_atoms`), a monovalent
    halogen substituent (P-35.2.1), or an atom matching one of
    `group_checks` are allowed; anything else raises `heteroatom_message`
    verbatim (module-specific wording, not reconstructed here).

    `group_checks`: an ordered list of `(atomic_num, allowed_idxs, message)`
    triples, one per non-carbon/non-halogen heteroatom kind the calling
    module allows -- an atom of that atomic number is accepted only if its
    index is in `allowed_idxs`, else `message` is raised. Two groups that
    share an atomic number (e.g. an ether oxygen and an aldehyde carbonyl
    oxygen, both O) are expressed as a single triple with the caller's own
    union of both groups' indices and one combined message, mirroring how
    the modules with this shape already merge that check themselves."""
    allowed_atomic_nums = {6, *HALOGEN_PREFIXES} | {atomic_num for atomic_num, _, _ in group_checks}
    checks_by_atomic_num = {atomic_num: (allowed_idxs, message) for atomic_num, allowed_idxs, message in group_checks}
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_nums:
            raise UnsupportedStructure(heteroatom_message)
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        else:
            allowed_idxs, message = checks_by_atomic_num[atomic_num]
            if atom.GetIdx() not in allowed_idxs:
                raise UnsupportedStructure(message)
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")


def non_single_bonds(mol):
    """List of (begin_atom_idx, end_atom_idx, bond_order) for every bond whose
    order isn't 1.0 (single). Used to classify a molecule's degree of
    unsaturation for dispatch (see `core.py`)."""
    return [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0
    ]


def kekulized_copy(mol):
    """Kekulizes `mol`'s aromatic bonds (P-31.1.4.2's "a benzene ring is
    treated as a cyclohexatriene" rule) only when doing so is unambiguous
    and within this project's current von Baeyer scope: exactly one
    fully-aromatic ring, and no other non-single bond anywhere else in the
    molecule. Two (or more) independently-aromatic rings (e.g. a
    cyclophane), or genuine unsaturation coexisting outside the aromatic
    ring (e.g. a dihydronaphthalene's own ring double bond, or a
    substituent on a bridged-aromatic shape that is not covered), are deliberately left alone -- returns `mol`
    unchanged for those (and for the ordinary all-saturated case), so
    every other caller's existing "not supported yet" behavior stays
    intact. Atom indices are preserved when a copy is made: `Chem.Kekulize`
    only rewrites bond orders/aromatic flags in place, never renumbers
    atoms -- so an already-computed atom-index-based ring topology (e.g. a
    von Baeyer `core`) stays valid against the result."""
    ring_info = mol.GetRingInfo()
    aromatic_rings = [
        ring for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
    ]
    if len(aromatic_rings) != 1:
        return mol
    if any(not bond.GetIsAromatic() and bond.GetBondTypeAsDouble() != 1.0 for bond in mol.GetBonds()):
        return mol
    copy = Chem.RWMol(mol)
    Chem.Kekulize(copy, clearAromaticFlags=True)
    return copy.GetMol()


def reject_unsaturated_substituents(mol, atoms):
    """Raise `UnsupportedStructure` if any bond in `mol` touching `atoms`
    (a characteristic-group module's own N-substituent atom indices) is
    unsaturated -- shared by every urea-family module (urea/thiourea/
    selenourea/tellurourea/guanidine), none of which support an
    unsaturated N-substituent yet."""
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")


def sanitize_probe(mol):
    # a failed trial sanitize is an expected "no", so keep RDKit's own log line off the caller's stderr
    with rdBase.BlockLogs():
        Chem.SanitizeMol(mol)


def reparsable_smiles(mol):
    """Canonical SMILES, or its Kekule form when the canonical one does not parse back, as for a ring whose aluminium,
    gallium, indium or thallium atom is written without an aromatic symbol beside aromatic carbons."""
    smiles = Chem.MolToSmiles(mol)
    with rdBase.BlockLogs():
        if Chem.MolFromSmiles(smiles) is not None:
            return smiles
    kekule = Chem.Mol(mol)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    return Chem.MolToSmiles(kekule, kekuleSmiles=True)


def adjacency(mol):
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        graph[a].append(b)
        graph[b].append(a)
    return graph


def carbon_adjacency(mol):
    """Like `adjacency`, but restricted to carbon atoms and the bonds directly
    between them. Chain/ring-skeleton search (P-44.3.2) must use this instead
    of `adjacency` so a terminal halogen substituent is never mistaken for a
    chain-extending skeletal atom; substituent-detection code should keep
    using the full `adjacency(mol)` so it can still find that halogen."""
    graph = {atom.GetIdx(): [] for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph[a].append(b)
            graph[b].append(a)
    return graph


def ring_chain_attachment(graph, ring_atoms, excluded):
    """If `ring_atoms` (a plain, otherwise-unsubstituted monocyclic ring)
    has exactly one exocyclic branch, hanging off exactly one ring atom,
    return (ring_atom, chain_root); else None. `excluded` are atom indices
    ignored as branch roots (e.g. a ring hydroxyl oxygen already accounted
    for separately) -- pass an empty set/frozenset if there are none."""
    attachments = []
    for atom in ring_atoms:
        branch_roots = [n for n in graph[atom] if n not in ring_atoms and n not in excluded]
        if branch_roots:
            attachments.append((atom, branch_roots))
    if len(attachments) != 1:
        return None
    ring_atom, branch_roots = attachments[0]
    if len(branch_roots) != 1:
        return None
    return ring_atom, branch_roots[0]


def _is_plain_alkyl_branch(mol, graph, root, ring_atom):
    seen = {root}
    stack = [(root, ring_atom)]
    while stack:
        node, previous = stack.pop()
        atom = mol.GetAtomWithIdx(node)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope():
            return False
        for neighbor in graph[node]:
            if neighbor != previous and neighbor not in seen:
                seen.add(neighbor)
                stack.append((neighbor, node))
    return True


def ring_branch_attachment(mol, graph, ring_atoms, excluded=frozenset(), known=()):
    """(ring_atom, chain_root) for the one branch of `ring_atoms` that is not a ring substituent (a halogen, an
    atom in `known`, or a plain acyclic all-carbon group, named by the general ring-group namer); None when there
    is no such branch or several."""
    halogens = halogen_substituents(mol)
    attachment = None
    for atom in ring_atoms:
        roots = [n for n in graph[atom] if n not in ring_atoms and n not in excluded]
        if not roots:
            continue
        if len(roots) != 1:
            return None
        (root,) = roots
        if root in halogens or root in known or _is_plain_alkyl_branch(mol, graph, root, atom):
            continue
        if attachment is not None:
            return None
        attachment = (atom, root)
    return attachment


def separate_aromatic_monocycles(mol, graph):
    """The disjoint plain aromatic monocycles (benzene or a one-heteroatom
    5/6-membered ring) when every ring in `mol` is one and no two rings
    share an atom or are bonded directly; else None. Lets a chain-parent
    module cite each as a substituent via `name_branch`."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings) < 2:
        return None
    for i, ring in enumerate(rings):
        if not (
            is_plain_benzene_ring(mol, ring)
            or heteroaromatic_monocycle_name(mol, ring_cycle(graph, list(ring))) is not None
        ):
            return None
        for other in rings[:i]:
            if ring & other or any(n in other for a in ring for n in graph[a]):
                return None
    return rings


def ring_hosting_anchors(mol, graph, rings, anchors):
    """The benzene ring directly bearing a functional group (one of `anchors`
    is bonded to one of its atoms) among several disjoint aromatic rings, or
    None when the group sits on a chain. Raises when two rings bear it or
    the hosting ring is not a plain benzene ring."""
    hosts = [r for r in rings if any(n in r for a in anchors for n in graph[a])]
    if not hosts:
        return None
    if len(hosts) > 1 or not is_plain_benzene_ring(mol, hosts[0]):
        raise UnsupportedStructure(
            "a functional group directly on a ring alongside another aromatic ring is not supported here yet"
        )
    return hosts[0]


def ring_branch_attachments(mol, graph, rings, excluded=frozenset(), known=()):
    """`ring_branch_attachment` for several rings: one (ring_atom, chain_root) per ring that carries a branch, or
    None if a ring atom has several exocyclic neighbors or a ring has several branches."""
    attachments = []
    for ring in rings:
        attachment = ring_branch_attachment(mol, graph, ring, excluded, known)
        if attachment is None:
            if any(len([n for n in graph[a] if n not in ring and n not in excluded]) > 1 for a in ring):
                return None
            halogens = halogen_substituents(mol)
            if any(
                n not in halogens and n not in known and not _is_plain_alkyl_branch(mol, graph, n, a)
                for a in ring
                for n in graph[a]
                if n not in ring and n not in excluded
            ):
                return None
            continue
        attachments.append(attachment)
    return attachments


def plain_saturated_ring_substituent_atoms(mol, graph, coming_from, root):
    """Ring atom set if `root` sits on a single plain, unsubstituted,
    saturated monocyclic ring whose only exocyclic bond is to
    `coming_from` -- else empty (a ring with a substituent/unsaturation, a
    polycyclic/spiro shape, or no ring at all). Used by any module naming
    a plain saturated ring as a substituent attached through one specific
    neighbor (e.g. `_carbamate.py`'s/`_ester.py`'s ring-on-the-alkoxy-side
    path) -- `name_branch`'s own `_simple_ring_substituent` already names
    the ring itself once it's let through the caller's own ring
    rejection; this only answers whether that shape applies here."""
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = set(ring)
        if root not in ring_atoms:
            continue
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            continue
        attachment = ring_chain_attachment(graph, ring_atoms, set())
        if attachment is None:
            continue
        ring_atom, chain_root = attachment
        if ring_atom == root and chain_root == coming_from:
            return ring_atoms
    return set()


def is_plain_benzene_ring(mol, ring_atoms):
    """True if `ring_atoms` is exactly a 6-membered ring of aromatic carbons
    (a plain, unsubstituted-shape benzene ring) -- used by any chain-parent
    module (`_carboxylic_acid.py`, `_ketone.py`, ...) that names such a ring
    as a 'phenyl' substituent prefix on the chain rather than as the parent
    itself."""
    return len(ring_atoms) == 6 and all(
        mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic()
        for idx in ring_atoms
    )


_HETEROAROMATIC_MONOCYCLE_NAMES = {
    (6, 7): "pyridine",
    (5, 8): "furan",
    (5, 16): "thiophene",
    (5, 7): "pyrrole",
}


def heteroaromatic_monocycle_name(mol, ring_order):
    """Parent hydride name if `ring_order` (an ordered ring walk, e.g. from
    `ring_cycle`) is a plain, fully aromatic 5- or 6-membered monocycle
    with exactly one heteroatom matching one of the four simple
    heteroaromatic monocycles this project names as a "-yl" substituent
    prefix -- pyridine, furan, thiophene, or pyrrole (P-29.3.4.1) -- else
    None. Composition only (no charges/isotopes, exactly one heteroatom,
    the rest aromatic carbon); doesn't check for a single exocyclic
    attachment point itself -- see `ring_chain_attachment` for that."""
    n = len(ring_order)
    if n not in (5, 6):
        return None
    if any(
        mol.GetAtomWithIdx(idx).GetFormalCharge() != 0 or mol.GetAtomWithIdx(idx).GetIsotope() != 0
        for idx in ring_order
    ):
        return None
    if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_order):
        return None
    heteroatoms = [idx for idx in ring_order if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    if len(heteroatoms) != 1:
        return None
    (heteroatom,) = heteroatoms
    return _HETEROAROMATIC_MONOCYCLE_NAMES.get((n, mol.GetAtomWithIdx(heteroatom).GetAtomicNum()))


def _heteroaromatic_monocycle_locant(mol, ring_order, attachment_atom):
    """(locant, indicated_hydrogen) for `attachment_atom` on a heteroaromatic
    monocycle whose heteroatom is fixed at locant 1 (its own established
    parent-hydride numbering) -- `attachment_atom` gets whichever of the two
    ring-walk directions gives it the lower locant (P-29.3.4.1's own worked
    example, "pyridin-2-yl"). `indicated_hydrogen` is "1H-" when pyrrole's
    own indicated hydrogen (P-25.7.1.3) still needs citing -- i.e.
    `attachment_atom` isn't the heteroatom itself, which would otherwise
    consume it (plain "pyrrol-1-yl"/"1-hydroperoxypyrrole", no citation
    needed);
    pyridine/furan/thiophene never need this, since their heteroatom
    carries no H to begin with. Used by `heteroaromatic_monocycle_prefix_name`
    (ring cited as the parent) -- the locant math is identical either way,
    only the surrounding name format differs."""
    n = len(ring_order)
    (heteroatom,) = [idx for idx in ring_order if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
    start = ring_order.index(heteroatom)
    target = ring_order.index(attachment_atom)
    forward = (target - start) % n
    backward = (start - target) % n
    locant = min(forward, backward) + 1
    indicated_hydrogen = ""
    hetero = mol.GetAtomWithIdx(heteroatom)
    if hetero.GetAtomicNum() == 7 and (hetero.GetTotalNumHs() > 0 or hetero.GetDegree() == 3):
        indicated_hydrogen = "1H-"
    return locant, indicated_hydrogen


def heteroaromatic_monocycle_prefix_name(mol, ring_order, attachment_atom, prefix):
    """Name for a heteroaromatic monocycle recognized by
    `heteroaromatic_monocycle_name`, cited as the *parent* hydride with
    `prefix` (e.g. "hydroperoxy") as a substituent prefix at
    `attachment_atom`'s locant -- e.g. "3-hydroperoxypyridine",
    "2-hydroperoxy-1H-pyrrole". Mirrors `_nitro.py`'s plain-benzene
    'nitrobenzene' shape (ring parent, no suffix form for this group, P-
    44.1.2.2) generalized to a heteroaromatic ring, where -- unlike
    benzene -- the substituent's position relative to the heteroatom needs
    an explicit locant. Returns None if `ring_order` isn't one of the four
    recognized rings."""
    name = heteroaromatic_monocycle_name(mol, ring_order)
    if name is None:
        return None
    locant, indicated_hydrogen = _heteroaromatic_monocycle_locant(mol, ring_order, attachment_atom)
    separator = "-" if indicated_hydrogen else ""
    return f"{locant}-{prefix}{separator}{indicated_hydrogen}{name}"


def two_separate_rings_with_plain_aromatic_substituent(mol, graph):
    """If `mol`'s RDKit SSSR reports exactly two rings sharing no atoms,
    joined by exactly one direct bond (not a longer chain -- a distinct,
    separately-handled shape), and exactly one of the two is a plain,
    otherwise-unsubstituted aromatic monocycle (a benzene ring or one of
    the four simple heteroaromatic monocycles, `is_plain_benzene_ring`/
    `heteroaromatic_monocycle_name`), return (other_ring_atoms,
    aromatic_ring_atoms, other_attachment_atom, aromatic_attachment_atom)
    -- else None (three or more rings, a genuinely fused/bridged/spiro
    system whose SSSR rings share atoms, more than one connecting bond,
    neither or both rings aromatic, or the aromatic ring carrying any
    exocyclic substituent of its own beyond the one connecting bond).
    Used by any parent-hydride module (`_ketone.py`, ...) generalizing its
    existing single-ring-plus-plain-substituent dispatch to also accept
    this shape, citing the aromatic ring as a substituent via `name_branch`
    the same way it already cites a plain alkyl ring substituent."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    if len(rings) != 2:
        return None
    ring_x, ring_y = rings
    if ring_x & ring_y:
        return None
    connections = [(a, b) for a in ring_x for b in graph[a] if b in ring_y]
    if len(connections) != 1:
        return None
    x_atom, y_atom = connections[0]

    def is_aromatic_monocycle(ring):
        if is_plain_benzene_ring(mol, ring):
            return True
        return heteroaromatic_monocycle_name(mol, ring_cycle(graph, list(ring))) is not None

    x_is_aromatic = is_aromatic_monocycle(ring_x)
    y_is_aromatic = is_aromatic_monocycle(ring_y)
    if x_is_aromatic == y_is_aromatic:
        return None
    if x_is_aromatic:
        aromatic_ring, aromatic_atom, other_ring, other_atom = ring_x, x_atom, ring_y, y_atom
    else:
        aromatic_ring, aromatic_atom, other_ring, other_atom = ring_y, y_atom, ring_x, x_atom

    for atom in aromatic_ring:
        for neighbor in graph[atom]:
            if neighbor not in aromatic_ring and neighbor != other_atom:
                return None

    return other_ring, aromatic_ring, other_atom, aromatic_atom


def plain_phenyl_substituent_atoms(mol, graph, roots):
    """Union of ring atoms for every plain, unsubstituted benzene ring in
    `mol` that hangs directly off one of `roots` with no other exocyclic
    attachment -- i.e. a lone 'phenyl' substituent directly on one of
    `roots`."""
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


def ordered_chain(graph, root, coming_from, excluded):
    """The chain of atoms starting at `root` and extending away from
    `coming_from`, ignoring `excluded` atoms (e.g. hydroxyl/carboxyl
    oxygens) the same way `_substituents._longest_chains_from_root`
    ignores halogens; None if it branches (more than one non-excluded,
    non-`previous` neighbor at any point)."""
    chain = [root]
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous and n not in excluded]
        if not neighbors:
            return chain
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        chain.append(current)


def longest_branched_chain(graph, source, ring_boundary, excluded=frozenset(), halogens=frozenset()):
    """Generalizes `ordered_chain`: the longest simple chain starting at
    `source` (typically a suffix's own principal-characteristic-group
    carbon, always a tree leaf once its own heteroatoms are excluded)
    and extending through the tree, absorbing any branch encountered
    along the way into the chain itself whenever doing so makes it
    longer (P-44.3.2: the parent chain is one of the longest chains
    containing the principal characteristic group) -- e.g. ibuprofen's
    alpha-methyl becomes part of the parent chain ('propanoic acid'),
    not a separate '2-methyl-...' prefix on a shorter 'ethanoic acid'
    (PubChem CID 3672, '2-[4-(2-methylpropyl)phenyl]propanoic acid').

    `ring_boundary`: atoms the chain itself may never enter (typically
    every ring atom -- this project's own "phenyl chain" modules always
    cite a ring hanging off the chain as a substituent prefix, never as
    part of the parent chain), but which still shows up as an ordinary
    branch wherever it's adjacent to a chosen chain atom.
    `excluded`: atoms that are neither part of the chain nor ever cited
    as a separate branch (typically the suffix's own carbonyl/hydroxyl
    oxygens, already accounted for by the suffix name itself).
    `halogens`: halogen atom indices (typically `halogen_substituents(mol)`,
    the same dict passed to `name_branch`) -- excluded from the chain
    itself (a halogen is always monovalent and would otherwise dead-end
    the BFS at a leaf, silently absorbing it into the chain as if it were
    carbon), but, like `ring_boundary`, still surfaced as an ordinary
    branch wherever one is adjacent to a chosen chain atom (P-35.2.1: a
    halogen is always cited as a substituent prefix, never part of the
    parent chain).

    On a tie for longest, prefers the chain giving the greater number of
    substituents cited as prefixes (P-44.3.2's own next tie-break after
    chain length), then the lowest set of locants among those (P-14.5.2)
    -- e.g. 'c1ccccc1CC(C)C(=O)O' (both a phenyl-bearing carbon and a
    methyl-bearing carbon are one bond from the acid's C2, tied for
    farthest) is '2-methyl-3-phenylpropanoic acid' (two substituents:
    'methyl' + 'phenyl'), not '2-(phenylmethyl)propanoic acid' (one
    compound 'benzyl'-shaped substituent) -- PubChem PUG REST verified.

    Returns (chain, branches): `chain` is the winning path as an
    atom-index list with `source` first (so 1-based `enumerate(chain,
    start=1)` locants match the suffix's own fixed-C1 convention);
    `branches` is {position -> [branch_root_atom, ...]} for every
    neighbor of a chain atom that isn't itself on the chain or in
    `excluded` -- name each via `name_branch`, same as any other
    substituent (a ring-atom branch root is named as a ring substituent
    automatically, since `name_branch` already recognizes one)."""
    blocked = set(ring_boundary) | set(excluded) | set(halogens)
    dist = {source: 0}
    parent = {source: None}
    queue = [source]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor in blocked or neighbor in dist:
                    continue
                dist[neighbor] = dist[node] + 1
                parent[neighbor] = node
                next_queue.append(neighbor)
        queue = next_queue
    farthest = max(dist.values())

    def branches_for(chain):
        chain_set = set(chain)
        branches = {}
        for position, atom in enumerate(chain, start=1):
            roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
            if roots:
                branches[position] = roots
        return branches

    best_chain, best_branches, best_key = None, None, None
    for node, d in dist.items():
        if d != farthest:
            continue
        chain = path_between(parent, source, node)
        branches = branches_for(chain)
        substituent_count = sum(len(roots) for roots in branches.values())
        locants = lowest_locant_set(pos for pos, roots in branches.items() for _ in roots)
        key = (-substituent_count, locants)
        if best_key is None or key < best_key:
            best_key, best_chain, best_branches = key, chain, branches
    return best_chain, best_branches


def longest_branched_chain_through(graph, required, ring_boundary, excluded=frozenset(), halogens=frozenset(), principal=frozenset()):
    """Like `longest_branched_chain`, but `required` need not be a chain
    terminus (e.g. a ketone's own carbonyl carbon, always internal once
    its aryl-ketone case is separately rejected) -- finds one of the
    longest chains in the tree that includes `required` somewhere along
    it, by combining `required`'s two longest 'arms' (P-44.3.2, the same
    "longest chain containing the principal characteristic group" rule,
    here allowing the group to sit anywhere on the chain instead of
    fixing it at C1).

    On a tie for longest, prefers the combination giving the greater
    number of substituents cited as prefixes, then the lowest set of
    locants among those -- the same P-44.3.2 tie-break `longest_branched_
    chain` already applies, generalized to the two-arm case (found via
    real-data testing: a branch point with a halogen-bearing 1-carbon arm
    tied in length against a plain methyl arm, e.g. 'C(F)(F)Br' vs 'C'
    off the same alcohol carbon, was previously resolved arbitrarily by
    BFS/adjacency insertion order rather than by this rule, e.g.
    'CC(O)(Cc1cccc(F)c1)C(F)(F)Br' wrongly named
    '2-(bromodifluoromethyl)-1-(3-fluorophenyl)propan-2-ol' instead of
    PubChem's '1-bromo-1,1-difluoro-3-(3-fluorophenyl)-2-methylpropan-2-ol',
    which absorbs the halogen-bearing carbon into the chain instead,
    citing 5 prefix substituents instead of 2). Every candidate two-arm
    combination (across every pair of distinct starting neighbors of
    `required`, and every farthest node tied within each) is enumerated
    and scored the same way, since `required` rarely has more than two or
    three non-excluded neighbors -- this stays cheap.

    `ring_boundary`/`excluded`/`halogens`: same meaning as
    `longest_branched_chain`. `principal`: carbons bearing a principal
    characteristic group -- a chain holding more of them wins before length
    does (P-44.1.1).

    Returns (chain, branches): `chain` is the winning path as an
    atom-index list, in an arbitrary direction -- the caller tries both
    ways (same as it already does for a plain unbranched chain) to give
    `required` its own lowest locant; `branches` is {position ->
    [branch_root_atom, ...]}, 1-based against this `chain`'s order, same
    shape as `longest_branched_chain`."""
    blocked = set(ring_boundary) | set(excluded) | set(halogens)
    neighbors = [n for n in graph[required] if n not in blocked]

    groups = []
    for start in neighbors:
        # `required` itself must stay off-limits here -- otherwise this
        # arm's search loops back through it into the *other* arm(s),
        # corrupting the two-arm split with a self-crossing path.
        arm_blocked = blocked | {required}
        dist = {start: 0}
        parent = {start: None}
        queue = [start]
        while queue:
            next_queue = []
            for node in queue:
                for neighbor in graph[node]:
                    if neighbor in arm_blocked or neighbor in dist:
                        continue
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
            queue = next_queue
        far = max(dist.values())
        groups.append([
            [required] + path_between(parent, start, node)
            for node, d in dist.items()
            if d == far or (principal and node in principal)
        ])

    def candidate_chains():
        if not groups:
            yield [required]
        elif len(groups) == 1:
            for arm in groups[0]:
                yield list(reversed(arm))
        else:
            for i in range(len(groups)):
                for j in range(len(groups)):
                    if i == j:
                        continue
                    for arm1 in groups[i]:
                        for arm2 in groups[j]:
                            yield list(reversed(arm1)) + arm2[1:]

    def branches_for(chain):
        chain_set = set(chain)
        branches = {}
        for position, atom in enumerate(chain, start=1):
            # Any neighbor not on the two chosen arms -- including a
            # third+ arm off `required` itself, e.g. a ketone carbon with
            # more than two carbon substituents (shouldn't normally
            # arise, but falls through safely here rather than being
            # silently dropped) -- is an ordinary branch, named via
            # `name_branch` like any other.
            roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
            if roots:
                branches[position] = roots
        return branches

    best_chain, best_branches, best_key = None, None, None
    for chain in candidate_chains():
        branches = branches_for(chain)
        substituent_count = sum(len(roots) for roots in branches.values())
        locants = lowest_locant_set(pos for pos, roots in branches.items() for _ in roots)
        key = (-sum(1 for atom in chain if atom in principal), -len(chain), -substituent_count, locants)
        if best_key is None or key < best_key:
            best_key, best_chain, best_branches = key, chain, branches
    return best_chain, best_branches


def is_nitro_nitrogen(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 7 or atom.GetDegree() != 3:
        return False
    oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
    if len(oxygens) != 2 or atom.GetFormalCharge() not in (0, 1):
        return False
    orders = sorted(mol.GetBondBetweenAtoms(idx, o.GetIdx()).GetBondTypeAsDouble() for o in oxygens)
    charges = sorted(o.GetFormalCharge() for o in oxygens)
    if atom.GetFormalCharge():
        return orders == [1.0, 2.0] and charges == [-1, 0]
    return orders == [2.0, 2.0] and charges == [0, 0]


YLO_MAP_NUMBER = 9001


def named_prefix(atom):
    """The prefix a placeholder atom stands for ('ylo' marks a radical centre cited as a prefix, P-71.5); the map
    number carries it through a SMILES round trip."""
    if atom.HasProp("_named_prefix"):
        return atom.GetProp("_named_prefix")
    return "ylo" if atom.GetAtomMapNum() == YLO_MAP_NUMBER else None


def halogen_prefix(atom):
    return named_prefix(atom) or HALOGEN_PREFIXES[atom.GetAtomicNum()]


def halogen_substituents(mol):
    """{atom_idx -> substituent prefix name} for every halogen atom in `mol`
    (P-35.2.1). Passed down into `name_branch` so it can name a halogen leaf
    and exclude halogens from a compound substituent's own internal chain
    search, the same way `carbon_adjacency` does for a parent hydride."""
    return {
        atom.GetIdx(): halogen_prefix(atom)
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetDegree() == 1
    }


def find_primary_amines(mol, exclude=frozenset()):
    """Atom indices of every primary amine nitrogen in `mol` (-NH2 bonded
    to exactly one carbon) not in `exclude` -- shared by every coexisting-
    functional-group pairwise module whose junior group is a primary
    amine (`_alcohol_amine.py`, `_carboxylic_acid_amine.py`, ...); a
    caller with an atom to exclude (e.g. an already-accounted-for amide/
    ester/ether oxygen's own attachment) passes it in `exclude`."""
    amines = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() in exclude or atom.GetAtomicNum() != 7 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 or atom.GetTotalNumHs() != 2:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            amines.add(atom.GetIdx())
    return amines


def find_ether_oxygens(mol, exclude=frozenset(), require_sole=False):
    """Atom(s) shaped like a plain ether oxygen (degree 2, both neighbors
    carbon) -- shared by every coexisting-functional-group pairwise
    module whose junior group is an ether. `require_sole=False`
    (default): every ether-shaped oxygen in the molecule not in
    `exclude` (a caller's own bridging oxygen with the same local shape,
    e.g. an ester's R-CO-O-R'). `require_sole=True`: the molecule's
    *sole* oxygen overall, only if it's ether-shaped -- a genuinely
    different, stricter shape (a molecule with a second, non-ether-shaped
    oxygen is rejected outright, not just filtered), used by a junior
    module whose senior group's own oxygen(s) aren't otherwise excluded
    at the call site."""
    if require_sole:
        oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
        if len(oxygens) != 1:
            return []
        (oxygen,) = oxygens
        if oxygen.GetDegree() == 2 and all(n.GetAtomicNum() == 6 for n in oxygen.GetNeighbors()):
            return [oxygen]
        return []
    return [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 8
        and atom.GetIdx() not in exclude
        and atom.GetDegree() == 2
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    ]


def unsubstituted_sulfonamide_sulfur_atoms(mol):
    """Sulfur atoms shaped like an *unsubstituted* sulfonamide group
    (-SO2NH2): bonded to exactly one carbon, two double-bonded (terminal)
    oxygens, and one single-bonded nitrogen that is itself terminal (two
    hydrogens, no other substituents). N-alkylated sulfonamides don't
    match -- deliberately narrower than `_sulfonamide.py`'s own sulfonamide-
    sulfur finder (which does accept N-alkylated forms), shared by the
    coexisting-functional-group pairwise modules that only support the
    unsubstituted case (`_carboxylic_acid_sulfonamide.py`,
    `_sulfonic_acid_sulfonamide.py`)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(oxygens) != 2 or len(nitrogens) != 1:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(double_os) != 2 or any(o.GetDegree() != 1 for o in double_os):
            continue
        (nitrogen,) = nitrogens
        if mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if nitrogen.GetDegree() != 1 or nitrogen.GetTotalNumHs() != 2:
            continue
        matches.append(atom)
    return matches


def bfs(graph, start):
    dist = {start: 0}
    parent = {start: None}
    queue = [start]
    while queue:
        next_queue = []
        for node in queue:
            for neighbor in graph[node]:
                if neighbor not in dist:
                    dist[neighbor] = dist[node] + 1
                    parent[neighbor] = node
                    next_queue.append(neighbor)
        queue = next_queue
    return dist, parent


def component_subgraph(graph, start):
    """`graph` restricted to the connected component containing `start`
    (e.g. one acyl branch of a symmetric multi-fragment structure once the
    bridging atom(s) are excluded from `graph` itself)."""
    dist, _ = bfs(graph, start)
    nodes = set(dist)
    return {node: [n for n in graph[node] if n in nodes] for node in nodes}


def path_between(parent, start, end):
    path = [end]
    while path[-1] != start:
        path.append(parent[path[-1]])
    return list(reversed(path))


def ring_cycle(graph, ring_atoms):
    """Order a monocyclic ring's atoms into a single walk around the ring,
    starting from `ring_atoms[0]` (an arbitrary RDKit ring-atom listing
    order) -- both traversal directions are still tried by the caller when
    picking the winning numbering."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def linear_branch(graph, root, coming_from):
    """Walk a branch outward; return its atom count, or None if it forks
    (a "compound" substituent, P-29.4, not yet supported)."""
    length = 1
    previous, current = coming_from, root
    while True:
        neighbors = [n for n in graph[current] if n != previous]
        if len(neighbors) == 0:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]
        length += 1


def unbranched_chain_length(mol, root_idx, exclude_idx):
    """Length of the straight, unbranched, saturated all-carbon chain
    starting at `root_idx` and walking away from `exclude_idx` -- or None
    if the chain branches, rings, or leaves carbon at any point."""
    length = 0
    previous = exclude_idx
    current = root_idx
    while True:
        atom = mol.GetAtomWithIdx(current)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
            return None
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != previous]
        length += 1
        if not neighbors:
            return length
        if len(neighbors) > 1:
            return None
        previous, current = current, neighbors[0]


def lowest_locant_set(locants):
    return tuple(sorted(locants))


def multiplied_word(count, base):
    """Multiplying-prefix word for `count` occurrences of a suffix like 'ol'/
    'one'/'al' (P-14.2.1): omitted for zero, bare for one, else prefixed with
    the basic numerical term ('di', 'tri', ...).

    P-16.3.3: a multiplying prefix's terminal 'a' (tetra, penta, hexa, ...)
    is elided before a suffix beginning with 'a' or 'o' -- e.g. 'tetra' +
    'ol' -> 'tetrol', not 'tetraol' (PubChem CID 8998, confirming
    "butane-1,2,3,4-tetrol"); 'tetra' + 'amine' -> 'tetramine' (CID
    6395580). 'di'/'tri' never end in 'a', so they're never affected."""
    if count == 0:
        return ""
    if count == 1:
        return base
    prefix = numerical_term(count)
    if prefix.endswith("a") and base[:1] in "ao":
        prefix = prefix[:-1]
    return prefix + base


def elides_before(word: str) -> bool:
    """True if a preceding word ending in 'e' (e.g. an 'ene' locant-suffix
    segment assembled by `multiplied_word`) should have that 'e' elided
    right before `word`. P-31.1.1.1's 'dien-...-yne' elision is specific to
    a following yne/diyne/triyne/... segment -- it always applies there
    regardless of the multiplying prefix, since 'yne' itself is 'y'-initial
    (found via real-data testing: 'deca-1,2,3-trien-5,7,9-triyne' and
    'methyl 18-bromooctadeca-9,17-dien-5,7-diynoate' were wrongly left
    unelided by the ~30 suffix modules that had inlined this check against
    the final word's own first letter). For every other following word
    (oic, diazonium, oate, ...) elision follows that word's own literal
    first letter, not its underlying un-prefixed form -- 'but-2-enedioic
    acid' and 'pent-4-ene-1-diazonium' both stay unelided even though
    'dioic'/'diazonium' start with a consonant only because of their own
    multiplying prefix."""
    if word.endswith("yne"):
        return True
    return word[:1] in "aeiouy"


def suffix_body(ene_locants, yne_locants, own_word, own_locants=None):
    """Locant-and-suffix string combining 'ene'/'yne' unsaturation segments
    with a module's own functional-group suffix (e.g. '4-en-1-ol',
    '2-enoate'), plus whether the chain stem's trailing 'e' should be
    elided at the stem/first-segment boundary (only relevant when there's
    no 'ene'/'yne' segment).

    `own_word` is the already-multiplied suffix word (`multiplied_word`
    applied by the caller with that module's own count -- what's being
    counted varies per functional group, so this stays caller-owned).
    `own_locants` selects citation shape: a non-empty iterable cites those
    locants as one more hyphenated segment just like ene/yne (e.g.
    '1-ol', '1-sulfonic acid'); `None`/empty glues `own_word` directly
    onto the end with no locant cited (P-14.3.3: this group's own locant
    is always 1 and never written), e.g. 'oate', 'oic', 'imidamide'."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    if own_locants:
        segments.append((sorted(own_locants), own_word))
        words = [word for _, word in segments]
        for i in range(len(words) - 1):
            if words[i].endswith("e") and elides_before(words[i + 1]):
                words[i] = words[i][:-1]
        parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words)
        ]
        body = "-".join(parts)
    else:
        words = [word for _, word in segments] + [own_word]
        for i in range(len(words) - 1):
            if words[i].endswith("e") and elides_before(words[i + 1]):
                words[i] = words[i][:-1]
        if segments:
            locant_parts = [
                f"{','.join(str(loc) for loc in locants)}-{word}"
                for (locants, _), word in zip(segments, words[:-1])
            ]
            body = "-".join(locant_parts) + words[-1]
        else:
            body = words[-1]

    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def should_omit_mononuclear_locants(chain_length, own_locants, has_unsaturation, substituted=False):
    """True when P-14.3.4.2(a) (`chain_length == 1`) or P-14.3.4.2(b) (a
    saturated two-carbon chain whose suffix has exactly one own locant)
    applies, so the caller's own `format_substituent_prefixes` call should
    pass `omit_locants=True`. At `chain_length == 2` the suffix's own
    locant is always forced to be C1 by P-44.4.1.8 (the suffix locant is
    minimized ahead of every other numbering criterion), so citing it is
    never informative regardless of how many other substituents are
    present or where they sit -- confirmed against PubChem across several
    functional-group families, e.g. 'CC(Cl)S(=O)(=O)O' ->
    '1-chloroethanesulfonic acid', 'ClCCNCC' ->
    '2-chloro-N-ethylethanamine', 'CC(=O)C1CCCCC1' ->
    '1-cyclohexylethanone' (all PubChem PUG REST-verified). A suffix with
    2+ own locants (e.g. a diol) is unaffected by this rule -- P-14.3.3's
    'no locant when it's the sole possible position' only applies when
    there's exactly one occurrence to place, so 'ethane-1,2-diol' still
    cites both locants ('OCCO', PubChem-verified).

    Safe even when non-numeric (e.g. 'N-') locants are also present in
    the caller's own prefix call, since `format_substituent_prefixes`
    never omits those regardless of `omit_locants` -- see its own
    docstring."""
    own_locants = own_locants or ()
    if chain_length == 1:
        return True
    return chain_length == 2 and not has_unsaturation and not substituted and len(own_locants) == 1


RETAINED_ACYL_STEMS = {
    (1, "amide"): "formamide",
    (2, "amide"): "acetamide",
    (1, "hydrazide"): "formohydrazide",
    (2, "hydrazide"): "acetohydrazide",
    (2, "diamide"): "oxamide",
    (2, "nitrile"): "acetonitrile",
    (2, "dinitrile"): "oxalonitrile",
    (1, "al"): "formaldehyde",
    (2, "al"): "acetaldehyde",
}


def name_from_substituents(
    chain_length, ene_locants, yne_locants, own_word, own_locants=None, force_own_locant=False, substituted=False
):
    """Assemble `<stem>[a]<separator><suffix body>` for an acyclic
    chain-parent suffix module (the caller still prepends its own
    `format_substituent_prefixes(grouped)` -- `_substituents.py` already
    imports from this module, so this one can't import back from it):
    pick the alkane stem (unsaturated form with `-ene`/`-yne` elided when
    `ene_locants`/`yne_locants` are present, inserting the linking `a` for
    2+ multiple bonds), delegate the ene/yne-plus-own-suffix segment to
    `suffix_body`, then glue everything together -- hyphenated whenever
    the suffix body starts with a cited locant (either because the chain
    is unsaturated, P-14.3.4.2(b)/(c), or because `own_locants` was
    passed, P-31.1.4.3.4), glued directly otherwise (P-14.3.3: a group
    whose own locant is always 1 is never cited).

    `force_own_locant=True` skips the P-14.3.4.2(b) two-carbon omission
    for a caller whose own suffix locant coincides with a cited
    substituent-prefix locant, where the Blue Book cites it anyway (e.g.
    P-72.2.2.1's 'acetyl anion' -> '1-oxoethan-1-ide (PIN)', not
    'ethanide' as it would be for the same suffix alone -- `_carbanide.py`
    is the only caller that passes this)."""
    has_unsaturation = bool(ene_locants or yne_locants)

    retained = None if has_unsaturation or force_own_locant else RETAINED_ACYL_STEMS.get((chain_length, own_word))
    if retained:
        return retained

    if not force_own_locant and should_omit_mononuclear_locants(
        chain_length, own_locants, has_unsaturation, substituted
    ):
        stem = alkane_name(chain_length)
        if own_word[0] in "aeiouy":
            stem = stem[:-1]
        return stem + own_word

    if (
        chain_length == 2
        and not force_own_locant
        and not substituted
        and len(ene_locants or ()) + len(yne_locants or ()) == 1
        and (own_locants is None or list(own_locants) == [1])
    ):
        stem = "ethyn" if yne_locants else "ethen"
        return stem + own_word if own_word[0] in "aeiouy" else stem + "e" + own_word

    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = suffix_body(ene_locants, yne_locants, own_word, own_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if (has_unsaturation or own_locants) else ""
    return stem + ("a" if needs_stem_a else "") + separator + body


def ring_name_from_substituents(ring_size, ene_locants, yne_locants, prefix, total_subs, own_word, own_locants):
    """Assemble a ring-parent suffix name from `"cyclo" + alkane_name(ring_size)`,
    the caller-supplied `prefix` (`format_substituent_prefixes(grouped)` --
    see `name_from_substituents` for why this stays caller-owned) and
    `total_subs` (`sum(len(info["locants"]) for info in grouped.values())`).

    Saturated ring (no `ene_locants`/`yne_locants`): elides the stem's
    trailing vowel-initial `own_word` (`cyclohexan` + `one`), and when
    `own_word`'s locants are the ring's only substituent (P-14.3.3) omits
    the locant entirely, e.g. 'cyclohexanone'; otherwise cites it as a
    hyphenated segment, e.g. '4-methylcyclohexan-1-one'.

    Unsaturated ring (P-31.1.3): the ring double/triple bond means the
    suffix locant is never omittable even when it's the sole substituent,
    e.g. 'cyclohex-2-en-1-one' -- delegates the ene/yne-plus-own-suffix
    segment to `suffix_body`."""
    has_unsaturation = bool(ene_locants or yne_locants)
    parent = "cyclo" + alkane_name(ring_size)

    if not has_unsaturation:
        elide = own_word[0] in "aeiouy"
        stem = parent[:-1] if elide else parent
        if total_subs == 0 and len(own_locants) == 1:
            return stem + own_word
        loc_str = ",".join(str(loc) for loc in sorted(own_locants))
        return f"{prefix}{stem}-{loc_str}-{own_word}"

    stem = parent[:-3]
    needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    body, _ = suffix_body(ene_locants, yne_locants, own_word, own_locants)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def group_substituents(substituents):
    """{position -> [(name, is_compound), ...]} -> {name -> {"locants": [...],
    "compound": bool}}, merging same-named substituents at different
    positions so they can be cited once with a multiplying prefix (P-16)."""
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def alpha_sort_key(name: str) -> str:
    """P-14.5.2: alphanumerical ordering ignores locants and italicized
    prefixes like 'tert-' -- only the rest of the name counts (so
    'tert-butyl' sorts under 'b', not 't')."""
    stripped = name
    previous = None
    while stripped != previous:
        previous = stripped
        stripped = _LEADING_ISOTOPE_RE.sub("", stripped)
        stripped = _LEADING_STEREO_RE.sub("", stripped)
        stripped = _LEADING_ANOMER_RE.sub("", stripped)
        stripped = _LEADING_LOCANTS_RE.sub("", stripped)
        if stripped[:1] in ("(", "[", "{") and not (_LEADING_STEREO_RE.match(stripped) or _LEADING_ISOTOPE_RE.match(stripped)):
            stripped = stripped[1:]
    stripped = _ITALIC_PREFIX_RE.sub("", stripped)
    return re.sub(r"[^a-z]", "", stripped.lower())


def citation_order_key(name: str):
    """Alphanumerical order, an isotopically modified substituent ahead of the unmodified one (P-82.2.2.1), then configuration (P-45.6.3)."""
    return alpha_sort_key(name), not _LEADING_ISOTOPE_RE.match(name), _descriptor_precedence(name)


_DESCRIPTOR_PRECEDENCE = {"Z": 0, "E": 1, "R": 2, "S": 3, "r": 4, "s": 5}


def _descriptor_precedence(name: str):
    """P-45.6.3: when names tie, Z precedes E, then R precedes S, compared in cited order."""
    match = _LEADING_STEREO_RE.match(_LEADING_ISOTOPE_RE.sub("", name))
    if match is None:
        return ()
    return tuple(_DESCRIPTOR_PRECEDENCE[c] for c in match.group() if c in _DESCRIPTOR_PRECEDENCE)


def substituent_locant_set_and_citation(grouped):
    """The P-45.2.1-2.3 "lowest-locant-set, then alphabetical citation
    order" computation shared by every suffix module's own `_candidate_key`/
    `_ring_candidate_key` sort key: `grouped` is `group_substituents`'s own
    output. Returns `(locant_set, total_count, citation_locants)` -- each
    module's own extra suffix/ene/yne locant-set fields, final sort-key
    tuple assembly, and name construction stay module-specific and are not
    part of this shared computation."""
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    return locant_set, total_count, citation_locants


def longest_chains(graph):
    """Every longest path (by atom count) through an undirected acyclic
    `graph` (P-44.3, the principal chain must be a longest chain candidate)."""
    nodes = list(graph)
    distances = {}
    parents = {}
    for node in nodes:
        dist, parent = bfs(graph, node)
        distances[node] = dist
        parents[node] = parent

    diameter = max(d for dist in distances.values() for d in dist.values())
    chains = []
    seen = set()
    for u in nodes:
        for v, d in distances[u].items():
            if d == diameter and (v, u) not in seen:
                seen.add((u, v))
                chains.append(path_between(parents[u], u, v))
    return chains


def _cycle_nodes(graph):
    """Nodes lying on a cycle of `graph`: what remains after repeatedly
    stripping degree-1 nodes is the cycle-bearing core, and a node of that
    core on no cycle is a bridge path between rings, so the test is a
    reachability one."""
    degree = {n: len(nbrs) for n, nbrs in graph.items()}
    stack = [n for n, d in degree.items() if d <= 1]
    removed = set()
    while stack:
        n = stack.pop()
        if n in removed:
            continue
        removed.add(n)
        for m in graph[n]:
            if m not in removed:
                degree[m] -= 1
                if degree[m] <= 1:
                    stack.append(m)
    core = set(graph) - removed
    on_cycle = set()
    for n in core:
        for m in graph[n]:
            if m in core:
                seen, queue = {n}, [m]
                found = False
                while queue and not found:
                    x = queue.pop()
                    if x == n:
                        found = True
                        break
                    if x in seen:
                        continue
                    seen.add(x)
                    queue.extend(y for y in graph[x] if y in core and not (x == m and y == n))
                if found:
                    on_cycle.add(n)
                    break
    return on_cycle


def carbon_on_ring(mol, graph, carbons):
    """True if any of the group-bearing `carbons` is bonded to a ring atom
    (a ring-attached group takes the 'carbo...' suffix instead of a chain)."""
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    return any(n in ring_atoms for c in carbons for n in graph[c])


def all_chains(graph):
    """Every simple path (one orientation each, including single atoms) through
    an undirected acyclic `graph` -- the candidates among which the principal
    chain is the longest one still carrying every principal characteristic
    group (P-44.1.1, then P-44.3)."""
    ring_atoms = _cycle_nodes(graph)
    graph = {n: [m for m in nbrs if m not in ring_atoms] for n, nbrs in graph.items() if n not in ring_atoms}
    chains = []
    for start in graph:
        dist, parent = bfs(graph, start)
        for end in dist:
            if end >= start:
                chains.append(path_between(parent, start, end))
    return chains


def bond_locant(chain, bond_atoms):
    """1-based position along `chain` of the bond between `bond_atoms` (a
    2-tuple of atom indices), or None if that bond doesn't lie on `chain`."""
    bond_set = set(bond_atoms)
    for i in range(len(chain) - 1):
        if {chain[i], chain[i + 1]} == bond_set:
            return i + 1
    return None


def bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for every (a, b, order) bond in `bonds` that
    lies on `chain`, or None if any bond doesn't (a shorter/wrong chain
    candidate)."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = bond_locant(chain, (a, b))
        if locant is None:
            return None
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return ene, yne


def unsaturation_suffix(ene_locants, yne_locants):
    """Locant-and-suffix string (e.g. '1,3-dien-5-yne') plus whether the
    parent stem needs its euphonic trailing 'a' (P-31.1.1.2), for a chain's
    full set of multiple bonds. 'ene' is always cited before 'yne'
    (P-31.1.1.1), with its final 'e' always elided when a 'yne' part
    follows, whether or not that 'yne' itself carries a multiplying
    prefix -- the elision is triggered by the underlying 'yne' word
    starting with a vowel sound, not by the final prefixed word's own
    first letter ('deca-1,2,3-trien-5,7,9-triyne', PubChem-verified: a
    'diyne'/'triyne' elides exactly like a plain 'yne' does, contrary to
    this function's own former assumption that a multiplying-prefixed
    'yne' word "begins with a consonant and elides nothing")."""
    ene_locants = sorted(ene_locants)
    yne_locants = sorted(yne_locants)
    ene_count, yne_count = len(ene_locants), len(yne_locants)
    ene_word = multiplied_word(ene_count, "ene")
    yne_word = multiplied_word(yne_count, "yne")

    if ene_count and yne_count:
        ene_part = ene_word[:-1]
        ene_loc_str = ",".join(str(loc) for loc in ene_locants)
        yne_loc_str = ",".join(str(loc) for loc in yne_locants)
        body = f"{ene_loc_str}-{ene_part}-{yne_loc_str}-{yne_word}"
    elif ene_count:
        body = f"{','.join(str(loc) for loc in ene_locants)}-{ene_word}"
    else:
        body = f"{','.join(str(loc) for loc in yne_locants)}-{yne_word}"

    needs_stem_a = (ene_count >= 2) if ene_count else (yne_count >= 2)
    return body, needs_stem_a


def chain_bond_locants(chain, bonds):
    """(ene_locants, yne_locants) for the multiple bonds lying on `chain`;
    any other multiple bond is cited inside a substituent prefix."""
    ene, yne = [], []
    for a, b, order in bonds:
        locant = bond_locant(chain, (a, b))
        if locant is not None:
            (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return ene, yne


def most_multiple_bonds(chains, bonds):
    """The `chains` carrying the greatest number of multiple bonds, then of
    double bonds (P-44.4.1.1, P-44.4.1.2) -- the choice among chains already
    tied on every senior criterion."""
    if not bonds or not chains:
        return chains

    def score(chain):
        ene, yne = chain_bond_locants(chain, bonds)
        return len(ene) + len(yne), len(ene)

    best = max(score(chain) for chain in chains)
    return [chain for chain in chains if score(chain) == best]


def von_baeyer_bond_citation(position, a, b):
    """(primary_locant, display, is_compound) for one double/triple bond
    of a von Baeyer bicyclic/polycyclic parent under a fixed `position`
    numbering (atom index -> 1-based locant) -- P-31.1.4.1's simple case
    when the two atoms have consecutive locants ('2', no parentheses),
    P-31.1.4.2(1)'s compound-locant case otherwise (the higher locant
    cited in parentheses after the lower one, e.g. '1(7)'). Unlike a
    chain's or a monocyclic ring's own bond locant, a von Baeyer parent's
    numbering has no ring-wraparound between its highest and lowest
    locants in general, so this never wraps."""
    pa, pb = position[a], position[b]
    lo, hi = min(pa, pb), max(pa, pb)
    if hi - lo == 1:
        return lo, str(lo), False
    return lo, f"{lo}({hi})", True


def von_baeyer_bond_stereo(mol, position, bond_stereo):
    """(sorted primary locants of the Z bonds, '(4Z,6E)-' prefix or '') of the specified double-bond stereo elements
    `bond_stereo` (from `specified_double_bond_stereo`) under a fixed von Baeyer `position` numbering (P-93.5.1.4)."""
    if not bond_stereo:
        return (), ""
    cited = []
    for bond_idx, code in bond_stereo:
        bond = mol.GetBondWithIdx(bond_idx)
        primary, display, _ = von_baeyer_bond_citation(position, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        cited.append((primary, f"{display}{code}", code))
    cited.sort()
    return tuple(p for p, _, c in cited if c == "Z"), "(" + ",".join(d for _, d, _ in cited) + ")-"


def von_baeyer_unsaturation_citations(position, bonds):
    """(ene_citations, yne_citations, compound_count, primary_locants,
    full_locants) for every (a, b, order) von Baeyer bond under a fixed
    `position` numbering -- each citation is a (primary_locant, display)
    pair from `von_baeyer_bond_citation`. `primary_locants` is every
    bond's own primary (never-parenthesized) locant, for P-31.1.4.2(2)'s
    "ignore parenthesized numbers" comparison; `full_locants` additionally
    includes each compound bond's own parenthesized locant, for
    P-31.1.4.2(3)'s full-locant-set tie-break."""
    ene_citations, yne_citations = [], []
    compound_count = 0
    primary_locants, full_locants = [], []
    for a, b, order in bonds:
        primary, display, is_compound = von_baeyer_bond_citation(position, a, b)
        primary_locants.append(primary)
        full_locants.append(primary)
        if is_compound:
            full_locants.append(max(position[a], position[b]))
            compound_count += 1
        (ene_citations if order == ENE_BOND_ORDER else yne_citations).append((primary, display))
    return ene_citations, yne_citations, compound_count, primary_locants, full_locants


def ring_bond_locant(position_of, bond_atoms, ring_size):
    """1-based ring locant of the lower-numbered atom of a ring bond, under
    `position_of` (atom index -> 1-based position), wrapping so the bond
    between positions 1 and `ring_size` reports as `ring_size`."""
    pa, pb = position_of[bond_atoms[0]], position_of[bond_atoms[1]]
    return ring_size if {pa, pb} == {1, ring_size} else min(pa, pb)


def ring_bond_locants(position_of, bonds, ring_size):
    """(ene_locants, yne_locants), both sorted, for every ring C=C/C#C bond
    under this ring numbering."""
    ene, yne = [], []
    for a, b, order in bonds:
        if a not in position_of or b not in position_of:
            continue
        locant = ring_bond_locant(position_of, (a, b), ring_size)
        (ene if order == ENE_BOND_ORDER else yne).append(locant)
    return sorted(ene), sorted(yne)


def specified_stereocenters(mol):
    """Scan `mol` for stereo elements (RDKit's `Chem.FindPotentialStereo`).
    If there are none at all, or every one present is left unspecified
    (no `@`/`@@`/E-Z bond marker anywhere in the input), return None --
    the caller should proceed exactly as if stereochemistry weren't a
    factor. This matches both this project's existing, long-standing
    behavior and the real IUPAC/PubChem convention for a name that doesn't
    specify configuration at all (e.g. this project's own pre-existing,
    PubChem-verified '2-fluorobutan-1-ol', whose C2 is a genuine but
    undrawn stereocenter) -- it is deliberately NOT treated as a new
    rejection case, unlike a *partially* specified molecule (see below).

    If one or more stereo elements are *specified*, no unspecified one
    alongside them, and every specified one is a tetrahedral atom
    stereocenter (no double-bond E/Z element), return a list of
    (atom_idx, "R" or "S") pairs, one per specified stereocenter, via
    `rdCIPLabeler` -- in the same order `Chem.FindPotentialStereo` reports
    them, not yet locant-sorted (the caller only learns each atom's locant
    once its own chain/ring numbering is fixed; see P-91.3's ascending-
    locant citation order in the module docstring). Otherwise -- a
    specified stereocenter mixed with an unspecified one, or any
    double-bond E/Z stereo element -- raise `UnsupportedStructure`
    explicitly (P-93 double-bond stereo combined with P-92 tetrahedral
    stereo, and any partially-specified molecule, are both out of
    scope)."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if any(e.type != Chem.StereoType.Atom_Tetrahedral for e in specified):
        raise UnsupportedStructure(
            "stereochemistry beyond one or more specified tetrahedral "
            "stereocenters (with no unspecified one alongside them) is not "
            "supported yet (a specified stereocenter mixed with an "
            "unspecified one, or any C=C/C#N double-bond E/Z stereo -- see "
            "P-92/P-93)"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        atom_idx = element.centeredOn
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp("_CIPCode"):
            raise UnsupportedStructure(
                "could not determine a CIP R/S label for this stereocenter"
            )
        code = atom.GetProp("_CIPCode")
        if code not in ("R", "S", "r", "s"):
            raise UnsupportedStructure("could not determine a CIP R/S/r/s label for this stereocenter")
        labels.append((atom_idx, code))
    return labels


def stereo_locants_prefix(stereo, position_of):
    """Format `specified_stereocenters`' `(atom_idx, "R"/"S")` pairs into a
    `"(<locant><R/S>,...)-"` prefix (P-91.3's ascending-locant citation
    order), once the caller's own numbering has fixed each atom's locant
    via `position_of` (`{atom_idx: locant}`). Raises `UnsupportedStructure`
    if any stereocenter isn't a key of `position_of` -- P-92: a
    stereocenter outside the caller's own principal chain/ring skeleton
    (e.g. on a substituent branch) is out of scope everywhere this is
    used. Shared by `_acyclic.py`, `_spiro.py`, and `_spiro_heteroatom.py`
    (each otherwise had its own copy of this exact formatting logic)."""
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain/ring skeleton is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    if len(position_of) == 1:
        # P-93.5: the one skeletal atom of a mononuclear parent needs no locant
        return f"({labels[0][1]})-"
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-"


def heteroatom_stereo_prefix(mol, heteroatom_idx):
    """None if `mol` has no specified stereocenter at all (the caller
    proceeds exactly as before). "(R)-"/"(S)-" if the molecule's sole
    specified stereocenter is `heteroatom_idx` itself -- P-93.3.3.2/
    P-93.3.4.1: a trigonal pyramidal center (a sulfoxide/sulfinyl/
    seleninyl sulfur or selenium, one lone pair standing in for a phantom
    low-priority ligand) is assigned an ordinary R/S descriptor "in the
    manner described for tetrahedral stereogenic centers," cited with no
    locant when it's the molecule's only stereocenter (Blue Book worked
    examples: "(S)-(methanesulfinyl)ethane", "ethyl (R)-4-nitrobenzene-
    1-sulfinate"). Raises `UnsupportedStructure` for every other
    specified-stereocenter shape (a chain/ring carbon, more than one
    specified center, ...) -- still out of scope, same rejection every
    caller already had before this existed."""
    stereo = specified_stereocenters(mol)
    if stereo is None:
        return None
    if len(stereo) == 1 and stereo[0][0] == heteroatom_idx:
        return f"({stereo[0][1]})-"
    raise UnsupportedStructure(
        "a specified stereocenter other than the sole heteroatom center "
        "is not supported yet (see P-92/P-93.3.4.1)"
    )


def specified_double_bond_stereo(mol):
    """`Bond_Double` analogue of `specified_stereocenters` above (see its
    docstring for the same unspecified-vs-rejected reasoning). None if
    there's no *specified* double-bond stereo element at all -- covering
    both a non-stereogenic double bond (e.g. `C=C(C)C`, where
    `FindPotentialStereo` doesn't report an element at all) and one left
    unspecified in the input (plain `C=C`, no `/`/`\\`) -- so the caller
    proceeds exactly as before (no E/Z prefix) in either case.

    If one or more stereo elements are *specified*, no unspecified one
    alongside them, and every specified one is a `Bond_Double` element (no
    tetrahedral atom stereocenter), return a list of
    (bond_idx, "E" or "Z") pairs via `rdCIPLabeler` -- in the same order
    `Chem.FindPotentialStereo` reports them, not yet locant-sorted (the
    caller only learns each bond's locant once its own chain numbering is
    fixed). Otherwise -- a specified element mixed with an unspecified one,
    or a tetrahedral stereocenter -- raise `UnsupportedStructure`
    explicitly."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if any(e.type != Chem.StereoType.Bond_Double for e in specified):
        raise UnsupportedStructure(
            "stereochemistry beyond one or more specified C=C double-bond "
            "E/Z elements (with no unspecified one alongside them) is not "
            "supported yet (a specified element mixed with an unspecified "
            "one, or a tetrahedral stereocenter -- see P-92/P-93)"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        bond_idx = element.centeredOn
        bond = mol.GetBondWithIdx(bond_idx)
        if not bond.HasProp("_CIPCode"):
            raise UnsupportedStructure(
                "could not determine a CIP E/Z label for this double bond"
            )
        labels.append((bond_idx, bond.GetProp("_CIPCode")))
    return labels


def stereo_element_atoms(mol, stereo):
    """Every atom a `specified_stereo_elements` result touches (a bond
    element contributes both ends), for requiring them all on the chain."""
    atoms = []
    for kind, idx, _ in stereo:
        if kind == "atom":
            atoms.append(idx)
        else:
            bond = mol.GetBondWithIdx(idx)
            atoms.extend((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    return atoms


def specified_stereo_elements(mol):
    """Like `specified_stereocenters`/`specified_double_bond_stereo`, but
    allows a specified tetrahedral stereocenter and a specified C=C
    double-bond E/Z element to coexist in the same molecule: P-91.3's own
    worked example, '(2Z,5R,7E)-nona-2,7-dien-5-ol (PIN)', cites both kinds
    together in one locant-ascending group, confirmed directly from the
    primary source text.

    Returns None if there is no specified stereo element at all (same
    unspecified-vs-rejected policy as the two functions above). Otherwise
    returns a list of ("atom" or "bond", atom_or_bond_idx, "R"/"S"/"E"/"Z")
    triples, one per specified element, not yet locant-sorted -- the
    caller only learns each element's locant once its own chain numbering
    is fixed. Raises `UnsupportedStructure` for an unspecified element
    mixed in, a pseudoasymmetric (lowercase r/s) stereocenter, or a CIP
    label RDKit couldn't determine."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    rdCIPLabeler.AssignCIPLabels(mol)
    labels = []
    for element in specified:
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            atom = mol.GetAtomWithIdx(element.centeredOn)
            if not atom.HasProp("_CIPCode"):
                raise UnsupportedStructure(
                    "could not determine a CIP R/S label for this stereocenter"
                )
            code = atom.GetProp("_CIPCode")
            if code not in ("R", "S", "r", "s"):
                raise UnsupportedStructure("could not determine a CIP R/S/r/s label for this stereocenter")
            labels.append(("atom", element.centeredOn, code))
        elif element.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(element.centeredOn)
            if not bond.HasProp("_CIPCode"):
                raise UnsupportedStructure(
                    "could not determine a CIP E/Z label for this double bond"
                )
            labels.append(("bond", element.centeredOn, bond.GetProp("_CIPCode")))
        else:
            raise UnsupportedStructure(
                "stereochemistry beyond a tetrahedral R/S stereocenter or a "
                "C=C double-bond E/Z element is not supported yet (see "
                "P-92/P-93)"
            )
    return labels


def stereo_locant_rank(mol, stereo, position_of):
    if not stereo:
        return ((), ())
    labels = []
    for item in stereo:
        kind, idx, code = item if len(item) == 3 else ("atom", *item)
        if kind == "atom":
            labels.append((position_of.get(idx, 0), code))
        else:
            bond = mol.GetBondWithIdx(idx)
            labels.append((min(position_of.get(bond.GetBeginAtomIdx(), 0), position_of.get(bond.GetEndAtomIdx(), 0)), code))
    return tuple(tuple(sorted(loc for loc, code in labels if code in wanted)) for wanted in (("Z",), ("R", "r")))


_SUPERSCRIPT_CHARS = str.maketrans("0123456789abcdefgh", "⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ")


def superscript_locant(primary, local):
    """Locant of atom `local` of the ring or amplificant numbered `primary`, the atom locant raised (P-26.4, P-28.3.1)."""
    return f"{primary}{str(local).translate(_SUPERSCRIPT_CHARS)}"


_STANDARD_BONDING = {
    "N": 3, "P": 3, "As": 3, "Sb": 3, "Bi": 3, "B": 3, "Al": 3, "Ga": 3, "In": 3, "Tl": 3,
    "O": 2, "S": 2, "Se": 2, "Te": 2, "Si": 4, "Ge": 4, "Sn": 4, "Pb": 4,
}


# Ion names cite the bonding numbers of their centres themselves, on a neutral parent named by a nested call.
CITE_SKELETAL_LAMBDA = ContextVar("cite_skeletal_lambda", default=True)


def nonstandard_bonding(atom):
    """The bonding number n of a neutral, non-aromatic skeletal heteroatom above its standard one (P-14.1), for the λn
    convention (P-15.4.1.3, P-22.2.7); charged atoms and oxo-bearing ones are named by the ion and heterone rules."""
    standard = _STANDARD_BONDING.get(atom.GetSymbol())
    if standard is None or atom.GetFormalCharge() or atom.GetIsAromatic() or not CITE_SKELETAL_LAMBDA.get():
        return None
    if any(atom.HasProp(prop) for prop in ("_anion", "_anion_word", "_anion_lambda", "_ring_cation_centre")):
        return None
    if any(b.GetBondTypeAsDouble() > 1.0 and b.GetOtherAtom(atom).GetDegree() == 1 for b in atom.GetBonds()):
        return None
    valence = int(atom.GetTotalValence())
    return valence if valence > standard else None


def lambda_cited(mol, atom_idx, locant):
    """`locant` with the λn mark of a nonstandard bonding number on skeletal atom `atom_idx` (P-14.1, P-15.4.1.3)."""
    bonding = nonstandard_bonding(mol.GetAtomWithIdx(atom_idx))
    return f"{locant}\u03bb{bonding}" if bonding else str(locant)


def assembly_join(prefix, core):
    """Prefixes directly precede an enclosing bracket of an assembly name and are hyphenated from a locant (P-16.5)."""
    if not prefix:
        return core
    return prefix + core if core[0] == "[" else f"{prefix}-{core}"


def alphanumerical_name_key(name):
    """P-14.5: whole names compare by their letters alone, locants and enclosing marks ignored."""
    return re.sub(r"[^a-z]", "", name.lower())
