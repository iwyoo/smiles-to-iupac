"""Naming of a common alpha-amino carboxylic acid using its P-103 retained
name and L/D stereodescriptor, for the side chains below, per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained names, keyed by a canonical
  fragment SMILES of the side chain itself (see `_SIDE_CHAIN_TABLE`)
  rather than a hand-coded graph walk per shape -- a growing per-shape
  `if`/`elif` chain (this module's own earlier form, #1057/#1068/#1070)
  is a narrow-enumeration smell even when each addition individually
  looks justified; a canonical-fragment table lookup is the general
  mechanism that actually covers the whole family in one place. New side
  chains are added by extending the table, not by writing a new
  structural walk.
- P-103.1.3.1: 'L' corresponds to the CIP 'S' configuration at the
  alpha-carbon for every common amino acid *except* cysteine, whose
  side-chain sulfur outranks the ring-ward carbon in CIP priority and so
  flips the correspondence to L='R'/D='S' -- confirmed against
  cysteine's real PubChem L-/D- structures (CIDs 92851, 5862). 'D'
  corresponds to 'R' (or 'S' for cysteine). Confirmed against real
  PubChem L-/D- structures for alanine/valine/leucine/serine/aspartic
  acid/glutamic acid during scoping (#1056 M1 steps 1-3) -- no exception
  for aspartic/glutamic acid, the ordinary S=L/R=D rule holds. Glycine's
  alpha-carbon bears two hydrogens, so it's never a stereocenter and gets
  no L-/D- prefix at all.

Scope, deliberately narrow (table entries only): exactly the side chains
in `_SIDE_CHAIN_TABLE`, each already correctly recognized as a plain
alpha-amino-acid backbone by `_carboxylic_acid_amine.py`'s shared
detection -- this module only intercepts those specific side chains
ahead of that module's own generic (CIP-only, no-retained-name) fallback
in `core.py`'s dispatch order, unchanged for every other amino acid or
amine/acid combination. The backbone-anchor check in `_match` (see its
own docstring) tolerates a side chain with its own extra nitrogen
(lysine, arginine) or its own ring (phenylalanine, tyrosine, tryptophan)
as long as the table recognizes it -- an *unrecognized* multi-amine or
ring-bearing shape still falls through unmatched, same as any other
unrecognized side chain. Isoleucine/threonine (a second side-chain
stereocenter and 'allo' complexity) remain out of scope regardless, since
this whole mechanism only ever looks up *one* alpha-stereocenter.
"""

import itertools
from contextvars import ContextVar

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure, adjacency, find_primary_amines

SYSTEMATIC_ACID_PROBE = ContextVar("SYSTEMATIC_ACID_PROBE", default=False)

_ALPHA_TO_LD = {"S": "L", "R": "D"}
_ALPHA_TO_LD_CYSTEINE = {"R": "L", "S": "D"}  # P-103.1.3.1's stated exception
_SULFUR_ON_C3 = {"cysteine", "cysteic acid"}

_SIDE_CHAIN_SMILES = {
    "alanine": "*C",
    "valine": "*C(C)C",
    "isoleucine": "*C(C)CC",
    "threonine": "*C(C)O",
    "leucine": "*CC(C)C",
    "serine": "*CO",
    "cysteine": "*CS",
    "aspartic acid": "*CC(=O)O",
    "glutamic acid": "*CCC(=O)O",
    "asparagine": "*CC(=O)N",
    "glutamine": "*CCC(=O)N",
    "methionine": "*CCSC",
    "lysine": "*CCCCN",
    "arginine": "*CCCNC(N)=N",
    "phenylalanine": "*Cc1ccccc1",
    "tyrosine": "*Cc1ccc(O)cc1",
    "tryptophan": "*Cc1c[nH]c2ccccc12",
    "ornithine": "*CCCN",
    "allysine": "*CCCC=O",
    "citrulline": "*CCCNC(N)=O",
    "cysteic acid": "*CS(=O)(=O)O",
    "homocysteine": "*CCS",
    "homoserine": "*CCO",
    "dopa": "*Cc1ccc(O)c(O)c1",
}
# Canonicalized at import time (rather than hardcoding the already-
# canonical strings above) so a future RDKit version's canonicalization
# doesn't silently desync the table from what `_side_chain_fragment`
# actually produces.
_SIDE_CHAIN_TABLE = {Chem.CanonSmiles(smiles): name for name, smiles in _SIDE_CHAIN_SMILES.items()}


def _is_plain_carbon(atom):
    return atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetIsotope() == 0


def _is_terminal_carboxylic_acid_carbon(mol, graph, atom_idx, coming_from):
    """True for a terminal -COOH carbon (bonded only to `coming_from` and
    its own carbonyl/hydroxyl oxygens) hanging off `coming_from` -- used
    to find the *main-chain* acid carbon among the alpha-carbon's own
    neighbors, not for side-chain recognition (aspartic/glutamic acid's
    second, side-chain-terminal carboxylic acid is matched via
    `_SIDE_CHAIN_TABLE` instead, like every other side chain)."""
    atom = mol.GetAtomWithIdx(atom_idx)
    if not _is_plain_carbon(atom):
        return False
    neighbors = [n for n in graph[atom_idx] if n != coming_from]
    if len(neighbors) != 2:
        return False
    oxygens = [n for n in neighbors if mol.GetAtomWithIdx(n).GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return False
    carbonyls = [
        o
        for o in oxygens
        if mol.GetAtomWithIdx(o).GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom_idx, o).GetBondTypeAsDouble() == 2.0
    ]
    hydroxyls = [
        o
        for o in oxygens
        if mol.GetAtomWithIdx(o).GetDegree() == 1
        and mol.GetBondBetweenAtoms(atom_idx, o).GetBondTypeAsDouble() == 1.0
        and mol.GetAtomWithIdx(o).GetTotalNumHs() == 1
    ]
    return len(carbonyls) == 1 and len(hydroxyls) == 1


def _side_chain_atom_indices(graph, alpha_carbon, side_root):
    """Every atom reachable from `side_root` without passing back through
    `alpha_carbon` -- the side chain's own atoms, in the original
    molecule's indexing (unlike `_side_chain_fragment`'s reindexed copy).
    Used only to scope the stereo check below to the side chain, not for
    side-chain identification itself (that's `_side_chain_fragment` +
    `_SIDE_CHAIN_TABLE`)."""
    seen = {side_root}
    stack = [side_root]
    while stack:
        current = stack.pop()
        for neighbor in graph[current]:
            if neighbor != alpha_carbon and neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return seen


def _alpha_stereo_label(mol, alpha_carbon, side_chain_atoms):
    """The alpha-carbon's own CIP label ('R'/'S'), or None if it has no
    specified stereocenter -- ignoring any potential stereo element
    entirely confined to `side_chain_atoms` (e.g. arginine's guanidino
    C=N, which RDKit reports as a genuine potential E/Z element even
    though it's never configurationally specified in practice -- the
    tautomeric/resonance exchange across the guanidino nitrogens makes it
    non-configurational, not something RDKit's graph-only stereo
    perception can infer on its own). A stereo element outside the
    recognized side chain (there is at most the alpha-carbon itself,
    given this module's own narrow scope) still applies the ordinary
    all-specified-or-none rule, same as `_common.specified_stereocenters`
    -- this is a narrower, side-chain-aware version of that shared
    helper, not a replacement for it elsewhere."""
    elements = Chem.FindPotentialStereo(mol)
    relevant = []
    for element in elements:
        if element.type == Chem.StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(element.centeredOn)
            if bond.GetBeginAtomIdx() in side_chain_atoms and bond.GetEndAtomIdx() in side_chain_atoms:
                continue
        elif element.type == Chem.StereoType.Atom_Tetrahedral and element.centeredOn in side_chain_atoms:
            continue
        relevant.append(element)
    specified = [e for e in relevant if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if len(specified) != len(relevant) or any(e.type != Chem.StereoType.Atom_Tetrahedral for e in specified):
        raise UnsupportedStructure(
            "stereochemistry beyond the alpha-carbon's own specified "
            "tetrahedral stereocenter (with no unspecified one alongside "
            "it, outside the recognized side chain) is not supported yet"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    for element in specified:
        if element.centeredOn == alpha_carbon:
            return mol.GetAtomWithIdx(alpha_carbon).GetPropsAsDict().get("_CIPCode")
    return None


def _side_chain_fragment(mol, alpha_carbon, other_alpha_neighbors):
    """The side-chain fragment hanging off `alpha_carbon`, as an RDKit Mol
    plus its canonical SMILES -- the cut bond to `alpha_carbon` is marked
    by a dummy atom (`*`) so the fragment's canonical form is stable
    regardless of what's on the other side of that bond (the amine/acid
    this module already excludes). `(None, None)` if `alpha_carbon` isn't
    itself a plain, unbranched attachment point (shouldn't happen given
    how callers use this, but keeps this function total)."""
    rw = Chem.RWMol(mol)
    for neighbor in other_alpha_neighbors:
        rw.RemoveBond(alpha_carbon, neighbor)
    dummy = rw.GetAtomWithIdx(alpha_carbon)
    dummy.SetAtomicNum(0)
    dummy.SetFormalCharge(0)
    dummy.SetNoImplicit(True)
    dummy.SetNumExplicitHs(0)
    frag_mol = rw.GetMol()
    Chem.SanitizeMol(frag_mol)
    for frag in Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True):
        if any(a.GetAtomicNum() == 0 for a in frag.GetAtoms()):
            Chem.RemoveStereochemistry(frag)
            return frag, Chem.MolToSmiles(frag)
    return None, None


def _backbone_candidate(mol, graph, amine_n):
    """(alpha_carbon, acid_carbon_idx) if `amine_n` directly anchors a
    plain alpha-amino-acid backbone (bonded to a plain carbon that itself
    has exactly one terminal-carboxylic-acid-carbon neighbor), else None.
    Doesn't care whether `mol` has other amines or rings elsewhere -- a
    side chain with its own extra nitrogen (lysine, arginine) or its own
    ring (phenylalanine, tyrosine, tryptophan) is fine, as long as
    exactly one amine in the whole molecule anchors a backbone this way;
    that side chain is handled entirely by `_side_chain_fragment` +
    `_SIDE_CHAIN_TABLE` below, not by this function."""
    neighbors = graph[amine_n]
    if len(neighbors) != 1:
        return None
    (alpha_carbon,) = neighbors
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if not _is_plain_carbon(alpha_atom):
        return None
    acid_candidates = [
        n
        for n in graph[alpha_carbon]
        if n != amine_n and _is_terminal_carboxylic_acid_carbon(mol, graph, n, alpha_carbon)
    ]
    if len(acid_candidates) != 1:
        return None
    return alpha_carbon, acid_candidates[0]


def _match_skeleton(mol):
    """(retained_name, alpha_carbon_idx_or_None, side_chain_atoms) if
    `mol` is a plain alpha-amino acid whose side chain matches a
    `_SIDE_CHAIN_TABLE` entry, else None. `alpha_carbon_idx` is None only
    for glycine (no stereocenter to look up an L/D descriptor for, and no
    side chain either -- `side_chain_atoms` is empty in that case).

    Anchors the backbone by finding which amine (of possibly several --
    lysine/arginine's side chains carry their own extra nitrogen) is
    directly bonded to a plain carbon that itself neighbors a terminal
    carboxylic acid carbon, rather than requiring the *whole molecule*
    have exactly one amine and no ring anywhere (confirmed this session:
    that whole-molecule check rejected lysine/arginine/asparagine/
    glutamine outright over their side chain's own nitrogen, and
    phenylalanine/tyrosine/tryptophan over their side chain's own ring,
    long before side-chain matching was ever reached)."""
    graph = adjacency(mol)
    amines = find_primary_amines(mol)
    backbones = []
    for amine_n in amines:
        found = _backbone_candidate(mol, graph, amine_n)
        if found is not None:
            backbones.append((amine_n, *found))
    if len(backbones) != 1:
        return None
    amine_n, alpha_carbon, acid_carbon_idx = backbones[0]
    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)

    side_neighbors = [n for n in graph[alpha_carbon] if n not in (amine_n, acid_carbon_idx)]

    if not side_neighbors:
        if alpha_atom.GetTotalNumHs() != 2:
            return None
        # amine N + acid C + its 2 oxygens + alpha C, nothing else anywhere
        # in the molecule (rejects e.g. a disconnected salt fragment).
        if mol.GetNumAtoms() != 5:
            return None
        return "glycine", None, frozenset(), None

    if len(side_neighbors) != 1 or alpha_atom.GetTotalNumHs() != 1:
        return None
    (side_root,) = side_neighbors

    frag, frag_smiles = _side_chain_fragment(mol, alpha_carbon, (amine_n, acid_carbon_idx))
    if frag_smiles is None:
        return None
    name = _SIDE_CHAIN_TABLE.get(frag_smiles)
    if name is None:
        phosphoryl = _serine_phosphoryl(mol, graph, alpha_carbon, side_root)
        if phosphoryl is None or mol.GetNumAtoms() != 5 + len(_side_chain_atom_indices(graph, alpha_carbon, side_root)):
            return None
        return "serine", alpha_carbon, _side_chain_atom_indices(graph, alpha_carbon, side_root), phosphoryl
    # amine N + acid C + its 2 oxygens + alpha C + the side chain's own
    # atoms (frag includes a dummy atom standing in for alpha_carbon, so
    # its real atom count is frag.GetNumAtoms() - 1) - same disconnected-
    # fragment guard as glycine's case above, generalized via the
    # fragment's own atom count instead of a per-name lookup table.
    if mol.GetNumAtoms() != 4 + frag.GetNumAtoms():
        return None
    return name, alpha_carbon, _side_chain_atom_indices(graph, alpha_carbon, side_root), None


# P-103.2.1: N, O, S locants, with the numbering of the alpha-amino acid where several nitrogens exist
_SIDE_CHAIN_SUBSTITUTION_SITE = {
    "serine": "O",
    "cysteine": "S",
    "tyrosine": "O",
    "lysine": "N",
    "asparagine": "N",
    "glutamine": "N",
    "arginine": "N",
    "ornithine": "N",
    "homoserine": "O",
    "homocysteine": "S",
}
# P-103.1.3.2.2: CIP label of C-3 in the L (alpha S) series; 'allo' inverts it, and the D series mirrors both centres
_BETA_NATURAL = {"isoleucine": "S", "threonine": "R"}
_ALPHA_NITROGEN_LOCANT = {"lysine": "N2", "asparagine": "N2", "glutamine": "N2", "arginine": "Nα", "ornithine": "N2", "citrulline": "N2"}
_SIDE_NITROGEN_LOCANT = {"lysine": "N6", "asparagine": "N4", "glutamine": "N5", "ornithine": "N5"}
_MAX_CUT_CANDIDATES = 8
_ALPHA_AMINO_ACID = Chem.MolFromSmarts("[NX3;+0][CX4][CX3](=O)[OX2H1]")


def _cut_candidates(mol):
    bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        for hetero, other in (
            (bond.GetBeginAtom(), bond.GetEndAtom()),
            (bond.GetEndAtom(), bond.GetBeginAtom()),
        ):
            if hetero.GetAtomicNum() in (7, 8, 16) and hetero.GetFormalCharge() == 0 and not hetero.IsInRing():
                bonds.append((hetero.GetIdx(), other.GetIdx()))
    return bonds


def _skeleton_after_cuts(mol, cuts):
    """The fragment left by cutting each (site, root) bond that is a table amino acid, with each atom's
    original index in `_orig`; else None."""
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetIntProp("_orig", atom.GetIdx())
    for site, root in cuts:
        rw.RemoveBond(site, root)
        atom = rw.GetAtomWithIdx(site)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    cut_mol = rw.GetMol()
    Chem.SanitizeMol(cut_mol)
    for frag in Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False):
        if _match_skeleton(frag) is not None:
            return frag
    return None


def _match(mol):
    """`_match_skeleton`'s tuple plus ((locant, site, root), ...) for the substituents on the amino, thiol or
    hydroxy group that must be cut to reach a table amino acid (P-103.2.3)."""
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    found = _match_skeleton(mol)
    if found is not None:
        return (*found, ())
    if not mol.HasSubstructMatch(_ALPHA_AMINO_ACID):
        return None
    candidates = _cut_candidates(mol)
    if not candidates or len(candidates) > _MAX_CUT_CANDIDATES:
        return None
    best = None
    for size in range(1, len(candidates) + 1):
        for cuts in itertools.combinations(candidates, size):
            if len({root for _, root in cuts}) != size:
                continue
            frag = _skeleton_after_cuts(mol, cuts)
            if frag is None:
                continue
            name, alpha, side_chain, _ = _match_skeleton(frag)
            if alpha is None:
                continue
            to_orig = {a.GetIdx(): a.GetIntProp("_orig") for a in frag.GetAtoms()}
            alpha = to_orig[alpha]
            side_chain = frozenset(to_orig[i] for i in side_chain)
            amine = next(n.GetIdx() for n in mol.GetAtomWithIdx(alpha).GetNeighbors() if n.GetAtomicNum() == 7)
            allowed = {amine: _ALPHA_NITROGEN_LOCANT.get(name, "N")}
            for site in side_chain:
                element = mol.GetAtomWithIdx(site).GetSymbol()
                if element == _SIDE_CHAIN_SUBSTITUTION_SITE.get(name):
                    allowed[site] = _SIDE_NITROGEN_LOCANT.get(name, element)
            if name == "arginine":
                allowed.update(_arginine_locants(mol, alpha, side_chain, {site for site, _ in cuts}))
            if any(site not in allowed for site, _ in cuts):
                continue
            if len(frag.GetSubstructMatches(_ALPHA_AMINO_ACID)) != len(mol.GetSubstructMatches(_ALPHA_AMINO_ACID)):
                continue
            found = (name, alpha, side_chain, None, tuple((allowed[site], site, root) for site, root in cuts))
            if best is None or frag.GetNumAtoms() > best[0]:
                best = (frag.GetNumAtoms(), found)
    return best[1] if best else None


def _arginine_locants(mol, alpha, side_chain, substituted):
    """Nδ for the nitrogen on the chain, Nω (then Nω′ for a second substituted one) for the terminal guanidine
    nitrogens (P-103.2.1), told apart by their distance from the alpha carbon (4 and 6 bonds)."""
    distance = Chem.GetDistanceMatrix(mol)
    locants = {}
    marks = iter(("Nω", "Nω′"))
    for a in sorted(a for a in side_chain if mol.GetAtomWithIdx(a).GetAtomicNum() == 7):
        if distance[alpha][a] == 4:
            locants[a] = "Nδ"
        elif distance[alpha][a] == 6:
            locants[a] = next(marks) if a in substituted else "Nω"
    return locants


def _serine_phosphoryl(mol, graph, alpha_carbon, side_root):
    """(oxygen, phosphorus) of a serine side chain CH2-O-P(=O)(O..)(O..), else None."""
    side = mol.GetAtomWithIdx(side_root)
    if not _is_plain_carbon(side) or side.GetTotalNumHs() != 2:
        return None
    oxygens = [n for n in graph[side_root] if n != alpha_carbon]
    if len(oxygens) != 1 or mol.GetAtomWithIdx(oxygens[0]).GetAtomicNum() != 8:
        return None
    phosphorus = [n for n in graph[oxygens[0]] if n != side_root]
    if len(phosphorus) != 1 or mol.GetAtomWithIdx(phosphorus[0]).GetAtomicNum() != 15:
        return None
    return oxygens[0], phosphorus[0]


def _phosphoryl_group_name(mol, graph, oxygen, phosphorus, side_chain_atoms):
    from ._common import halogen_substituents
    from ._functional_prefixes import _enclose, functional_names
    from ._substituents import BRANCH_STEREO

    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    atoms = {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in probe.GetAtoms()
        if a.GetIdx() in side_chain_atoms and a.HasProp("_CIPCode")
    }
    bonds = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
        for b in probe.GetBonds()
        if b.HasProp("_CIPCode") and b.GetBeginAtomIdx() in side_chain_atoms
    }
    context = {"atoms": atoms, "bonds": bonds, "used": set()}
    token = BRANCH_STEREO.set(context)
    try:
        named, _, _ = functional_names(mol, graph, [(oxygen, phosphorus)], {oxygen}, halogen_substituents(mol))
    finally:
        BRANCH_STEREO.reset(token)
    if phosphorus not in named or any(("atom", a) not in context["used"] for a in atoms):
        raise UnsupportedStructure("a stereocentre of the phosphoryl substituent is not cited by any supported name")
    name = named[phosphorus][0]
    if name == "phosphono":
        return name
    return _enclose(name)


def _beta_atom(mol, alpha_carbon, side_chain_atoms):
    return next(n.GetIdx() for n in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors() if n.GetIdx() in side_chain_atoms)


def _beta_specified(mol, beta):
    return any(
        e.type == Chem.StereoType.Atom_Tetrahedral and e.centeredOn == beta and e.specified == Chem.StereoSpecified.Specified
        for e in Chem.FindPotentialStereo(mol)
    )


def has_amino_acid_shape(mol) -> bool:
    found = _match(mol)
    if found is None:
        return False
    name, alpha_carbon, side_chain_atoms = found[:3]
    if name in _BETA_NATURAL and alpha_carbon is not None:
        beta = _beta_specified(mol, _beta_atom(mol, alpha_carbon, side_chain_atoms))
        try:
            alpha = _alpha_stereo_label(mol, alpha_carbon, side_chain_atoms) is not None
        except UnsupportedStructure:
            return True
        return beta == alpha
    return True


def _substituent_prefixes(mol, graph, substituents):
    from ._common import halogen_substituents
    from ._substituents import BRANCH_STEREO, format_substituent_prefixes, name_branch

    halogens = halogen_substituents(mol)
    aromatic = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    inside = set()
    for _, site, root in substituents:
        stack = [root]
        while stack:
            node = stack.pop()
            if node not in inside:
                inside.add(node)
                stack.extend(n for n in graph[node] if n != site and n not in inside)
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    atoms = {a.GetIdx(): a.GetProp("_CIPCode") for a in probe.GetAtoms() if a.GetIdx() in inside and a.HasProp("_CIPCode")}
    bonds = {
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()): b.GetProp("_CIPCode")
        for b in probe.GetBonds()
        if b.HasProp("_CIPCode") and b.GetBeginAtomIdx() in inside and b.GetEndAtomIdx() in inside
    }
    context = {"atoms": atoms, "bonds": bonds, "used": set()}
    grouped = {}
    token = BRANCH_STEREO.set(context)
    try:
        for locant, site, root in substituents:
            name, compound = name_branch(graph, root, site, halogens, aromatic, mol)
            grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    finally:
        BRANCH_STEREO.reset(token)
    if any(("atom", a) not in context["used"] for a in atoms) or any(("bond", b) not in context["used"] for b in bonds):
        raise UnsupportedStructure("a stereo element of an N-substituent is not cited by any supported name")
    return format_substituent_prefixes(grouped)


def name_amino_acid(mol) -> str:
    name, alpha_carbon, side_chain_atoms, phosphoryl, substituents = _match(mol)
    if alpha_carbon is None:
        return name
    label = _alpha_stereo_label(mol, alpha_carbon, side_chain_atoms)
    mapping = _ALPHA_TO_LD_CYSTEINE if name in _SULFUR_ON_C3 else _ALPHA_TO_LD
    if name in _BETA_NATURAL and label:
        beta = mol.GetAtomWithIdx(_beta_atom(mol, alpha_carbon, side_chain_atoms)).GetProp("_CIPCode")
        natural = _BETA_NATURAL[name] if label == "S" else {"S": "R", "R": "S"}[_BETA_NATURAL[name]]
        if beta != natural:
            name = "allo" + name
    prefix = ""
    if phosphoryl is not None:
        prefix = f"O-{_phosphoryl_group_name(mol, adjacency(mol), *phosphoryl, side_chain_atoms)}-"
    elif substituents:
        prefix = _substituent_prefixes(mol, adjacency(mol), substituents) + "-"
    return f"{prefix}{mapping[label]}-{name}" if label else f"{prefix}{name}"
