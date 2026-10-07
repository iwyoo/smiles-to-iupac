"""Skeletal replacement ('a' prefix) naming of a saturated linear polyspiro
ring system (`_polyspiro.py`'s `find_linear_polyspiro_chain` shape)
containing exactly one ring heteroatom, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-24.2.4 / P-24.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf)
  and the spiro-nomenclature companion document's SP-1.8.1
  (https://iupac.qmul.ac.uk/spiro/sp16p.html: "Low locants are allocated
  for heteroatoms indicated by replacement ['a'] terms as a set"): the
  same rule `_spiro_heteroatom.py` already applies to a monospiro system
  extends here unchanged -- the von Baeyer spiro descriptor and the
  numbering it drives (`_polyspiro.py`'s spiro-locant-set, then
  descriptor-citation-order tie-breaks) are never modified by introducing
  a heteroatom; the heteroatom's own locant is minimized only among
  whichever numbering choices those structural rules still leave open
  (chain direction, arc choice, ring-walk direction), ahead of
  substituent locants. Confirmed against PubChem's own computed
  IUPACName for a from-scratch-built one-oxygen dispiro skeleton with no
  such remaining choice (the heteroatom sits alone on the descriptor's
  length-1 internal arc, directly bonded to both spiro atoms, so its
  locant is structurally forced either way): '6-oxadispiro[4.1.4^7.2^5]
  tridecane' (PubChem CID 45098582, built from the same all-carbon
  skeleton as this project's own dispiro[4.1.4^7.2^5]tridecane test).
  For shapes where a genuine numbering choice *does* remain open (e.g. a
  terminal ring whose two non-spiro atoms are both bonded straight to the
  spiro atom, so either can be visited first), PubChem's own generated
  names for two from-scratch structures of that shape (ConnectivitySMILES
  matching CID 114820437 and CID 134451841) instead give the heteroatom
  the *higher* of the two available locants, contradicting SP-1.8.1's own
  quoted text above -- treated here as a PubChem generator limitation for
  this uncommon shape rather than followed, so this module always
  minimizes the heteroatom locant per SP-1.8.1 regardless.
- The spiro atoms themselves must stay carbon, same reasoning as
  `_spiro_heteroatom.py`: P-24.2.2's descriptor numbers each spiro atom
  using all four of its ring bonds, a valence no neutral O/N/S atom can
  supply.

Explicitly out of scope (raise `UnsupportedStructure`):
- More than one ring heteroatom (of any kind), a heteroatom at a spiro
  atom itself, or a heteroatom outside the ring system entirely.
- Any heteroatom other than O, N, or S.
- The branched polyspiro shape (`_polyspiro.find_branched_polyspiro_hub`) --
  this module's 1st-pass scope is the linear chain only.
- Unsaturation, charged/isotopic atoms, or anything else
  `_polyspiro.name_linear_polyspiro`'s own validation already rejects for
  the all-carbon case.
"""

from itertools import product

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    lambda_cited,
    non_single_bonds,
    specified_stereocenters,
    stereo_locants_prefix,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name, numerical_term
from ._polyspiro import _arc_choice_options, _build_sequence, _chain_direction_candidates
from ._substituents import format_substituent_prefixes, substituents_for_ring

_HETEROATOM_PREFIXES = {8: "oxa", 7: "aza", 16: "thia"}
_ALLOWED_ATOMIC_NUMS = {6, *_HETEROATOM_PREFIXES, *HALOGEN_PREFIXES}


def _chain_ring_atoms(mol, chain):
    ring_order, _ = chain
    atom_rings = mol.GetRingInfo().AtomRings()
    ring_atoms = set()
    for idx in ring_order:
        ring_atoms |= set(atom_rings[idx])
    return ring_atoms


def has_single_ring_heteroatom_shape(mol, chain) -> bool:
    """True iff `mol` has exactly one O/N/S atom, and it's a ring atom of
    `chain` (a `_polyspiro.find_linear_polyspiro_chain` result) other than
    one of its spiro atoms -- the shape this module accepts."""
    _, spiro_atoms = chain
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    if len(heteroatoms) != 1:
        return False
    idx = heteroatoms[0].GetIdx()
    return idx in _chain_ring_atoms(mol, chain) and idx not in spiro_atoms


def _candidate_key(parent, spiro_locants, descriptor, heteroatom_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = parent if not grouped else format_substituent_prefixes(grouped) + "-" + parent
    return (spiro_locants, descriptor, heteroatom_locant, locant_set, citation_locants, name)


def name_linear_polyspiro_heteroatom(mol, chain) -> str:
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure(
            "multi-fragment structures are not supported yet (see P-13.6, multiplicative nomenclature)"
        )
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a single ring O/N/S (P-24.2.4) and "
                "halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num in HALOGEN_PREFIXES and atom.GetDegree() != 1:
            raise UnsupportedStructure(
                "a halogen atom must be a monovalent substituent (P-35.2.1); "
                "polyvalent or bridging halogen structures are not supported"
            )
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated spiro ring systems are not supported yet (see "
            "P-31.1.5, unsaturated alicyclic spiro ring systems)"
        )

    ring_order, spiro_atoms = chain
    heteroatoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _HETEROATOM_PREFIXES]
    ring_atoms = _chain_ring_atoms(mol, chain)
    if len(heteroatoms) != 1 or heteroatoms[0] not in ring_atoms or heteroatoms[0] in spiro_atoms:
        raise UnsupportedStructure(
            "exactly one ring heteroatom, at a non-spiro ring position, is "
            "supported here (two or more ring heteroatoms, a heteroatom at "
            "a spiro atom itself, and an exocyclic heteroatom substituent "
            "are out of scope)"
        )
    (heteroatom_idx,) = heteroatoms
    a_prefix = _HETEROATOM_PREFIXES[mol.GetAtomWithIdx(heteroatom_idx).GetAtomicNum()]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    atom_rings = mol.GetRingInfo().AtomRings()
    spiro_count = len(spiro_atoms)
    spiro_prefix = numerical_term(spiro_count) + "spiro"

    best_key = None
    best_locants = None
    for order, spiros in _chain_direction_candidates(ring_order, spiro_atoms, atom_rings):
        arc_options = _arc_choice_options(graph, atom_rings, order, spiros)
        first_ring_non_spiro = set(atom_rings[order[0]]) - {spiros[0]}
        last_ring_non_spiro = set(atom_rings[order[-1]]) - {spiros[-1]}
        start_dirs = [a for a in graph[spiros[0]] if a in first_ring_non_spiro]
        end_dirs = [a for a in graph[spiros[-1]] if a in last_ring_non_spiro]
        for start_dir in start_dirs:
            for end_dir in end_dirs:
                for arc_choices in product(*arc_options):
                    seq, descriptor, superscripts = _build_sequence(
                        graph, atom_rings, order, spiros, start_dir, end_dir, arc_choices
                    )
                    locants = {atom: pos for pos, atom in enumerate(seq, start=1)}
                    descriptor_str = ".".join(
                        str(num) if sup is None else f"{num}^{locants[sup]}"
                        for num, sup in zip(descriptor, superscripts)
                    )
                    heteroatom_locant = locants[heteroatom_idx]
                    parent = f"{lambda_cited(mol, heteroatom_idx, heteroatom_locant)}-{a_prefix}{spiro_prefix}[{descriptor_str}]{alkane_name(len(seq))}"
                    spiro_locants = tuple(sorted(locants[s] for s in spiros))
                    substituents = substituents_for_ring(graph, seq, halogens)
                    key = _candidate_key(parent, spiro_locants, tuple(descriptor), heteroatom_locant, substituents)
                    if best_key is None or key < best_key:
                        best_key, best_locants = key, locants

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return best_key[-1]
    return stereo_locants_prefix(stereo, best_locants) + best_key[-1]
