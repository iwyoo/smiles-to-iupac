"""Naming of proline (pyrrolidine-2-carboxylic acid) using its P-103
retained name and L/D stereodescriptor, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-103.1.1.1 (Table 10.4): the retained name 'proline'.
- P-103.1.2: "The atoms in proline are numbered as in pyrrolidine, the
  nitrogen atom being numbered '1', and the carbon atom bonded to the
  carboxy group is numbered '2'." Unlike every other common amino acid,
  proline's own amino group is the ring nitrogen itself (a secondary
  amine), not an exocyclic primary amine -- `_carboxylic_acid_amine.py`
  explicitly excludes rings and requires a primary amine, so it can never
  reach this shape; a separate module is needed, not an extension of that
  one (confirmed empirically during scoping, #1056 M2 step 1 -- both
  L-/D-proline currently misroute to `_ketone.py`'s hetero-ring-ketone
  path, an unrelated module that claims any saturated single-heteroatom
  5/6/7-membered ring purely by shape).
- P-103.1.3.1: 'L' corresponds to CIP 'S' at the alpha-carbon of the
  unsubstituted proline.

Scope: a pyrrolidine ring whose ring nitrogen carries only hydrogen and whose C-2 bears the exocyclic -COOH.
Substituents on C-3, C-4 and C-5 are cited as prefixes, and their centres by CIP descriptors before the prefixes while
the alpha centre keeps D or L (P-103.1.3.2.1, P-103.2.3).
"""

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._carboxylic_acid_amine import _find_carboxylic_acid_carbon
from ._common import adjacency, specified_stereocenters

_ALPHA_TO_LD = {"S": "L", "R": "D"}


def _match(mol):
    """alpha_carbon_idx if `mol` is plain proline's exact shape, else
    None."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring = ring_info.AtomRings()[0]
    if len(ring) != 5 or any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return None
    ring_set = set(ring)

    nitrogens = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 7]
    if len(nitrogens) != 1:
        return None
    (ring_n,) = nitrogens
    if any(mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in ring if a != ring_n):
        return None

    n_atom = mol.GetAtomWithIdx(ring_n)
    if n_atom.GetFormalCharge() != 0 or n_atom.GetDegree() != 2 or n_atom.GetTotalNumHs() != 1:
        return None

    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return None
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    acid_carbon_idx = acid_carbon.GetIdx()
    if acid_carbon_idx in ring_set:
        return None

    graph = adjacency(mol)
    alpha_candidates = [n for n in graph[acid_carbon_idx] if n in ring_set]
    if len(alpha_candidates) != 1:
        return None
    (alpha_carbon,) = alpha_candidates
    if ring_n not in graph[alpha_carbon]:
        return None

    alpha_atom = mol.GetAtomWithIdx(alpha_carbon)
    if alpha_atom.GetFormalCharge() != 0 or alpha_atom.GetDegree() != 3 or alpha_atom.GetTotalNumHs() != 1:
        return None

    path = [alpha_carbon]
    while len(path) < 4:
        path.append(next(n for n in graph[path[-1]] if n in ring_set and n not in path and n != ring_n))
    locant_of = {path[0]: "2", path[1]: "3", path[2]: "4", path[3]: "5"}
    acid = {acid_carbon_idx, carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    substituents = []
    for ring_carbon in ring_set - {ring_n, alpha_carbon}:
        atom = mol.GetAtomWithIdx(ring_carbon)
        if atom.GetFormalCharge() != 0 or atom.GetDegree() not in (2, 3, 4):
            return None
        for n in graph[ring_carbon]:
            if n not in ring_set:
                substituents.append((locant_of[ring_carbon], ring_carbon, n))
    covered = set(ring_set) | acid
    for _, site, root in substituents:
        stack = [root]
        while stack:
            node = stack.pop()
            if node in covered:
                if node in ring_set:
                    return None
                continue
            covered.add(node)
            stack.extend(n for n in graph[node] if n != site)
    if len(covered) != mol.GetNumAtoms():
        return None
    return alpha_carbon, substituents, locant_of


def has_proline_shape(mol) -> bool:
    return _match(mol) is not None


def _alpha_descriptor(mol, alpha_carbon, substituents):
    """D or L from the alpha carbon of the proline with every ring substituent replaced by hydrogen: a substituent such
    as a sulfur on C-3 can outrank the carboxy group and reverse the CIP label (P-103.1.3.1)."""
    graph = adjacency(mol)
    removed = set()
    for _, site, root in substituents:
        stack = [root]
        while stack:
            node = stack.pop()
            if node not in removed:
                removed.add(node)
                stack.extend(n for n in graph[node] if n != site)
    editable = Chem.RWMol(mol)
    editable.GetAtomWithIdx(alpha_carbon).SetIntProp("_alpha", 1)
    for index in sorted(removed, reverse=True):
        editable.RemoveAtom(index)
    parent = editable.GetMol()
    Chem.SanitizeMol(parent)
    rdCIPLabeler.AssignCIPLabels(parent)
    atom = next(a for a in parent.GetAtoms() if a.HasProp("_alpha"))
    return _ALPHA_TO_LD[atom.GetProp("_CIPCode")]


def name_proline(mol) -> str:
    from ._amino_acid import _substituent_prefixes

    alpha_carbon, substituents, locant_of = _match(mol)
    prefixes = _substituent_prefixes(mol, adjacency(mol), substituents) if substituents else ""
    labels = dict(specified_stereocenters(mol) or ())
    other = sorted(
        (int(locant_of[atom]), f"{locant_of[atom]}{label}") for atom, label in labels.items() if atom != alpha_carbon
    )
    descriptors = f"({','.join(text for _, text in other)})-" if other else ""
    stem = f"{_alpha_descriptor(mol, alpha_carbon, substituents)}-proline" if alpha_carbon in labels else "proline"
    return descriptors + prefixes + ("-" if prefixes and alpha_carbon in labels else "") + stem
