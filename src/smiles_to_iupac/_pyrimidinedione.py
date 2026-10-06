"""Naming of the pyrimidine-2,4(1H,3H)-dione tautomer (uracil/thymine
diketo shape) -- a 6-membered mancude 1,3-diazine ring bearing two
ring-carbon oxo groups (at positions 2 and 4, flanking each nitrogen) and
two ring-nitrogen indicated hydrogens, named via the Blue Book's
parenthesized multi-locant "added indicated hydrogen" convention
(P-14.7.2/P-58.2.2), per the IUPAC 2013 Recommendations:

- Confirmed directly against two independent primary-source passages (not
  extrapolated from a single example): the Blue Book
  ('pyrimidine-4,6(1H,5H)-dione (PIN)', an explicitly PIN-marked dione on
  this exact ring skeleton needing two parenthesized added-hydrogen
  locants) and the Blue Book (barbituric acid's
  mancude-based alternative name 'pyrimidine-2,4,6(1H,3H,5H)-trione',
  confirming the same bracket convention scales to a third oxo/added-H
  pair -- though the PIN there, line 3722, switches to the fully
  saturated '1,3-diazinane-2,4,6-trione' once every ring carbon is oxo and
  no ring double bond survives, which is why that all-oxo case is out of
  scope here).
- Real structures (PubChem): uracil CID 1174 and thymine (5-methyluracil)
  CID 1135 -- PubChem's own auto-generated IUPACName field for both uses a
  looser front-form convention ('1H-pyrimidine-2,4-dione'), not the strict
  parenthesized PIN form built here.

`_pyrimidinone.py` already ships the single-oxo, single added-hydrogen
case (P-31.1.4.3.4's simpler front-form 'H-' prefix); this module is its
two-oxo extension, needing the harder parenthesized construction since two
indicated-hydrogen positions must now be cited together.

Scope, deliberately narrow: exactly one 6-membered, mancude 1,3-diazine
ring, two ring-carbon oxo groups at the two positions flanking each
nitrogen (leaving exactly one ring C=C double bond, at the two remaining
ring carbons), and at most one plain hydrocarbon substituent on one of
those two remaining carbons. The fully saturated triketo case (barbituric
acid, no remaining ring double bond), any other oxo pattern, any non-1,3
nitrogen relationship, or any non-hydrocarbon substituent fall through to
another module's existing rejection.
"""

from ._common import adjacency, ring_cycle
from ._substituents import format_substituent_prefixes, name_branch


def _match_pyrimidinedione(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = list(ring_info.AtomRings()[0])
    if len(ring_atoms) != 6:
        return None
    ring_set = set(ring_atoms)
    if any(not mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring_atoms):
        return None

    nitrogens = [atom for atom in ring_atoms if mol.GetAtomWithIdx(atom).GetAtomicNum() == 7]
    if len(nitrogens) != 2:
        return None
    for n_atom in nitrogens:
        n_data = mol.GetAtomWithIdx(n_atom)
        if n_data.GetTotalNumHs() != 1 or n_data.GetFormalCharge() != 0:
            return None
        if any(n.GetIdx() not in ring_set for n in n_data.GetNeighbors()):
            return None

    carbons = [atom for atom in ring_atoms if atom not in nitrogens]
    if any(mol.GetAtomWithIdx(atom).GetAtomicNum() != 6 for atom in carbons):
        return None

    oxo_atoms = []
    oxo_oxygens = set()
    other_carbons = []
    for atom in carbons:
        data = mol.GetAtomWithIdx(atom)
        exocyclic = [n for n in data.GetNeighbors() if n.GetIdx() not in ring_set]
        if not exocyclic:
            if data.GetTotalNumHs() != 1:
                return None
            other_carbons.append((atom, None))
            continue
        if len(exocyclic) != 1 or data.GetTotalNumHs() != 0:
            return None
        (nb,) = exocyclic
        bond = mol.GetBondBetweenAtoms(atom, nb.GetIdx())
        if nb.GetAtomicNum() == 8 and nb.GetFormalCharge() == 0 and nb.GetDegree() == 1 and bond.GetBondTypeAsDouble() == 2.0:
            oxo_atoms.append(atom)
            oxo_oxygens.add(nb.GetIdx())
        elif nb.GetAtomicNum() == 6 and bond.GetBondTypeAsDouble() == 1.0:
            other_carbons.append((atom, nb.GetIdx()))
        else:
            return None
    if len(oxo_atoms) != 2 or len(other_carbons) != 2:
        return None

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_set or idx in oxo_oxygens:
            continue
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0:
            return None

    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ring_set and b in ring_set:
            continue
        if a in ring_set or b in ring_set:
            continue
        if bond.GetBondTypeAsDouble() != 1.0:
            return None

    graph = adjacency(mol)
    ring_order = ring_cycle(graph, ring_atoms)
    plain_atoms = {atom for atom, _ in other_carbons}

    best = None
    for start_atom in nitrogens:
        start = ring_order.index(start_atom)
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
            if mol.GetAtomWithIdx(candidate[0]).GetAtomicNum() != 7:
                continue
            if mol.GetAtomWithIdx(candidate[2]).GetAtomicNum() != 7:
                continue
            if candidate[1] not in oxo_atoms or candidate[3] not in oxo_atoms:
                continue
            if candidate[4] not in plain_atoms or candidate[5] not in plain_atoms:
                continue

            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            sub_locant, sub_name = None, None
            for atom, nb in other_carbons:
                if nb is not None:
                    sub_locant = position_of[atom]
                    sub_name = name_branch(graph, nb, atom, {}, mol=mol)
            if best is None or (sub_locant is not None and (best[0] is None or sub_locant < best[0])):
                best = (sub_locant, sub_name)
    return best


def has_pyrimidinedione_shape(mol) -> bool:
    return _match_pyrimidinedione(mol) is not None


def name_pyrimidinedione(mol) -> str:
    sub_locant, sub_name = _match_pyrimidinedione(mol)
    prefix = ""
    if sub_locant is not None:
        name, is_compound = sub_name
        grouped = {name: {"locants": [sub_locant], "compound": is_compound}}
        prefix = format_substituent_prefixes(grouped)
    return f"{prefix}pyrimidine-2,4(1H,3H)-dione"
