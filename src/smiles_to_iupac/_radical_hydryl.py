"""A radical made by adding H. to a mancude ring system is named with the suffix 'hydryl' at the locant of the atom
that received the hydrogen, every other group a prefix (P-71.1, P-71.2.1.3)."""

import networkx as nx
from rdkit import Chem

from ._common import adjacency
from ._diester_ring_diyl import _system_of, evaluate_skeleton

_LONE_PAIR_RING_ATOMS = {8, 16, 34, 52}


def _classify(kek, centre, atoms):
    """(pi atoms, the one sp3 atom that received the hydrogen) of the ring system, else None."""
    pi, received = {centre}, []
    for index in atoms:
        if index == centre:
            continue
        atom = kek.GetAtomWithIdx(index)
        doubles = [b for b in atom.GetBonds() if b.GetBondTypeAsDouble() > 1.0]
        if any(b.GetOtherAtomIdx(index) not in atoms or b.GetBondTypeAsDouble() > 2.0 for b in doubles):
            return None
        if doubles:
            if len(doubles) > 1 or atom.GetFormalCharge():
                return None
            pi.add(index)
        elif atom.GetAtomicNum() in _LONE_PAIR_RING_ATOMS and atom.GetDegree() == 2 and not atom.GetTotalNumHs():
            continue
        elif atom.GetAtomicNum() in (6, 7) and atom.GetTotalNumHs() and not atom.GetFormalCharge() and (
            atom.GetAtomicNum() == 6 or atom.GetDegree() == 2
        ):
            received.append(index)
        else:
            return None
    return (pi | set(received), received[0]) if len(received) == 1 else None


def _mancude_parent(kek, centre, atoms):
    """The mancude ring system plus substituents whose H. adduct is the radical, and the atom that received the
    hydrogen; None when the radical is not such an adduct."""
    found = _classify(kek, centre, atoms)
    if found is None or any(b.GetBondTypeAsDouble() != 1.0 for b in kek.GetAtomWithIdx(centre).GetBonds()):
        return None
    pi, received = found
    editable = Chem.RWMol(kek)
    for index in (centre, received):
        atom = editable.GetAtomWithIdx(index)
        hydrogens = atom.GetTotalNumHs() - (index == received)
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(hydrogens)
    network = nx.Graph()
    network.add_nodes_from(pi)
    for bond in kek.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in pi and b in pi:
            editable.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.SINGLE)
            network.add_edge(a, b)
    matching = nx.max_weight_matching(network, maxcardinality=True)
    if 2 * len(matching) != len(pi):
        return None
    for a, b in matching:
        editable.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.DOUBLE)
    parent = editable.GetMol()
    try:
        Chem.SanitizeMol(parent)
    except Chem.rdchem.MolSanitizeException:
        return None
    return parent, received


def hydryl_radical_name(mol, centre):
    """Name of the radical `centre` of `mol` as a hydryl radical of a mancude ring system, else None."""
    if (
        centre.GetAtomicNum() not in (6, 7)
        or centre.GetNumRadicalElectrons() != 1
        or not centre.IsInRing()
        or any(a.GetFormalCharge() or a.GetIsotope() or a.HasProp("_named_prefix") for a in mol.GetAtoms())
    ):
        return None
    kek = Chem.Mol(mol)
    Chem.Kekulize(kek, clearAromaticFlags=True)
    system = _system_of(kek, centre.GetIdx())
    if system is None:
        return None
    rings, atoms = system
    built = _mancude_parent(kek, centre.GetIdx(), atoms)
    if built is None:
        return None
    parent, received = built
    ranks = Chem.CanonicalRankAtoms(parent, breakTies=False)
    if len({ranks[a] for a in atoms}) == 1:
        return None
    found = evaluate_skeleton(parent, adjacency(parent), "ring", rings, atoms, [received], set(), "hydryl")
    return found[1] if found is not None else None
