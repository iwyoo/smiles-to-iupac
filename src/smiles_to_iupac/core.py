from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._alcohol import name_alcohol
from ._aldehyde import name_aldehyde
from ._amide import has_amide_shape, name_amide
from ._amine import name_amine
from ._aromatic import find_aromatic_fused_core, name_aromatic_fused
from ._bicyclic import find_bicyclic_core, name_bicycloalkane
from ._branched_fused_aromatic import has_retained_branched_fused_name, name_retained_branched_fused
from ._carboxylic_acid import has_carboxylic_acid_shape, name_carboxylic_acid
from ._common import UnsupportedStructure, non_single_bonds
from ._cyclic import name_cycloalkane
from ._ester import has_ester_shape, name_ester
from ._ether import has_ether_shape, name_ether
from ._fullerene import has_fullerene_name, name_fullerene
from ._heteroaromatic_fused import has_retained_heteroaromatic_fused_name, name_retained_heteroaromatic_fused
from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic
from ._ketone import name_ketone
from ._nitrile import has_nitrile_shape, name_nitrile
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._peri_fused_aromatic import has_retained_peri_fused_name, name_retained_peri_fused
from ._polyspiro import find_linear_polyspiro_chain, name_linear_polyspiro
from ._ring_assembly import find_ring_assembly_core, name_ring_assembly
from ._silane_chain import has_silane_chain_shape, name_silane_chain
from ._spiro import find_monospiro_atom, name_monospiro
from ._spiro_heteroatom import (
    has_single_ring_heteroatom_shape as has_single_spiro_heteroatom_shape,
    name_spiro_heteroatom,
)
from ._sulfide import has_sulfide_shape, name_sulfide
from ._thiol import has_thiol_shape, name_thiol
from ._tricyclic import find_propellane_core, name_propellane
from ._unsaturated import name_acyclic_unsaturated
from ._von_baeyer_heteroatom import (
    has_single_ring_heteroatom_shape as has_single_bicyclic_heteroatom_shape,
    name_von_baeyer_heteroatom,
)


def _is_aldehyde_shaped(carbonyl_oxygen):
    (carbon,) = carbonyl_oxygen.GetNeighbors()
    return carbon.GetAtomicNum() == 6 and sum(1 for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6) == 1


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    # An all-silicon skeleton (P-21.2.1's silane chain naming) has no
    # carbon at all, so it must be routed here before every other branch
    # below, all of which assume at least one carbon atom.
    if has_silane_chain_shape(mol):
        return name_silane_chain(mol)

    # buckminsterfullerene (P-27's '[60]fullerene', a fixed 12-pentagon/
    # 20-hexagon cage) is recognized by exact whole-molecule match --
    # see _fullerene.py's module docstring; none of the ring modules
    # below understand a cage shape at all.
    if has_fullerene_name(mol):
        return name_fullerene(mol)

    # pyrene/acenaphthylene (P-25.1.2's peri-fused retained names) are
    # recognized by exact whole-molecule match, independent of every other
    # branch below -- see _peri_fused_aromatic.py's module docstring for
    # why (acenaphthylene in particular has a non-6-membered, non-aromatic
    # ring that none of the other dispatch branches expect).
    if has_retained_peri_fused_name(mol):
        return name_retained_peri_fused(mol)

    # triphenylene (P-25.1.2's branched-fusion retained name) is recognized
    # the same way -- see _branched_fused_aromatic.py's module docstring;
    # _aromatic.py's chain-only algorithm explicitly rejects this shape.
    if has_retained_branched_fused_name(mol):
        return name_retained_branched_fused(mol)

    # quinoline/1H-indole (P-25.2.1's heteroaromatic retained names) are
    # recognized the same way -- see _heteroaromatic_fused.py's module
    # docstring; neither has an all-carbon skeleton, so _aromatic.py's
    # dispatch would never even consider them.
    if has_retained_heteroaromatic_fused_name(mol):
        return name_retained_heteroaromatic_fused(mol)

    # oxirane/thiane/piperidine etc. (P-22.2.1's Hantzsch-Widman
    # saturated-monocyclic retained names) are recognized the same way --
    # see _hetero_monocyclic.py's module docstring; none of the O/N
    # branches below understand a plain heteroatom ring at all.
    if has_hetero_monocyclic_name(mol):
        return name_hetero_monocyclic(mol)

    # A single O/N/S skeletal atom in an otherwise-carbon von Baeyer
    # bicyclic ring (P-23.2.1's 'a'-prefix skeletal replacement) or
    # monospiro ring system (P-24.2.1's) must be routed here before any of
    # the O/N-triggered branches below, none of which understand a ring
    # heteroatom at all.
    bicyclic_core = find_bicyclic_core(mol)
    if bicyclic_core is not None and has_single_bicyclic_heteroatom_shape(mol, bicyclic_core):
        return name_von_baeyer_heteroatom(mol, bicyclic_core)
    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None and has_single_spiro_heteroatom_shape(mol, spiro_atom):
        return name_spiro_heteroatom(mol, spiro_atom)

    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        # A plain -O- ether (P-63.2.1) has no suffix, so it must be routed
        # here before the carbonyl/alcohol checks below, none of which
        # accept a degree-2 oxygen at all.
        if has_ether_shape(mol):
            return name_ether(mol)
        # A carbon bearing both a carbonyl oxygen and a second, carbon-bonded
        # oxygen is an ester (-COO-), which must be routed before the
        # carboxylic-acid/aldehyde/ketone checks below: its carbonyl half
        # would otherwise look aldehyde/ketone-shaped, and (for a rejected,
        # out-of-scope case) its non-carbonyl oxygen would never satisfy the
        # carboxylic acid module's hydroxyl (O-H) requirement anyway.
        if has_ester_shape(mol):
            return name_ester(mol)
        # A carbon bearing both a carbonyl and a hydroxyl oxygen is a -COOH
        # group (Table 3.3's most senior suffix here) and must be routed
        # before the aldehyde/ketone/alcohol checks below, which would
        # otherwise misread its carbonyl or hydroxyl half in isolation.
        if has_carboxylic_acid_shape(mol):
            return name_carboxylic_acid(mol)
        # A carbon bearing both a carbonyl oxygen and a primary-amide
        # nitrogen (-CONH2) is an amide (junior only to the acid/ester
        # suffixes above in Table 3.3) and must be routed before the
        # aldehyde/ketone checks below: an amide carbon looks
        # aldehyde-shaped to `_is_aldehyde_shaped` (it counts only carbon
        # neighbors, ignoring the nitrogen).
        if has_amide_shape(mol):
            return name_amide(mol)
        # A doubly-bonded, monovalent oxygen is carbonyl-shaped (aldehyde or
        # ketone, depending on how many carbon neighbors its carbon has);
        # anything else falls to the alcohol module, which itself rejects a
        # coexisting carbonyl oxygen it finds among otherwise hydroxyl-only
        # atoms (Table 3.3 seniority between 'ol'/'one'/'al' isn't handled).
        carbonyl_oxygens = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetBonds()[0].GetBondTypeAsDouble() == 2.0
        ]
        if carbonyl_oxygens:
            if any(_is_aldehyde_shaped(o) for o in carbonyl_oxygens):
                return name_aldehyde(mol)
            return name_ketone(mol)
        return name_alcohol(mol)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        # A nitrile nitrogen (-C#N, P-66.5) has no oxygen, so it reaches this
        # branch alongside plain amines; it must be routed here before
        # name_amine, which doesn't recognize a triple-bonded nitrogen at all.
        if has_nitrile_shape(mol):
            return name_nitrile(mol)
        return name_amine(mol)
    if has_sulfide_shape(mol):
        # A plain -S- sulfide (P-63.2.1) has no suffix, so it must be routed
        # here before has_thiol_shape below: _thiol.py's validation rejects
        # a degree-2 sulfur outright (not a monovalent -SH), so a sulfide
        # would otherwise raise the wrong error there instead of being named.
        return name_sulfide(mol)
    if has_thiol_shape(mol):
        # A thiol (-SH, P-63.1.1) has neither O nor N, so it only reaches
        # this branch once both are ruled out above -- this module doesn't
        # yet handle Table 3.3's alcohol/thiol/amine seniority coexistence,
        # so a molecule with O or N never reaches here at all (see
        # _thiol.py's module docstring).
        return name_thiol(mol)

    num_rings = mol.GetRingInfo().NumRings()
    # Two disjoint (unfused) benzene rings joined by a single bond -- e.g.
    # biphenyl -- must be routed here before find_aromatic_fused_core: that
    # function only checks each SSSR ring is a 6-membered aromatic carbocycle
    # and doesn't require the rings to be fused, so it would otherwise claim
    # this shape too and then fail in name_aromatic_fused's ring-fusion-graph
    # validation (no shared bond means no fusion edge at all).
    if num_rings == 2:
        ring_assembly_core = find_ring_assembly_core(mol)
        if ring_assembly_core is not None:
            return name_ring_assembly(mol, ring_assembly_core)
    # Aromatic rings carry non-single (order 1.5) bonds, which every other
    # ring module's non_single_bonds check rejects; an aromatic ring
    # system's carbon skeleton can also be graph-isomorphic to a *saturated*
    # bicyclic through pentacyclic core (e.g. naphthalene <-> decahydro-
    # naphthalene), so this check must run, and must succeed for any
    # in-scope aromatic shape, before num_rings==1 or any saturated
    # find_*_core below gets a chance to misdetect it and raise the wrong
    # ("unsaturated ... not supported yet") error (see _aromatic.py).
    if num_rings >= 1:
        aromatic_core = find_aromatic_fused_core(mol)
        if aromatic_core is not None:
            return name_aromatic_fused(mol, aromatic_core)
    if num_rings == 0:
        bonds = non_single_bonds(mol)
        if not bonds:
            return name_acyclic_alkane(mol)
        if all(order in (2.0, 3.0) for _, _, order in bonds):
            return name_acyclic_unsaturated(mol)
        raise UnsupportedStructure(
            "a bond order other than double or triple is not supported yet "
            "(see P-31.1.1.1)"
        )
    if num_rings == 1:
        return name_cycloalkane(mol)

    # num_rings >= 2 from here on. RDKit's SSSR can overcount rings for
    # symmetric bridged bicyclics (see _bicyclic.py's find_bicyclic_core
    # docstring), so bicyclic detection isn't gated on num_rings == 2 either.
    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro(mol, spiro_atom)
    polyspiro_chain = find_linear_polyspiro_chain(mol)
    if polyspiro_chain is not None:
        return name_linear_polyspiro(mol, polyspiro_chain)
    bicyclic_core = find_bicyclic_core(mol)
    if bicyclic_core is not None:
        return name_bicycloalkane(mol, bicyclic_core)
    for ring_count in (3, 4, 5):
        core = find_polycyclic_core(mol, ring_count)
        if core is not None:
            return name_polycycloalkane(mol, core, ring_count)
    propellane_core = find_propellane_core(mol)
    if propellane_core is not None:
        return name_propellane(mol, propellane_core)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
