from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._acyl_halide import has_acyl_halide_shape, name_acyl_halide
from ._anhydride import has_anhydride_shape, name_anhydride
from ._carbamate import has_carbamate_shape, name_carbamate
from ._alcohol import name_alcohol
from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._aldehyde import name_aldehyde
from ._carboxylic_acid_amine import has_carboxylic_acid_amine_shape, name_carboxylic_acid_amine
from ._aldehyde_carboxylic_acid import (
    has_aldehyde_carboxylic_acid_shape,
    name_aldehyde_carboxylic_acid,
)
from ._aldehyde_ketone import has_aldehyde_ketone_shape, name_aldehyde_ketone
from ._acetal import has_acetal_shape, name_acetal
from ._amide import has_amide_shape, name_amide
from ._hidden_amide_ketone import has_hidden_amide_shape, name_hidden_amide_ketone
from ._thiourea import has_thiourea_shape, name_thiourea
from ._selenourea import has_selenourea_shape, name_selenourea
from ._tellurourea import has_tellurourea_shape, name_tellurourea
from ._urea import has_urea_shape, name_urea
from ._amidine import has_amidine_shape, name_amidine
from ._guanidine import has_guanidine_shape, name_guanidine
from ._hydrazide import has_hydrazide_shape, name_hydrazide
from ._imide import has_imide_shape, name_imide
from ._amine import name_amine
from ._ring_amine import has_ring_amine_shape, name_ring_amine
from ._amine_oxide import has_amine_oxide_shape, name_amine_oxide
from ._aminide import has_aminide_shape, name_aminide
from ._ammonium import has_ammonium_shape, name_ammonium
from ._phosphonium import has_phosphonium_shape, name_phosphonium
from ._oxonium import has_oxonium_shape, name_oxonium
from ._carbenium import has_carbenium_shape, name_carbenium
from ._sulfonium import has_sulfonium_shape, name_sulfonium
from ._diazonium import has_diazonium_shape, name_diazonium
from ._radical import has_radical_shape, name_radical
from ._aromatic import find_aromatic_fused_core, name_aromatic_fused
from ._bicyclic import find_bicyclic_core, name_bicycloalkane
from ._bridged_aromatic import (
    find_bridged_anthracene_core,
    find_bridged_naphthalene_core,
    name_bridged_anthracene,
    name_bridged_naphthalene,
)
from ._borane import has_simple_borane_shape, name_simple_borane
from ._branched_fused_aromatic import has_retained_branched_fused_name, name_retained_branched_fused
from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._selenoate import has_selenoate_shape, name_selenoate
from ._thioate import has_thioate_shape, name_thioate
from ._carboxylic_acid import has_carboxylic_acid_shape, name_carboxylic_acid
from ._carboxylic_acid_sulfonamide import (
    has_carboxylic_acid_sulfonamide_shape,
    name_carboxylic_acid_sulfonamide,
)
from ._carboxylic_acid_seleninic_acid import (
    has_carboxylic_acid_seleninic_acid_shape,
    name_carboxylic_acid_seleninic_acid,
)
from ._carboxylic_acid_sulfinic_acid import (
    has_carboxylic_acid_sulfinic_acid_shape,
    name_carboxylic_acid_sulfinic_acid,
)
from ._carboxylic_acid_sulfonic_acid import (
    has_carboxylic_acid_sulfonic_acid_shape,
    name_carboxylic_acid_sulfonic_acid,
)
from ._common import UnsupportedStructure, non_single_bonds
from ._cyclic import name_cycloalkane
from ._cyclic_unsaturated import find_cyclic_unsaturated_core, name_cyclic_unsaturated
from ._dihydro_aromatic import find_dihydronaphthalene_core, name_dihydronaphthalene
from ._diester_acyloxy import has_diester_shape, name_diester_acyloxy
from ._ester import has_ester_shape, name_ester
from ._cyanate import has_cyanate_shape, name_cyanate
from ._ether import has_ether_shape, name_ether
from ._fullerene import has_fullerene_name, name_fullerene
from ._androstane import has_androstane_name, name_androstane
from ._gonane import has_gonane_name, name_gonane
from ._heteroaromatic_fused import has_retained_heteroaromatic_fused_name, name_retained_heteroaromatic_fused
from ._polycyclic_component_fusion import (
    has_polycyclic_component_fusion_name,
    name_polycyclic_component_fusion,
)
from ._hetero_monocyclic import (
    has_hetero_monocyclic_name,
    has_hetero_monocyclic_substituent_name,
    name_hetero_monocyclic,
    name_hetero_monocyclic_substituent,
)
from ._hydroxylamine import has_hydroxylamine_shape, name_hydroxylamine
from ._imine import has_simple_imine_shape, name_imine
from ._isotope import has_isotope_shape, name_isotope
from ._two_component_heterocycle_fusion import (
    has_two_component_heterocycle_fusion_name,
    name_two_component_heterocycle_fusion,
)
from ._ketone import (
    has_five_membered_1_2_ring_ketone_shape,
    has_five_membered_1_3_ring_ketone_shape,
    has_hetero_ring_ketone_shape,
    has_seven_membered_1_2_ring_ketone_shape,
    has_seven_membered_1_3_ring_ketone_shape,
    name_ketone,
)
from ._thione import has_thione_shape, name_thione
from ._selone import has_selone_shape, name_selone
from ._tellone import has_tellone_shape, name_tellone
from ._ketone_amide import has_ketone_amide_shape, name_ketone_amide
from ._ketone_ester import has_ketone_ester_shape, name_ketone_ester
from ._nitrile import has_nitrile_shape, name_nitrile
from ._selenocyanate import has_selenocyanate_shape, name_selenocyanate
from ._tellurocyanate import has_tellurocyanate_shape, name_tellurocyanate
from ._thiocyanate import has_thiocyanate_shape, name_thiocyanate
from ._azide import has_azide_shape, name_azide
from ._diazene import has_diazene_shape, name_diazene
from ._azine import has_azine_shape, name_azine
from ._hydrazine import has_hydrazine_shape, name_hydrazine
from ._hydrazone import has_hydrazone_shape, name_hydrazone
from ._diazo import has_diazo_shape, name_diazo
from ._isocyanate import has_isocyanate_shape, name_isocyanate
from ._isocyanide import has_isocyanide_shape, name_isocyanide
from ._isoselenocyanate import has_isoselenocyanate_shape, name_isoselenocyanate
from ._isotellurocyanate import has_isotellurocyanate_shape, name_isotellurocyanate
from ._isothiocyanate import has_isothiocyanate_shape, name_isothiocyanate
from ._nitro import has_nitro_shape, name_nitro
from ._nitroso import has_nitroso_shape, name_nitroso
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._cyclophane import has_cyclophane_name, name_cyclophane
from ._naphthalene_benzene_phane import (
    has_naphthalene_benzene_phane_name,
    name_naphthalene_benzene_phane,
)
from ._phosphane import has_simple_phosphane_shape, name_simple_phosphane
from ._phosphanone import has_phosphanone_shape, name_phosphanone
from ._phosphane_chain import has_phosphane_chain_shape, name_phosphane_chain
from ._peri_fused_aromatic import has_retained_peri_fused_name, name_retained_peri_fused
from ._hydroperoxide import has_hydroperoxide_shape, name_hydroperoxide
from ._peroxide import has_peroxide_shape, name_peroxide
from ._polyspiro import (
    find_branched_polyspiro_hub,
    find_linear_polyspiro_chain,
    name_branched_polyspiro,
    name_linear_polyspiro,
)
from ._polyspiro_heteroatom import (
    has_single_ring_heteroatom_shape as has_single_polyspiro_heteroatom_shape,
    name_linear_polyspiro_heteroatom,
)
from ._ring_assembly import find_ring_assembly_core, name_ring_assembly
from ._silane_chain import has_silane_chain_shape, name_silane_chain
from ._spiro import find_monospiro_atom, name_monospiro
from ._spiro_heteroatom import (
    has_single_ring_heteroatom_shape as has_single_spiro_heteroatom_shape,
    name_spiro_heteroatom,
)
from ._disulfide import has_disulfide_shape, name_disulfide
from ._sulfide import has_sulfide_shape, name_sulfide
from ._sulfinamide import has_sulfinamide_shape, name_sulfinamide
from ._sulfinic_acid import has_sulfinic_acid_shape, name_sulfinic_acid
from ._sulfonamide import has_sulfonamide_shape, name_sulfonamide
from ._seleninic_acid import has_seleninic_acid_shape, name_seleninic_acid
from ._selenonic_acid import has_selenonic_acid_shape, name_selenonic_acid
from ._tellurinic_acid import has_tellurinic_acid_shape, name_tellurinic_acid
from ._telluronic_acid import has_telluronic_acid_shape, name_telluronic_acid
from ._sulfonic_acid import has_sulfonic_acid_shape, name_sulfonic_acid
from ._sulfonic_acid_sulfinic_acid import (
    has_sulfonic_acid_sulfinic_acid_shape,
    name_sulfonic_acid_sulfinic_acid,
)
from ._sulfonic_acid_sulfonamide import (
    has_sulfonic_acid_sulfonamide_shape,
    name_sulfonic_acid_sulfonamide,
)
from ._sulfonic_acid_thiol import has_sulfonic_acid_thiol_shape, name_sulfonic_acid_thiol
from ._sulfone import has_sulfone_shape, name_sulfone
from ._sulfoxide import has_sulfoxide_shape, name_sulfoxide
from ._selenone import has_selenone_shape, name_selenone
from ._selenoxide import has_selenoxide_shape, name_selenoxide
from ._tellurone import has_tellurone_shape, name_tellurone
from ._telluroxide import has_telluroxide_shape, name_telluroxide
from ._diselenide import has_diselenide_shape, name_diselenide
from ._selenide import has_selenide_shape, name_selenide
from ._selenol import has_selenol_shape, name_selenol
from ._ditelluride import has_ditelluride_shape, name_ditelluride
from ._telluride import has_telluride_shape, name_telluride
from ._tellurol import has_tellurol_shape, name_tellurol
from ._selenoic_acid import has_selenoic_acid_shape, name_selenoic_acid
from ._telluroic_acid import has_telluroic_acid_shape, name_telluroic_acid
from ._thioic_acid import has_thioic_acid_shape, name_thioic_acid
from ._thiol import has_thiol_shape, name_thiol
from ._tricyclic import find_propellane_core, name_propellane
from ._unsaturated import name_acyclic_unsaturated
from ._von_baeyer_heteroatom import (
    has_mixed_element_heteroatom_shape as has_mixed_bicyclic_heteroatom_shape,
    has_multi_ring_heteroatom_shape as has_multi_bicyclic_heteroatom_shape,
    has_single_ring_heteroatom_shape as has_single_bicyclic_heteroatom_shape,
    has_single_ring_heteroatom_shape_polycyclic,
    name_von_baeyer_heteroatom,
    name_von_baeyer_heteroatom_mixed,
    name_von_baeyer_heteroatom_multi,
    name_von_baeyer_heteroatom_polycyclic,
)


def _is_aldehyde_shaped(carbonyl_oxygen):
    (carbon,) = carbonyl_oxygen.GetNeighbors()
    return carbon.GetAtomicNum() == 6 and sum(1 for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6) == 1


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    # An isotopically labeled atom (P-82's isotope descriptor nomenclature)
    # must be routed here before every other branch below: RDKit represents
    # an isotopically substituted hydrogen (e.g. 2H) as its own explicit
    # atom (atomic number 1), which none of the other branches recognize at
    # all -- every one of them would reject it outright as an unsupported
    # heteroatom.
    if has_isotope_shape(mol):
        return name_isotope(mol)

    # A radical center (P-71.2.1.1's 'yl' radical naming) must be routed
    # here before every other branch below: none of them recognize a
    # nonzero radical electron count at all -- an unbranched-chain or
    # monocyclic-ring radical would otherwise fall straight through to the
    # plain alkane/cycloalkane dispatch further down, which doesn't know a
    # hydrogen is missing.
    if has_radical_shape(mol):
        return name_radical(mol)

    # A charged ammonium nitrogen (P-73.1.1.2's hydron-addition cation
    # naming) must be routed here before every other branch below: none of
    # them recognize a charged atom at all -- `_amine.py` in particular
    # rejects any charged atom outright rather than attempting to name it.
    if has_ammonium_shape(mol):
        return name_ammonium(mol)

    # A secondary/tertiary amine N-oxide (P-62.5's zwitterionic N+-O-) has a
    # charged nitrogen too, for the same reason as ammonium above -- and
    # `has_ammonium_shape` itself doesn't match it (an oxide-bearing
    # nitrogen has 2-3 carbon neighbors plus the oxide oxygen, never the
    # single-carbon/three-H shape ammonium requires), so it needs its own
    # explicit routing here.
    if has_amine_oxide_shape(mol):
        return name_amine_oxide(mol)

    # A diazonium cation (R-N#N+, P-73) has a charged nitrogen too, for the
    # same reason as ammonium above -- routed here, unconditionally,
    # before every other branch.
    if has_diazonium_shape(mol):
        return name_diazonium(mol)

    # A phosphonium cation (P-73.1.1.2) has a charged phosphorus too, for
    # the same reason as ammonium above -- routed here, unconditionally,
    # before every other branch.
    if has_phosphonium_shape(mol):
        return name_phosphonium(mol)

    # A sulfonium cation (P-73.1.1.2) has a charged sulfur too, for the
    # same reason as ammonium above -- routed here, unconditionally,
    # before every other branch.
    if has_sulfonium_shape(mol):
        return name_sulfonium(mol)

    # An oxonium cation (P-73.1.1.2) has a charged oxygen too, for the
    # same reason as ammonium above -- routed here, unconditionally,
    # before every other branch.
    if has_oxonium_shape(mol):
        return name_oxonium(mol)

    # A carbenium cation (P-73.2.2.1.1's 'ylium' suffix naming) has a
    # charged carbon too, for the same reason as ammonium above -- routed
    # here, unconditionally, before every other branch.
    if has_carbenium_shape(mol):
        return name_carbenium(mol)

    # An all-silicon skeleton (P-21.2.1's silane chain naming) has no
    # carbon at all, so it must be routed here before every other branch
    # below, all of which assume at least one carbon atom.
    if has_silane_chain_shape(mol):
        return name_silane_chain(mol)

    # An all-phosphorus, multi-atom skeleton (diphosphane, triphosphane,
    # ...) has no carbon at all and must be routed before
    # has_simple_phosphane_shape below, which would otherwise misname it
    # via `_phosphane.py`'s own explicit "more than one phosphorus atom"
    # rejection.
    if has_phosphane_chain_shape(mol) and mol.GetNumAtoms() > 1:
        return name_phosphane_chain(mol)

    # A phosphine oxide (P-68.3.2.3.1's '-phosphanone' suffix, R-P(=O)<)
    # has its own phosphorus-bonded oxygen that `_phosphane.py` doesn't
    # expect at all (that module rejects any heteroatom besides its own
    # phosphorus outright) -- must be routed here first, before
    # has_simple_phosphane_shape below, for the same reason as
    # has_phosphane_chain_shape above.
    if has_phosphanone_shape(mol):
        return name_phosphanone(mol)

    # A phosphorus atom (P-68's phosphane substitutive nomenclature) must
    # be routed here before every other branch below: none of them
    # recognize phosphorus at all, and a phosphane carbon substituent would
    # otherwise reach the plain acyclic-alkane/amine dispatch further down
    # with no phosphorus handling.
    if has_simple_phosphane_shape(mol):
        return name_simple_phosphane(mol)

    # A boron atom (P-68's borane substitutive nomenclature, the same shape
    # as phosphane above with boron in place of phosphorus) must be routed
    # here for the same reason -- none of the branches below recognize
    # boron at all.
    if has_simple_borane_shape(mol):
        return name_simple_borane(mol)

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

    # [2.2]paracyclophane/[2.2]metacyclophane (P-26's phane nomenclature
    # retained-name-style recognition, see module docstring) are recognized
    # the same way, independent of every other branch below -- their two
    # -CH2CH2- bridges make their carbon skeletons look like bridged
    # aromatic ring systems to every other dispatch branch, none of which
    # understand phane nomenclature at all.
    if has_cyclophane_name(mol):
        return name_cyclophane(mol)

    # A naphthalene superatom + a benzene superatom joined by two bridges
    # (P-26's "different ring kinds" phane case, see module docstring) is
    # recognized the same way -- same reasoning as the plain cyclophane
    # check above.
    if has_naphthalene_benzene_phane_name(mol):
        return name_naphthalene_benzene_phane(mol)

    # gonane (the 1989 IUPAC steroid nomenclature's fundamental tetracyclic
    # parent, Rule 2.1 -- see module docstring) is recognized the same way,
    # independent of every other branch below: `_polycyclic.py`'s general
    # von Baeyer engine already names this exact skeleton (confirmed by
    # direct testing), so this check must come first or gonane would never
    # be reached.
    if has_gonane_name(mol):
        return name_gonane(mol)

    # androstane (gonane + the two angular C18/C19 methyls, Rule 3S-2.3 --
    # see module docstring) is recognized the same way, for the same
    # reason: its own skeleton would otherwise fall through to
    # `_polycyclic.py`'s general von Baeyer engine instead.
    if has_androstane_name(mol):
        return name_androstane(mol)

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

    # benzo[g]indole/benzo[e][1]benzofuran/benzo[g][1]benzofuran (P-25.3.1.3's
    # computed fusion-locant-letter mechanism, this time for a plain benzo
    # ring fused onto an already-bicyclic retained-name base component)
    # must be routed here before `_aromatic.py`'s own tricyclic dispatch
    # further below, which doesn't recognize a heteroatom at all.
    if has_polycyclic_component_fusion_name(mol):
        return name_polycyclic_component_fusion(mol)

    # thieno[2,3-b]thiophene/furo[2,3-b]furan/thieno[2,3-b]furan etc.
    # (P-25.3.1.3's computed fusion-locant-letter mechanism, for two
    # five-membered heteromonocycles -- identical or a mixed O/S pair --
    # self-fused) must be routed here before `_hetero_monocyclic.py` below,
    # which only understands a single ring.
    if has_two_component_heterocycle_fusion_name(mol):
        return name_two_component_heterocycle_fusion(mol)

    # oxirane/thiane/piperidine etc. (P-22.2.1's Hantzsch-Widman
    # saturated-monocyclic retained names) are recognized the same way --
    # see _hetero_monocyclic.py's module docstring; none of the O/N
    # branches below understand a plain heteroatom ring at all.
    if has_hetero_monocyclic_name(mol):
        return name_hetero_monocyclic(mol)

    # A single substituent on one of the same mancude parents above (P-22.2.1
    # heteroatom locants stay fixed; only one ring atom's H is replaced) --
    # see _hetero_monocyclic.py's module docstring for the role-sequence
    # matching this uses instead of the exact-match table above.
    if has_hetero_monocyclic_substituent_name(mol):
        return name_hetero_monocyclic_substituent(mol)

    # A single O/N/S skeletal atom in an otherwise-carbon von Baeyer
    # bicyclic ring (P-23.2.1's 'a'-prefix skeletal replacement) or
    # monospiro ring system (P-24.2.1's) must be routed here before any of
    # the O/N-triggered branches below, none of which understand a ring
    # heteroatom at all.
    bicyclic_core = find_bicyclic_core(mol)
    if bicyclic_core is not None and has_single_bicyclic_heteroatom_shape(mol, bicyclic_core):
        return name_von_baeyer_heteroatom(mol, bicyclic_core)
    if bicyclic_core is not None and has_multi_bicyclic_heteroatom_shape(mol, bicyclic_core):
        return name_von_baeyer_heteroatom_multi(mol, bicyclic_core)
    if bicyclic_core is not None and has_mixed_bicyclic_heteroatom_shape(mol, bicyclic_core):
        return name_von_baeyer_heteroatom_mixed(mol, bicyclic_core)
    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None and has_single_spiro_heteroatom_shape(mol, spiro_atom):
        return name_spiro_heteroatom(mol, spiro_atom)
    # Same reasoning, one step up the spiro chain: a single O/N/S ring atom
    # (never a spiro atom itself) in an otherwise-carbon linear polyspiro
    # ring system (P-24.2.4) must be routed here for the same reason --
    # `find_linear_polyspiro_chain` itself is topology-only and doesn't
    # care about atom identity, so it would otherwise match first further
    # down and hand the heteroatom-bearing molecule to the all-carbon
    # `name_linear_polyspiro`, which explicitly rejects any non-carbon atom.
    polyspiro_chain = find_linear_polyspiro_chain(mol)
    if polyspiro_chain is not None and has_single_polyspiro_heteroatom_shape(mol, polyspiro_chain):
        return name_linear_polyspiro_heteroatom(mol, polyspiro_chain)

    # A naphthalene skeleton with a single -CH2- or -O- bridge across one
    # ring's 1,4-positions (1,4-dihydro-1,4-methano-/epoxynaphthalene) must
    # be routed here before the plain "any O atom" branch below (the -O-
    # bridge variant would otherwise be misdetected as a plain ether) -- see
    # _bridged_aromatic.py's module docstring for why RDKit's own ring
    # perception can't be trusted for this shape either.
    bridged_core = find_bridged_naphthalene_core(mol)
    if bridged_core is not None:
        return name_bridged_naphthalene(mol, bridged_core)

    # Same shape, one ring larger: an anthracene skeleton bridged across
    # its own 9,10 meso positions (see _bridged_aromatic.py's module
    # docstring) must be routed here for the same reason.
    bridged_anthracene_core = find_bridged_anthracene_core(mol)
    if bridged_anthracene_core is not None:
        return name_bridged_anthracene(mol, bridged_anthracene_core)

    # A single O/N/S skeletal atom in an otherwise-carbon von Baeyer
    # polycyclic (ring_count>=3) ring (P-23.2.1's 'a'-prefix skeletal
    # replacement, the tricyclic+ generalization of the bicyclic case
    # routed near the top of this function) must be routed here, after the
    # more specific aromatic/bridged-ring routes just above (pyrene,
    # cyclophane, bridged naphthalene/anthracene) but before the plain
    # "any O atom" branch below (a degree-2 ring oxygen would otherwise be
    # misdetected as a plain ether, which explicitly rejects rings) -- see
    # _von_baeyer_heteroatom.py's module docstring.
    for ring_count in (3, 4, 5):
        polycyclic_hetero_core = find_polycyclic_core(mol, ring_count)
        if polycyclic_hetero_core is not None and has_single_ring_heteroatom_shape_polycyclic(
            mol, polycyclic_hetero_core
        ):
            return name_von_baeyer_heteroatom_polycyclic(mol, polycyclic_hetero_core, ring_count)

    # A sulfonic acid coexisting with a thiol (P-41/P-43, see
    # `_seniority.py`) has the same three-oxygen sulfonic sulfur as plain
    # sulfonic acid below, plus an extra thiol sulfur that `_sulfonic_
    # acid.py`'s own validation would otherwise reject outright -- must be
    # routed here first.
    if has_sulfonic_acid_thiol_shape(mol):
        return name_sulfonic_acid_thiol(mol)
    # A sulfonic acid coexisting with a sulfinic acid (P-41/P-43, see
    # `_seniority.py`) has the same four-oxygen sulfonic sulfur as plain
    # sulfonic acid below, plus an extra sulfinic sulfur that neither
    # `_sulfonic_acid.py` nor `_sulfinic_acid.py`'s own validation would
    # accept -- must be routed here first.
    if has_sulfonic_acid_sulfinic_acid_shape(mol):
        return name_sulfonic_acid_sulfinic_acid(mol)
    # A sulfonic acid coexisting with an unsubstituted sulfonamide
    # (P-41/P-43, see `_seniority.py`) has the same four-oxygen sulfonic
    # sulfur as plain sulfonic acid below, plus an extra sulfonamide
    # sulfur that neither `_sulfonic_acid.py` nor `_sulfonamide.py`'s own
    # validation would accept -- must be routed here first.
    if has_sulfonic_acid_sulfonamide_shape(mol):
        return name_sulfonic_acid_sulfonamide(mol)
    # A carboxylic acid coexisting with a sulfonic acid (P-41/P-43, see
    # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
    # would otherwise reject on sight of the extra sulfonic sulfur, and the
    # same four-oxygen sulfonic sulfur `_sulfonic_acid.py` would otherwise
    # reject on sight of the extra -COOH oxygens -- must be routed here
    # first, before either single-group module.
    if has_carboxylic_acid_sulfonic_acid_shape(mol):
        return name_carboxylic_acid_sulfonic_acid(mol)
    # A carboxylic acid coexisting with a sulfinic acid (P-41/P-43, see
    # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
    # would otherwise reject on sight of the extra sulfinic sulfur, and the
    # same three-oxygen sulfinic sulfur `_sulfinic_acid.py` would otherwise
    # reject on sight of the extra -COOH oxygens -- must be routed here
    # first, before either single-group module.
    if has_carboxylic_acid_sulfinic_acid_shape(mol):
        return name_carboxylic_acid_sulfinic_acid(mol)
    # A carboxylic acid coexisting with a seleninic acid (P-41/P-43, see
    # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
    # would otherwise reject on sight of the extra seleninic selenium, and
    # the same three-oxygen-cluster seleninic selenium
    # `_seleninic_acid.py` would otherwise reject on sight of the extra
    # -COOH oxygens -- must be routed here first, before either
    # single-group module.
    if has_carboxylic_acid_seleninic_acid_shape(mol):
        return name_carboxylic_acid_seleninic_acid(mol)
    # A carboxylic acid coexisting with an unsubstituted sulfonamide
    # (P-41/P-43, see `_seniority.py`) has the same -COOH carbon shape
    # `_carboxylic_acid.py` would otherwise reject on sight of the extra
    # sulfonamide sulfur, and the same four-oxygen-cluster sulfonamide
    # sulfur `_sulfonamide.py` would otherwise reject on sight of the
    # extra -COOH oxygens -- must be routed here first, before either
    # single-group module.
    if has_carboxylic_acid_sulfonamide_shape(mol):
        return name_carboxylic_acid_sulfonamide(mol)
    # A sulfonic acid (-SO3H, P-65.3.1) has three oxygens on its own sulfur,
    # so it must be routed here before the plain "any O atom" branch below --
    # none of the ether/ester/carboxylic-acid/aldehyde/ketone/alcohol checks
    # in that branch understand a sulfur-centered oxygen cluster at all.
    if has_sulfonic_acid_shape(mol):
        return name_sulfonic_acid(mol)
    # A selenonic acid (-Se(=O)(=O)OH, P-65.3.1) has the same oxygen-cluster
    # shape as sulfonic acid above, just on selenium instead of sulfur, so
    # it too must be routed before the plain "any O atom" branch.
    if has_selenonic_acid_shape(mol):
        return name_selenonic_acid(mol)
    # A telluronic acid (-Te(=O)(=O)OH, P-65.3.1) has the same oxygen-
    # cluster shape as sulfonic/selenonic acid above, just on tellurium, so
    # it too must be routed before the plain "any O atom" branch.
    if has_telluronic_acid_shape(mol):
        return name_telluronic_acid(mol)
    # A tellurinic acid (-Te(=O)OH, P-65.3.1) has the same oxygen-cluster
    # shape as sulfinic/seleninic acid, just on tellurium, so it too must
    # be routed before the plain "any O atom" branch.
    if has_tellurinic_acid_shape(mol):
        return name_tellurinic_acid(mol)
    # A seleninic acid (-Se(=O)OH, P-65.3.1) has the same oxygen-cluster
    # shape as sulfinic acid, just on selenium instead of sulfur, so it too
    # must be routed before the plain "any O atom" branch.
    if has_seleninic_acid_shape(mol):
        return name_seleninic_acid(mol)
    # A sulfonamide (-SO2NH2, P-65.3.1) has two oxygens on its own sulfur,
    # the same reasoning as sulfonic acid above, plus a nitrogen that would
    # otherwise be mistaken for a plain amine -- so it too must be routed
    # before both the "any O atom" and "any N atom" branches below.
    if has_sulfonamide_shape(mol):
        return name_sulfonamide(mol)
    # A sulfinic acid (-SO2H, P-65.3.1) has two oxygens on its own sulfur --
    # the same reasoning as sulfonic acid above -- so it too must be routed
    # before the plain "any O atom" branch.
    if has_sulfinic_acid_shape(mol):
        return name_sulfinic_acid(mol)
    # A sulfinamide (-S(=O)NH2, P-65.3.1) has one oxygen and one nitrogen on
    # its own sulfur -- the same reasoning as sulfonamide above -- so it too
    # must be routed before both the "any O atom" and "any N atom" branches.
    if has_sulfinamide_shape(mol):
        return name_sulfinamide(mol)
    # A sulfone (-SO2-, P-63.6) has two oxygens on its own sulfur, just like
    # a sulfinic/sulfonic acid's cluster above, so it must be routed here for
    # the same reason -- before it, since a sulfone's sulfur has two carbon
    # neighbors instead of the acid's hydroxyl, which would otherwise never
    # match `_sulfonic_acid.py`'s own shape check anyway, but routing it
    # alongside its acid relatives keeps this family together.
    if has_sulfone_shape(mol):
        return name_sulfone(mol)
    # A sulfoxide (-S(=O)-, P-63.6) has one oxygen on its own sulfur, same
    # reasoning as the sulfinic acid check above.
    if has_sulfoxide_shape(mol):
        return name_sulfoxide(mol)
    # A selenone (-Se(=O)(=O)-, P-63.6) is the selenium analogue of a
    # sulfone -- same reasoning, checked before the selenoxide below since
    # its selenium has two oxygens instead of one.
    if has_selenone_shape(mol):
        return name_selenone(mol)
    # A selenoxide (-Se(=O)-, P-63.6) is the selenium analogue of a
    # sulfoxide.
    if has_selenoxide_shape(mol):
        return name_selenoxide(mol)
    # A tellurone (-Te(=O)(=O)-, P-63.6) is the tellurium analogue of a
    # sulfone/selenone -- same reasoning, checked before the telluroxide
    # below since its tellurium has two oxygens instead of one.
    if has_tellurone_shape(mol):
        return name_tellurone(mol)
    # A telluroxide (-Te(=O)-, P-63.6) is the tellurium analogue of a
    # sulfoxide/selenoxide.
    if has_telluroxide_shape(mol):
        return name_telluroxide(mol)
    # A nitro group (-NO2, P-61.5.1) has its own nitrogen and two oxygens
    # neither the ether/carbonyl/alcohol checks below nor the plain-amine
    # branch further down expect, so it must be routed before both -- a
    # nitro-bearing molecule always has an oxygen atom, so it would
    # otherwise be swallowed by the "any O atom" branch's unconditional
    # `name_alcohol` fallback and never even reach the nitrogen branch.
    if has_nitro_shape(mol):
        return name_nitro(mol)
    # A nitroso group (-N=O, P-61.5.1's sibling prefix) has the same "own
    # oxygen" issue as nitro above, so it must be routed here for the same
    # reason.
    if has_nitroso_shape(mol):
        return name_nitroso(mol)
    # An isocyanate group (-N=C=O, P-61.8) has the same "own oxygen" issue
    # as nitro above, so it must be routed here for the same reason.
    if has_isocyanate_shape(mol):
        return name_isocyanate(mol)
    # An oxime (=N-OH/=N-O-R, P-68.3.1.1.2) is imine-shaped (C=N) but has
    # its own oxygen the ether/carbonyl/alcohol checks below don't expect
    # at all, so it must be routed before the "any O atom" branch for the
    # same reason as nitro/isocyanate above -- it would otherwise be
    # swallowed by that branch's unconditional `name_alcohol` fallback and
    # never reach the nitrogen branch further down. Gated on an oxygen
    # actually being present so this doesn't preempt a diazo group's own
    # C=N bond (checked later, in the nitrogen branch) -- diazo has no
    # oxygen at all, so it's unaffected either way, but this keeps the
    # check's intent explicit. A non-oxime imine (no oxygen) is still
    # correctly caught by has_simple_imine_shape's other call site further
    # down.
    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()) and has_simple_imine_shape(mol):
        return name_imine(mol)

    # Thiourea (H2N-C(=S)-NH2) has no oxygen at all, so it would otherwise
    # fall straight through the oxygen-gated block below (and every other
    # check in it) to the plain-amine fallback at the very end of this
    # function -- it must be checked here, unconditionally, before that
    # gate.
    if has_thiourea_shape(mol):
        return name_thiourea(mol)

    # Selenourea/tellurourea (H2N-C(=Se/Te)-NH2) have no oxygen either, for
    # the same reason as thiourea above.
    if has_selenourea_shape(mol):
        return name_selenourea(mol)

    if has_tellurourea_shape(mol):
        return name_tellurourea(mol)

    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        # An acyl group bonded directly to the nitrogen of an otherwise-
        # plain saturated monocyclic ring (e.g. 1-acetylpiperidine) is a
        # P-64.1.2.1(b) 'pseudoketone'/'hidden amide' -- it must be routed
        # before `has_hetero_ring_ketone_shape` below, which claims any
        # single-heteroatom saturated ring alongside any oxygen in the
        # molecule without checking the ring itself actually contains a
        # ketone (it would otherwise misname this shape's acyl branch as
        # an "N-alkyl substituent"), and before `has_amide_shape` further
        # down, whose `name_amide` unconditionally rejects any non-
        # benzene ring.
        if has_hidden_amide_shape(mol):
            return name_hidden_amide_ketone(mol)
        # A plain 1,4-N,O saturated six-membered ring whose nitrogen
        # carries one substituent (e.g. 4-methylmorpholine) looks
        # ether-shaped to `has_ether_shape` below (its ring oxygen is a
        # perfectly ordinary degree-2 ether oxygen) -- must be routed
        # before that check, mirroring the same "claim the specific ring
        # shape before the generic one" pattern used throughout this
        # branch. The N,N/N,S sibling shapes (piperazine/thiomorpholine)
        # have no oxygen at all and are instead routed in the
        # oxygen-free "any nitrogen" branch further down.
        if has_ring_amine_shape(mol):
            return name_ring_amine(mol)
        # A ketone carbonyl sitting directly between two ring heteroatoms
        # in a five-membered 1,3-related saturated ring (e.g.
        # 1,3-dioxolan-2-one, imidazolidin-2-one) looks ether-, acetal-,
        # or carbamate-shaped (depending on the element pair) to several
        # checks below and would otherwise be misnamed by their
        # acyclic-only constructions -- `_ketone.py` already names this
        # narrow ring shape correctly, so claim it here first, ahead of
        # every other oxygen-containing check in this branch. Unlike
        # `has_hetero_ring_ketone_shape` below, this shape can never
        # collide with a cyclic anhydride (see that function's own
        # docstring), so routing it this early is safe.
        if has_five_membered_1_3_ring_ketone_shape(mol):
            return name_ketone(mol)
        # Same reasoning for the 7-membered 1,3-related shape (e.g.
        # 1,3-diazepan-2-one) -- an O+N pair here looks carbamate-shaped
        # (N-C(=O)-O) instead, but the collision-safety argument is
        # identical.
        if has_seven_membered_1_3_ring_ketone_shape(mol):
            return name_ketone(mol)
        # Same reasoning for the five-membered 1,2-related shape (e.g.
        # pyrazolidin-3-one) -- it looks aldehyde- or ether-shaped
        # depending on the element pair instead, but the collision-safety
        # argument is identical.
        if has_five_membered_1_2_ring_ketone_shape(mol):
            return name_ketone(mol)
        # Same reasoning for the 7-membered 1,2-related shape (e.g.
        # 1,2-diazepan-3-one) -- the collision-safety argument is
        # identical, and empirically this shape falls all the way through
        # to the generic aldehyde/ketone fallback further below (misnamed
        # as an aldehyde, or rejected by its heteroatom check) if not
        # claimed here first.
        if has_seven_membered_1_2_ring_ketone_shape(mol):
            return name_ketone(mol)
        # An alkoxide anion (R-O(-), P-72.2.2.2.2) has a formal-charge -1
        # oxygen none of the neutral-oxygen checks below (or `name_alcohol`'s
        # own fallback) expect, so it must be routed first in this branch.
        if has_alkoxide_shape(mol):
            return name_alkoxide(mol)
        # Hydroxylamine (H2N-OH, P-68.3.1.1.1) and its O-substituted
        # derivatives (H2N-O-R) have a nitrogen the ether/carbonyl/alcohol
        # checks below don't expect at all, so it must be routed before all
        # of them. N-substituted forms (R-NH-OH) don't match this shape
        # (see _hydroxylamine.py's module docstring) and fall through to
        # name_amine below instead.
        if has_hydroxylamine_shape(mol):
            return name_hydroxylamine(mol)
        # A cyanate (R-O-C#N, P-6) has a degree-2 oxygen bonded to two
        # carbons, the same shape `has_ether_shape` looks for -- it must be
        # routed here first to let `_cyanate.py` claim it before
        # `_ether.py` would otherwise misname it as a plain ether.
        if has_cyanate_shape(mol):
            return name_cyanate(mol)
        # A plain -O- ether (P-63.2.1) has no suffix, so it must be routed
        # here before the carbonyl/alcohol checks below, none of which
        # accept a degree-2 oxygen at all.
        if has_ether_shape(mol):
            return name_ether(mol)
        # An acetal/ketal carbon (two alkoxy oxygens on the same carbon,
        # P-66.6.5.1) has two ether-type oxygens, not the single one
        # has_ether_shape requires, so it must be routed here before the
        # peroxide/carbonyl/alcohol checks below, none of which accept two
        # degree-2 oxygens on one carbon.
        if has_acetal_shape(mol):
            return name_acetal(mol)
        # A plain -O-O- peroxide (P-63.2.5) also has no suffix and no
        # carbonyl, so with no other check to intercept it, it would
        # otherwise fall all the way through to the alcohol module below
        # (which doesn't accept a degree-2 oxygen either).
        if has_peroxide_shape(mol):
            return name_peroxide(mol)
        # A hydroperoxide (R-O-O-H, P-56.1) has a 'peroxol' suffix rather
        # than no suffix at all, but its degree-1 terminal oxygen doesn't
        # match any check above or `name_alcohol`'s own fallback below, so
        # it must be routed here too, right alongside its R-O-O-R' cousin.
        if has_hydroperoxide_shape(mol):
            return name_hydroperoxide(mol)
        # An anhydride's bridging oxygen (-C(=O)-O-C(=O)-) is also
        # ester-shaped from either acyl carbon's point of view (a carbonyl
        # oxygen plus a second, carbon-bonded oxygen), so it must be routed
        # before has_ester_shape below, which would otherwise misname it.
        if has_anhydride_shape(mol):
            return name_anhydride(mol)
        # A carbamate's carbon (R-O-C(=O)-NH2) is simultaneously
        # ester-shaped (carbonyl + a second, carbon-bonded oxygen) and
        # amide-shaped (carbonyl + a primary -NH2) on the very same carbon,
        # so it must be routed before both has_ester_shape and
        # has_amide_shape below, either of which would otherwise misname it.
        if has_carbamate_shape(mol):
            return name_carbamate(mol)
        # A carboxylate anion's carbon (R-COO-, P-72.2.2.2.1.1) bears a
        # neutral carbonyl oxygen and a formal-charge -1 oxygen, so it must
        # be routed before every check below: its carbonyl half looks
        # aldehyde/ketone-shaped to the generic carbonyl fallback further
        # down, and its charged oxygen would otherwise be rejected outright
        # by every other module here, none of which expect a charged atom.
        if has_carboxylate_shape(mol):
            return name_carboxylate(mol)
        # A thioate anion's carbon (R-CO-S(-)/R-CS-O(-), P-72.2.2.2.1.1) is
        # the chalcogen analogue of a carboxylate anion -- same reasoning,
        # routed here before the thioic-acid/carboxylic-acid checks below,
        # neither of which expect a charged chalcogen.
        if has_thioate_shape(mol):
            return name_thioate(mol)
        # A selenoate anion (R-CO-Se(-)/R-CSe-O(-), P-72.2.2.2.1.1) is the
        # selenium analogue of a thioate anion -- same reasoning, routed
        # here for the same reason.
        if has_selenoate_shape(mol):
            return name_selenoate(mol)
        # A ketone carbonyl directly bonded to a saturated single- or
        # 1,4-two-heteroatom ring's own O/S heteroatom (a lactone, e.g.
        # oxan-2-one/1,4-dioxan-2-one) looks ester-shaped to
        # `has_ester_shape` below (carbonyl + a second, carbon-bonded
        # oxygen) and would otherwise be misnamed by `_ester.py`'s
        # acyclic-only construction -- `_ketone.py` already names this
        # narrow ring shape correctly (heteroatom always locant 1, ketone
        # locant set minimized), so claim it here first, mirroring the
        # same routing already used ahead of `has_amide_shape`.
        if has_hetero_ring_ketone_shape(mol):
            return name_ketone(mol)
        # A carbon bearing both a carbonyl oxygen and a second, carbon-bonded
        # oxygen is an ester (-COO-), which must be routed before the
        # carboxylic-acid/aldehyde/ketone checks below: its carbonyl half
        # would otherwise look aldehyde/ketone-shaped, and (for a rejected,
        # out-of-scope case) its non-carbonyl oxygen would never satisfy the
        # carboxylic acid module's hydroxyl (O-H) requirement anyway.
        if has_diester_shape(mol):
            # P-65.6.3: a diester on a shared diol chain is not a simple
            # multiplied 'oate' suffix -- one ester stays the suffix parent
            # and the other is demoted to an 'acyloxy' prefix, a different
            # construction from `_ester.py`'s single-ester "exactly one"
            # rejection below, so it must be routed first.
            return name_diester_acyloxy(mol)
        if has_ester_shape(mol):
            # P-41/Table 3.3: 'oate' outranks 'one', so an ester whose acyl
            # chain also carries one or more ketones names the ester as the
            # suffix and demotes each ketone to an 'oxo' prefix instead of
            # `_ester.py`'s own "coexisting oxygen" rejection.
            if has_ketone_ester_shape(mol):
                return name_ketone_ester(mol)
            return name_ester(mol)
        # A thioic acid (-CO-SH/-CS-OH, P-65.1.5) has a carbonyl-shaped
        # chalcogen cluster that would otherwise look like a plain
        # carboxylic acid's carbonyl (if the double-bonded atom is O) or
        # trip the aromatic/heteroatom checks in other O-only modules (if
        # it's S) -- route it here, before all of those.
        if has_thioic_acid_shape(mol):
            return name_thioic_acid(mol)
        # A selenoic acid (-CO-SeH/-CSe-OH, P-65.1.5) is the selenium
        # analogue of a thioic acid -- same reasoning, routed here for the
        # same reason.
        if has_selenoic_acid_shape(mol):
            return name_selenoic_acid(mol)
        # A telluroic acid (-CO-TeH/-CTe-OH, P-65.1.5) is the tellurium
        # analogue of a thioic/selenoic acid -- same reasoning, routed here
        # for the same reason.
        if has_telluroic_acid_shape(mol):
            return name_telluroic_acid(mol)
        # A carbon bearing both a carbonyl and a hydroxyl oxygen is a -COOH
        # group (Table 3.3's most senior suffix here) and must be routed
        # before the aldehyde/ketone/alcohol checks below, which would
        # otherwise misread its carbonyl or hydroxyl half in isolation.
        if has_carboxylic_acid_shape(mol):
            # P-41/Table 3.3: 'oic acid' outranks 'al', so a carboxylic acid
            # that also carries one or more aldehydes names the acid as the
            # suffix and demotes each aldehyde to an 'oxo' prefix instead of
            # `_carboxylic_acid.py`'s own "coexisting oxygen" rejection.
            if has_aldehyde_carboxylic_acid_shape(mol):
                return name_aldehyde_carboxylic_acid(mol)
            # P-41/Table 3.3: 'oic acid' also far outranks 'amine', so a
            # carboxylic acid that also carries a primary amine names the
            # acid as the suffix and demotes the amine to an 'amino' prefix
            # instead of `_carboxylic_acid.py`'s own "coexisting nitrogen"
            # rejection.
            if has_carboxylic_acid_amine_shape(mol):
                return name_carboxylic_acid_amine(mol)
            return name_carboxylic_acid(mol)
        # A carbon bearing both a carbonyl oxygen and a primary-amide
        # nitrogen (-CONH2) is an amide (junior only to the acid/ester
        # suffixes above in Table 3.3) and must be routed before the
        # aldehyde/ketone checks below: an amide carbon looks
        # aldehyde-shaped to `_is_aldehyde_shaped` (it counts only carbon
        # neighbors, ignoring the nitrogen).
        # Urea (H2N-C(=O)-NH2) has a carbonyl carbon with two qualifying
        # nitrogens, the same shape `has_amide_shape` looks for -- routing
        # it first here lets `_urea.py` claim it before `_amide.py` would
        # otherwise (correctly, but unhelpfully) reject it as "not exactly
        # one nitrogen". It must also be checked before the hydrazide check
        # below: semicarbazide (H2N-NH-C(=O)-NH2) is urea-shaped but its
        # amino nitrogen also happens to match hydrazide's -CO-NH-NH2
        # pattern (hydrazide's shape check doesn't look at the carbonyl
        # carbon's other substituents), so urea must claim it first.
        if has_urea_shape(mol):
            return name_urea(mol)
        # A hydrazide carbon (-CO-NH-NH2, P-66.3.1.1) has a carbonyl plus
        # a two-nitrogen chain that has_amide_shape's own single-nitrogen
        # check doesn't match (its first nitrogen has degree 2, not 1), so
        # it wouldn't collide with the amide check below either way -- but
        # routing it first here keeps the two suffix-shaped carbonyl
        # groups together.
        if has_hydrazide_shape(mol):
            return name_hydrazide(mol)
        if has_amide_shape(mol):
            # P-41/Table 3.3: 'amide' outranks 'one', so an amide that also
            # carries one or more ketones names the amide as the suffix and
            # demotes each ketone to an 'oxo' prefix instead of
            # `_amide.py`'s own "coexisting carbonyl" rejection.
            if has_ketone_amide_shape(mol):
                return name_ketone_amide(mol)
            return name_amide(mol)
        # A one-H nitrogen bridging two carbonyl carbons (-C(=O)-NH-C(=O)-)
        # is an imide, junior only to the acid/ester/amide suffixes above
        # in Table 3.3; it must be routed before the generic carbonyl
        # fallback below, which would otherwise misname each acyl carbon as
        # a plain aldehyde.
        if has_imide_shape(mol):
            return name_imide(mol)
        # A carbon bearing a carbonyl oxygen and a halogen (-C(=O)X) is an
        # acyl halide, junior only to the acid/ester/amide suffixes above in
        # Table 3.3; it must be routed before the generic carbonyl fallback
        # below, which would otherwise misname it as a plain aldehyde/ketone
        # and leave the halogen as an ordinary substituent prefix instead.
        if has_acyl_halide_shape(mol):
            return name_acyl_halide(mol)
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
            # P-41/Table 3.3: 'al' outranks 'one', so a molecule combining a
            # terminal aldehyde with one or more ketones names the aldehyde
            # as the suffix and demotes each ketone to an 'oxo' prefix
            # instead of raising the "coexisting carbonyl" rejection either
            # single-shape module would otherwise hit on its own.
            if has_aldehyde_ketone_shape(mol):
                return name_aldehyde_ketone(mol)
            if any(_is_aldehyde_shaped(o) for o in carbonyl_oxygens):
                return name_aldehyde(mol)
            return name_ketone(mol)
        return name_alcohol(mol)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        # An aminide anion (-NH(-), P-72.2.2.2.3) has a formal-charge -1
        # nitrogen none of the neutral-nitrogen checks below (or
        # `name_amine`'s own fallback) expect, so it must be routed first.
        if has_aminide_shape(mol):
            return name_aminide(mol)
        # An azide group (-N3, P-61.7) has no oxygen and three nitrogens
        # name_amine doesn't recognize at all, so it must be routed here
        # before name_amine for the same reason as the nitrile/imine checks
        # below.
        if has_azide_shape(mol):
            return name_azide(mol)
        # An isocyanide group (-NC, P-61.9) similarly has no oxygen and a
        # nitrogen shape (triple-bonded to a carbon) name_amine doesn't
        # recognize -- must be routed here for the same reason.
        if has_isocyanide_shape(mol):
            return name_isocyanide(mol)
        # An isothiocyanate group (-N=C=S, P-61.8) has a C=N double bond
        # on its own skeletal carbon that would otherwise confuse the
        # imine check further below (a different, unrelated C=N shape) --
        # must be routed here first, same reason as diazo below.
        if has_isothiocyanate_shape(mol):
            return name_isothiocyanate(mol)
        # An isoselenocyanate group (-N=C=Se, P-61.8) has the same C=N
        # issue as isothiocyanate above -- must be routed here first too.
        if has_isoselenocyanate_shape(mol):
            return name_isoselenocyanate(mol)
        # An isotellurocyanate group (-N=C=Te, P-61.8) has the same C=N
        # issue as the other chalcogen analogues above -- must be routed
        # here first too.
        if has_isotellurocyanate_shape(mol):
            return name_isotellurocyanate(mol)
        # A diazo group (=N2, P-61.4) has a C=N double bond on its own
        # skeletal carbon that would otherwise confuse the imine check
        # below (a different, unrelated C=N shape) -- must be routed here
        # first.
        if has_diazo_shape(mol):
            return name_diazo(mol)
        # A diazene/azo skeleton (R-N=N-R', P-68.3.1.3.2) has its own two
        # skeletal nitrogens with no carbon parent chain at all -- must be
        # routed before name_amine for the same reason as the checks above.
        if has_diazene_shape(mol):
            return name_diazene(mol)
        # A symmetric azine (R2C=N-N=CR2, P-68.3.1.2.3) has two C=N double
        # bonds, which would otherwise look like a polyimine to the plain
        # imine check further below (`_imine.py` rejects more than one
        # C=N bond outright) -- must be routed here first.
        if has_azine_shape(mol):
            return name_azine(mol)
        # A hydrazone (R2C=N-NH2, P-68.3.1.2.2) also has two skeletal
        # nitrogens joined by a single bond -- it would otherwise look
        # hydrazine-shaped to the check below (`_hydrazine.py`'s own shape
        # check doesn't look for a C=N bond on the other nitrogen) -- must
        # be routed here first.
        if has_hydrazone_shape(mol):
            return name_hydrazone(mol)
        # A hydrazine skeleton (H2N-NH2, P-68.3.1.2.1) has its own two
        # skeletal nitrogens (single-bonded, not double-bonded like
        # diazene above) with no carbon parent chain at all -- must be
        # routed before name_amine for the same reason as the checks
        # above.
        if has_hydrazine_shape(mol):
            return name_hydrazine(mol)
        # A thiocyanate (R-S-C#N, P-6) has a nitrile-shaped -C#N group with
        # a sulfur instead of a carbon on its other side -- `_nitrile.py`'s
        # own validation already correctly (but unhelpfully) rejects the
        # stray sulfur, so this must be routed first to let
        # `_thiocyanate.py` claim it instead.
        if has_thiocyanate_shape(mol):
            return name_thiocyanate(mol)
        # A selenocyanate (R-Se-C#N, P-6) is the selenium analogue of
        # thiocyanate above -- same reasoning, routed here for the same
        # reason.
        if has_selenocyanate_shape(mol):
            return name_selenocyanate(mol)
        # A tellurocyanate (R-Te-C#N, P-6) is the tellurium analogue of
        # thiocyanate/selenocyanate above -- same reasoning, routed here
        # for the same reason.
        if has_tellurocyanate_shape(mol):
            return name_tellurocyanate(mol)
        # A nitrile nitrogen (-C#N, P-66.5) has no oxygen, so it reaches this
        # branch alongside plain amines; it must be routed here before
        # name_amine, which doesn't recognize a triple-bonded nitrogen at all.
        if has_nitrile_shape(mol):
            return name_nitrile(mol)
        # An amidine carbon (-C(=NH)NH2, P-66.4.1.1) has a C=N double bond
        # that would otherwise look imine-shaped to the check below (a
        # different, unrelated interpretation of the same C=N bond) -- must
        # be routed here first.
        # Guanidine (HN=C(NH2)2, P-66.4.1) has a carbon with the same
        # imino+amino nitrogen shape `has_amidine_shape` looks for (just
        # with two amino nitrogens instead of one) -- `_amidine.py`'s own
        # docstring explicitly anticipates this "geminal diamidine" shape
        # as one its validation rejects, so it must be routed here first
        # to let `_guanidine.py` claim it before that (correct, but
        # unhelpful) rejection.
        if has_guanidine_shape(mol):
            return name_guanidine(mol)
        if has_amidine_shape(mol):
            return name_amidine(mol)
        # A plain (non-oxime) imine (C=N, P-62.3) has no oxygen, so it
        # reaches this branch alongside plain amines -- an oxime (which
        # does have an oxygen) is already caught by the has_simple_imine_shape
        # check much earlier, before the "any O atom" branch.
        if has_simple_imine_shape(mol):
            return name_imine(mol)
        # A plain saturated monocyclic amine whose sole ring nitrogen
        # carries one substituent (e.g. 1-methylpiperidine, or the N,N/
        # N,S piperazine/thiomorpholine siblings -- the N,O morpholine
        # sibling has an oxygen and is instead routed in the oxygen-gated
        # branch above) has the ring itself as the parent hydride, unlike
        # every other shape reaching `name_amine` below (which always
        # treats an acyclic chain as the parent) -- `name_amine` itself
        # explicitly defers this shape, so it must be routed here first.
        # These oxygen-free variants never collide with
        # `_hidden_amide_ketone.py`'s acyl-on-ring-nitrogen path, which
        # requires the substituent's own carbonyl oxygen and is routed
        # separately in the oxygen-gated branch above.
        if has_ring_amine_shape(mol):
            return name_ring_amine(mol)
        return name_amine(mol)
    if has_thione_shape(mol):
        # A thione (C=S, P-64.6.1) has no oxygen or nitrogen, so it only
        # reaches this branch once both are ruled out above. Its own shape
        # check is precise (a real C=S double bond), unlike thiol's/
        # sulfide's own loose "any sulfur atom" checks, so it's safe to
        # check here regardless of order relative to them.
        return name_thione(mol)
    if has_disulfide_shape(mol):
        # A disulfide (R-S-S-R') has two sulfurs -- it would otherwise
        # look thiol-shaped to the check below (that check just looks for
        # the presence of any sulfur atom) -- must be routed here first.
        return name_disulfide(mol)
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
    if has_selone_shape(mol):
        # A selone (C=Se, P-64.6.1) has the same precise-shape reasoning
        # as thione above (a real C=Se double bond), so it's safe to check
        # here regardless of order relative to the selenide/selenol chain.
        return name_selone(mol)
    if has_diselenide_shape(mol):
        # A diselenide (R-Se-Se-R') has two seleniums -- it would
        # otherwise look selenol-shaped to the check below (that check
        # just looks for the presence of any selenium atom) -- must be
        # routed here first.
        return name_diselenide(mol)
    if has_selenide_shape(mol):
        # A plain -Se- selenide (P-63.2.1) has no suffix, so it must be
        # routed here before has_selenol_shape below for the same reason
        # as has_sulfide_shape above (that check doesn't look at degree
        # at all, so a degree-2 selenide would otherwise raise the wrong
        # error inside `name_selenol`'s degree-1 validation).
        return name_selenide(mol)
    if has_selenol_shape(mol):
        # A selenol (-SeH, P-63.1.1) is the next chalcogen analogue after
        # a thiol -- has neither O, N, nor S, so it only reaches this
        # branch once all three are ruled out above.
        return name_selenol(mol)
    if has_tellone_shape(mol):
        # A tellone (C=Te, P-64.6.1) has the same precise-shape reasoning
        # as thione/selone above, so it's safe to check here regardless of
        # order relative to the telluride/tellurol chain.
        return name_tellone(mol)
    if has_ditelluride_shape(mol):
        # A ditelluride (R-Te-Te-R') has two telluriums -- it would
        # otherwise look tellurol-shaped to the check below (that check
        # just looks for the presence of any tellurium atom) -- must be
        # routed here first.
        return name_ditelluride(mol)
    if has_telluride_shape(mol):
        # A plain -Te- telluride (P-63.2.1) has no suffix, so it must be
        # routed here before has_tellurol_shape below for the same reason
        # as has_selenide_shape above (that check doesn't look at degree
        # at all, so a degree-2 telluride would otherwise raise the wrong
        # error inside `name_tellurol`'s degree-1 validation).
        return name_telluride(mol)
    if has_tellurol_shape(mol):
        # A tellurol (-TeH, P-63.1.1) is the next chalcogen analogue after
        # a selenol -- has neither O, N, S, nor Se, so it only reaches
        # this branch once all four are ruled out above.
        return name_tellurol(mol)

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
    # A naphthalene skeleton with one adjacent ring-atom pair saturated
    # (1,2- or 1,4-dihydronaphthalene) has one fully aromatic ring and one
    # partially reduced ring, so find_aromatic_fused_core above correctly
    # rejects it (not every ring atom is aromatic); it must be routed here
    # before find_bicyclic_core below, which would otherwise misdetect its
    # carbon skeleton (graph-isomorphic to a saturated bicyclic) and reject
    # it with a confusing "unsaturated bicyclics ... not supported" error
    # instead of this module's own name.
    if num_rings == 2:
        dihydro_core = find_dihydronaphthalene_core(mol)
        if dihydro_core is not None:
            return name_dihydronaphthalene(mol, dihydro_core)
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
        unsaturated_ring = find_cyclic_unsaturated_core(mol)
        if unsaturated_ring is not None:
            return name_cyclic_unsaturated(mol, unsaturated_ring)
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
    # A hub ring with three distinct spiro atoms, each fused to its own
    # terminal ring (P-24.2.3/SP-1.5's minimal branched-polyspiro shape),
    # must be routed here before find_bicyclic_core/find_polycyclic_core
    # below: those only understand shared-bridgehead-atom bridging, not a
    # spiro-atom hub, so they would otherwise fail to recognize this shape
    # and fall through to the generic "not supported" error.
    branched_hub = find_branched_polyspiro_hub(mol)
    if branched_hub is not None:
        return name_branched_polyspiro(mol, branched_hub)
    bicyclic_core = find_bicyclic_core(mol)
    if bicyclic_core is not None:
        return name_bicycloalkane(mol, bicyclic_core)
    for ring_count in (3, 4, 5, 6):
        core = find_polycyclic_core(mol, ring_count)
        if core is not None:
            return name_polycycloalkane(mol, core, ring_count)
    propellane_core = find_propellane_core(mol)
    if propellane_core is not None:
        return name_propellane(mol, propellane_core)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
