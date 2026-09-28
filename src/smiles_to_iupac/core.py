from rdkit import Chem

from ._zwitterion import has_zwitterion_shape, name_zwitterion
from ._salt import has_salt_shape, name_salt
from ._hydrohalide_salt import has_hydrohalide_salt_shape, name_hydrohalide_salt
from ._hydrate_adduct import has_hydrate_adduct_shape, name_hydrate_adduct
from ._acyclic import name_acyclic_alkane
from ._acyl_halide import has_acyl_halide_shape, name_acyl_halide
from ._anhydride import has_anhydride_shape, name_anhydride
from ._carbamate import has_carbamate_shape, name_carbamate
from ._alcohol import name_alcohol
from ._alcohol_amine import has_alcohol_amine_shape, name_alcohol_amine
from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._aldehyde import name_aldehyde
from ._aldehyde_amine import has_aldehyde_amine_shape, name_aldehyde_amine
from ._ketone_amine import has_ketone_amine_shape, name_ketone_amine
from ._amino_acid import has_amino_acid_shape, name_amino_acid
from ._histidine import has_histidine_shape, name_histidine
from ._proline import has_proline_shape, name_proline
from ._carboxylic_acid_amine import has_carboxylic_acid_amine_shape, name_carboxylic_acid_amine
from ._aldehyde_carboxylic_acid import (
    has_aldehyde_carboxylic_acid_shape,
    name_aldehyde_carboxylic_acid,
)
from ._carboxylic_acid_amide import has_carboxylic_acid_amide_shape, name_carboxylic_acid_amide
from ._aldehyde_ketone import has_aldehyde_ketone_shape, name_aldehyde_ketone
from ._carbohydrate import (
    has_cyclic_aldofuranose_shape,
    has_cyclic_aldopyranose_shape,
    has_cyclic_ketohexopyranose_shape,
    has_open_chain_2_ketose_shape,
    has_open_chain_aldose_shape,
    name_cyclic_aldofuranose,
    name_cyclic_aldopyranose,
    name_cyclic_ketohexopyranose,
    name_open_chain_2_ketose,
    name_open_chain_aldose,
)
from ._inositol import has_inositol_shape, name_inositol
from ._acetal import has_acetal_shape, name_acetal
from ._amide import has_amide_shape, name_amide
from ._amide_amine import has_amide_amine_shape, name_amide_amine
from ._ester_amine import has_ester_amine_shape, name_ester_amine
from ._ether_ester import has_ether_ester_shape, name_ether_ester
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
from ._hetero_ring_amine import has_hetero_ring_amine_shape, name_hetero_ring_amine
from ._ring_amine import has_ring_amine_shape, has_ring_amine_sulfonyl_shape, name_ring_amine
from ._amine_oxide import has_amine_oxide_shape, name_amine_oxide
from ._aminide import has_aminide_shape, name_aminide
from ._ammonium import has_ammonium_shape, name_ammonium
from ._ylide import has_nitrogen_ylide_shape, has_pos_ylide_shape, name_nitrogen_ylide, name_pos_ylide
from ._amine_imide import has_amine_imide_shape, name_amine_imide
from ._phosphonium import has_phosphonium_shape, name_phosphonium
from ._oxonium import has_oxonium_shape, name_oxonium
from ._carbenium import has_acylium_shape, has_carbenium_shape, name_acylium, name_carbenium
from ._carbanide import has_carbanide_shape, name_carbanide
from ._benzenide import has_benzenide_shape, name_benzenide
from ._cyclopentadienide import has_cyclopentadienide_shape, name_cyclopentadienide
from ._sulfonium import has_sulfonium_shape, name_sulfonium
from ._diazonium import has_diazonium_shape, name_diazonium
from ._radical import has_radical_shape, name_radical
from ._radical_ion import has_radical_ion_shape, name_radical_ion
from ._aromatic import find_aromatic_fused_core, name_aromatic_fused
from ._bicyclic import find_bicyclic_core, name_bicycloalkane
from ._bridged_aromatic import (
    find_bridged_anthracene_core,
    find_bridged_aromatic_core,
    name_bridged_anthracene,
    name_bridged_aromatic,
)
from ._alkaloid_parent_hydrides import has_alkaloid_morphinan_name, name_alkaloid_morphinan
from ._bridged_alicyclic_parent import has_bridged_steroid_name, name_bridged_steroid_parent
from ._borane import has_simple_borane_shape, name_simple_borane
from ._boronic_acid import has_boronic_acid_shape, name_boronic_acid
from ._borinic_acid import has_borinic_acid_shape, name_borinic_acid
from ._group1_2_organometallic import has_group1_2_organometallic_shape, name_group1_2_organometallic
from ._group13_hydride import (
    has_group13_hydride_shape,
    has_group14_hydride_shape,
    has_group15_hydride_shape,
    name_group13_hydride,
    name_group14_hydride,
    name_group15_hydride,
)
from ._branched_fused_aromatic import has_retained_branched_fused_name, name_retained_branched_fused
from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._selenoate import has_selenoate_shape, name_selenoate
from ._sulfonate import has_sulfonate_shape, name_sulfonate
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
from ._cyclic import find_exocyclic_ylidene_core, name_cycloalkane, name_exocyclic_ylidene
from ._disjoint_ring_substituents import find_disjoint_ring_pair_core, name_disjoint_ring_pair
from ._cyclic_unsaturated import find_cyclic_unsaturated_core, name_cyclic_unsaturated
from ._dihydro_aromatic import (
    find_decahydronaphthalene_core,
    find_dihydronaphthalene_core,
    find_partially_unsaturated_naphthalene_core,
    name_decahydronaphthalene,
    name_dihydronaphthalene,
    name_partially_unsaturated_naphthalene,
)
from ._diester_acyloxy import has_diester_shape, name_diester_acyloxy
from ._ester import has_ester_shape, name_ester
from ._cyanate import has_cyanate_shape, name_cyanate
from ._ether import has_ether_shape, name_ether
from ._ether_amine import has_ether_amine_shape, name_ether_amine
from ._ether_aldehyde import has_ether_aldehyde_shape, name_ether_aldehyde
from ._ether_amide import has_ether_amide_shape, name_ether_amide
from ._ether_hydroperoxide import has_ether_hydroperoxide_shape, name_ether_hydroperoxide
from ._ether_ketone import has_ether_ketone_shape, name_ether_ketone
from ._ether_thiol import has_ether_thiol_shape, name_ether_thiol
from ._fullerene import has_fullerene_name, name_fullerene
from ._nucleoside import has_nucleoside_name, name_nucleoside
from ._nucleotide import has_nucleotide_name, name_nucleotide
from ._metallacycle import has_metallacycle_shape, name_metallacycle
from ._metallocene import has_metallocene_name, name_metallocene
from ._fused_hetero_ring_oxide import has_fused_hetero_ring_oxide_shape, name_fused_hetero_ring_oxide
from ._hetero_ring_oxide import has_hetero_ring_oxide_shape, name_hetero_ring_oxide
from ._pyridinone import has_pyridinone_shape, name_pyridinone
from ._pyrimidinedione import has_pyrimidinedione_shape, name_pyrimidinedione
from ._pyrimidinone import has_pyrimidinone_shape, name_pyrimidinone
from ._homo_steroid import has_homo_steroid_shape, name_homo_steroid
from ._cyclo_steroid import has_cyclo_steroid_shape, name_cyclo_steroid
from ._dinor_steroid import has_dinor_steroid_shape, name_dinor_steroid
from ._nor_steroid import has_nor_steroid_shape, name_nor_steroid
from ._seco_steroid import has_seco_steroid_shape, name_seco_steroid
from ._steroid_parent_hydrides import (
    has_steroid_aromatic_a_ring_name,
    has_steroid_parent_hydride_name,
    has_steroid_unsaturated_name,
    name_steroid_aromatic_a_ring,
    name_steroid_parent_hydride,
    name_steroid_unsaturated,
)
from ._heteroaromatic_fused import has_retained_heteroaromatic_fused_name, name_retained_heteroaromatic_fused
from ._phenanthroline_naphthyridine import (
    find_phenanthroline_naphthyridine_core,
    name_phenanthroline_naphthyridine,
)
from ._polycyclic_component_fusion import (
    has_polycyclic_component_fusion_name,
    name_polycyclic_component_fusion,
)
from ._pyridine_bicyclic_fusion import has_pyridine_bicyclic_fusion_name, name_pyridine_bicyclic_fusion
from ._anthracene_fusion import has_anthracene_fusion_name, name_anthracene_fusion
from ._tetracene_fusion import has_tetracene_fusion_name, name_tetracene_fusion
from ._pentacene_fusion import has_pentacene_fusion_name, name_pentacene_fusion
from ._hexacene_fusion import has_hexacene_fusion_name, name_hexacene_fusion
from ._heptacene_fusion import has_heptacene_fusion_name, name_heptacene_fusion
from ._phenanthrene_fusion import has_phenanthrene_fusion_name, name_phenanthrene_fusion
from ._pyrene_fusion import has_pyrene_fusion_name, name_pyrene_fusion
from ._chrysene_fusion import has_chrysene_fusion_name, name_chrysene_fusion
from ._picene_fusion import has_picene_fusion_name, name_picene_fusion
from ._pentaphene_fusion import has_pentaphene_fusion_name, name_pentaphene_fusion
from ._triphenylene_fusion import has_triphenylene_fusion_name, name_triphenylene_fusion
from ._fluoranthene_fusion import has_fluoranthene_fusion_name, name_fluoranthene_fusion
from ._aceanthrylene_fusion import has_aceanthrylene_fusion_name, name_aceanthrylene_fusion
from ._acephenanthrylene_fusion import (
    has_acephenanthrylene_fusion_name,
    name_acephenanthrylene_fusion,
)
from ._hetero_monocyclic import (
    has_hetero_monocyclic_name,
    has_hetero_monocyclic_substituent_name,
    has_pyran_indicated_hydrogen_name,
    name_hetero_monocyclic,
    name_hetero_monocyclic_substituent,
    name_pyran_indicated_hydrogen,
)
from ._didehydro_ring import has_didehydro_ring_name, name_didehydro_ring
from ._benzo_bis_heterocycle_fusion import (
    has_benzo_bis_heterocycle_fusion_name,
    name_benzo_bis_heterocycle_fusion,
)
from ._bridgehead_heteroatom_fusion import (
    has_bridgehead_heteroatom_fusion_name,
    name_bridgehead_heteroatom_fusion,
)
from ._hydroxylamine import has_hydroxylamine_shape, name_hydroxylamine
from ._imine import has_simple_imine_shape, name_imine
from ._dipole_oxide import (
    has_nitrile_oxide_shape,
    has_nitrone_shape,
    name_nitrile_oxide,
    name_nitrone,
)
from ._isotope import has_isotope_shape, name_isotope
from ._isotope_alcohol import has_isotope_alcohol_shape, name_isotope_alcohol
from ._isotope_carboxylic_acid import has_isotope_carboxylic_acid_shape, name_isotope_carboxylic_acid
from ._isotope_ketone import has_isotope_ketone_shape, name_isotope_ketone
from ._pyridine_bis_heterocycle_fusion import (
    has_pyridine_bis_heterocycle_fusion_name,
    name_pyridine_bis_heterocycle_fusion,
)
from ._pyridine_heterocycle_fusion import (
    has_pyridine_heterocycle_fusion_name,
    name_pyridine_heterocycle_fusion,
)
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
from ._nitrate_ester import has_nitrate_ester_shape, name_nitrate_ester
from ._carbonic_acid import has_carbonic_acid_shape, name_carbonic_acid
from ._nitrite_ester import has_nitrite_ester_shape, name_nitrite_ester
from ._nitro import has_nitro_shape, name_nitro
from ._nitroso import has_nitroso_shape, name_nitroso
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._cyclophane import has_cyclophane_name, name_cyclophane
from ._naphthalene_benzene_phane import (
    has_naphthalene_benzene_phane_name,
    name_naphthalene_benzene_phane,
)
from ._phosphane import has_simple_phosphane_shape, name_simple_phosphane
from ._functional_replacement_oxoacid import (
    has_functional_replacement_oxoacid_shape,
    name_functional_replacement_oxoacid,
)
from ._phosphanone import has_phosphanone_shape, name_phosphanone
from ._phosphate import has_phosphate_shape, name_phosphate
from ._dinuclear_oxoacid import has_dinuclear_oxoacid_shape, name_dinuclear_oxoacid
from ._phosphite import has_phosphite_shape, name_phosphite
from ._sulfate import has_sulfate_shape, name_sulfate
from ._sulfite import has_sulfite_shape, name_sulfite
from ._phosphonic_acid import has_phosphonic_acid_shape, name_phosphonic_acid
from ._phosphinic_acid import has_phosphinic_acid_shape, name_phosphinic_acid
from ._phosphane_chain import has_phosphane_chain_shape, name_phosphane_chain
from ._peri_fused_aromatic import has_retained_peri_fused_name, name_retained_peri_fused
from ._fluorene_parent import has_fluorene_parent_name, name_fluorene_parent
from ._cyclopenta_naphthalene import has_cyclopenta_naphthalene_name, name_cyclopenta_naphthalene
from ._hydroperoxide import has_hydroperoxide_shape, name_hydroperoxide
from ._hydroperoxide_amine import has_hydroperoxide_amine_shape, name_hydroperoxide_amine
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
from ._ring_assembly_chain import find_ring_assembly_chain_core, name_ring_assembly_chain
from ._ring_assembly_ylidene import find_ring_assembly_ylidene_core, name_ring_assembly_ylidene
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
from ._thiol_amine import has_thiol_amine_shape, name_thiol_amine
from ._tricyclic import find_propellane_core, name_propellane
from ._unsaturated import name_acyclic_unsaturated
from ._von_baeyer_heteroatom import (
    has_mixed_element_heteroatom_shape as has_mixed_bicyclic_heteroatom_shape,
    has_mixed_element_heteroatom_shape_polycyclic,
    has_multi_ring_heteroatom_shape as has_multi_bicyclic_heteroatom_shape,
    has_multi_ring_heteroatom_shape_polycyclic,
    has_single_ring_heteroatom_shape as has_single_bicyclic_heteroatom_shape,
    has_single_ring_heteroatom_shape_polycyclic,
    name_von_baeyer_heteroatom,
    name_von_baeyer_heteroatom_mixed,
    name_von_baeyer_heteroatom_mixed_polycyclic,
    name_von_baeyer_heteroatom_multi,
    name_von_baeyer_heteroatom_multi_polycyclic,
    name_von_baeyer_heteroatom_polycyclic,
)


def _is_aldehyde_shaped(carbonyl_oxygen):
    (carbon,) = carbonyl_oxygen.GetNeighbors()
    return carbon.GetAtomicNum() == 6 and sum(1 for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6) == 1


def smiles_to_iupac(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")

    # The 7 retained nucleoside names (P-105.1) are recognized by exact
    # whole-molecule match, so they must be routed before every other
    # branch below: a purine/pyrimidine base's fused-ring nitrogen pattern
    # crashes the ring-assembly-chain detector (thymidine/uridine, an
    # unrelated pre-existing bug) or gets misrouted as an amino-acid/
    # aromatic-ring shape long before any ring-count-specific check runs.
    if has_nucleoside_name(mol):
        return name_nucleoside(mol)

    # The 7 retained nucleotide names (P-106.1) are likewise recognized by
    # exact whole-molecule match and must be routed right after the
    # nucleoside case above, for the same dispatch-ordering reason: the
    # phosphate ester's own oxygens would otherwise reach `_phosphate.py`'s
    # generic dispatch, which has no path for a nucleoside-shaped R group.
    if has_nucleotide_name(mol):
        return name_nucleotide(mol)

    # The 7 retained metallocene names (P-69.2.7) are likewise recognized
    # by exact whole-molecule match and must be routed here for the same
    # reason: a bare metal atom plus two disjoint cyclopentadienyl/
    # cyclopentadienide fragments has no shape any other branch below
    # expects, and reaches a radical/carbanide dispatch that rejects the
    # multi-fragment SMILES outright long before any ring-count check.
    if has_metallocene_name(mol):
        return name_metallocene(mol)

    # A bare metallacyclic parent hydride (P-69.4's skeletal-replacement
    # ring, e.g. '1-titanacyclobutane') has a transition-metal ring atom
    # too, for the same reason as the metallocene case above: no other
    # branch below recognizes a non-carbon ring atom, so it must be routed
    # here before the generic cycloalkane dispatch's own heteroatom
    # rejection would otherwise claim it.
    if has_metallacycle_shape(mol):
        return name_metallacycle(mol)

    # A chalcogen ring-oxide (P-62.5's functional-class "oxide" pattern,
    # not limited to acyclic amines) breaks the ring's own aromaticity as
    # RDKit perceives it, so it must be routed here before any ring-shape
    # or aromatic dispatch below ever gets a chance to reject it outright.
    if has_hetero_ring_oxide_shape(mol):
        return name_hetero_ring_oxide(mol)

    # The fused-bicyclic analogue of the chalcogen ring-oxide above (e.g.
    # benzothiophene 1-oxide) has the same aromaticity-breaking shape, and
    # would otherwise be misrouted into the von Baeyer bicyclic-heteroatom
    # dispatch further below -- routed here, right alongside its
    # single-ring sibling.
    if has_fused_hetero_ring_oxide_shape(mol):
        return name_fused_hetero_ring_oxide(mol)

    # The pyridinone tautomer (P-31.1.4.3.4's indicated-hydrogen oxo form)
    # keeps its ring-carbon aromatic despite the exocyclic oxo, so it must
    # be routed here before `_ketone.py`'s own generic aryl-ketone
    # rejection below ever gets a chance to claim it.
    if has_pyridinone_shape(mol):
        return name_pyridinone(mol)

    # The pyrimidinone tautomer (a second ring nitrogen alongside the
    # pyridinone shape above) is checked right after it, for the same
    # aromatic-aryl-ketone dispatch-ordering reason.
    if has_pyrimidinone_shape(mol):
        return name_pyrimidinone(mol)

    # The uracil/thymine diketo tautomer (both ring nitrogens carrying
    # their own indicated hydrogen, P-58.2.2's parenthesized multi-locant
    # convention) is checked right after the single-oxo pyrimidinone case
    # above, for the same aromatic-aryl-ketone dispatch-ordering reason.
    if has_pyrimidinedione_shape(mol):
        return name_pyrimidinedione(mol)

    # An amino-acid/betaine-type zwitterion (P-74.1.3's ammonium-nitrogen-
    # prefix-on-a-carboxylate-parent citation order) must be routed here
    # before `has_salt_shape` below: it's a single connected fragment that
    # nonetheless has an ammonium-shaped nitrogen, which `_salt.py`'s own
    # cation loop would otherwise try (and fail) to name as if the entire
    # molecule were a bare ammonium cation, crashing rather than falling
    # through.
    if has_zwitterion_shape(mol):
        return name_zwitterion(mol)

    # A multi-fragment SMILES (P-77 salts) must be routed here before every
    # other branch below: those all assume one connected molecule and would
    # reject a foreign atom like sodium outright, never getting a chance to
    # recognize the two fragments as a cation/anion pair.
    if has_salt_shape(mol):
        return name_salt(mol)

    # A bare hydrogen halide fragment (P-77.1.3(3)'s 'hydrochloride'-style
    # general nomenclature, see _hydrohalide_salt.py) must likewise be
    # routed here before every other branch below, for the same reason.
    if has_hydrohalide_salt_shape(mol):
        return name_hydrohalide_salt(mol, smiles_to_iupac)

    # A hydrate adduct (P-14.8's em-dash notation, see
    # _hydrate_adduct.py) -- one organic fragment plus one or more
    # separate water molecules -- must likewise be routed here before
    # every other branch below, for the same reason.
    if has_hydrate_adduct_shape(mol):
        return name_hydrate_adduct(mol, smiles_to_iupac)

    # P-103.1.1.1: a common amino acid's retained name + L/D descriptor
    # must be routed here, before every ring-count/functional-group
    # dispatch branch below: a side chain recognized by `_amino_acid.py`'s
    # table can carry its own extra nitrogen (lysine, arginine) or its
    # own ring (phenylalanine, tyrosine, tryptophan), which would
    # otherwise be misrouted first -- confirmed empirically: arginine's
    # guanidino C=N was caught by `_imine.py`'s dispatch and tryptophan's
    # indole ring by the bicyclic-heteroatom dispatch, both well before
    # this check's original position further down ever ran.
    if has_amino_acid_shape(mol):
        return name_amino_acid(mol)

    # An unbranched chain of 3-6 disjoint mancude rings (aromatic benzo,
    # saturated cycloalkane, or pyridine -- P-28.3) -- e.g. terphenyl,
    # tercyclopropane, terpyridine -- must be routed here before every
    # single-heteroatom dispatch branch below: a pyridine-ring assembly's
    # own nitrogen atoms would otherwise reach `_amine.py`'s own
    # nitrogen-presence branch (which has no ring-assembly awareness) long
    # before this shape's own num_rings==3..6 check would run if it stayed
    # down with the other ring modules further below.
    if 3 <= mol.GetRingInfo().NumRings() <= 6:
        ring_assembly_chain_core = find_ring_assembly_chain_core(mol)
        if ring_assembly_chain_core is not None:
            return name_ring_assembly_chain(mol, ring_assembly_chain_core)

    # Two disjoint (unfused) identical rings joined by a single bond -- e.g.
    # biphenyl, bipyridine, bifuran -- must likewise be routed here before
    # every single-heteroatom dispatch branch below, for the same reason as
    # the 3-6-ring case just above (and before find_aromatic_fused_core
    # further down: that function only checks each SSSR ring is a
    # 6-membered aromatic carbocycle and doesn't require the rings to be
    # fused, so it would otherwise claim the benzo case too and then fail in
    # name_aromatic_fused's ring-fusion-graph validation).
    if mol.GetRingInfo().NumRings() == 2:
        ring_assembly_core = find_ring_assembly_core(mol)
        if ring_assembly_core is not None:
            return name_ring_assembly(mol, ring_assembly_core)

    # An isotopically labeled hydroxyl oxygen and/or skeletal carbon
    # combined with the '-ol' suffix (P-82.5.1/P-82.5.2) must be routed
    # here before `_isotope.py`'s own plain chain/methane path just below,
    # which rejects any oxygen outright.
    if has_isotope_alcohol_shape(mol):
        return name_isotope_alcohol(mol)

    # An isotopically labeled carbonyl oxygen and/or skeletal carbon
    # combined with the '-one' suffix (P-82.5.1/P-82.5.2) must likewise be
    # routed here before `_isotope.py`'s own plain path.
    if has_isotope_ketone_shape(mol):
        return name_isotope_ketone(mol)

    # An isotopically labeled carboxyl oxygen and/or skeletal carbon
    # combined with the '-oic acid' suffix (P-82.5.1/P-82.5.2) must
    # likewise be routed here before `_isotope.py`'s own plain path.
    if has_isotope_carboxylic_acid_shape(mol):
        return name_isotope_carboxylic_acid(mol)

    # An isotopically labeled atom (P-82's isotope descriptor nomenclature)
    # must be routed here before every other branch below: RDKit represents
    # an isotopically substituted hydrogen (e.g. 2H) as its own explicit
    # atom (atomic number 1), which none of the other branches recognize at
    # all -- every one of them would reject it outright as an unsupported
    # heteroatom.
    if has_isotope_shape(mol):
        return name_isotope(mol)

    # A radical ion on an ionic suffix group (P-75.3.1's 'aminiumyl'
    # radical cation) has both a charge and a radical electron on the
    # same nitrogen -- checked ahead of `has_radical_shape` below, whose
    # own broader "any nonzero radical electron count" check would
    # otherwise claim it first and misroute it into the plain-radical
    # dispatch, which rejects any charged atom outright.
    if has_radical_ion_shape(mol):
        return name_radical_ion(mol)

    # A radical center (P-71.2.1.1's 'yl' radical naming) must be routed
    # here before every other branch below: none of them recognize a
    # nonzero radical electron count at all -- an unbranched-chain or
    # monocyclic-ring radical would otherwise fall straight through to the
    # plain alkane/cycloalkane dispatch further down, which doesn't know a
    # hydrogen is missing.
    if has_radical_shape(mol):
        return name_radical(mol)

    # A nitrogen ylide (P-74.2.1.1.1's zwitterionic anion-carbon-parent
    # naming) has its own charged nitrogen too -- checked before
    # `has_ammonium_shape` below, whose own check already explicitly
    # excludes this shape (a second charged atom elsewhere) rather than
    # naming it, per that module's own docstring.
    if has_nitrogen_ylide_shape(mol):
        return name_nitrogen_ylide(mol)

    # A phosphorus/oxygen/sulfur ylide (P-74.2.1.1.2/.3/.4) has its own
    # charged cation atom too -- checked before `has_phosphonium_shape`/
    # `has_oxonium_shape`/`has_sulfonium_shape` below for the same reason
    # as the nitrogen ylide above.
    if has_pos_ylide_shape(mol):
        return name_pos_ylide(mol)

    # An amine imide (P-74.2.1.3's zwitterionic hydrazinium-ide naming)
    # has two charged nitrogens (one +1, one -1) -- checked before
    # `has_ammonium_shape` below for the same reason as the nitrogen
    # ylide above: its own +1 nitrogen would otherwise match ammonium's
    # shape check and get misnamed as a plain quaternary ammonium,
    # silently dropping the -1 nitrogen fragment.
    if has_amine_imide_shape(mol):
        return name_amine_imide(mol)

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

    # An acylium cation (P-73.2.3.1's 'oylium'/'ylium' suffix naming) has
    # a charged carbon too, checked ahead of the plain carbenium case
    # below since a C=O double bond gives the cation carbon degree 2, not
    # `has_carbenium_shape`'s own required degree 3 -- the two shapes
    # never overlap, so order between them doesn't otherwise matter.
    if has_acylium_shape(mol):
        return name_acylium(mol)

    # A carbenium cation (P-73.2.2.1.1's 'ylium' suffix naming) has a
    # charged carbon too, for the same reason as ammonium above -- routed
    # here, unconditionally, before every other branch.
    if has_carbenium_shape(mol):
        return name_carbenium(mol)

    # The cyclopentadienide anion (P-72.2.2.1's ring worked example) has a
    # charged ring carbon too, but `_carbanide.py` below is explicitly
    # acyclic-only and would reject any ring outright -- so this shape must
    # be routed here first.
    if has_cyclopentadienide_shape(mol):
        return name_cyclopentadienide(mol)

    # The benzenide anion (phenyl anion, P-72.2.2.1's other ring worked
    # example) has a charged ring carbon too, for the same reason as the
    # cyclopentadienide case above -- routed here, right alongside it.
    if has_benzenide_shape(mol):
        return name_benzenide(mol)

    # A carbanion center (P-72.2.2.1's '-ide' suffix naming) has a charged
    # carbon too, the anionic mirror of carbenium above -- routed here,
    # unconditionally, before every other branch.
    if has_carbanide_shape(mol):
        return name_carbanide(mol)

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

    # A phosphonic acid (P-67.1.1's R-P(=O)(OH)2) has two hydroxyl oxygens
    # on phosphorus that `_phosphanone.py`'s own phosphine-oxide shape
    # doesn't expect -- `has_phosphanone_shape` would otherwise still
    # match (it only checks for a single P=O) and then fail inside
    # `name_phosphanone`'s validation, so this must be routed first.
    if has_phosphonic_acid_shape(mol):
        return name_phosphonic_acid(mol)

    # A phosphinic acid (P-67.1.1's R2-P(=O)-OH) has the same
    # phosphanone-shape collision as phosphonic acid above (its own
    # hydroxyl oxygen isn't expected by `_phosphanone.py`), so it must be
    # routed here for the same reason.
    if has_phosphinic_acid_shape(mol):
        return name_phosphinic_acid(mol)

    # Diphosphoric acid (P-67.2.1's own preselected dinuclear-acid name)
    # has each phosphorus individually shaped like a phosphate ester (the
    # other phosphorus group standing in as the "R" of a P-O-R ester
    # oxygen), so `has_phosphate_shape` would otherwise also match it and
    # then fail inside `name_phosphate`'s own single-phosphorus-only
    # validation -- must be routed here first.
    if has_dinuclear_oxoacid_shape(mol):
        return name_dinuclear_oxoacid(mol)

    # A phosphate ester (P-67.1.3.2's P(=O)(OR)3) has three P-O-R ester
    # oxygens that `_phosphanone.py`'s own phosphine-oxide shape doesn't
    # expect -- `has_phosphanone_shape` would otherwise still match (it
    # only checks for a single P=O) and then fail inside
    # `name_phosphanone`'s validation, so this must be routed first, same
    # reason as phosphonic/phosphinic acid above.
    if has_phosphate_shape(mol):
        return name_phosphate(mol)

    # A phosphite ester (P-67.1.3.2's P(OR)3, no P=O) has three P-O-R
    # ester oxygens that `_phosphane.py`'s own plain-phosphane shape
    # doesn't expect (it rejects any heteroatom besides its own
    # phosphorus outright) -- must be routed here first, same reason as
    # phosphate above.
    if has_phosphite_shape(mol):
        return name_phosphite(mol)

    # A sulfate ester (P-67.1.3.2's S(=O)(=O)(OR)2) has two S-O-R ester
    # oxygens no other sulfur module expects (they all assume a direct
    # S-C bond) -- must be routed before any of them for the same
    # ether-oxygen-rejection reason phosphate/phosphite are routed early.
    if has_sulfate_shape(mol):
        return name_sulfate(mol)

    # A sulfite ester (P-67.1.3.2's S(=O)(OR)2, one fewer double-bonded
    # oxygen than sulfate) needs the same early routing, for the same
    # ether-oxygen-rejection reason as sulfate above.
    if has_sulfite_shape(mol):
        return name_sulfite(mol)

    # A nitrate ester (P-67.1.3.2's O-NO2) has an N-O-R ester oxygen
    # `_nitro.py`'s own nitrogen shape doesn't expect (that module
    # requires a direct N-C bond) -- routed here for the same early-ester
    # reasoning as sulfate/sulfite above.
    if has_nitrate_ester_shape(mol):
        return name_nitrate_ester(mol)

    # A nitrite ester (P-67.1.3.2's O-N=O) needs the same early routing
    # as nitrate above, for the same N-O-R ester-oxygen reason.
    if has_nitrite_ester_shape(mol):
        return name_nitrite_ester(mol)

    # Carbonic acid or one of its esters (P-65.2.1's O=C(OR)(OR')) has a
    # central carbon with two -O-R/-OH oxygens neither `_ether.py` (which
    # rejects an oxygen bonded to more than one heavy atom outright) nor
    # `_carboxylic_acid.py` (which expects exactly one -OH, not two)
    # expects -- routed here for the same early-ester reasoning as
    # sulfate/nitrate above.
    if has_carbonic_acid_shape(mol):
        return name_carbonic_acid(mol)

    # A phosphine oxide (P-68.3.2.3.1's '-phosphanone' suffix, R-P(=O)<)
    # has its own phosphorus-bonded oxygen that `_phosphane.py` doesn't
    # expect at all (that module rejects any heteroatom besides its own
    # phosphorus outright) -- must be routed here first, before
    # has_simple_phosphane_shape below, for the same reason as
    # has_phosphane_chain_shape above.
    if has_phosphanone_shape(mol):
        return name_phosphanone(mol)

    # Thiophosphoric acid (P-67.1.2's own preselected infix-modified
    # oxoacid name) has a phosphorus with 4 substituents (=S plus three
    # -OH), which `_phosphane.py` rejects outright (more than three
    # substituents) -- must be routed here first.
    if has_functional_replacement_oxoacid_shape(mol):
        return name_functional_replacement_oxoacid(mol)

    # A phosphorus atom (P-68's phosphane substitutive nomenclature) must
    # be routed here before every other branch below: none of them
    # recognize phosphorus at all, and a phosphane carbon substituent would
    # otherwise reach the plain acyclic-alkane/amine dispatch further down
    # with no phosphorus handling.
    if has_simple_phosphane_shape(mol):
        return name_simple_phosphane(mol)

    # A boronic acid (P-68.1.4.1's R-B(OH)2) has two hydroxyl oxygens on
    # boron that `_borane.py`'s own plain-borane shape doesn't expect --
    # `has_simple_borane_shape` matches any molecule with a boron atom at
    # all, so this must be routed first, before it misfires on the two
    # -OH oxygens as unsupported heteroatoms.
    if has_boronic_acid_shape(mol):
        return name_boronic_acid(mol)

    # A borinic acid (P-68.1.4.1's R2-B-OH) has the same borane-shape
    # collision as boronic acid above, so it must be routed here for the
    # same reason.
    if has_borinic_acid_shape(mol):
        return name_borinic_acid(mol)

    # A boron atom (P-68's borane substitutive nomenclature, the same shape
    # as phosphane above with boron in place of phosphorus) must be routed
    # here for the same reason -- none of the branches below recognize
    # boron at all.
    if has_simple_borane_shape(mol):
        return name_simple_borane(mol)

    # A Group 13 metal (Al/Ga/In/Tl, P-69.1) is the same substitutive-
    # naming shape as boron/phosphorus above, generalized as one shared
    # mechanism -- must be routed here for the same reason: none of the
    # branches below recognize any of these elements at all.
    if has_group13_hydride_shape(mol):
        return name_group13_hydride(mol)
    if has_group14_hydride_shape(mol):
        return name_group14_hydride(mol)
    if has_group15_hydride_shape(mol):
        return name_group15_hydride(mol)

    # A Group 1/2 metal (Li/Na/K/Mg/Ca, P-69.3) is a different additive
    # naming mechanism from Group 13 above, but the same reasoning for
    # dispatch order applies: none of the branches below recognize any of
    # these elements at all.
    if has_group1_2_organometallic_shape(mol):
        return name_group1_2_organometallic(mol)

    # buckminsterfullerene (P-27's '[60]fullerene', a fixed 12-pentagon/
    # 20-hexagon cage) is recognized by exact whole-molecule match --
    # see _fullerene.py's module docstring; none of the ring modules
    # below understand a cage shape at all.
    if has_fullerene_name(mol):
        return name_fullerene(mol)

    # pyrene/acenaphthylene/fluoranthene/aceanthrylene/acephenanthrylene
    # (P-25.1.2's peri-fused retained names) are recognized by exact
    # whole-molecule match, independent of every other branch below -- see
    # _peri_fused_aromatic.py's module docstring for why (each has a
    # non-6-membered, non-aromatic ring that none of the other dispatch
    # branches expect).
    if has_retained_peri_fused_name(mol):
        return name_retained_peri_fused(mol)

    # 9H-fluorene (P-25.1.2) is recognized the same way -- its central
    # ring's sp3 CH2 fails the aromatic-fused dispatch's precondition too.
    if has_fluorene_parent_name(mol):
        return name_fluorene_parent(mol)

    # cyclopenta[a/b]naphthalene (P-25.3.1.3) needs the same early dispatch
    # for the same reason -- its 5-ring is never fully aromatic either.
    if has_cyclopenta_naphthalene_name(mol):
        return name_cyclopenta_naphthalene(mol)

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

    # The seven 1989 IUPAC steroid parent ring hydrides (gonane through
    # ergostane, Rule 2.1/3S-2.2/2.3/2.4 -- see module docstring) are
    # recognized the same way, independent of every other branch below:
    # `_polycyclic.py`'s general von Baeyer engine already names the bare
    # gonane skeleton (confirmed by direct testing), so this check must
    # come first or gonane would never be reached.
    if has_steroid_parent_hydride_name(mol):
        return name_steroid_parent_hydride(mol)

    # A steroid parent hydride with exactly one ring C=C double bond at a
    # standard, non-ring-fusion locant (e.g. 'androst-5-ene') must be
    # checked right alongside the bare-skeleton case above, for the same
    # von-Baeyer-engine-would-otherwise-claim-it reason.
    if has_steroid_unsaturated_name(mol):
        return name_steroid_unsaturated(mol)

    # A steroid parent hydride whose A-ring is aromatic (the mancude
    # 1,3,5(10)-triene, e.g. 'estra-1,3,5(10)-triene') is checked right
    # alongside the single-double-bond case above, for the same reason.
    if has_steroid_aromatic_a_ring_name(mol):
        return name_steroid_aromatic_a_ring(mol)

    # A steroid parent hydride plus one O/S/N one-atom bridge across an
    # already-adjacent ring bond (e.g. '5,6-epoxycholestane') is checked
    # right alongside the other steroid-skeleton-plus-one-modification
    # cases above, for the same reason.
    if has_bridged_steroid_name(mol):
        return name_bridged_steroid_parent(mol)

    # The morphinan retained parent hydride (Appendix 3 / P-101), with
    # morphine/codeine's exact substituent shape (N-methyl, one or two
    # O-substituents, a transannular epoxy bridge, one extra ring
    # double bond) is checked right alongside the other skeleton-dict
    # cases above, for the same reason.
    if has_alkaloid_morphinan_name(mol):
        return name_alkaloid_morphinan(mol)

    # A steroid parent hydride missing one non-fusion ring atom or angular
    # methyl (P-101.3.1's 'nor' prefix) is checked right after the exact
    # parent match above, so a structure that happens to reproduce a
    # different retained-name parent (e.g. androstane minus its C19 methyl
    # is exactly estrane) is already claimed by that check first.
    if has_nor_steroid_shape(mol):
        return name_nor_steroid(mol)

    # A steroid parent hydride missing two non-fusion ring atoms/angular
    # methyls together (P-101.3.1.1's 'dinor' prefix) is checked right
    # after the single-atom 'nor' check above, for the same dispatch-
    # ordering reason.
    if has_dinor_steroid_shape(mol):
        return name_dinor_steroid(mol)

    # A steroid parent hydride with one extra methylene inserted into a
    # ring bond or angular methyl (P-101.3.2's 'homo' prefix) is checked
    # right after 'nor' for the same reason -- an insertion happening to
    # reproduce a different retained-name parent would already be claimed
    # above (none found empirically, but the ordering stays defensive).
    if has_homo_steroid_shape(mol):
        return name_homo_steroid(mol)

    # A steroid parent hydride missing one ring bond (P-101.3.4.1's
    # 'seco' prefix) is checked right after 'homo' for the same
    # dispatch-ordering reason -- checked before the general von Baeyer
    # engine below, which either misnames a plain-ring-bond cleavage as
    # an unrelated tricyclic system or outright rejects a ring-fusion
    # cleavage as an unsupported polycyclic shape.
    if has_seco_steroid_shape(mol):
        return name_seco_steroid(mol)

    # A steroid parent hydride with one extra ring bond formed between
    # two non-adjacent ring atoms (P-101.3.3's 'cyclo' prefix) is checked
    # right after 'seco' for the same dispatch-ordering reason.
    if has_cyclo_steroid_shape(mol):
        return name_cyclo_steroid(mol)

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

    # phenanthroline/naphthyridine (P-2 Table 2.8's locanted diaza retained
    # names) must be routed here too, before the fusion-letter dispatches
    # below and _aromatic.py's own all-carbon dispatch, which would reject
    # the nitrogen atoms outright.
    phenanthroline_naphthyridine_core = find_phenanthroline_naphthyridine_core(mol)
    if phenanthroline_naphthyridine_core is not None:
        return name_phenanthroline_naphthyridine(mol, phenanthroline_naphthyridine_core)

    # benzo[g]indole/benzo[e][1]benzofuran/benzo[g][1]benzofuran (P-25.3.1.3's
    # computed fusion-locant-letter mechanism, this time for a plain benzo
    # ring fused onto an already-bicyclic retained-name base component)
    # must be routed here before `_aromatic.py`'s own tricyclic dispatch
    # further below, which doesn't recognize a heteroatom at all.
    if has_polycyclic_component_fusion_name(mol):
        return name_polycyclic_component_fusion(mol)

    # benzo[g]quinoline/benzo[h]isoquinoline (P-25.3.1.3's same computed
    # fusion-locant-letter mechanism, this time for quinoline/isoquinoline
    # as the base) must be routed here for the same reason as the check
    # just above.
    if has_pyridine_bicyclic_fusion_name(mol):
        return name_pyridine_bicyclic_fusion(mol)

    # benzo[a]anthracene (P-25.3.1.3's computed fusion-locant-letter
    # mechanism again, this time for a plain benzo ring fused onto
    # anthracene itself -- the Blue Book's own worked example for this
    # section) must be routed here for the same reason as the check just
    # above: before `_aromatic.py`'s own tetracyclic dispatch further
    # below, which only recognizes the *linear* fusion (tetracene) and
    # would otherwise reject the angular one outright. The linear shape
    # is deliberately excluded from `has_anthracene_fusion_name` itself
    # so it still falls through to that tetracene recognition unchanged.
    if has_anthracene_fusion_name(mol):
        return name_anthracene_fusion(mol)

    # benzo[a]tetracene (same mechanism again, for tetracene as the base
    # component -- the second plain catacondensed-chain base, after
    # anthracene) must be routed here for the same reason as the check
    # just above: before `_aromatic.py`'s own pentacyclic dispatch, which
    # only recognizes the *linear* fusion (pentacene) and would otherwise
    # reject the angular one outright. The linear shape is deliberately
    # excluded from `has_tetracene_fusion_name` itself so it still falls
    # through to that pentacene recognition unchanged.
    if has_tetracene_fusion_name(mol):
        return name_tetracene_fusion(mol)

    # benzo[a]pentacene (same mechanism again, for pentacene as the base
    # component -- the third plain catacondensed-chain base, after
    # anthracene/tetracene) must be routed here for the same reason as
    # the checks just above. The linear ('hexacene') shape is deliberately
    # excluded from `has_pentacene_fusion_name` itself so it still falls
    # through to `_aromatic.py`'s own hexacene recognition unchanged.
    if has_pentacene_fusion_name(mol):
        return name_pentacene_fusion(mol)

    # benzo[a]hexacene (same mechanism again, for hexacene as the base
    # component -- the fourth plain catacondensed-chain base, after
    # anthracene/tetracene/pentacene) must be routed here for the same
    # reason as the checks just above. The linear ('heptacene') shape is
    # deliberately excluded from `has_hexacene_fusion_name` itself so it
    # still falls through to `_aromatic.py`'s own heptacene recognition
    # unchanged.
    if has_hexacene_fusion_name(mol):
        return name_hexacene_fusion(mol)

    # benzo[a]heptacene (same mechanism again, for heptacene as the base
    # component -- the fifth plain catacondensed-chain base, after
    # anthracene/tetracene/pentacene/hexacene) must be routed here for
    # the same reason as the checks just above. The linear ('octacene')
    # shape is deliberately excluded from `has_heptacene_fusion_name`
    # itself so it still falls through to `_aromatic.py`'s own octacene
    # recognition unchanged.
    if has_heptacene_fusion_name(mol):
        return name_heptacene_fusion(mol)

    # chrysene/benzo[c]phenanthrene (same mechanism once more, for
    # phenanthrene as the base component) -- routed here for the same
    # reason as the two checks just above. The third structurally
    # possible letter on phenanthrene is deliberately excluded from
    # `has_phenanthrene_fusion_name` itself: 'b' is benzo[a]anthracene,
    # the exact same compound the anthracene check just above already
    # names via a different (senior) base component.
    if has_phenanthrene_fusion_name(mol):
        return name_phenanthrene_fusion(mol)

    # benzo[a]pyrene/benzo[e]pyrene (same mechanism again, for pyrene as
    # the base component -- the first *peri*-fused, not simply
    # catacondensed, base for this algorithm) -- routed here for the same
    # reason as the checks above.
    if has_pyrene_fusion_name(mol):
        return name_pyrene_fusion(mol)

    # benzo[b]chrysene/benzo[c]chrysene/benzo[g]chrysene/picene (same
    # mechanism again, for chrysene as the base component -- chrysene
    # itself is `_phenanthrene_fusion.py`'s letter 'a'; picene is this
    # module's own letter 'a') -- routed here for the same reason as the
    # checks above.
    if has_chrysene_fusion_name(mol):
        return name_chrysene_fusion(mol)

    # benzo[b]picene/benzo[c]picene (same mechanism again, for picene as
    # the base component -- picene itself is `_chrysene_fusion.py`'s
    # letter 'a') -- routed here for the same reason as the checks above.
    if has_picene_fusion_name(mol):
        return name_picene_fusion(mol)

    # benzo[a]pentaphene/hexaphene/benzo[c]pentaphene (same mechanism
    # again, for pentaphene as the base component) -- routed here for the
    # same reason as the checks above.
    if has_pentaphene_fusion_name(mol):
        return name_pentaphene_fusion(mol)

    # benzo[b]triphenylene (same mechanism again, for triphenylene as the
    # base component) -- routed here for the same reason as the checks
    # above. The other structurally possible letter is the exact same
    # compound as `_chrysene_fusion.py`'s letter 'g' and is excluded from
    # `has_triphenylene_fusion_name` itself to avoid a second route to it.
    if has_triphenylene_fusion_name(mol):
        return name_triphenylene_fusion(mol)

    # benzo[a]fluoranthene/benzo[b]fluoranthene/benzo[j]fluoranthene/
    # benzo[k]fluoranthene (same mechanism again, for fluoranthene as the
    # base component -- the second *peri*-fused base, after pyrene)
    # -- routed here for the same reason as the checks above.
    if has_fluoranthene_fusion_name(mol):
        return name_fluoranthene_fusion(mol)

    # benzo[d]aceanthrylene/benzo[e]aceanthrylene/benzo[j]aceanthrylene/
    # benzo[k]aceanthrylene/benzo[l]aceanthrylene (same mechanism again,
    # for aceanthrylene as the base component -- the third peri-fused
    # base) -- routed here for the same reason as the checks above. The
    # sixth structurally possible letter is the exact same compound as
    # `_fluoranthene_fusion.py`'s own letter 'a'/'f' and is excluded from
    # `has_aceanthrylene_fusion_name` itself to avoid a second route to
    # it.
    if has_aceanthrylene_fusion_name(mol):
        return name_aceanthrylene_fusion(mol)

    # benzo[a]acephenanthrylene/benzo[j]acephenanthrylene/
    # benzo[k]acephenanthrylene/benzo[l]acephenanthrylene (same mechanism
    # again, for acephenanthrylene as the base component -- the fourth
    # peri-fused base) -- routed here for the same reason as the checks
    # above. The other two structurally possible letters are the exact
    # same compounds as `_aceanthrylene_fusion.py`'s own letter 'e' and
    # `_fluoranthene_fusion.py`'s own letter 'b'/'e' respectively, and are
    # excluded from `has_acephenanthrylene_fusion_name` itself to avoid a
    # second route to them.
    if has_acephenanthrylene_fusion_name(mol):
        return name_acephenanthrylene_fusion(mol)

    # thieno[2,3-b]thiophene/furo[2,3-b]furan/thieno[2,3-b]furan etc.
    # (P-25.3.1.3's computed fusion-locant-letter mechanism, for two
    # five-membered heteromonocycles -- identical or a mixed O/S pair --
    # self-fused) must be routed here before `_hetero_monocyclic.py` below,
    # which only understands a single ring.
    if has_two_component_heterocycle_fusion_name(mol):
        return name_two_component_heterocycle_fusion(mol)

    # furo[3,2-b]pyridine/thieno[2,3-b]pyridine etc. (P-25.3.1.3's computed
    # fusion-locant-letter mechanism, pyridine as the base component with
    # a named five-membered O/S heteromonocycle attached) -- also routed
    # here before `_hetero_monocyclic.py` below, same reason as above.
    if has_pyridine_heterocycle_fusion_name(mol):
        return name_pyridine_heterocycle_fusion(mol)

    # imidazo[1,2-a]pyridine/imidazo[2,1-b]thiazole etc. (P-25.3.2.5.1's
    # bridgehead-heteroatom fusion -- the shared fusion atom is itself a
    # nitrogen common to both components) -- also routed here before
    # `_hetero_monocyclic.py` below, same reason as above.
    if has_bridgehead_heteroatom_fusion_name(mol):
        return name_bridgehead_heteroatom_fusion(mol)

    # benzo[1,2-b:4,5-b']dithiophene etc. (P-25.3.4.1.3's multiparent
    # fusion -- one benzo ring bridging two identical five-membered
    # heteromonocycles) -- also routed here before `_hetero_monocyclic.py`
    # below, same reason as above.
    if has_benzo_bis_heterocycle_fusion_name(mol):
        return name_benzo_bis_heterocycle_fusion(mol)

    # dithieno[2,3-b:3',2'-e]pyridine etc. (P-25.3.6.1's "identical
    # attached components" -- pyridine as the sole senior parent, with
    # two identical five-membered heteromonocycles as first-order
    # attached components) -- also routed here before
    # `_hetero_monocyclic.py` below, same reason as above.
    if has_pyridine_bis_heterocycle_fusion_name(mol):
        return name_pyridine_bis_heterocycle_fusion(mol)

    # 2,3-didehydrooxepane etc. (P-31.2.2/P-31.2.4.1's 'didehydro' prefix,
    # adding one ring double bond to a saturated Hantzsch-Widman/retained
    # parent) -- routed here before the exact-match check below, since a
    # didehydro ring's extra double bond means it never matches that
    # check's fully-saturated canonical SMILES anyway, but grouped here
    # for the shared `saturated_ring_name` dependency.
    if has_didehydro_ring_name(mol):
        return name_didehydro_ring(mol)

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

    # 2H-pyran/4H-pyran (P-25.7.1.3.1's indicated-hydrogen case) -- unlike
    # furan/thiophene above, RDKit doesn't treat this ring as aromatic at
    # all, so it needs its own recognition shape rather than an extension
    # of `_ROLE_SEQUENCES`; must be routed here before `_ether.py` below,
    # which otherwise rejects any ring outright.
    if has_pyran_indicated_hydrogen_name(mol):
        return name_pyran_indicated_hydrogen(mol)

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

    # A fused aromatic parent (naphthalene, phenanthrene, tetracene, ...)
    # with a single -CH2- or -O- bridge across one ring's 1,4-positions
    # (1,4-dihydro-1,4-methano-/epoxy-<parent>) must be routed here before
    # the plain "any O atom" branch below (the -O- bridge variant would
    # otherwise be misdetected as a plain ether) -- see
    # _bridged_aromatic.py's module docstring for why RDKit's own ring
    # perception can't be trusted for this shape either.
    bridged_core = find_bridged_aromatic_core(mol)
    if bridged_core is not None:
        return name_bridged_aromatic(mol, bridged_core)

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
    for ring_count in (3, 4, 5, 6):
        polycyclic_hetero_core = find_polycyclic_core(mol, ring_count)
        if polycyclic_hetero_core is None:
            continue
        if has_single_ring_heteroatom_shape_polycyclic(mol, polycyclic_hetero_core):
            return name_von_baeyer_heteroatom_polycyclic(mol, polycyclic_hetero_core, ring_count)
        if has_multi_ring_heteroatom_shape_polycyclic(mol, polycyclic_hetero_core):
            return name_von_baeyer_heteroatom_multi_polycyclic(mol, polycyclic_hetero_core, ring_count)
        if has_mixed_element_heteroatom_shape_polycyclic(mol, polycyclic_hetero_core):
            return name_von_baeyer_heteroatom_mixed_polycyclic(mol, polycyclic_hetero_core, ring_count)

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
    # A sulfonyl group on a plain saturated ring nitrogen (e.g.
    # 1-methylsulfonylpiperidine) looks sulfonamide-shaped to
    # `has_sulfonamide_shape` below, which doesn't know about this
    # ring-as-parent construction and would misclaim/reject it -- must be
    # routed here first (narrower than `has_ring_amine_shape` alone, so it
    # doesn't also preempt `_hidden_amide_ketone.py`'s unrelated
    # acyl-on-ring-nitrogen shape further down).
    if has_ring_amine_sulfonyl_shape(mol):
        return name_ring_amine(mol)
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
    # A nitrone (imine N-oxide, P-74.2.1.2) has its own N+/O- dipole
    # pair the plain imine/oxime checks below don't expect, and it's
    # C=N-bonded (like an ordinary imine) so it would otherwise be
    # swallowed by the oxime-gated `has_simple_imine_shape` branch further
    # down and misrouted into `_imine.py`'s own rejection there -- must be
    # routed before it.
    if has_nitrone_shape(mol):
        return name_nitrone(mol)
    # A nitrile oxide (P-74.2.2.2.1.2) has the same dipole-pair issue as
    # nitrone above, checked here for the same reason (before it would
    # otherwise fall through to a general heteroatom-allowlist rejection
    # further down, none of which know about this shape).
    if has_nitrile_oxide_shape(mol):
        return name_nitrile_oxide(mol)
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
        # An ether coexisting with a separate primary amine (P-41: ether
        # has no suffix at all, class 41, so it's always the 'alkoxy'
        # prefix, never competing for parent-hood) must be routed before
        # `has_ether_shape` below, which would otherwise misname/reject it
        # via `_ether.py`'s own "coexisting nitrogen" rejection.
        if has_ether_amine_shape(mol):
            return name_ether_amine(mol)
        # An ether coexisting with a separate thiol (same P-41 reasoning
        # as the amine case immediately above) must likewise be routed
        # before `has_ether_shape` below.
        if has_ether_thiol_shape(mol):
            return name_ether_thiol(mol)
        # An ether coexisting with a separate ketone (same P-41 reasoning
        # as the amine/thiol cases above) must likewise be routed before
        # `has_ether_shape` below.
        if has_ether_ketone_shape(mol):
            return name_ether_ketone(mol)
        # An ether coexisting with a separate aldehyde (same P-41
        # reasoning as the amine/thiol/ketone cases above) must likewise
        # be routed before `has_ether_shape` below.
        if has_ether_aldehyde_shape(mol):
            return name_ether_aldehyde(mol)
        # An ether coexisting with a separate unsubstituted primary amide
        # (same P-41 reasoning as the amine/thiol/ketone/aldehyde cases
        # above) must likewise be routed before `has_ether_shape` below.
        if has_ether_amide_shape(mol):
            return name_ether_amide(mol)
        # An ether coexisting with a separate hydroperoxide (same P-41
        # reasoning as the amine/thiol/ketone/aldehyde/amide cases above)
        # must likewise be routed before `has_ether_shape` below.
        if has_ether_hydroperoxide_shape(mol):
            return name_ether_hydroperoxide(mol)
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
            # P-41/Table 4.1: 'hydroperoxide' (class 18) also outranks
            # 'amine' (class 19), so a hydroperoxide that also carries a
            # separate primary amine names the hydroperoxide as the
            # suffix and demotes the amine to an 'amino' prefix instead of
            # `_hydroperoxide.py`'s own "coexisting nitrogen" rejection.
            # (Note: this follows Table 4.1's text directly, diverging
            # from PubChem's own auto-namer for this specific pair -- see
            # `_hydroperoxide_amine.py`'s module docstring.)
            if has_hydroperoxide_amine_shape(mol):
                return name_hydroperoxide_amine(mol)
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
        # A sulfonate anion's sulfur (R-SO3-, P-72.2.2.2.1.1) bears an
        # anionic oxygen (formal charge -1, no H) that every other module
        # here -- including `_sulfonic_acid.py`'s own neutral -SO3H check
        # further down -- would reject outright, so it must be routed
        # before those.
        if has_sulfonate_shape(mol):
            return name_sulfonate(mol)
        # A D-aldohexopyranose ring (a hemiacetal, no actual ring C=O
        # anywhere) still matches `has_hetero_ring_ketone_shape` below --
        # that check claims any saturated single-heteroatom 5/6/7-membered
        # ring purely by shape, regardless of whether a ketone is actually
        # present, so a plain pyranose would otherwise fall into
        # `_ketone.py`'s hetero-ring-ketone path and have its hydroxyls
        # rejected as unrecognized substituents. Must be routed first.
        if has_cyclic_aldopyranose_shape(mol):
            return name_cyclic_aldopyranose(mol)
        # Same reasoning for a furanose (5-membered) ring -- one of
        # `has_hetero_ring_ketone_shape`'s own accepted ring sizes too.
        if has_cyclic_aldofuranose_shape(mol):
            return name_cyclic_aldofuranose(mol)
        # A D-2-ketohexopyranose ring (e.g. fructopyranose) is a hemiketal,
        # not an actual ring ketone either -- same `has_hetero_ring_ketone_shape`
        # false-claim reasoning as the aldopyranose/aldofuranose checks above.
        if has_cyclic_ketohexopyranose_shape(mol):
            return name_cyclic_ketohexopyranose(mol)
        # Proline's pyrrolidine ring (P-103.1.2) is a plain secondary
        # cyclic amine, not an actual ring ketone either -- same
        # `has_hetero_ring_ketone_shape` false-claim reasoning as the
        # sugar-ring checks above.
        if has_proline_shape(mol):
            return name_proline(mol)
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
            # P-41 Table 4.1: an ether has no suffix at all (class 41), so
            # an ester whose acyl chain also carries a separate ether
            # names the ester as the suffix and demotes the ether to an
            # 'alkoxy' prefix instead of `_ester.py`'s own "coexisting
            # oxygen" rejection.
            if has_ether_ester_shape(mol):
                return name_ether_ester(mol)
            # P-41/Table 4.1: 'ester' (class 9) also outranks 'amine'
            # (class 19), so an ester that also carries a separate primary
            # amine on its acyl chain names the ester as the suffix and
            # demotes the amine to an 'amino' prefix instead of
            # `_ester.py`'s own "coexisting nitrogen" rejection.
            if has_ester_amine_shape(mol):
                return name_ester_amine(mol)
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
            # P-66.1.1.3.3: 'oic acid' also outranks 'amide', so a
            # carboxylic acid that also carries a coexisting primary amide
            # (an amic acid) names the acid as the suffix and demotes the
            # amide to 'amino'+'oxo' prefixes instead of
            # `_carboxylic_acid.py`'s own "coexisting nitrogen" rejection.
            if has_carboxylic_acid_amide_shape(mol):
                return name_carboxylic_acid_amide(mol)
            # P-103.1.1.1: histidine's imidazol-4-ylmethyl side chain also
            # gets its own retained name + L/D descriptor -- must be
            # checked before `has_carboxylic_acid_amine_shape` below,
            # whose own ring dispatch only recognizes a plain benzene ring
            # and would otherwise reject this shape outright.
            if has_histidine_shape(mol):
                return name_histidine(mol)
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
            # P-41/Table 3.3: 'amide' also outranks 'amine', so a primary
            # amide that also carries a separate primary amine names the
            # amide as the suffix and demotes the amine to an 'amino'
            # prefix instead of `_amide.py`'s own "coexisting nitrogen"
            # rejection.
            if has_amide_amine_shape(mol):
                return name_amide_amine(mol)
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
            # P-41/Table 3.3: 'one'/'al' also outrank 'amine', so a ketone
            # or aldehyde coexisting with a primary amine names the
            # carbonyl as the suffix and demotes the amine to an 'amino'
            # prefix instead of `_ketone.py`'s/`_aldehyde.py`'s own
            # "coexisting nitrogen" rejection.
            if has_ketone_amine_shape(mol):
                return name_ketone_amine(mol)
            if has_aldehyde_amine_shape(mol):
                return name_aldehyde_amine(mol)
            # A plain open-chain aldose (P-102.5.1/P-102.5.2.2, #1039 M1
            # step 1) is itself aldehyde-shaped (a terminal -CHO), so it
            # must be routed before the generic aldehyde fallback below,
            # which would otherwise name it as a plain polyhydroxy-
            # aldehyde substitutive name instead of its carbohydrate name.
            if has_open_chain_aldose_shape(mol):
                return name_open_chain_aldose(mol)
            if any(_is_aldehyde_shaped(o) for o in carbonyl_oxygens):
                return name_aldehyde(mol)
            # A plain open-chain 2-ketose (P-102.5.2.1/P-102.5.2.2, #1039
            # M1 step 2) is itself ketone-shaped (a non-terminal
            # carbonyl), so it must be routed before the generic ketone
            # fallback below, which would otherwise name it as a plain
            # polyhydroxy-ketone substitutive name instead of its
            # carbohydrate name.
            if has_open_chain_2_ketose_shape(mol):
                return name_open_chain_2_ketose(mol)
            return name_ketone(mol)
        # P-41/Table 3.3: '-ol' also outranks 'amine', so one or more
        # hydroxyls coexisting with a primary amine names the alcohol as
        # the suffix and demotes the amine to an 'amino' prefix instead of
        # `_alcohol.py`'s own "coexisting nitrogen" rejection.
        if has_alcohol_amine_shape(mol):
            return name_alcohol_amine(mol)
        # A real nitrile/thiol/selenol coexisting with an oxygen that isn't
        # a real alcohol at all (e.g. a furan ring substituent's own
        # non-functional ring oxygen, now recognized by these modules
        # themselves -- P-616 M4) has no coexisting-alcohol shape above to
        # catch it, so it would otherwise fall through to `name_alcohol`'s
        # own unconditional fallback -- must be routed here first, mirroring
        # the nitrogen-gated branch's own thiol/selenol checks before its
        # `name_amine` fallback below.
        if has_nitrile_shape(mol):
            return name_nitrile(mol)
        if has_thiol_shape(mol):
            return name_thiol(mol)
        if has_selenol_shape(mol):
            return name_selenol(mol)
        if has_tellurol_shape(mol):
            return name_tellurol(mol)
        # One of the 9 retained-name inositol stereoisomers (P-104.2.1) is
        # itself a plain cyclohexane-1,2,3,4,5,6-hexol shape, so it must be
        # routed before the generic systematic fallback below, which would
        # otherwise name it with CIP descriptors instead of its retained
        # name -- mirrors the open-chain-aldose/2-ketose carbohydrate-name-
        # before-generic-fallback pattern above.
        if has_inositol_shape(mol):
            return name_inositol(mol)
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
        # The reverse shape: the ring nitrogen itself is plain (no
        # substituent), but a ring carbon bears an exocyclic primary
        # amine (e.g. piperidin-4-amine) -- `name_amine` also explicitly
        # defers any ring bearing a nitrogen heteroatom, so this must be
        # routed here first too, same reasoning as `has_ring_amine_shape`
        # just above.
        if has_hetero_ring_amine_shape(mol):
            return name_hetero_ring_amine(mol)
        # P-41/Table 3.3: '-thiol' (sharing alcohol's rank as a chalcogen
        # analogue) outranks 'amine', so one or more thiols coexisting with
        # a primary amine names the thiol as the suffix and demotes the
        # amine to an 'amino' prefix instead of `name_amine`'s own
        # "coexisting sulfur" rejection.
        if has_thiol_amine_shape(mol):
            return name_thiol_amine(mol)
        # A thiol coexisting with a nitrogen that isn't a real amine at all
        # (e.g. a heteroaromatic ring substituent's own pyridine/pyrrole
        # nitrogen, now recognized by `_thiol.py` itself -- P-25 M2 step 1,
        # #628) has no coexisting-amine shape above to catch it (those all
        # require an actual primary amine), so it would otherwise fall
        # through to `name_amine`'s own unconditional rejection of a
        # nitrogen-free molecule it doesn't recognize as amine-shaped at
        # all -- must be routed here first, mirroring `_ketone.py`'s/
        # `_alcohol.py`'s own unconditional fallback to their suffix
        # module regardless of a coexisting non-functional nitrogen.
        if has_thiol_shape(mol):
            return name_thiol(mol)
        # A selenol coexisting with a nitrogen that isn't a real amine at
        # all (e.g. a heteroaromatic ring substituent's own pyridine/
        # pyrrole nitrogen) has no coexisting-amine shape above to catch
        # it either -- same reasoning as the thiol check just above, one
        # chalcogen row down.
        if has_selenol_shape(mol):
            return name_selenol(mol)
        # A tellurol coexisting with a nitrogen that isn't a real amine at
        # all has no coexisting-amine shape above to catch it either --
        # same reasoning as the thiol/selenol checks just above, the next
        # chalcogen row down.
        if has_tellurol_shape(mol):
            return name_tellurol(mol)
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
        # this branch once both are ruled out above (a thiol coexisting
        # with an amine is instead routed inside the nitrogen-gated branch
        # above, before its own `name_amine` fallback).
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
        # A fully-saturated naphthalene skeleton (decalin) is a mancude
        # ring system's hydro derivative (P-31.2.3.3.2), not a von Baeyer
        # system, even though its carbon skeleton is graph-isomorphic to
        # one -- must be routed here before find_bicyclic_core below for
        # the same reason the partial-hydro case above already is.
        decahydro_core = find_decahydronaphthalene_core(mol)
        if decahydro_core is not None:
            return name_decahydronaphthalene(mol, decahydro_core)
        # A naphthalene skeleton with 1-4 plain (non-aromatic) Kekule ring
        # double bonds -- anywhere between the dihydro case above and full
        # saturation -- is likewise a mancude ring system's hydro
        # derivative (P-31.2.3.3.2), not a von Baeyer system (#828, M2
        # step 1); must be routed here for the same reason as the two
        # cases above.
        partial_core = find_partially_unsaturated_naphthalene_core(mol)
        if partial_core is not None:
            return name_partially_unsaturated_naphthalene(mol, partial_core)
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
        # A plain exocyclic double bond (e.g. '=CH2'/'=CHR') on an
        # otherwise saturated, unsubstituted ring -- e.g.
        # methylidenecyclohexane -- is a different shape from both the
        # ring-internal-unsaturation case above and `name_cycloalkane`'s
        # own plain saturated ring below, so it must be routed here first;
        # `name_cycloalkane` itself now rejects any exocyclic non-single
        # bond that reaches it unclaimed (see its own docstring comment).
        exocyclic_ylidene_core = find_exocyclic_ylidene_core(mol)
        if exocyclic_ylidene_core is not None:
            return name_exocyclic_ylidene(mol, exocyclic_ylidene_core)
        return name_cycloalkane(mol)

    # Two disjoint (unfused) identical rings or ring systems joined by a C=C
    # double bond -- e.g. bi(cyclopentylidene), bi(bicyclo[2.2.1]heptan-
    # ylidene) -- must be routed here before find_bicyclic_core/
    # find_polycyclic_core below, for the same reason the aromatic
    # num_rings == 2 ring-assembly case above is. Not gated on a specific
    # num_rings value, same reasoning as the bicyclic detection below it:
    # a von Baeyer bicyclic side's own SSSR ring count can overcount for
    # symmetric bridging, and a bicyclic-sided assembly has two rings per
    # side to begin with (4 total, not 2) -- find_ring_assembly_ylidene_core
    # itself does the real, cheap shape check (exactly one non-aromatic C=C
    # bond, degree 3 on both ends).
    ring_assembly_ylidene_core = find_ring_assembly_ylidene_core(mol)
    if ring_assembly_ylidene_core is not None:
        return name_ring_assembly_ylidene(mol, ring_assembly_ylidene_core)

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
    # Two disjoint plain rings joined only through an acyclic bridge (e.g.
    # dicyclohexylmethane) have no shared atom and no ring-to-ring bond, so
    # none of the fused/spiro/bridged/assembly checks above ever claim
    # them -- must be routed last, right before the catch-all rejection.
    disjoint_ring_pair_core = find_disjoint_ring_pair_core(mol)
    if disjoint_ring_pair_core is not None:
        return name_disjoint_ring_pair(mol, disjoint_ring_pair_core)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
