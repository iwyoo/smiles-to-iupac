import re

from rdkit import Chem

from ._zwitterion import has_zwitterion_shape, name_zwitterion
from ._adduct import has_adduct_shape, name_adduct
from ._acyclic import name_acyclic_alkane
from ._acid_derivatives import name_acid_derivative, name_phosphorous_acid
from ._acid_salts import name_acid_salt
from ._hetero_carboxylic import name_hetero_parent_acid
from ._polycarbonic import name_polycarbonic
from ._carbonic_family import name_carbonic_family
from ._acyl_halide import has_acyl_halide_shape, name_acyl_halide
from ._anhydride import has_anhydride_shape, name_anhydride
from ._carbamate import has_carbamate_shape, name_carbamate
from ._alcohol import name_alcohol
from ._alcohol_amine import has_alcohol_amine_shape, name_alcohol_amine
from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._acetyl_names import acetyl_names
from ._anion import name_anion
from ._aldehyde import name_aldehyde
from ._aldehyde_amine import has_aldehyde_amine_shape, name_aldehyde_amine
from ._ketone_amine import has_ketone_amine_shape, name_ketone_amine
from ._amino_acid_derivative import has_amino_acid_shape, name_amino_acid
from ._peptide import has_peptide_shape, name_peptide
from ._mixed_onium import has_mixed_onium_shape, name_mixed_onium
from ._axial_stereo import cite_axial_stereo
from ._chalcone import has_chalcone_shape, name_chalcone
from ._hydrogen_cation import hydrogen_salt_name
from ._silicic_cyanate import silicic_cyanate_name
from ._borane_silane_amide import borane_silane_amide_name
from ._polyborane import lewis_adduct_mol, polyborane_name
from ._diacylamine import diacylamine_name
from ._carbonic_hydrazide import carbonic_hydrazide_name, name_carbonic_hydrazide
from ._chalcogen_hydrazide import name_chalcogen_hydrazide
from ._diacylhydrazine import diacylhydrazine_name
from ._chain_amide import chain_amide_name, chain_ketone_name
from ._ring_nitrogen_hydrazide import name_ring_nitrogen_hydrazide
from ._alternating_cage import has_alternating_cage_shape, name_alternating_cage
from ._dipolar import has_dipolar_shape, name_dipolar
from ._chalcogen_aldehyde import name_chalcogen_aldehyde
from ._condensed_guanidine import name_condensed_guanidine
from ._ring_heteroatom_nitrile import name_ring_heteroatom_nitrile
from ._chain_onium import has_chain_onium_shape, name_chain_onium
from ._hetero_acylium import has_hetero_acylium_shape, name_hetero_acylium
from ._chain_ylium import has_chain_ylium_shape, name_chain_ylium
from ._group_polycation import has_group_polycation_shape, name_group_polycation
from ._hydride_ylium import (
    has_hydride_onium_shape,
    has_hydride_ylium_shape,
    has_poly_ylium_shape,
    name_hydride_onium,
    name_hydride_ylium,
    name_poly_ylium,
)
from ._spiro_hub_atom import has_spiro_hub_atom_shape, name_spiro_hub_atom
from ._substituents import FORCED_BRANCH_NAMES
from ._glycoside import has_glycoside_shape, name_glycoside
from ._sugar_acid import has_sugar_alcohol_acid_shape, has_sugar_lactone_shape, name_sugar_lactone, name_sugar_alcohol_acid, sugar_acid_derivative_name
from ._sugar_substituted import has_substituted_sugar_shape, name_substituted_sugar
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
    has_cyclic_ketohexofuranose_shape,
    has_cyclic_ketohexopyranose_shape,
    has_open_chain_2_ketose_shape,
    has_open_chain_aldose_shape,
    name_cyclic_aldofuranose,
    name_cyclic_aldopyranose,
    name_cyclic_ketohexofuranose,
    name_cyclic_ketohexopyranose,
    name_open_chain_2_ketose,
    name_open_chain_aldose,
)
from ._inositol import has_inositol_shape, name_inositol
from ._inositol_derivative import has_inositol_derivative_shape, name_inositol_derivative
from ._sphingoid import has_sphingoid_shape, name_sphingoid
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
from ._ammonium import has_ammonium_shape, has_polyammonium_shape, name_ammonium, name_polyammonium
from ._uronium import has_uronium_shape, name_uronium
from ._polycation import has_polycation_shape, has_ring_nitrenium_shape, name_polycation
from ._polyspiro_union import has_spiro_union_shape, name_spiro_union
from ._ylium_ring import has_ylium_ring_shape, name_ylium_ring
from ._anisole import has_anisole_shape, name_anisole
from ._polynuclear_oxoacid import name_polynuclear_oxoacid
from ._common_hydride import has_common_hydride_shape, name_common_hydride
from ._chain_cation import has_chain_cation_shape, name_chain_cation
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
from ._appendix3_skeletons import name_appendix3_skeleton
from ._borane import has_simple_borane_shape, name_simple_borane
from ._boronic_acid import has_boronic_acid_shape, name_boronic_acid
from ._borinic_acid import has_borinic_acid_shape, name_borinic_acid
from ._metal_pair import has_metal_pair_shape, name_metal_pair
from ._coordination import has_coordination_shape, name_coordination
from ._group1_2_organometallic import has_group1_2_organometallic_shape, name_group1_2_organometallic
from ._nonstandard_hydride import has_nonstandard_hydride_shape, name_nonstandard_hydride
from ._group13_hydride import (
    has_group13_hydride_shape,
    has_group14_hydride_shape,
    has_group15_hydride_shape,
    name_group13_hydride,
    name_group14_hydride,
    name_group15_hydride,
)
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
from ._common import CITE_SKELETAL_LAMBDA, UnsupportedStructure, alphanumerical_name_key, non_single_bonds
from ._cyclic import name_cycloalkane
from ._disjoint_ring_substituents import find_disjoint_ring_pair_core, name_disjoint_ring_pair
from ._cyclic_unsaturated import find_cyclic_unsaturated_core, name_cyclic_unsaturated
from ._diester_acyloxy import has_diester_shape, has_polyester_of_one_polyol_shape, name_diester_acyloxy
from ._hydride_polyester import has_hydride_polyester_shape, name_hydride_polyester
from ._ester import has_ester_shape, name_ester
from ._ester_by_parts import name_ester_by_parts
from ._heteroacyclic import name_heteroacyclic, has_carbonless_acyl
from ._nitrogen_methylene_multiplicative import name_nitrogen_methylene_multiplicative
from ._chain_multiplicative import has_chain_multiplicative_shape
from ._np import PREFERRED_OPERATIONS, name_natural_product_ranked
from ._steroid_named import name_steroid
from ._polyfunctional import name_polyfunctional
from ._cyanate import has_cyanate_shape, name_cyanate
from ._ether import has_ether_shape, name_ether
from ._ether_amine import has_ether_amine_shape, name_ether_amine
from ._ether_aldehyde import has_ether_aldehyde_shape, name_ether_aldehyde
from ._ether_amide import has_ether_amide_shape, name_ether_amide
from ._ether_hydroperoxide import has_ether_hydroperoxide_shape, name_ether_hydroperoxide
from ._ether_ketone import has_ether_ketone_shape, name_ether_ketone
from ._ether_thiol import has_ether_thiol_shape, name_ether_thiol
from ._fusion_name import FUSION_NAME_REQUIRED, PREFER_VON_BAEYER, fused_ring_system_name
from ._hetero_prefixes import _has_senior_principal_group
from ._fullerene import (
    has_fullerene_name,
    has_substituted_fullerene_cage,
    name_fullerene,
    require_defined_fullerene_numbering,
)
from ._fullerene_numbering import name_cage_parent
from ._multiplicative import name_if_multiplicative
from ._nucleoside import has_nucleoside_name, name_nucleoside
from ._nucleoside_substituted import has_substituted_nucleoside_name, name_substituted_nucleoside
from ._oligonucleotide import oligonucleotide_name
from ._metallacycle import has_metallacycle_shape, name_metallacycle
from ._metallacycle_group import name_metallacycle_as_group
from ._metallafused import has_metallafused_shape, name_metallafused
from ._metallapolycycle import has_metallapolycycle_shape, name_metallapolycycle
from ._ocene import has_ocene_shape, name_ocene
from ._pin import enter, leave, mark, nested, outermost, reason_count, reasons_since, replay
from ._hydride_carbo_suffix import has_hydride_carbo_suffix_shape, name_hydride_carbo_suffix
from ._lambda_ring import has_lambda_ring_shape, name_lambda_ring
from ._ring_lambda_heterone import has_ring_lambda_heterone_shape, name_ring_lambda_heterone
from ._hetero_ring_oxide import has_hetero_ring_oxide_shape, name_hetero_ring_oxide
from ._pyridinone import has_pyridinone_shape, name_pyridinone
from ._pyrimidinedione import has_pyrimidinedione_shape, name_pyrimidinedione
from ._pyrimidinone import has_pyrimidinone_shape, name_pyrimidinone
from ._steroid_parent_hydrides import (
    has_steroid_aromatic_a_ring_name,
    has_steroid_parent_hydride_name,
    has_steroid_unsaturated_name,
    name_steroid_aromatic_a_ring,
    name_steroid_parent_hydride,
    name_steroid_unsaturated,
)
from ._hetero_monocyclic import (
    has_hetero_monocyclic_name,
    has_hetero_monocyclic_substituent_name,
    has_pyran_indicated_hydrogen_name,
    name_hetero_monocyclic,
    name_hetero_monocyclic_substituent,
    name_pyran_indicated_hydrogen,
)
from ._carbene_amine import has_carbene_amine_shape, name_carbene_amine
from ._carbon_monoxide import has_carbon_monoxide_shape, name_carbon_monoxide
from ._heteroaryne import has_heteroaryne_shape, name_heteroaryne
from ._didehydro_ring import has_didehydro_ring_name, name_didehydro_ring
from ._chalcogen_chain_heterone import name_chalcogen_chain_heterone
from ._halogen_acid_ester import name_halogen_acid_ester
from ._halogen_amide import name_halogen_amide
from ._inorganic_acid_derivative import has_inorganic_acid_derivative_shape, name_inorganic_acid_derivative
from ._halogen_oxo import name_halogen_oxo
from ._hydroxylamine import has_hydroxylamine_shape, name_hydroxylamine
from ._formazan import has_formazan_shape, name_formazan
from ._hydroxylamine_acid import has_hydroxylamine_acid_shape, name_hydroxylamine_acid
from ._hydroxylamine_general import has_o_substituted_hydroxylamine_shape, name_o_substituted_hydroxylamine
from ._imine import has_simple_imine_shape, name_imine
from ._dipole_oxide import (
    has_nitrile_oxide_prefix_shape,
    has_nitrile_oxide_shape,
    has_nitrone_shape,
    name_nitrile_oxide_prefix,
    name_nitrile_oxide,
    name_nitrone,
)
from ._isotope import has_isotope_shape, name_isotope
from ._radical_ion_skeleton import has_skeleton_radical_ion_shape, name_skeleton_radical_ion
from ._radical_group import has_group_cation_shape, has_radical_group_shape, name_group_cation, name_radical_group
from ._isotope_alcohol import has_isotope_alcohol_shape, name_isotope_alcohol
from ._isotope_carboxylic_acid import has_isotope_carboxylic_acid_shape, name_isotope_carboxylic_acid
from ._isotope_ketone import has_isotope_ketone_shape, name_isotope_ketone
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
from ._hydrazine import has_hydrazine_aminooxy_shape, has_hydrazine_shape, name_hydrazine
from ._hydrazine_multiplicative import has_hydrazine_multiplicative_shape, name_hydrazine_multiplicative
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
from ._cyclophane import has_cyclophane_name, name_cyclophane, name_nonpreferred_cyclophane
from ._linear_phane import has_linear_phane_shape, linear_phane_pin, name_linear_phane
from ._phosphane import has_simple_phosphane_shape, name_simple_phosphane
from ._polyphosphane import has_polyphosphane_shape, name_polyphosphane
from ._noncarbon_oxoacid import has_noncarbon_oxoacid_shape, name_noncarbon_oxoacid
from ._functional_replacement_oxoacid import (
    has_functional_replacement_oxoacid_shape,
    name_functional_replacement_oxoacid,
)
from ._ring_imine import has_ring_imine_shape, name_ring_imine
from ._methanediimine import has_methanediimine_shape, name_methanediimine
from ._phosphanone import has_phosphanimine_shape, has_phosphanone_shape, name_phosphanimine, name_phosphanone
from ._mononuclear_oxoacid import has_mononuclear_oxoacid_shape, name_mononuclear_oxoacid
from ._sulfuric_amide import has_sulfuric_amide_shape, name_sulfuric_amide
from ._phosphate import has_phosphate_shape, name_phosphate
from ._phosphorus_thioester import has_phosphorus_thioester_shape, name_boron_peroxy_ester, name_phosphorus_thioester
from ._dinuclear_oxoacid import has_dinuclear_oxoacid_shape, name_dinuclear_oxoacid
from ._phosphite import has_phosphite_shape, name_phosphite
from ._sulfate import has_sulfate_shape, name_sulfate
from ._sulfite import has_sulfite_shape, name_sulfite
from ._phosphonic_acid import has_phosphonic_acid_shape, name_phosphonic_acid
from ._phosphinic_acid import has_phosphinic_acid_shape, name_phosphinic_acid
from ._phosphorus_acid_derivative import has_phosphorus_acid_derivative_shape, name_phosphorus_acid_derivative
from ._phosphane_chain import has_phosphane_chain_shape, name_phosphane_chain
from ._ring_diyl_numbering import bridged_ring_system_name, is_hydro_fusion_system
from ._hydroperoxide import has_chalcogen_peroxol_shape, has_hydroperoxide_shape, name_hydroperoxide
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
from ._vb_ring_assembly import find_vb_ring_assembly_core, name_vb_ring_assembly
from ._ring_assembly_chain import find_ring_assembly_chain_core, name_ring_assembly_chain
from ._ring_assembly_ylidene import find_ring_assembly_ylidene_core, name_ring_assembly_ylidene
from ._silane_chain import has_silane_chain_shape, name_silane_chain
from ._spiro_tree import has_spiro_tree_shape, name_spiro_tree
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
from ._hetero_chain import contract_hetero_groups_candidates, name_hetero_macrocycle, name_skeletal_chain
from ._phosphanyl_group import contract_phosphanyl_groups
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
    if carbon.GetAtomicNum() != 6 or sum(1 for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6) != 1:
        return False
    return not any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(carbon).GetAtomicNum() == 6 for b in carbon.GetBonds())


_NO_PIN_ADDUCT = "the Blue Book assigns no PIN to Lewis adducts, whose preferred names are coordination names (P-68.1.6.2)"
_NO_PIN_ORGANOMETALLIC ="the Blue Book defines no PIN for this class of organometallic compound (P-69.0)"

_FALLBACKS_RUNNING = set()


# P-22.1.3: toluene and the xylenes are preferred names of the unsubstituted hydrocarbons only
_METHYLBENZENES = {"Cc1ccccc1": "toluene", "Cc1ccccc1C": "1,2-xylene", "Cc1cccc(C)c1": "1,3-xylene", "Cc1ccc(C)cc1": "1,4-xylene"}


_ADAMANTANE = re.compile(r"(?<!bi)(?<!ter)(?<!quater)(?<!yclo)tricyclo\[3\.3\.1\.1\^3,7\]decan(?=e|-)")


_CUBANE = re.compile(r"(?<!bi)(?<!ter)(?<!quater)(?<!yclo)pentacyclo\[4\.2\.0\.0\^2,5\.0\^3,8\.0\^4,7\]octan(?=e|-)")


_INDACENE_PREFIX = re.compile(r"([a-z\]\)])(as-indacen|(?<!a)s-indacen)")


def _retained_polycycle_names(name):
    """P-23.7: 'adamantane' replaces tricyclo[3.3.1.1^3,7]decane and 'cubane' pentacyclo[4.2.0.0^2,5.0^3,8.0^4,7]octane; the
    numbering is the same."""
    return _INDACENE_PREFIX.sub(r"\1-\2", _CUBANE.sub("cuban", _ADAMANTANE.sub("adamantan", name)))


def _is_nonbenzene_monocyclic_annulene(mol):
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) != 1 or len(rings[0]) == 6:
        return False
    return all(
        (a := mol.GetAtomWithIdx(i)).GetIsAromatic()
        and a.GetAtomicNum() in (6, 7, 8, 16, 34, 52)
        and (a.GetAtomicNum() == 6 or len(rings[0]) > 8)
        and not a.GetFormalCharge()
        for i in rings[0]
    )


def _is_aromatic_ring_without_double_bonds(mol):
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) != 1 or not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in rings[0]):
        return False
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return False
    ring = set(rings[0])
    return not any(
        b.GetBondTypeAsDouble() == 2.0 for b in kekule.GetBonds() if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring
    )


def _is_aromatic_ring_with_triple_bond(mol):
    return any(b.GetBondType() == Chem.BondType.TRIPLE and b.IsInRing() and b.GetIsAromatic() for b in mol.GetBonds())


def _neutral_group15_oxides(mol):
    """P+-O-, As+-O- and Sb+-O- written as zwitterions are the doubly bonded oxides of the lambda-convention
    (P-74.2.1.4); RDKit writes some neutral P=O groups this way."""
    pairs = [
        (atom.GetIdx(), n.GetIdx())
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() in (15, 33, 51) and atom.GetFormalCharge() == 1
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() in (8, 16, 34) and n.GetFormalCharge() == -1 and n.GetDegree() == 1
    ]
    if not pairs:
        return mol
    editable = Chem.RWMol(mol)
    seen = set()
    for centre, anion in pairs:
        if centre in seen:
            continue
        seen.add(centre)
        editable.GetAtomWithIdx(centre).SetFormalCharge(0)
        editable.GetAtomWithIdx(anion).SetFormalCharge(0)
        editable.GetBondBetweenAtoms(centre, anion).SetBondType(Chem.BondType.DOUBLE)
    converted = editable.GetMol()
    try:
        Chem.SanitizeMol(converted, Chem.SANITIZE_ALL ^ Chem.SANITIZE_CLEANUP)
    except Exception:
        return mol
    return converted


def _parse_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and (
        _is_nonbenzene_monocyclic_annulene(mol)
        or _is_aromatic_ring_without_double_bonds(mol)
        or _is_aromatic_ring_with_triple_bond(mol)
    ):
        # P-54.2: only benzene is named as an aromatic ring; larger annulenes take ene/yne endings, and
        # RDKit's aromatic perception would drop their E/Z bond stereo. A ring of NH-type atoms with no double
        # bond is saturated although RDKit counts its lone pairs as an aromatic sextet.
        kekule = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.SanitizeMol(kekule, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
        Chem.AssignStereochemistry(kekule, cleanIt=True, force=True)
        return kekule
    if mol is not None:
        for atom in mol.GetAtoms():
            # RDKit knows only the thallium(I) valence, so the standard TlH3 of P-68.1.1.1 looks like a radical
            if atom.GetAtomicNum() == 81 and atom.GetNumRadicalElectrons() and atom.GetTotalValence() == 3:
                atom.SetNumRadicalElectrons(0)
        return _neutral_group15_oxides(mol)
    # hypervalent anionic centers (lambda-convention parents) fail RDKit's valence check only
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None or not any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
        return None
    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    except Exception:
        return None
    mol.SetProp("_hypervalent_anion", "1")
    return mol


def _has_free_anion(mol) -> bool:
    if any(
        a.GetFormalCharge() > 0 and a.GetDegree() and not any(n.GetFormalCharge() < 0 for n in a.GetNeighbors())
        for a in mol.GetAtoms()
    ):
        return False
    return any(
        a.GetFormalCharge() < 0 and not any(n.GetFormalCharge() > 0 for n in a.GetNeighbors()) for a in mol.GetAtoms()
    )


def smiles_to_iupac(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError(f"smiles must be a str, not {type(smiles).__name__}")
    if not smiles.strip():
        raise ValueError(f"invalid SMILES: {smiles!r}")
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(False)
    try:
        if "[H+]" in smiles:
            hydrogen_salt = hydrogen_salt_name(smiles, smiles_to_iupac)
            if hydrogen_salt is not None:
                return hydrogen_salt
        name = cite_axial_stereo(smiles, _retained_polycycle_names(_smiles_to_iupac_unabridged(smiles)))
        mol = _parse_smiles(smiles)
        if mol is not None:
            _require_isotopes_cited(mol, name)
            _require_radicals_cited(mol, name)
        if mol is not None and _has_free_anion(mol):
            name = acetyl_names(name)
        return name
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)


_RADICAL_ENDINGS = ("yl", "ylidene", "ylidyne", "yne", "ylium", "yloxy", "yliumyl")


_HYDRIDE_RADICAL_ELEMENTS = {5, 13, 14, 15, 31, 32, 33, 49, 50, 51, 81, 82, 83}


def _require_radicals_cited(mol, name):
    if has_carbon_monoxide_shape(mol) or has_carbene_amine_shape(mol):
        return
    if any(
        a.GetNumRadicalElectrons() and a.GetAtomicNum() in _HYDRIDE_RADICAL_ELEMENTS and (a.GetDegree() or not a.GetFormalCharge())
        for a in mol.GetAtoms()
    ) and "λ" not in name and not name.rstrip(")]} ").endswith(_RADICAL_ENDINGS):
        raise UnsupportedStructure("the name does not cite the radical centre of a skeletal heteroatom")
    if (
        not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        or any(a.GetFormalCharge() for a in mol.GetAtoms())
        or len(Chem.GetMolFrags(mol)) > 1
        or any(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds())
        or has_nonstandard_hydride_shape(mol)
    ):
        return
    if not name.rstrip(")]} ").endswith(_RADICAL_ENDINGS):
        raise UnsupportedStructure("the name does not cite the radical centre")


def _has_nameable_radical_group(mol) -> bool:
    if not has_radical_group_shape(mol):
        return False
    try:
        name_radical_group(mol)
    except UnsupportedStructure:
        return False
    return True


def _has_isotope_label(mol) -> bool:
    return any(a.GetIsotope() for a in mol.GetAtoms())


def _name_isotope_label(mol) -> str:
    try:
        return name_polyfunctional(mol)
    except UnsupportedStructure:
        pass
    for has_shape, namer in (
        (has_isotope_alcohol_shape, name_isotope_alcohol),
        (has_isotope_ketone_shape, name_isotope_ketone),
        (has_isotope_carboxylic_acid_shape, name_isotope_carboxylic_acid),
        (has_isotope_shape, name_isotope),
    ):
        if has_shape(mol):
            try:
                return namer(mol)
            except UnsupportedStructure as error:
                try:
                    return _name_labelled_substituents(mol)
                except UnsupportedStructure:
                    raise error from None
    return _name_labelled_substituents(mol)


def _name_labelled_substituents(mol):
    """Any parent named without its nuclides whose labelled atoms all sit in substituent groups that cite them."""
    from ._isotope_labels import split_isotopes
    from ._substituents import ISOTOPE_LABELS

    split = split_isotopes(mol)
    if split is None or any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in split[0].GetAtoms()):
        raise UnsupportedStructure("this isotopically modified structure is not supported yet")
    clean, labels, _ = split
    context = {"labels": labels, "consumed": set(), "mol": clean}
    token = ISOTOPE_LABELS.set(context)
    try:
        name = _name_mol(clean)
    finally:
        ISOTOPE_LABELS.reset(token)
    if set(labels) - context["consumed"]:
        raise UnsupportedStructure("a labelled atom of the parent is not cited by this name")
    _require_isotopes_cited(mol, name)
    return name


_MULTIPLIER_VALUE = {"di": 2, "bis": 2, "tri": 3, "tris": 3, "tetra": 4, "tetrakis": 4, "penta": 5, "hexa": 6}
_MULTIPLIED_GROUP = re.compile(r"(di|bis|tri|tris|tetra|tetrakis|penta|hexa)[\[{(]$")
_NUCLIDE_ITEM = re.compile(r"^(\d+)([A-Z][a-z]?)(\d*)$")


def _cited_nuclides(name):
    table = Chem.GetPeriodicTable()
    cited = {}
    for found in re.finditer(r"\(([^()]*)\)", name):
        group = found.group(1)
        multiplier = _MULTIPLIED_GROUP.search(name[: found.start()])
        factor = _MULTIPLIER_VALUE[multiplier.group(1)] if multiplier else 1
        items = group.split(",")
        located, pending = [], []
        for item in items:
            locant, _, tail = item.rpartition("-")
            token = _NUCLIDE_ITEM.match(tail or locant)
            if token is None:
                pending.append(item)
                continue
            mass, symbol, count = int(token.group(1)), token.group(2), token.group(3)
            try:
                number = table.GetAtomicNumber(symbol)
            except Exception:
                pending.append(item)
                continue
            if mass < number:
                pending.append(item)
                continue
            locants = pending + ([locant] if tail else [])
            pending = []
            amount = (int(count) if count else max(len(locants), 1)) * factor
            key = f"{mass}{symbol}"
            cited[key] = cited.get(key, 0) + amount
    return cited


def _require_isotopes_cited(mol, name):
    expected = {}
    for atom in mol.GetAtoms():
        if atom.GetIsotope():
            key = f"{atom.GetIsotope()}{atom.GetSymbol()}"
            expected[key] = expected.get(key, 0) + 1
    if not expected:
        return
    cited = _cited_nuclides(name)
    for key, count in expected.items():
        if key not in name:
            raise UnsupportedStructure(f"the name does not cite the {key} nuclide label")
        if cited.get(key, 0) != count:
            raise UnsupportedStructure(f"the name does not cite every {key} nuclide label")


_NESTED_NAMES: dict = {}
_NESTED_NAMES_MAX = 4096


def _smiles_to_iupac_unabridged(smiles: str) -> str:
    # Acyl/substituent namers re-name the same fragment hundreds of times while ranking candidates.
    if not nested() or FORCED_BRANCH_NAMES.get():
        return _name_unabridged(smiles)
    entry = _NESTED_NAMES.get(smiles)
    if entry is None:
        start = reason_count()
        try:
            outcome = (_name_unabridged(smiles), None)
        except Exception as error:
            outcome = (None, error)
        if len(_NESTED_NAMES) >= _NESTED_NAMES_MAX:
            _NESTED_NAMES.clear()
        entry = _NESTED_NAMES[smiles] = (outcome, tuple(reasons_since(start)))
    else:
        replay(entry[1])
    name, error = entry[0]
    if error is not None:
        raise error
    return name


def _name_unabridged(smiles: str) -> str:
    enter()
    result = None
    try:
        result = _name_unabridged_body(smiles)
        return result
    finally:
        leave(result)


def _name_unabridged_body(smiles: str) -> str:
    name = None
    lambda_token = None
    try:
        parsed = _parse_smiles(smiles)
        if parsed is not None and outermost() and Chem.MolToSmiles(parsed) in _METHYLBENZENES:
            return _METHYLBENZENES[Chem.MolToSmiles(parsed)]
        if parsed is not None and any(a.GetFormalCharge() for a in parsed.GetAtoms()):
            lambda_token = CITE_SKELETAL_LAMBDA.set(False)
        if parsed is not None and has_carbon_monoxide_shape(parsed):
            return name_carbon_monoxide(parsed)
        if parsed is not None and has_carbene_amine_shape(parsed):
            return name_carbene_amine(parsed)
        if parsed is not None and has_hydride_ylium_shape(parsed):
            return name_hydride_ylium(parsed)
        if parsed is not None and has_poly_ylium_shape(parsed):
            return name_poly_ylium(parsed)
        if parsed is not None and has_group_polycation_shape(parsed):
            return name_group_polycation(parsed)
        if parsed is not None and has_ring_nitrenium_shape(parsed):
            return name_polycation(parsed)
        if parsed is not None and has_skeleton_radical_ion_shape(parsed):
            return name_skeleton_radical_ion(parsed)
        if parsed is not None and parsed.HasProp("_hypervalent_anion"):
            name = name_anion(parsed)
            return name
        if parsed is not None and has_inositol_derivative_shape(parsed):
            return name_inositol_derivative(parsed)
        if parsed is not None and has_common_hydride_shape(parsed):
            name = name_common_hydride(parsed)
            return name
        if parsed is not None and has_anisole_shape(parsed):
            return name_anisole(parsed)
        if parsed is not None and has_noncarbon_oxoacid_shape(parsed):
            return name_noncarbon_oxoacid(parsed)
        if parsed is not None and has_functional_replacement_oxoacid_shape(parsed):
            return name_functional_replacement_oxoacid(parsed)
        if parsed is not None and parsed.GetNumAtoms() > 4:
            try:
                return name_polynuclear_oxoacid(parsed, priority=True)
            except UnsupportedStructure:
                pass
        if parsed is not None and has_spiro_hub_atom_shape(parsed):
            return name_spiro_hub_atom(parsed)
        if parsed is not None and has_alternating_cage_shape(parsed):
            cage = name_alternating_cage(parsed)
            if cage is not None:
                return cage
        if parsed is not None and has_nonstandard_hydride_shape(parsed):
            return name_nonstandard_hydride(parsed)
        if parsed is not None and has_sphingoid_shape(parsed):
            return name_sphingoid(parsed)
        if parsed is not None and has_nucleoside_name(parsed):
            return name_nucleoside(parsed)
        if parsed is not None and oligonucleotide_name(parsed) is not None:
            return oligonucleotide_name(parsed)
        if parsed is not None and has_substituted_nucleoside_name(parsed):
            name = name_substituted_nucleoside(parsed)
            return name
        if parsed is not None and has_chalcone_shape(parsed):
            return name_chalcone(parsed)
        if parsed is not None and silicic_cyanate_name(parsed) is not None:
            return silicic_cyanate_name(parsed)
        if parsed is not None and borane_silane_amide_name(parsed) is not None:
            return borane_silane_amide_name(parsed)
        if parsed is not None:
            polyborane = polyborane_name(parsed)
            if polyborane is not None:
                return polyborane
            parts = lewis_adduct_mol(parsed)
            if parts is not None and has_adduct_shape(parts):
                return mark(name_adduct(parts, smiles_to_iupac), _NO_PIN_ADDUCT)
        if parsed is not None:
            diacylamine = diacylamine_name(parsed)
            if diacylamine is not None:
                return diacylamine
            diacylhydrazine = diacylhydrazine_name(parsed)
            if diacylhydrazine is not None:
                return diacylhydrazine
            chain_amide = chain_amide_name(parsed)
            if chain_amide is not None:
                return chain_amide
            chain_ketone = chain_ketone_name(parsed)
            if chain_ketone is not None:
                return chain_ketone
            carbonic_hydrazidine = carbonic_hydrazide_name(parsed)
            if carbonic_hydrazidine is not None:
                return carbonic_hydrazidine
        if parsed is not None and has_dipolar_shape(parsed):
            return name_dipolar(parsed)
        if parsed is not None and has_sugar_lactone_shape(parsed):
            return name_sugar_lactone(parsed)
        if parsed is not None and has_glycoside_shape(parsed):
            return name_glycoside(parsed)
        if parsed is not None and has_sugar_alcohol_acid_shape(parsed):
            return name_sugar_alcohol_acid(parsed)
        if parsed is not None and sugar_acid_derivative_name(parsed) is not None:
            return sugar_acid_derivative_name(parsed)
        if parsed is not None and has_substituted_sugar_shape(parsed):
            return name_substituted_sugar(parsed)
        if parsed is not None and has_peptide_shape(parsed):
            return name_peptide(parsed)
        if parsed is not None and has_amino_acid_shape(parsed):
            return name_amino_acid(parsed)
        if parsed is not None and not has_sphingoid_shape(parsed):
            if has_carbonless_acyl(parsed):
                try:
                    skeletal = name_heteroacyclic(parsed)
                except UnsupportedStructure:
                    skeletal = None
                if skeletal is not None:
                    return skeletal
            stereo_specified = _has_specified_stereo(parsed)
            for namer in (
                name_acid_salt,
                name_polycarbonic,
                name_carbonic_family,
                name_acid_derivative,
                name_boron_peroxy_ester,
                name_hetero_parent_acid,
                name_phosphorous_acid,
            ):
                try:
                    candidate = namer(parsed)
                except UnsupportedStructure:
                    continue
                # these namers do not cite stereodescriptors, so a flat name would silently drop the stereo
                if stereo_specified and namer is not name_acid_salt and not _cites_every_stereo_element(parsed, candidate):
                    continue
                return candidate
        if parsed is not None and has_spiro_union_shape(parsed):
            try:
                name = name_spiro_union(parsed)
                return name
            except UnsupportedStructure:
                pass
        beyond_preferred = None
        if parsed is not None:
            name = name_heteroacyclic(parsed)
            if name is not None:
                return name
            name = name_nitrogen_methylene_multiplicative(parsed)
            if name is not None:
                return name
            name = name_appendix3_skeleton(parsed)
            if name is not None:
                return name
            natural, operations = name_natural_product_ranked(parsed)
            if natural is not None:
                phane_name = linear_phane_pin(parsed)
                if phane_name is not None:
                    return phane_name
            if natural is not None and operations <= PREFERRED_OPERATIONS:
                return natural
            beyond_preferred = natural
            steroid = name_steroid(parsed) if parsed.GetRingInfo().NumRings() == 4 else None
            if steroid is not None:
                return steroid
            if has_chain_multiplicative_shape(parsed) and not has_phosphate_shape(parsed) and not has_borinic_acid_shape(parsed) and not has_boronic_acid_shape(parsed):
                try:
                    return name_polyfunctional(parsed)
                except UnsupportedStructure:
                    pass
            if _has_aromatic_oxo(parsed):
                try:
                    return name_polyfunctional(parsed)
                except UnsupportedStructure:
                    pass
        hydro_fusion = parsed is not None and _has_hydro_fusion_system(parsed)
        try:
            name = _smiles_to_iupac_dispatch(smiles)
            if parsed is not None and _drops_anionic_charge(parsed, name):
                raise UnsupportedStructure("the negative charge of this structure is not cited by any supported name")
        except UnsupportedStructure as original:
            try:
                name = _run_fallbacks(smiles, original)
            except UnsupportedStructure:
                if beyond_preferred is None:
                    raise
                name = beyond_preferred
        except Exception:
            name = _hydro_fusion_name(parsed) if hydro_fusion else None
            if name is None:
                raise
        if hydro_fusion and re.search(r"cyclo\[", name):
            name = _hydro_fusion_name(parsed) or name
        if parsed is not None and _has_specified_stereo(parsed):
            tokens = bool(_STEREO_TOKENS.search(name) or _STEREO_IN_RETAINED_NAME.search(name))
            if tokens and not _cites_every_double_bond(parsed, name):
                raise UnsupportedStructure("a stereodefined double bond of this structure is not cited by any supported name")
            if tokens and name.startswith("("):
                cited = _engine_name(parsed)
                if cited is not None and cited.startswith("(") and cited != name:
                    strip = lambda text: re.sub(r"^\([^()]*\)-", "", text)
                    if strip(cited) == strip(name):
                        name = cited
            if not tokens or parsed.GetRingInfo().NumRings() == 0:
                cited = _engine_name(parsed)
                strip = (
                    (lambda text: _DESCRIPTOR_GROUP.sub("", text))
                    if not tokens
                    else (lambda text: re.sub(r"^\([^()]*\)-", "", text))
                )
                if cited is not None and _STEREO_TOKENS.search(cited) and (
                    strip(cited) == strip(name) or (not tokens and parsed.GetRingInfo().NumRings() == 0)
                ):
                    name = cited
                elif not tokens:
                    raise UnsupportedStructure("the stereochemistry of this structure is not cited by any supported name")
        return _sulfinyl_descriptor_in_front(name)
    finally:
        if lambda_token is not None:
            CITE_SKELETAL_LAMBDA.reset(lambda_token)


_SULFINYL_DESCRIPTOR = re.compile(r"\[\(([RS])\)-([a-z]+(?:sulfinyl|seleninyl|tellurinyl))\]([a-z]+)")


def _sulfinyl_descriptor_in_front(name):
    """P-93.3.4.1: with a single sulfinyl-type group on a simple parent the descriptor of its stereogenic atom leads the
    name, '(S)-(methanesulfinyl)ethane'."""
    match = _SULFINYL_DESCRIPTOR.fullmatch(name)
    return f"({match.group(1)})-({match.group(2)}){match.group(3)}" if match else name


_ANION_NAME_ENDING = re.compile(r"(?:ide|uide|ate|ite|ato|ido|elide)\b|(?:ide|uide|ate|ite)-|id(?:yl|ylidene|ylidyne)\b|-id-\d")


def _drops_anionic_charge(mol, name) -> bool:
    if len(Chem.GetMolFrags(mol)) != 1 or sum(a.GetFormalCharge() for a in mol.GetAtoms()) >= 0:
        return False
    return _ANION_NAME_ENDING.search(name) is None


_HYDRO_FUSION_RUNNING = set()
_STEREO_TOKENS = re.compile(
    r"(?<=[\d'a-z⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])[RSEZ](?=[,)])|(?<=[\d'a⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])[rs](?=[,)])|\((?:R|S|E|Z)\)|\b(?:[DL]|alpha|beta)-|\((?:T|SP|SS|TBPY|OC|SPY|TPR|PBPY|CU|SAPR|TPRS)-|cis-|trans-|rel-|rac-"
)


_DESCRIPTOR_GROUP = re.compile(r"\((?:\d+[a-z]?[\u2032']*)?[RSEZrs](?:,(?:\d+[a-z]?[\u2032']*)?[RSEZrs])*\)-")
_STEREO_IN_RETAINED_NAME = re.compile(r"inositol|(?:adenos|guanos|inos|xanthos|cytid|urid|thymid)in")


def _cites_every_stereo_element(mol, name) -> bool:
    cited = len(_STEREO_TOKENS.findall(name))
    if re.search(r"\b(?:bis|tris|tetrakis)\b|\b(?:di|tri|tetra)\(", name):
        return cited > 0
    return cited >= sum(1 for s in Chem.FindPotentialStereo(mol) if s.specified == Chem.StereoSpecified.Specified)


def _cites_every_double_bond(mol, name) -> bool:
    if re.search(r"\b(?:bis|tris|tetrakis)\b|\b(?:di|tri|tetra)\(", name):
        return True
    bonds = sum(
        1
        for s in Chem.FindPotentialStereo(mol)
        if s.type == Chem.StereoType.Bond_Double and s.specified == Chem.StereoSpecified.Specified
    )
    return len(re.findall(r"(?<=[\d'a-z⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])[EZ](?=[,)])|\([EZ]\)", name)) >= bonds


def _has_specified_stereo(mol) -> bool:
    return any(s.specified == Chem.StereoSpecified.Specified for s in Chem.FindPotentialStereo(mol))


def _engine_name(mol):
    return _hydro_fusion_name(mol)


def _hydro_fusion_name(mol):
    key = Chem.MolToSmiles(mol)
    if key in _HYDRO_FUSION_RUNNING:
        return None
    _HYDRO_FUSION_RUNNING.add(key)
    try:
        return name_polyfunctional(mol)
    except UnsupportedStructure:
        return None
    finally:
        _HYDRO_FUSION_RUNNING.discard(key)


def _has_hydro_fusion_system(mol) -> bool:
    from ._diester_ring_diyl import _system_of

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if any(ring_info.NumAtomRings(a) > 1 for a in ring):
            _, atoms = _system_of(mol, ring[0])
            if is_hydro_fusion_system(mol, atoms):
                return True
    return False


def _has_aromatic_oxo(mol) -> bool:
    return any(
        b.GetBondTypeAsDouble() == 2.0
        and not b.GetIsAromatic()
        and ((b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetAtomicNum() == 8) or (b.GetEndAtom().GetIsAromatic() and b.GetBeginAtom().GetAtomicNum() == 8))
        for b in mol.GetBonds()
    )


def _name_o_substituted_hydroxylamine(mol):
    if not has_o_substituted_hydroxylamine_shape(mol):
        raise UnsupportedStructure("not an O-substituted hydroxylamine")
    return name_o_substituted_hydroxylamine(mol)


def _name_chain_onium(mol):
    if not has_chain_onium_shape(mol):
        raise UnsupportedStructure("not an onium cation of a heteroatom chain")
    return name_chain_onium(mol)


def _name_hydride_onium(mol):
    if not has_hydride_onium_shape(mol):
        raise UnsupportedStructure("not an onium cation of a mononuclear hydride")
    return name_hydride_onium(mol)


def _run_fallbacks(smiles, original):
    mol = _parse_smiles(smiles)
    key = Chem.MolToSmiles(mol)
    if key in _FALLBACKS_RUNNING:
        raise original
    _FALLBACKS_RUNNING.add(key)
    try:
        for skeletal in (name_skeletal_chain, name_hetero_macrocycle):
            try:
                name = skeletal(mol)
            except UnsupportedStructure:
                continue
            if name is not None:
                return name
        for fallback in (
            name_chalcogen_aldehyde,
            name_carbonic_hydrazide,
            name_chalcogen_hydrazide,
            name_ring_nitrogen_hydrazide,
            name_condensed_guanidine,
            name_ring_heteroatom_nitrile,
            name_polynuclear_oxoacid,
            name_halogen_amide,
            name_halogen_acid_ester,
            name_halogen_oxo,
            name_polyfunctional,
            _name_o_substituted_hydroxylamine,
            _name_hydride_onium,
            _name_chain_onium,
            name_anion,
            name_ester_by_parts,
            name_chalcogen_chain_heterone,
        ):
            try:
                return fallback(mol)
            except UnsupportedStructure:
                continue
        name = _name_via_fallbacks(mol)
        if name is not None:
            return name
        try:
            return name_nonpreferred_cyclophane(mol)
        except UnsupportedStructure:
            pass
    finally:
        _FALLBACKS_RUNNING.discard(key)
    raise original


def _smiles_to_iupac_dispatch(smiles: str) -> str:
    mol = _parse_smiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    return _name_mol(mol)


def _name_via_fallbacks(mol):
    for contract in (contract_phosphanyl_groups, contract_hetero_groups_candidates):
        try:
            contracted = contract(mol)
            candidates = contracted if isinstance(contracted, list) else [contracted]
            names = []
            for candidate in candidates:
                if candidate is None:
                    continue
                try:
                    names.append(_name_mol(candidate))
                except UnsupportedStructure:
                    continue
            if not names:
                continue
        except UnsupportedStructure:
            continue
        names = [n for n in names if not ("iodo" in n and not any(a.GetAtomicNum() == 53 for a in mol.GetAtoms()))]
        if names:
            return min(names, key=alphanumerical_name_key)
    return None


def _kekule_forms_without_fusion_name(mol):
    """P-52.2.4.1: a bicyclic system with an aromatic ring but without two rings of five or more members has no
    preferred fusion name, so it is named as an unsaturated von Baeyer system from its Kekule structures."""
    if FUSION_NAME_REQUIRED.get():
        return []
    info = mol.GetRingInfo()
    if info.NumRings() != 2 or sum(len(ring) >= 5 for ring in info.AtomRings()) >= 2:
        return []
    if not any(a.GetIsAromatic() for a in mol.GetAtoms()) or find_bicyclic_core(mol) is None:
        return []
    base = Chem.Mol(mol)
    try:
        Chem.Kekulize(base, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return []
    forms = [base]
    (benzene,) = [ring for ring in info.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)] or [None]
    if benzene is not None and len(benzene) == 6:
        ring = set(benzene)
        ring_bonds = [b.GetIdx() for b in base.GetBonds() if b.GetBeginAtomIdx() in ring and b.GetEndAtomIdx() in ring]
        if sum(base.GetBondWithIdx(i).GetBondType() == Chem.BondType.DOUBLE for i in ring_bonds) == 3:
            flipped = Chem.RWMol(base)
            for i in ring_bonds:
                bond = flipped.GetBondWithIdx(i)
                bond.SetBondType(
                    Chem.BondType.SINGLE if bond.GetBondType() == Chem.BondType.DOUBLE else Chem.BondType.DOUBLE
                )
            other = flipped.GetMol()
            Chem.SanitizeMol(other, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
            forms.append(other)
    return forms


def _name_mol(mol) -> str:
    candidates = []
    token = PREFER_VON_BAEYER.set(True)
    try:
        for form in _kekule_forms_without_fusion_name(mol):
            try:
                candidates.append(_name_mol(form))
            except UnsupportedStructure:
                continue
    finally:
        PREFER_VON_BAEYER.reset(token)
    if candidates:
        # a double bond between non-consecutive atoms is cited as 1(6); the unparenthesised locants are cited first
        return min(candidates, key=lambda name: (len(re.findall(r"\d\(\d+\)", name)), name))
    fused = fused_ring_system_name(mol)
    if fused is not None:
        return fused

    # The 7 retained nucleoside names (P-105.1) are recognized by exact
    # whole-molecule match, so they must be routed before every other
    # branch below: a purine/pyrimidine base's fused-ring nitrogen pattern
    # crashes the ring-assembly-chain detector (thymidine/uridine, an
    # unrelated pre-existing bug) or gets misrouted as an amino-acid/
    # aromatic-ring shape long before any ring-count-specific check runs.
    if has_nucleoside_name(mol):
        return name_nucleoside(mol)

    # P-105.2 substituted nucleosides and P-106 nucleotides are recognized
    # here for dispatch-ordering reasons: the base's fused-ring nitrogens,
    # the sugar hydroxyls and the phosphate esters would otherwise be
    # claimed by unrelated generic branches below.
    oligonucleotide = oligonucleotide_name(mol)
    if oligonucleotide is not None:
        return oligonucleotide
    if has_substituted_nucleoside_name(mol):
        return name_substituted_nucleoside(mol)

    # The 7 retained metallocene names (P-69.2.7) are likewise recognized
    # by exact whole-molecule match and must be routed here for the same
    # reason: a bare metal atom plus two disjoint cyclopentadienyl/
    # cyclopentadienide fragments has no shape any other branch below
    # expects, and reaches a radical/carbanide dispatch that rejects the
    # multi-fragment SMILES outright long before any ring-count check.
    if has_ocene_shape(mol):
        return name_ocene(mol)

    # A bare metallacyclic parent hydride (P-69.4's skeletal-replacement
    # ring, e.g. '1-titanacyclobutane') has a transition-metal ring atom
    # too, for the same reason as the metallocene case above: no other
    # branch below recognizes a non-carbon ring atom, so it must be routed
    # here before the generic cycloalkane dispatch's own heteroatom
    # rejection would otherwise claim it.
    for has_shape, namer in (
        (has_metallacycle_shape, name_metallacycle),
        (has_metallafused_shape, name_metallafused),
        (has_metallapolycycle_shape, name_metallapolycycle),
    ):
        if has_shape(mol):
            try:
                return mark(namer(mol), _NO_PIN_ORGANOMETALLIC)
            except UnsupportedStructure as first:
                try:
                    return mark(name_metallacycle_as_group(mol), _NO_PIN_ORGANOMETALLIC)
                except UnsupportedStructure:
                    raise first

    # A polyspiro tree of monocycles that the linear and hub namers do not take (P-24.2.3) precedes every
    # heteroatom-parent dispatch below.
    if (
        has_spiro_tree_shape(mol)
        and find_linear_polyspiro_chain(mol) is None
        and find_branched_polyspiro_hub(mol) is None
    ):
        return name_spiro_tree(mol)

    # Group 3-12 metal complexes (P-69.2 coordination naming) must precede
    # every heteroatom-parent dispatch below, which would otherwise claim
    # a metal-bound phosphane/amine/ether ligand's donor atom.
    if has_coordination_shape(mol):
        if has_group1_2_organometallic_shape(mol):
            try:
                return mark(name_group1_2_organometallic(mol), _NO_PIN_ORGANOMETALLIC)
            except UnsupportedStructure:
                pass
        return mark(name_coordination(mol), _NO_PIN_ORGANOMETALLIC)

    # Two or more Group 13-15 metals (P-69.5.3) must precede the
    # single-metal hydride dispatches below, which reject a second metal.
    if has_metal_pair_shape(mol) and not _has_senior_principal_group(mol):
        return name_metal_pair(mol)
    if has_inorganic_acid_derivative_shape(mol):
        try:
            return name_inorganic_acid_derivative(mol)
        except UnsupportedStructure:
            pass
    for has_shape, namer in (
        # P-103.1.1.1: a common amino acid's retained name + L/D descriptor
        # must be routed here, before every ring-count/functional-group
        # dispatch branch below: a side chain recognized by `_amino_acid.py`'s
        # table can carry its own extra nitrogen (lysine, arginine) or its
        # own ring (phenylalanine, tyrosine, tryptophan), which would
        # otherwise be misrouted first -- confirmed empirically: arginine's
        # guanidino C=N was caught by `_imine.py`'s dispatch and tryptophan's
        # indole ring by the bicyclic-heteroatom dispatch, both well before
        # this check's original position further down ever ran.
        (has_peptide_shape, name_peptide),
        (has_amino_acid_shape, name_amino_acid),
        # A ring-system diester of one polyol (P-65.6.3.3.3) is claimed before every
        # ring/functional-group shape check below, which would misread its esters.
        (has_hydride_polyester_shape, name_hydride_polyester),
        (has_polyester_of_one_polyol_shape, name_diester_acyloxy),
        (has_hydride_carbo_suffix_shape, name_hydride_carbo_suffix),
        (has_lambda_ring_shape, name_lambda_ring),
        (has_ring_lambda_heterone_shape, name_ring_lambda_heterone),
        # A chalcogen ring-oxide (P-62.5's functional-class "oxide" pattern,
        # not limited to acyclic amines) breaks the ring's own aromaticity as
        # RDKit perceives it, so it must be routed here before any ring-shape
        # or aromatic dispatch below ever gets a chance to reject it outright.
        (has_hetero_ring_oxide_shape, name_hetero_ring_oxide),
        # The pyridinone tautomer (P-31.1.4.3.4's indicated-hydrogen oxo form)
        # keeps its ring-carbon aromatic despite the exocyclic oxo, so it must
        # be routed here before `_ketone.py`'s own generic aryl-ketone
        # rejection below ever gets a chance to claim it.
        (has_pyridinone_shape, name_pyridinone),
        # The pyrimidinone tautomer (a second ring nitrogen alongside the
        # pyridinone shape above) is checked right after it, for the same
        # aromatic-aryl-ketone dispatch-ordering reason.
        (has_pyrimidinone_shape, name_pyrimidinone),
        # The uracil/thymine diketo tautomer (both ring nitrogens carrying
        # their own indicated hydrogen, P-58.2.2's parenthesized multi-locant
        # convention) is checked right after the single-oxo pyrimidinone case
        # above, for the same aromatic-aryl-ketone dispatch-ordering reason.
        (has_pyrimidinedione_shape, name_pyrimidinedione),
        # An amino-acid/betaine-type zwitterion (P-74.1.3's ammonium-nitrogen-
        # prefix-on-a-carboxylate-parent citation order) must be routed here
        # before `has_salt_shape` below: it's a single connected fragment that
        # nonetheless has an ammonium-shaped nitrogen, which `_salt.py`'s own
        # cation loop would otherwise try (and fail) to name as if the entire
        # molecule were a bare ammonium cation, crashing rather than falling
        # through.
        (has_zwitterion_shape, name_zwitterion),
    ):
        if has_shape(mol):
            return namer(mol)

    # Neutral adducts and solvates (P-14.8, see _adduct.py) must likewise be
    # routed here before every other branch below, for the same reason.
    if has_adduct_shape(mol):
        return name_adduct(mol, smiles_to_iupac)

    vb_assembly_core = find_vb_ring_assembly_core(mol)
    if vb_assembly_core is not None:
        return name_vb_ring_assembly(mol, vb_assembly_core)

    # P-54.3: an assembly of three or more otherwise identical rings that mixes mancude and saturated rings takes
    # hydro prefixes, ahead of any substitutive or multiplicative name built on the saturated ring.
    if 3 <= mol.GetRingInfo().NumRings() <= 6:
        hydro_chain_core = find_ring_assembly_chain_core(mol)
        if hydro_chain_core is not None and hydro_chain_core[3]:
            try:
                return name_ring_assembly_chain(mol, hydro_chain_core)
            except UnsupportedStructure:
                pass

    multiplicative_name = name_if_multiplicative(mol)
    if multiplicative_name is not None:
        return multiplicative_name

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
        ylidene_core = find_ring_assembly_ylidene_core(mol)
        if ylidene_core is not None and any(a.GetAtomicNum() != 6 and a.IsInRing() for a in mol.GetAtoms()):
            try:
                return name_ring_assembly_ylidene(mol, ylidene_core)
            except UnsupportedStructure:
                pass
    for has_shape, namer in (
        (_has_isotope_label, _name_isotope_label),
        (has_o_substituted_hydroxylamine_shape, name_o_substituted_hydroxylamine),
        (has_uronium_shape, name_uronium),
        # An isotopically labeled hydroxyl oxygen and/or skeletal carbon
        # combined with the '-ol' suffix (P-82.5.1/P-82.5.2) must be routed
        # here before `_isotope.py`'s own plain chain/methane path just below,
        # which rejects any oxygen outright.
        (has_isotope_alcohol_shape, name_isotope_alcohol),
        # An isotopically labeled carbonyl oxygen and/or skeletal carbon
        # combined with the '-one' suffix (P-82.5.1/P-82.5.2) must likewise be
        # routed here before `_isotope.py`'s own plain path.
        (has_isotope_ketone_shape, name_isotope_ketone),
        # An isotopically labeled carboxyl oxygen and/or skeletal carbon
        # combined with the '-oic acid' suffix (P-82.5.1/P-82.5.2) must
        # likewise be routed here before `_isotope.py`'s own plain path.
        (has_isotope_carboxylic_acid_shape, name_isotope_carboxylic_acid),
        # An isotopically labeled atom (P-82's isotope descriptor nomenclature)
        # must be routed here before every other branch below: RDKit represents
        # an isotopically substituted hydrogen (e.g. 2H) as its own explicit
        # atom (atomic number 1), which none of the other branches recognize at
        # all -- every one of them would reject it outright as an unsupported
        # heteroatom.
        (has_isotope_shape, name_isotope),
        # A radical ion on an ionic suffix group (P-75.3.1's 'aminiumyl'
        # radical cation) has both a charge and a radical electron on the
        # same nitrogen -- checked ahead of `has_radical_shape` below, whose
        # own broader "any nonzero radical electron count" check would
        # otherwise claim it first and misroute it into the plain-radical
        # dispatch, which rejects any charged atom outright.
        (has_chain_ylium_shape, name_chain_ylium),
        (has_hetero_acylium_shape, name_hetero_acylium),
        (has_skeleton_radical_ion_shape, name_skeleton_radical_ion),
        (has_radical_ion_shape, name_radical_ion),
        # A radical center (P-71.2.1.1's 'yl' radical naming) must be routed
        # here before every other branch below: none of them recognize a
        # nonzero radical electron count at all -- an unbranched-chain or
        # monocyclic-ring radical would otherwise fall straight through to the
        # plain alkane/cycloalkane dispatch further down, which doesn't know a
        # hydrogen is missing.
        (_has_nameable_radical_group, name_radical_group),
        (has_radical_shape, name_radical),
        # A nitrogen ylide (P-74.2.1.1.1's zwitterionic anion-carbon-parent
        # naming) has its own charged nitrogen too -- checked before
        # `has_ammonium_shape` below, whose own check already explicitly
        # excludes this shape (a second charged atom elsewhere) rather than
        # naming it, per that module's own docstring.
        (has_nitrogen_ylide_shape, name_nitrogen_ylide),
        # A phosphorus/oxygen/sulfur ylide (P-74.2.1.1.2/.3/.4) has its own
        # charged cation atom too -- checked before `has_phosphonium_shape`/
        # `has_oxonium_shape`/`has_sulfonium_shape` below for the same reason
        # as the nitrogen ylide above.
        (has_pos_ylide_shape, name_pos_ylide),
        # An amine imide (P-74.2.1.3's zwitterionic hydrazinium-ide naming)
        # has two charged nitrogens (one +1, one -1) -- checked before
        # `has_ammonium_shape` below for the same reason as the nitrogen
        # ylide above: its own +1 nitrogen would otherwise match ammonium's
        # shape check and get misnamed as a plain quaternary ammonium,
        # silently dropping the -1 nitrogen fragment.
        (has_amine_imide_shape, name_amine_imide),
    ):
        if has_shape(mol):
            return namer(mol)

    # A charged ammonium nitrogen (P-73.1.1.2's hydron-addition cation
    # naming) must be routed here before every other branch below: none of
    # them recognize a charged atom at all -- `_amine.py` in particular
    # rejects any charged atom outright rather than attempting to name it.
    if has_chain_cation_shape(mol):
        return name_chain_cation(mol)
    if has_ylium_ring_shape(mol):
        try:
            return name_ylium_ring(mol)
        except UnsupportedStructure:
            pass
    if has_mixed_onium_shape(mol):
        return name_mixed_onium(mol)
    if has_polyammonium_shape(mol):
        return name_polyammonium(mol)
    if has_polycation_shape(mol):
        return name_polycation(mol)
    if has_ammonium_shape(mol):
        if has_phosphate_shape(mol) and any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
            return name_phosphate(mol)
        return name_ammonium(mol)
    for has_shape, namer in (
        # A secondary/tertiary amine N-oxide (P-62.5's zwitterionic N+-O-) has a
        # charged nitrogen too, for the same reason as ammonium above -- and
        # `has_ammonium_shape` itself doesn't match it (an oxide-bearing
        # nitrogen has 2-3 carbon neighbors plus the oxide oxygen, never the
        # single-carbon/three-H shape ammonium requires), so it needs its own
        # explicit routing here.
        (has_amine_oxide_shape, name_amine_oxide),
        # A diazonium cation (R-N#N+, P-73) has a charged nitrogen too, for the
        # same reason as ammonium above -- routed here, unconditionally,
        # before every other branch.
        (has_diazonium_shape, name_diazonium),
        # A phosphonium cation (P-73.1.1.2) has a charged phosphorus too, for
        # the same reason as ammonium above -- routed here, unconditionally,
        # before every other branch.
        (has_phosphonium_shape, name_phosphonium),
        # A sulfonium cation (P-73.1.1.2) has a charged sulfur too, for the
        # same reason as ammonium above -- routed here, unconditionally,
        # before every other branch.
        (has_sulfonium_shape, name_sulfonium),
        # An oxonium cation (P-73.1.1.2) has a charged oxygen too, for the
        # same reason as ammonium above -- routed here, unconditionally,
        # before every other branch.
        (has_oxonium_shape, name_oxonium),
        # An acylium cation (P-73.2.3.1's 'oylium'/'ylium' suffix naming) has
        # a charged carbon too, checked ahead of the plain carbenium case
        # below since a C=O double bond gives the cation carbon degree 2, not
        # `has_carbenium_shape`'s own required degree 3 -- the two shapes
        # never overlap, so order between them doesn't otherwise matter.
        (has_acylium_shape, name_acylium),
        # A carbenium cation (P-73.2.2.1.1's 'ylium' suffix naming) has a
        # charged carbon too, for the same reason as ammonium above -- routed
        # here, unconditionally, before every other branch.
        (has_group_cation_shape, name_group_cation),
        (has_carbenium_shape, name_carbenium),
        # The cyclopentadienide anion (P-72.2.2.1's ring worked example) has a
        # charged ring carbon too, but `_carbanide.py` below is explicitly
        # acyclic-only and would reject any ring outright -- so this shape must
        # be routed here first.
        (has_cyclopentadienide_shape, name_cyclopentadienide),
        # The benzenide anion (phenyl anion, P-72.2.2.1's other ring worked
        # example) has a charged ring carbon too, for the same reason as the
        # cyclopentadienide case above -- routed here, right alongside it.
        (has_benzenide_shape, name_benzenide),
        # A carbanion center (P-72.2.2.1's '-ide' suffix naming) has a charged
        # carbon too, the anionic mirror of carbenium above -- routed here,
        # unconditionally, before every other branch.
        (has_carbanide_shape, name_carbanide),
        # An all-silicon skeleton (P-21.2.1's silane chain naming) has no
        # carbon at all, so it must be routed here before every other branch
        # below, all of which assume at least one carbon atom.
        (has_silane_chain_shape, name_silane_chain),
        (has_carbon_monoxide_shape, name_carbon_monoxide),
        (has_methanediimine_shape, name_methanediimine),
        (has_ring_imine_shape, name_ring_imine),
    ):
        if has_shape(mol):
            return namer(mol)

    # An all-phosphorus, multi-atom skeleton (diphosphane, triphosphane,
    # ...) has no carbon at all and must be routed before
    # has_simple_phosphane_shape below, which would otherwise misname it
    # via `_phosphane.py`'s own explicit "more than one phosphorus atom"
    # rejection.
    if has_phosphane_chain_shape(mol) and mol.GetNumAtoms() > 1:
        return name_phosphane_chain(mol)
    for has_shape, namer in (
        # A phosphonic acid (P-67.1.1's R-P(=O)(OH)2) has two hydroxyl oxygens
        # on phosphorus that `_phosphanone.py`'s own phosphine-oxide shape
        # doesn't expect -- `has_phosphanone_shape` would otherwise still
        # match (it only checks for a single P=O) and then fail inside
        # `name_phosphanone`'s validation, so this must be routed first.
        (has_phosphonic_acid_shape, name_phosphonic_acid),
        # A phosphinic acid (P-67.1.1's R2-P(=O)-OH) has the same
        # phosphanone-shape collision as phosphonic acid above (its own
        # hydroxyl oxygen isn't expected by `_phosphanone.py`), so it must be
        # routed here for the same reason.
        (has_phosphinic_acid_shape, name_phosphinic_acid),
        (has_phosphorus_acid_derivative_shape, name_phosphorus_acid_derivative),
        # Diphosphoric acid (P-67.2.1's own preselected dinuclear-acid name)
        # has each phosphorus individually shaped like a phosphate ester (the
        # other phosphorus group standing in as the "R" of a P-O-R ester
        # oxygen), so `has_phosphate_shape` would otherwise also match it and
        # then fail inside `name_phosphate`'s own single-phosphorus-only
        # validation -- must be routed here first.
        (has_dinuclear_oxoacid_shape, name_dinuclear_oxoacid),
        # A phosphate ester (P-67.1.3.2's P(=O)(OR)3) has three P-O-R ester
        # oxygens that `_phosphanone.py`'s own phosphine-oxide shape doesn't
        # expect -- `has_phosphanone_shape` would otherwise still match (it
        # only checks for a single P=O) and then fail inside
        # `name_phosphanone`'s validation, so this must be routed first, same
        # reason as phosphonic/phosphinic acid above.
        (has_mononuclear_oxoacid_shape, name_mononuclear_oxoacid),
        (has_sulfuric_amide_shape, name_sulfuric_amide),
        (has_phosphate_shape, name_phosphate),
        # Thio analogues of the phosphate, phosphonate and phosphinate esters carry sulfur on the phosphorus.
        (has_phosphorus_thioester_shape, name_phosphorus_thioester),
        (has_hydroxylamine_acid_shape, name_hydroxylamine_acid),
        (has_formazan_shape, name_formazan),
        # A phosphite ester (P-67.1.3.2's P(OR)3, no P=O) has three P-O-R
        # ester oxygens that `_phosphane.py`'s own plain-phosphane shape
        # doesn't expect (it rejects any heteroatom besides its own
        # phosphorus outright) -- must be routed here first, same reason as
        # phosphate above.
        (has_phosphite_shape, name_phosphite),
        # A sulfate ester (P-67.1.3.2's S(=O)(=O)(OR)2) has two S-O-R ester
        # oxygens no other sulfur module expects (they all assume a direct
        # S-C bond) -- must be routed before any of them for the same
        # ether-oxygen-rejection reason phosphate/phosphite are routed early.
        (has_sulfate_shape, name_sulfate),
        # A sulfite ester (P-67.1.3.2's S(=O)(OR)2, one fewer double-bonded
        # oxygen than sulfate) needs the same early routing, for the same
        # ether-oxygen-rejection reason as sulfate above.
        (has_sulfite_shape, name_sulfite),
        # A nitrate ester (P-67.1.3.2's O-NO2) has an N-O-R ester oxygen
        # `_nitro.py`'s own nitrogen shape doesn't expect (that module
        # requires a direct N-C bond) -- routed here for the same early-ester
        # reasoning as sulfate/sulfite above.
        (has_nitrate_ester_shape, name_nitrate_ester),
        # A nitrite ester (P-67.1.3.2's O-N=O) needs the same early routing
        # as nitrate above, for the same N-O-R ester-oxygen reason.
        (has_nitrite_ester_shape, name_nitrite_ester),
        # Carbonic acid or one of its esters (P-65.2.1's O=C(OR)(OR')) has a
        # central carbon with two -O-R/-OH oxygens neither `_ether.py` (which
        # rejects an oxygen bonded to more than one heavy atom outright) nor
        # `_carboxylic_acid.py` (which expects exactly one -OH, not two)
        # expects -- routed here for the same early-ester reasoning as
        # sulfate/nitrate above.
        (has_carbonic_acid_shape, name_carbonic_acid),
        # A phosphine oxide (P-68.3.2.3.1's '-phosphanone' suffix, R-P(=O)<)
        # has its own phosphorus-bonded oxygen that `_phosphane.py` doesn't
        # expect at all (that module rejects any heteroatom besides its own
        # phosphorus outright) -- must be routed here first, before
        # has_simple_phosphane_shape below, for the same reason as
        # has_phosphane_chain_shape above.
        (has_phosphanone_shape, name_phosphanone),
        (has_phosphanimine_shape, name_phosphanimine),
        # Thiophosphoric acid (P-67.1.2's own preselected infix-modified
        # oxoacid name) has a phosphorus with 4 substituents (=S plus three
        # -OH), which `_phosphane.py` rejects outright (more than three
        # substituents) -- must be routed here first.
        (has_noncarbon_oxoacid_shape, name_noncarbon_oxoacid),
        (has_functional_replacement_oxoacid_shape, name_functional_replacement_oxoacid),
    ):
        if has_shape(mol):
            return namer(mol)

    # Phosphinine/phosphinoline/isophosphinoline (P-25's own P-ring
    # counterparts to pyridine/quinoline/isoquinoline) have P in a ring,
    # which `_phosphane.py` never expects -- routed here first.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        if has_hetero_monocyclic_name(mol):
            return name_hetero_monocyclic(mol)
        if has_hetero_monocyclic_substituent_name(mol):
            return name_hetero_monocyclic_substituent(mol)
    for has_shape, namer in (
        # A phosphorus atom (P-68's phosphane substitutive nomenclature) must
        # be routed here before every other branch below: none of them
        # recognize phosphorus at all, and a phosphane carbon substituent would
        # otherwise reach the plain acyclic-alkane/amine dispatch further down
        # with no phosphorus handling.
        (has_polyphosphane_shape, name_polyphosphane),
        (has_simple_phosphane_shape, name_simple_phosphane),
        # A boronic acid (P-68.1.4.1's R-B(OH)2) has two hydroxyl oxygens on
        # boron that `_borane.py`'s own plain-borane shape doesn't expect --
        # `has_simple_borane_shape` matches any molecule with a boron atom at
        # all, so this must be routed first, before it misfires on the two
        # -OH oxygens as unsupported heteroatoms.
        (has_boronic_acid_shape, name_boronic_acid),
        # A borinic acid (P-68.1.4.1's R2-B-OH) has the same borane-shape
        # collision as boronic acid above, so it must be routed here for the
        # same reason.
        (has_borinic_acid_shape, name_borinic_acid),
        # A boron atom (P-68's borane substitutive nomenclature, the same shape
        # as phosphane above with boron in place of phosphorus) must be routed
        # here for the same reason -- none of the branches below recognize
        # boron at all.
        (has_simple_borane_shape, name_simple_borane),
        # A Group 13 metal (Al/Ga/In/Tl, P-69.1) is the same substitutive-
        # naming shape as boron/phosphorus above, generalized as one shared
        # mechanism -- must be routed here for the same reason: none of the
        # branches below recognize any of these elements at all.
        (has_group13_hydride_shape, name_group13_hydride),
        (has_group14_hydride_shape, name_group14_hydride),
        (has_group15_hydride_shape, name_group15_hydride),
    ):
        if has_shape(mol):
            return namer(mol)

    # A Group 1/2 metal (Li/Na/K/Mg/Ca, P-69.3) is a different additive
    # naming mechanism from Group 13 above, but the same reasoning for
    # dispatch order applies: none of the branches below recognize any of
    # these elements at all.
    if has_group1_2_organometallic_shape(mol):
        return mark(name_group1_2_organometallic(mol), _NO_PIN_ORGANOMETALLIC)

    # buckminsterfullerene (P-27's '[60]fullerene', a fixed 12-pentagon/
    # 20-hexagon cage) is recognized by exact whole-molecule match --
    # see _fullerene.py's module docstring; none of the ring modules
    # below understand a cage shape at all.
    if has_fullerene_name(mol):
        return name_fullerene(mol)
    if has_substituted_fullerene_cage(mol):
        require_defined_fullerene_numbering(mol, {a for r in mol.GetRingInfo().AtomRings() for a in r})
    for has_shape, namer in (
        # [2.2]paracyclophane/[2.2]metacyclophane (P-26's phane nomenclature
        # retained-name-style recognition, see module docstring) are recognized
        # the same way, independent of every other branch below -- their two
        # -CH2CH2- bridges make their carbon skeletons look like bridged
        # aromatic ring systems to every other dispatch branch, none of which
        # understand phane nomenclature at all.
        (has_cyclophane_name, name_cyclophane),
        (has_linear_phane_shape, name_linear_phane),
        # The seven 1989 IUPAC steroid parent ring hydrides (gonane through
        # ergostane, Rule 2.1/3S-2.2/2.3/2.4 -- see module docstring) are
        # recognized the same way, independent of every other branch below:
        # `_polycyclic.py`'s general von Baeyer engine already names the bare
        # gonane skeleton (confirmed by direct testing), so this check must
        # come first or gonane would never be reached.
        (has_steroid_parent_hydride_name, name_steroid_parent_hydride),
        # A steroid parent hydride with exactly one ring C=C double bond at a
        # standard, non-ring-fusion locant (e.g. 'androst-5-ene') must be
        # checked right alongside the bare-skeleton case above, for the same
        # von-Baeyer-engine-would-otherwise-claim-it reason.
        (has_steroid_unsaturated_name, name_steroid_unsaturated),
        # A steroid parent hydride whose A-ring is aromatic (the mancude
        # 1,3,5(10)-triene, e.g. 'estra-1,3,5(10)-triene') is checked right
        # alongside the single-double-bond case above, for the same reason.
        (has_steroid_aromatic_a_ring_name, name_steroid_aromatic_a_ring),
        # 2,3-didehydrooxepane etc. (P-31.2.2/P-31.2.4.1's 'didehydro' prefix,
        # adding one ring double bond to a saturated Hantzsch-Widman/retained
        # parent) -- routed here before the exact-match check below, since a
        # didehydro ring's extra double bond means it never matches that
        # check's fully-saturated canonical SMILES anyway, but grouped here
        # for the shared `saturated_ring_name` dependency.
        (has_didehydro_ring_name, name_didehydro_ring),
        (has_heteroaryne_shape, name_heteroaryne),
        # oxirane/thiane/piperidine etc. (P-22.2.1's Hantzsch-Widman
        # saturated-monocyclic retained names) are recognized the same way --
        # see _hetero_monocyclic.py's module docstring; none of the O/N
        # branches below understand a plain heteroatom ring at all.
        (has_hetero_monocyclic_name, name_hetero_monocyclic),
        # A single substituent on one of the same mancude parents above (P-22.2.1
        # heteroatom locants stay fixed; only one ring atom's H is replaced) --
        # see _hetero_monocyclic.py's module docstring for the role-sequence
        # matching this uses instead of the exact-match table above.
        (has_hetero_monocyclic_substituent_name, name_hetero_monocyclic_substituent),
        # 2H-pyran/4H-pyran (P-25.7.1.3.1's indicated-hydrogen case) -- unlike
        # furan/thiophene above, RDKit doesn't treat this ring as aromatic at
        # all, so it needs its own recognition shape rather than an extension
        # of `_ROLE_SEQUENCES`; must be routed here before `_ether.py` below,
        # which otherwise rejects any ring outright.
        (has_pyran_indicated_hydrogen_name, name_pyran_indicated_hydrogen),
    ):
        if has_shape(mol):
            return namer(mol)

    # A fused ring system with bridges (P-25.4) is preferred to a von Baeyer name (P-52.2.5.2.1); it is checked here, after
    # the retained steroid parents and their cyclo/seco/nor modifications, which keep their own names.
    bridged = bridged_ring_system_name(mol)
    if bridged is not None:
        return bridged

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
    for has_shape, namer in (
        # A sulfonic acid coexisting with a thiol (P-41/P-43, see
        # `_seniority.py`) has the same three-oxygen sulfonic sulfur as plain
        # sulfonic acid below, plus an extra thiol sulfur that `_sulfonic_
        # acid.py`'s own validation would otherwise reject outright -- must be
        # routed here first.
        (has_sulfonic_acid_thiol_shape, name_sulfonic_acid_thiol),
        # A sulfonic acid coexisting with a sulfinic acid (P-41/P-43, see
        # `_seniority.py`) has the same four-oxygen sulfonic sulfur as plain
        # sulfonic acid below, plus an extra sulfinic sulfur that neither
        # `_sulfonic_acid.py` nor `_sulfinic_acid.py`'s own validation would
        # accept -- must be routed here first.
        (has_sulfonic_acid_sulfinic_acid_shape, name_sulfonic_acid_sulfinic_acid),
        # A sulfonic acid coexisting with an unsubstituted sulfonamide
        # (P-41/P-43, see `_seniority.py`) has the same four-oxygen sulfonic
        # sulfur as plain sulfonic acid below, plus an extra sulfonamide
        # sulfur that neither `_sulfonic_acid.py` nor `_sulfonamide.py`'s own
        # validation would accept -- must be routed here first.
        (has_sulfonic_acid_sulfonamide_shape, name_sulfonic_acid_sulfonamide),
        # A carboxylic acid coexisting with a sulfonic acid (P-41/P-43, see
        # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
        # would otherwise reject on sight of the extra sulfonic sulfur, and the
        # same four-oxygen sulfonic sulfur `_sulfonic_acid.py` would otherwise
        # reject on sight of the extra -COOH oxygens -- must be routed here
        # first, before either single-group module.
        (has_carboxylic_acid_sulfonic_acid_shape, name_carboxylic_acid_sulfonic_acid),
        # A carboxylic acid coexisting with a sulfinic acid (P-41/P-43, see
        # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
        # would otherwise reject on sight of the extra sulfinic sulfur, and the
        # same three-oxygen sulfinic sulfur `_sulfinic_acid.py` would otherwise
        # reject on sight of the extra -COOH oxygens -- must be routed here
        # first, before either single-group module.
        (has_carboxylic_acid_sulfinic_acid_shape, name_carboxylic_acid_sulfinic_acid),
        # A carboxylic acid coexisting with a seleninic acid (P-41/P-43, see
        # `_seniority.py`) has the same -COOH carbon shape `_carboxylic_acid.py`
        # would otherwise reject on sight of the extra seleninic selenium, and
        # the same three-oxygen-cluster seleninic selenium
        # `_seleninic_acid.py` would otherwise reject on sight of the extra
        # -COOH oxygens -- must be routed here first, before either
        # single-group module.
        (has_carboxylic_acid_seleninic_acid_shape, name_carboxylic_acid_seleninic_acid),
        # A carboxylic acid coexisting with an unsubstituted sulfonamide
        # (P-41/P-43, see `_seniority.py`) has the same -COOH carbon shape
        # `_carboxylic_acid.py` would otherwise reject on sight of the extra
        # sulfonamide sulfur, and the same four-oxygen-cluster sulfonamide
        # sulfur `_sulfonamide.py` would otherwise reject on sight of the
        # extra -COOH oxygens -- must be routed here first, before either
        # single-group module.
        (has_carboxylic_acid_sulfonamide_shape, name_carboxylic_acid_sulfonamide),
        # A sulfonic acid (-SO3H, P-65.3.1) has three oxygens on its own sulfur,
        # so it must be routed here before the plain "any O atom" branch below --
        # none of the ether/ester/carboxylic-acid/aldehyde/ketone/alcohol checks
        # in that branch understand a sulfur-centered oxygen cluster at all.
        (has_sulfonic_acid_shape, name_sulfonic_acid),
        # A selenonic acid (-Se(=O)(=O)OH, P-65.3.1) has the same oxygen-cluster
        # shape as sulfonic acid above, just on selenium instead of sulfur, so
        # it too must be routed before the plain "any O atom" branch.
        (has_selenonic_acid_shape, name_selenonic_acid),
        # A telluronic acid (-Te(=O)(=O)OH, P-65.3.1) has the same oxygen-
        # cluster shape as sulfonic/selenonic acid above, just on tellurium, so
        # it too must be routed before the plain "any O atom" branch.
        (has_telluronic_acid_shape, name_telluronic_acid),
        # A tellurinic acid (-Te(=O)OH, P-65.3.1) has the same oxygen-cluster
        # shape as sulfinic/seleninic acid, just on tellurium, so it too must
        # be routed before the plain "any O atom" branch.
        (has_tellurinic_acid_shape, name_tellurinic_acid),
        # A seleninic acid (-Se(=O)OH, P-65.3.1) has the same oxygen-cluster
        # shape as sulfinic acid, just on selenium instead of sulfur, so it too
        # must be routed before the plain "any O atom" branch.
        (has_seleninic_acid_shape, name_seleninic_acid),
        # A sulfonyl group on a plain saturated ring nitrogen (e.g.
        # 1-methylsulfonylpiperidine) looks sulfonamide-shaped to
        # `has_sulfonamide_shape` below, which doesn't know about this
        # ring-as-parent construction and would misclaim/reject it -- must be
        # routed here first (narrower than `has_ring_amine_shape` alone, so it
        # doesn't also preempt `_hidden_amide_ketone.py`'s unrelated
        # acyl-on-ring-nitrogen shape further down).
        (has_ring_amine_sulfonyl_shape, name_ring_amine),
        # A sulfonamide (-SO2NH2, P-65.3.1) has two oxygens on its own sulfur,
        # the same reasoning as sulfonic acid above, plus a nitrogen that would
        # otherwise be mistaken for a plain amine -- so it too must be routed
        # before both the "any O atom" and "any N atom" branches below.
        (has_sulfonamide_shape, name_sulfonamide),
        # A sulfinic acid (-SO2H, P-65.3.1) has two oxygens on its own sulfur --
        # the same reasoning as sulfonic acid above -- so it too must be routed
        # before the plain "any O atom" branch.
        (has_sulfinic_acid_shape, name_sulfinic_acid),
        # A sulfinamide (-S(=O)NH2, P-65.3.1) has one oxygen and one nitrogen on
        # its own sulfur -- the same reasoning as sulfonamide above -- so it too
        # must be routed before both the "any O atom" and "any N atom" branches.
        (has_sulfinamide_shape, name_sulfinamide),
        # A sulfone (-SO2-, P-63.6) has two oxygens on its own sulfur, just like
        # a sulfinic/sulfonic acid's cluster above, so it must be routed here for
        # the same reason -- before it, since a sulfone's sulfur has two carbon
        # neighbors instead of the acid's hydroxyl, which would otherwise never
        # match `_sulfonic_acid.py`'s own shape check anyway, but routing it
        # alongside its acid relatives keeps this family together.
        (has_sulfone_shape, name_sulfone),
        # A sulfoxide (-S(=O)-, P-63.6) has one oxygen on its own sulfur, same
        # reasoning as the sulfinic acid check above.
        (has_sulfoxide_shape, name_sulfoxide),
        # A selenone (-Se(=O)(=O)-, P-63.6) is the selenium analogue of a
        # sulfone -- same reasoning, checked before the selenoxide below since
        # its selenium has two oxygens instead of one.
        (has_selenone_shape, name_selenone),
        # A selenoxide (-Se(=O)-, P-63.6) is the selenium analogue of a
        # sulfoxide.
        (has_selenoxide_shape, name_selenoxide),
        # A tellurone (-Te(=O)(=O)-, P-63.6) is the tellurium analogue of a
        # sulfone/selenone -- same reasoning, checked before the telluroxide
        # below since its tellurium has two oxygens instead of one.
        (has_tellurone_shape, name_tellurone),
        # A telluroxide (-Te(=O)-, P-63.6) is the tellurium analogue of a
        # sulfoxide/selenoxide.
        (has_telluroxide_shape, name_telluroxide),
        # A nitrone (imine N-oxide, P-74.2.1.2) has its own N+/O- dipole
        # pair the plain imine/oxime checks below don't expect, and it's
        # C=N-bonded (like an ordinary imine) so it would otherwise be
        # swallowed by the oxime-gated `has_simple_imine_shape` branch further
        # down and misrouted into `_imine.py`'s own rejection there -- must be
        # routed before it.
        (has_nitrone_shape, name_nitrone),
        # A nitrile oxide (P-74.2.2.2.1.2) has the same dipole-pair issue as
        # nitrone above, checked here for the same reason (before it would
        # otherwise fall through to a general heteroatom-allowlist rejection
        # further down, none of which know about this shape).
        (has_nitrile_oxide_shape, name_nitrile_oxide),
        (has_nitrile_oxide_prefix_shape, name_nitrile_oxide_prefix),
        # A nitro group (-NO2, P-61.5.1) has its own nitrogen and two oxygens
        # neither the ether/carbonyl/alcohol checks below nor the plain-amine
        # branch further down expect, so it must be routed before both -- a
        # nitro-bearing molecule always has an oxygen atom, so it would
        # otherwise be swallowed by the "any O atom" branch's unconditional
        # `name_alcohol` fallback and never even reach the nitrogen branch.
        (has_nitro_shape, name_nitro),
        # A nitroso group (-N=O, P-61.5.1's sibling prefix) has the same "own
        # oxygen" issue as nitro above, so it must be routed here for the same
        # reason.
        (has_nitroso_shape, name_nitroso),
        # An isocyanate group (-N=C=O, P-61.8) has the same "own oxygen" issue
        # as nitro above, so it must be routed here for the same reason.
        (has_isocyanate_shape, name_isocyanate),
    ):
        if has_shape(mol):
            return namer(mol)
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
    for has_shape, namer in (
        # Thiourea (H2N-C(=S)-NH2) has no oxygen at all, so it would otherwise
        # fall straight through the oxygen-gated block below (and every other
        # check in it) to the plain-amine fallback at the very end of this
        # function -- it must be checked here, unconditionally, before that
        # gate.
        (has_thiourea_shape, name_thiourea),
        # Selenourea/tellurourea (H2N-C(=Se/Te)-NH2) have no oxygen either, for
        # the same reason as thiourea above.
        (has_selenourea_shape, name_selenourea),
        (has_tellurourea_shape, name_tellurourea),
    ):
        if has_shape(mol):
            return namer(mol)

    if has_hydrazine_aminooxy_shape(mol):
        return name_hydrazine(mol)
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
        if has_cyclic_ketohexofuranose_shape(mol):
            return name_cyclic_ketohexofuranose(mol)
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
        if has_guanidine_shape(mol):
            return name_guanidine(mol)
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
            try:
                return name_hydrazone(mol)
            except UnsupportedStructure:
                pass
        # A hydrazine skeleton (H2N-NH2, P-68.3.1.2.1) has its own two
        # skeletal nitrogens (single-bonded, not double-bonded like
        # diazene above) with no carbon parent chain at all -- must be
        # routed before name_amine for the same reason as the checks
        # above.
        if has_hydrazine_multiplicative_shape(mol):
            return name_hydrazine_multiplicative(mol)
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
    for has_shape, namer in (
        (has_chalcogen_peroxol_shape, name_hydroperoxide),
        (has_thione_shape, name_thione),  # A thione (C=S, P-64.6.1) has no oxygen or nitrogen, so it only # reaches this branch once both are ruled out above. Its own shape # check is precise (a real C=S double bond), unlike thiol's/ # sulfide's own loose "any sulfur atom" checks, so it's safe to # check here regardless of order relative to them.
        (has_disulfide_shape, name_disulfide),  # A disulfide (R-S-S-R') has two sulfurs -- it would otherwise # look thiol-shaped to the check below (that check just looks for # the presence of any sulfur atom) -- must be routed here first.
        (has_sulfide_shape, name_sulfide),  # A plain -S- sulfide (P-63.2.1) has no suffix, so it must be routed # here before has_thiol_shape below: _thiol.py's validation rejects # a degree-2 sulfur outright (not a monovalent -SH), so a sulfide # would otherwise raise the wrong error there instead of being named.
        (has_thiol_shape, name_thiol),  # A thiol (-SH, P-63.1.1) has neither O nor N, so it only reaches # this branch once both are ruled out above (a thiol coexisting # with an amine is instead routed inside the nitrogen-gated branch # above, before its own `name_amine` fallback).
        (has_selone_shape, name_selone),  # A selone (C=Se, P-64.6.1) has the same precise-shape reasoning # as thione above (a real C=Se double bond), so it's safe to check # here regardless of order relative to the selenide/selenol chain.
        (has_diselenide_shape, name_diselenide),  # A diselenide (R-Se-Se-R') has two seleniums -- it would # otherwise look selenol-shaped to the check below (that check # just looks for the presence of any selenium atom) -- must be # routed here first.
        (has_selenide_shape, name_selenide),  # A plain -Se- selenide (P-63.2.1) has no suffix, so it must be # routed here before has_selenol_shape below for the same reason # as has_sulfide_shape above (that check doesn't look at degree # at all, so a degree-2 selenide would otherwise raise the wrong # error inside `name_selenol`'s degree-1 validation).
        (has_selenol_shape, name_selenol),  # A selenol (-SeH, P-63.1.1) is the next chalcogen analogue after # a thiol -- has neither O, N, nor S, so it only reaches this # branch once all three are ruled out above.
        (has_tellone_shape, name_tellone),  # A tellone (C=Te, P-64.6.1) has the same precise-shape reasoning # as thione/selone above, so it's safe to check here regardless of # order relative to the telluride/tellurol chain.
        (has_ditelluride_shape, name_ditelluride),  # A ditelluride (R-Te-Te-R') has two telluriums -- it would # otherwise look tellurol-shaped to the check below (that check # just looks for the presence of any tellurium atom) -- must be # routed here first.
        (has_telluride_shape, name_telluride),  # A plain -Te- telluride (P-63.2.1) has no suffix, so it must be # routed here before has_tellurol_shape below for the same reason # as has_selenide_shape above (that check doesn't look at degree # at all, so a degree-2 telluride would otherwise raise the wrong # error inside `name_tellurol`'s degree-1 validation).
        (has_tellurol_shape, name_tellurol),  # A tellurol (-TeH, P-63.1.1) is the next chalcogen analogue after # a selenol -- has neither O, N, S, nor Se, so it only reaches # this branch once all four are ruled out above.
    ):
        if has_shape(mol):
            return namer(mol)

    num_rings = mol.GetRingInfo().NumRings()
    # Aromatic rings carry non-single (order 1.5) bonds, which every other
    # ring module's non_single_bonds check rejects; an aromatic ring
    # system's carbon skeleton can also be graph-isomorphic to a *saturated*
    # bicyclic through pentacyclic core (e.g. naphthalene <-> decahydro-
    # naphthalene), so this check must run, and must succeed for any
    # in-scope aromatic shape, before num_rings==1 or any saturated
    # find_*_core below gets a chance to misdetect it and raise the wrong
    # ("unsaturated ... not supported yet") error (see _aromatic.py).
    if num_rings == 1:
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
        unsaturated_ring = find_cyclic_unsaturated_core(mol)
        if unsaturated_ring is not None:
            return name_cyclic_unsaturated(mol, unsaturated_ring)
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
    if has_substituted_fullerene_cage(mol):
        return name_cage_parent(mol)
    raise UnsupportedStructure(
        "polycyclic ring systems are not supported yet (see P-23/P-24/P-25)"
    )
