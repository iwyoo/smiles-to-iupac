import contextlib
import contextvars
import re

from rdkit import Chem, rdBase

from ._zwitterion import has_zwitterion_shape, name_zwitterion
from ._adduct import has_adduct_shape, has_bare_inorganic_component, has_inorganic_component, name_adduct
from ._acyclic import name_acyclic_alkane
from ._acid_derivatives import (
    has_isotope_acyl_halide_shape,
    in_acid_derivative_scope,
    name_acid_derivative,
    name_isotope_acyl_halide,
    name_labelled_anhydride,
    name_phosphorous_acid,
    name_sulfur_acid_ester,
)
from ._acid_salts import name_acid_salt
from ._hetero_carboxylic import name_hetero_parent_acid
from ._polycarbonic import name_polycarbonic
from ._carbonic_ester_anion import name_carbonic_ester_anion
from ._carbonic_family import name_carbonic_family
from ._carbamate import has_carbamate_shape, name_carbamate
from ._carbamate_ester import name_carbamate_ester
from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._acetyl_names import acetyl_names
from ._anion import name_anion
from ._amino_acid import SYSTEMATIC_ACID_PROBE
from ._locant_omission import CITE_SUBSTITUENT_LOCANTS
from ._amino_acid_derivative import has_amino_acid_shape, name_amino_acid
from ._peptide import has_peptide_shape, name_peptide
from ._mixed_onium import has_mixed_onium_shape, name_mixed_onium
from ._axial_stereo import annotate_ring_allene_stereo, cite_axial_stereo
from ._ring_diyl_numbering import SYSTEMATIC_FUSION
from ._spiro_stereo import cite_spiro_stereo
from ._chalcone import has_chalcone_shape, name_chalcone
from ._hydrogen_cation import hydrogen_salt_name
from ._ion_adduct import has_ion_adduct_shape, name_ion_adduct
from ._silicic_cyanate import silicic_cyanate_name
from ._borane_silane_amide import borane_silane_amide_name
from ._polyborane import lewis_adduct_mol, polyborane_name
from ._diacylamine import diacylamine_name
from ._chalcogen_chain_parent import name_chalcogen_chain_multiplicative, name_chalcogen_chain_parent
from ._chalcogen_polychain import name_polychalcogen_hydride
from ._carbonic_amide import (
    carbonohydrazonic_diamide_name,
    cyanamide_name,
    name_labelled_cyanamide,
    diazenecarbohydrazide_name,
    hydrazine_dicarboxamide_name,
    name_polycarbonic_amide,
    polycarbonic_amide_name,
)
from ._carbonic_hydrazide import carbonic_hydrazide_name, name_carbonic_hydrazide
from ._chalcogen_hydrazide import name_chalcogen_hydrazide
from ._diacylhydrazine import diacylhydrazine_name
from ._chain_amide import chain_amide_name, chain_ketone_name
from ._ring_nitrogen_hydrazide import name_ring_nitrogen_hydrazide
from ._alternating_cage import has_alternating_cage_shape, name_alternating_cage
from ._dipolar import azoxy_prefix_name, has_dipolar_shape, name_dipolar
from ._chalcogen_aldehyde import name_chalcogen_aldehyde
from ._condensed_guanidine import has_condensed_guanidine_shape, name_condensed_guanidine
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
from ._substituents import BRANCH_STEREO, FORCED_BRANCH_NAMES, ISOTOPE_LABELS
from ._glycoside import has_glycoside_shape, name_glycoside
from ._sugar_acid import has_ketoaldonic_chain_shape, name_ketoaldonic_chain, has_sugar_ring_acid_shape, name_sugar_ring_acid, has_sugar_lactam_shape, name_sugar_lactam, has_substituted_chain_shape, name_substituted_chain, has_sugar_alcohol_acid_shape, has_sugar_lactone_shape, sugar_carboxylate_name, name_sugar_lactone, name_sugar_alcohol_acid, sugar_acid_derivative_name
from ._sugar_substituted import has_substituted_sugar_shape, name_substituted_sugar, has_sugar_dianhydride_shape, name_sugar_dianhydride
from ._histidine import has_histidine_shape, name_histidine
from ._proline import has_proline_shape, name_proline
from ._aldehyde_carboxylic_acid import (
    has_aldehyde_carboxylic_acid_shape,
    name_aldehyde_carboxylic_acid,
)
from ._carboxylic_acid_amide import has_carboxylic_acid_amide_shape, name_carboxylic_acid_amide
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
from ._hidden_amide_ketone import has_hidden_amide_shape, name_hidden_amide_ketone
from ._thiourea import has_thiourea_shape, name_thiourea
from ._selenourea import has_selenourea_shape, name_selenourea
from ._tellurourea import has_tellurourea_shape, name_tellurourea
from ._urea import has_urea_shape, name_urea
from ._chalcogenourea import carries_only_element_hydride_groups
from ._guanidine import has_carbonic_amide_parent, has_guanidine_shape, name_guanidine
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
from ._diazonium import has_diazonium_shape, has_isotope_diazonium_shape, name_diazonium, name_isotope_diazonium
from ._radical_hetero import hetero_radical_name
from ._radical_ion import has_radical_ion_shape, name_radical_ion
from ._aromatic import find_aromatic_fused_core, name_aromatic_fused
from ._bicyclic import find_bicyclic_core, name_bicycloalkane
from ._isotope_acid_hydrogen import has_isotope_acid_hydrogen_shape, name_isotope_acid_hydrogen
from ._appendix3_skeletons import name_appendix3_skeleton
from ._borane import has_simple_borane_shape, name_simple_borane
from ._boronic_acid import boron_acid_is_parent, has_boronic_acid_shape, name_boron_acid, name_boronic_acid
from ._borinic_acid import has_borinic_acid_shape, name_borinic_acid
from ._metal_pair import has_metal_pair_shape, name_metal_pair
from ._coordination import has_coordination_shape, name_coordination
from ._group1_2_organometallic import has_group1_2_organometallic_shape, name_group1_2_organometallic
from ._chalcogen_imide import has_chalcogen_imide_shape, name_chalcogen_imide
from ._nonstandard_hydride import has_nonstandard_hydride_shape, name_nonstandard_hydride
from ._sulfonate import has_sulfonate_shape, name_sulfonate
from ._acid_groups import acid_group_at
from ._common import CITE_SKELETAL_LAMBDA, InternalInvariantError, InvalidSmilesError, has_hydrazone_shape, has_simple_imine_shape, has_amide_shape, has_amidine_shape, has_imide_shape, find_carboxylic_acid_carbon, HALOGEN_PREFIXES, adjacency, ISOTOPIC_CIP, STEREO_OVERRIDE, UnsupportedStructure, phantom_capped, genuine_stereo_centres, alphanumerical_name_key, assign_cip_labels, non_single_bonds, spec_chirality_as_parsed
from ._cyclic import name_cycloalkane
from ._disjoint_ring_substituents import find_disjoint_ring_pair_core, name_disjoint_ring_pair
from ._cyclic_unsaturated import find_cyclic_unsaturated_core, name_cyclic_unsaturated
from ._diester_acyloxy import has_diester_shape, has_polyester_of_one_polyol_shape, name_diester_acyloxy
from ._inorganic_polyester import has_inorganic_polyester_shape, name_inorganic_polyester
from ._hydride_polyester import has_chain_polyester_shape, has_hydride_polyester_shape, name_hydride_polyester
from ._ester_by_parts import name_ester_by_parts
from ._heteroacyclic import amine_suffix_capacity, name_heteroacyclic, has_carbonless_acyl
from ._nitrogen_methylene_multiplicative import name_nitrogen_methylene_multiplicative
from ._chain_multiplicative import has_chain_multiplicative_shape
from ._np import PREFERRED_OPERATIONS, name_natural_product_ranked
from ._steroid_named import name_steroid
from ._polyfunctional import AMINIUM, CATION_RING_ATOM, FORCED_PRINCIPAL, RING_CENTER, name_polyfunctional
from ._cyanate import has_cyanate_shape, name_cyanate
from ._ether import has_ether_shape, name_ether
from ._fusion_name import FUSION_NAME_REQUIRED, PREFER_VON_BAEYER, fused_ring_system_name
from ._hetero_prefixes import CATION_PARENT, _has_senior_principal_group
from ._fullerene import (
    has_fullerene_name,
    fullerene_cage_atoms,
    has_substituted_fullerene_cage,
    name_fullerene,
    require_defined_fullerene_numbering,
)
from ._fullerene_attached import name_fullerene_derivative
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
from ._fused_ring_yne import has_fused_ring_yne_shape, name_fused_ring_yne
from ._heteroaryne import has_heteroaryne_shape, name_heteroaryne
from ._didehydro_ring import has_didehydro_ring_name, name_didehydro_ring
from ._chalcogen_chain_heterone import name_chalcogen_chain_heterone
from ._halogen_acid_ester import name_halogen_acid_ester
from ._halogen_amide import has_halogen_amide_shape, name_halogen_amide
from ._inorganic_acid_derivative import has_inorganic_acid_derivative_shape, name_inorganic_acid_derivative
from ._halogen_oxo import name_halogen_oxo
from ._hydroxylamine import has_hydroxylamine_shape, name_hydroxylamine
from ._formazan import has_formazan_shape, name_formazan
from ._hydroxylamine_acid import has_hydroxylamine_acid_shape, name_hydroxylamine_acid
from ._hydroxylamine_general import has_o_substituted_hydroxylamine_shape, name_o_substituted_hydroxylamine
from ._dipole_oxide import (
    has_nitrile_oxide_prefix_shape,
    has_nitrile_oxide_shape,
    has_nitrone_shape,
    name_nitrile_oxide_prefix,
    name_nitrile_oxide,
    name_nitrone,
)
from ._isotope import (
    has_isotope_hydride_shape,
    has_isotope_substituted_hydride_shape,
    name_isotope_hydride,
    name_isotope_substituted_hydride,
)
from ._isotope_steroid import has_isotope_steroid_shape, name_isotope_steroid
from ._radical_ion_skeleton import has_skeleton_radical_ion_shape, name_skeleton_radical_ion
from ._radical_group import (
    has_group_cation_shape,
    has_radical_group_shape,
    has_radical_shape,
    name_group_cation,
    name_radical_group,
)
from ._isotope_one_carbon import has_isotope_one_carbon_shape, name_isotope_one_carbon
from ._selenocyanate import has_selenocyanate_shape, name_selenocyanate
from ._tellurocyanate import has_tellurocyanate_shape, name_tellurocyanate
from ._thiocyanate import has_thiocyanate_shape, name_thiocyanate
from ._azide import has_azide_shape, name_azide
from ._diazene import has_diazene_shape, name_diazene
from ._azine import has_azine_shape, name_azine
from ._hydrazine import has_hydrazine_aminooxy_shape, has_hydrazine_shape, name_hydrazine, name_hydrazine_with_heteroatom_groups
from ._hydrazine_multiplicative import has_hydrazine_multiplicative_shape, name_hydrazine_multiplicative
from ._diazo import has_diazo_shape, name_diazo
from ._isocyanate import has_isocyanate_shape, name_isocyanate
from ._isocyanide import has_isocyanide_shape, name_isocyanide
from ._isoselenocyanate import has_isoselenocyanate_shape, name_isoselenocyanate
from ._isotellurocyanate import has_isotellurocyanate_shape, name_isotellurocyanate
from ._isothiocyanate import has_isothiocyanate_shape, name_isothiocyanate
from ._inorganic_ester_parent import name_inorganic_ester_parent
from ._nitrate_ester import has_nitrate_ester_shape, name_nitrate_ester
from ._carbonic_acid import has_carbonic_acid_shape, name_carbonic_acid
from ._nitrite_ester import has_nitrite_ester_shape, name_nitrite_ester
from ._nitro import has_nitro_shape, name_nitro
from ._nitroso import has_nitroso_shape, name_nitroso
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._cyclophane import has_cyclophane_name, name_cyclophane, name_nonpreferred_cyclophane
from ._linear_phane import has_linear_phane_shape, linear_phane_pin, name_linear_phane
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
from ._phosphate import has_phosphate_shape, name_phosphate
from ._phosphorus_thioester import has_phosphorus_thioester_shape, name_boron_peroxy_ester, name_phosphorus_thioester
from ._dinuclear_oxoacid import has_dinuclear_oxoacid_shape, name_dinuclear_oxoacid
from ._phosphite import has_phosphite_shape, name_phosphite
from ._sulfate import has_sulfate_shape, name_sulfate
from ._sulfite import has_sulfite_shape, name_sulfite
from ._phosphonic_acid import has_phosphorus_acid_shape, name_phosphorus_acid
from ._phosphorus_acid_derivative import has_phosphorus_acid_derivative_shape, name_phosphorus_acid_derivative
from ._azonic_acid import has_azonic_acid_shape, name_azonic_acid
from ._phosphane_chain import has_phosphane_chain_shape, name_phosphane_chain
from ._ring_diyl_numbering import bridged_ring_system_name, is_hydro_fusion_system
from ._ring_system_seniority import system_type
from ._hydroperoxide import has_chalcogen_peroxol_shape, has_hydroperoxide_shape, name_hydroperoxide
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
from ._sulfone import has_sulfone_shape, name_sulfone
from ._sulfoxide import has_sulfoxide_shape, name_sulfoxide
from ._selenone import has_selenone_shape, name_selenone
from ._selenoxide import has_selenoxide_shape, name_selenoxide
from ._tellurone import has_tellurone_shape, name_tellurone
from ._telluroxide import has_telluroxide_shape, name_telluroxide
from ._diselenide import has_diselenide_shape, name_diselenide
from ._selenide import has_selenide_shape, name_selenide
from ._ditelluride import has_ditelluride_shape, name_ditelluride
from ._telluride import has_telluride_shape, name_telluride
from ._selenoic_acid import has_selenoic_acid_shape, name_selenoic_acid
from ._telluroic_acid import has_telluroic_acid_shape, name_telluroic_acid
from ._thioic_acid import has_thioic_acid_shape, name_thioic_acid
from ._hetero_chain import contract_hetero_groups_candidates, name_hetero_macrocycle, name_skeletal_chain
from ._nitro_group import contract_nitro_groups
from ._phosphanyl_group import PREFIX_PROP, contract_phosphanyl_groups
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


_NO_PIN_ADDUCT = "the Blue Book assigns no PIN to Lewis adducts, whose preferred names are coordination names (P-68.1.6.2)"
_NO_PIN_MIXED_ADDUCT = "the Blue Book assigns no PIN to adducts with inorganic components (P-14.8.2)"
_NO_PIN_ION_ADDUCT = "the Blue Book assigns no PIN to an adduct of ions with other components or of ions that do not balance (P-14.8.2)"
_NO_PIN_SKELETAL_AMINE = "the preferred name keeps the amine suffixes of the substitutive name, which a skeletal replacement name does not express (P-51.4.1.3)"
_NO_PIN_ORGANOMETALLIC ="the Blue Book defines no PIN for this class of organometallic compound (P-69.0)"

_FALLBACKS_RUNNING = set()


# P-22.1.3: toluene and the xylenes are preferred names of the unsubstituted hydrocarbons only
_METHYLBENZENES = {"Cc1ccccc1": "toluene", "Cc1ccccc1C": "1,2-xylene", "Cc1cccc(C)c1": "1,3-xylene", "Cc1ccc(C)cc1": "1,4-xylene"}


_ADAMANTANE = re.compile(r"(?<!bi)(?<!ter)(?<!quater)(?<!yclo)tricyclo\[3\.3\.1\.1\^3,7\]decan(?=e|-)")


_CUBANE = re.compile(r"(?<!bi)(?<!ter)(?<!quater)(?<!yclo)pentacyclo\[4\.2\.0\.0\^2,5\.0\^3,8\.0\^4,7\]octan(?=[a-z-])")


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


def _is_carbocycle_with_exocyclic_double_bond_perceived_aromatic(mol):
    """RDKit counts one electron for a ring carbon with an exocyclic C=C, so a quinodimethane looks like benzene; a
    benzene ring has no exocyclic double bond (P-29.3.4.2)."""
    rings = [r for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(a).GetIsAromatic() and mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in r)]
    return any(
        b.GetBondType() == Chem.BondType.DOUBLE
        and not b.IsInRing()
        and any(b.GetBeginAtomIdx() in r or b.GetEndAtomIdx() in r for r in rings)
        and b.GetBeginAtom().GetAtomicNum() == 6
        and b.GetEndAtom().GetAtomicNum() == 6
        for b in mol.GetBonds()
    )


def _is_aromatic_ring_with_triple_bond(mol):
    return any(b.GetBondType() == Chem.BondType.TRIPLE and b.IsInRing() and b.GetIsAromatic() for b in mol.GetBonds())


def _neutral_group15_oxides(mol):
    """P+-O-, As+-O- and Sb+-O- written as zwitterions, and the S+-O- of a sulfoxide with three bonds, are the doubly
    bonded oxides of the lambda-convention (P-74.2.1.4, P-74.2.2.1.8); RDKit writes some neutral P=O groups this way."""
    pairs = [
        (atom.GetIdx(), n.GetIdx())
        for atom in mol.GetAtoms()
        if atom.GetFormalCharge() == 1 and (atom.GetAtomicNum() in (15, 33, 51) or (atom.GetAtomicNum() in (16, 34, 52) and atom.GetDegree() == 3))
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


_C_GLYCOSYL = Chem.MolFromSmarts("[c]-[C]1O[C]([CH2][OX2])[C]([OX2])[C]([OX2])[C]1[OX2]")


def _parse_smiles(smiles):
    mol = _parse_smiles_unannotated(smiles)
    if mol is not None:
        annotate_ring_allene_stereo(mol, smiles)
    return mol


def _keep_radical_chirality(smiles, mol):
    """Perception of stereo ignores an unpaired electron, so a centre whose substituents differ only by a radical centre
    loses its tag at parse; the tags are restored for the namer that cites the radical as a phantom ligand (P-92.1.4.1)."""
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)
    try:
        tagged = Chem.MolFromSmiles(smiles)
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)
    if tagged is None or tagged.GetNumAtoms() != mol.GetNumAtoms():
        return mol
    for atom, source in zip(mol.GetAtoms(), tagged.GetAtoms()):
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED and not atom.GetNumRadicalElectrons():
            atom.SetChiralTag(source.GetChiralTag())
    return mol


def _parse_smiles_unannotated(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and (
        _is_nonbenzene_monocyclic_annulene(mol)
        or _is_aromatic_ring_without_double_bonds(mol)
        or _is_aromatic_ring_with_triple_bond(mol)
        or _is_carbocycle_with_exocyclic_double_bond_perceived_aromatic(mol)
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
            # RDKit knows only the thallium(I) valence, so the standard TlH3 of P-68.1.1.1 looks like a radical and a
            # thallium atom with one or two bonds looks saturated
            if atom.GetAtomicNum() == 81 and atom.GetNumRadicalElectrons() and atom.GetTotalValence() == 3:
                atom.SetNumRadicalElectrons(0)
            elif atom.GetAtomicNum() == 81 and atom.GetDegree() and not atom.GetFormalCharge() and atom.GetTotalValence() < 3:
                atom.SetNumRadicalElectrons(3 - atom.GetTotalValence())
        if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
            return _keep_radical_chirality(smiles, mol)
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
    if isinstance(smiles, str) and not _NAMING.get():
        with rdBase.BlockLogs():
            smiles = spec_chirality_as_parsed(smiles)
    return smiles_to_iupac_as_parsed(smiles)


def smiles_to_iupac_as_parsed(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError(f"smiles must be a str, not {type(smiles).__name__}")
    if not smiles.strip():
        raise InvalidSmilesError(f"invalid SMILES: {smiles!r}")
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(False)
    nested = _NAMING.get()
    token = _NAMING.set(True)
    labels_token = ISOTOPE_LABELS.set(None)
    try:
        with rdBase.BlockLogs():
            name = _name_checked(smiles)
        return name if nested else _group_primes(name)
    finally:
        ISOTOPE_LABELS.reset(labels_token)
        _NAMING.reset(token)
        Chem.SetUseLegacyStereoPerception(legacy)


_NAMING = contextvars.ContextVar("naming", default=False)
_PRIME_RUN = re.compile("[\u2032']{4,}")


def _group_primes(name):
    return _PRIME_RUN.sub(lambda m: " ".join(m.group()[i : i + 3] for i in range(0, len(m.group()), 3)), name)


def _name_checked(smiles: str) -> str:
    if "[H+]" in smiles:
        hydrogen_salt = hydrogen_salt_name(smiles, smiles_to_iupac)
        if hydrogen_salt is not None:
            return hydrogen_salt
    name = cite_spiro_stereo(smiles, cite_axial_stereo(smiles, _retained_polycycle_names(_smiles_to_iupac_unabridged(smiles))))
    mol = _parse_smiles(smiles)
    if mol is not None:
        if not nuclides_cited(mol, name):
            raise InternalInvariantError(f"a name without every nuclide of {smiles!r} left the naming pipeline: {name!r}")
        _require_radicals_cited(mol, name)
        _require_ring_triple_bonds_cited(mol, name)
        _require_cumulated_ring_atoms_cited(mol, name)
    if mol is not None and _has_free_anion(mol):
        name = acetyl_names(name)
    return name


_MULTIPLE_ENE = re.compile(r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)en(?:e|-|$)")


def _require_cumulated_ring_atoms_cited(mol, name):
    """A ring atom with two ring double bonds in a ring system named as a mancude parent is described only by the
    delta convention (P-25.7.2); a mancude name without it is the name of another compound."""
    ring_info = mol.GetRingInfo()
    rings = [set(r) for r in ring_info.AtomRings()]
    for atom in mol.GetAtoms():
        if sum(1 for b in atom.GetBonds() if b.IsInRing() and b.GetBondTypeAsDouble() == 2.0) < 2:
            continue
        system = {atom.GetIdx()}
        grown = True
        while grown:
            grown = False
            for ring in rings:
                if ring & system and not ring <= system:
                    system |= ring
                    grown = True
        members = [r for r in rings if r <= system]
        ortho_fused = any(len(a & b) == 2 for i, a in enumerate(members) for b in members[i + 1 :])
        mancude_named = ortho_fused or any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).GetIsAromatic() for i in system)
        if mancude_named and "δ" not in name and not _MULTIPLE_ENE.search(name):
            raise UnsupportedStructure("a ring atom with two ring double bonds needs the delta convention (P-25.7.2)")


def _require_ring_triple_bonds_cited(mol, name):
    ring_triple = any(
        b.GetBondTypeAsDouble() == 3.0 and b.IsInRing() and b.GetBeginAtom().GetAtomicNum() == 6 and b.GetEndAtom().GetAtomicNum() == 6
        for b in mol.GetBonds()
    )
    if ring_triple and "yn" not in name and "dehydro" not in name:
        raise UnsupportedStructure("the name does not cite the ring triple bond")


_RADICAL_ENDINGS = ("yl", "ylidene", "ylidyne", "yne", "ylium", "yloxy", "yliumyl")
# the parent hydride of an assembly of benzene rings ends in "yl" without being a radical
_PHENYL_ASSEMBLY_PARENT = re.compile(r"(?:bi|ter|quater|quinque|sexi|septi|octi|novi|deci)phenyl$")


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
        or any(a.GetFormalCharge() > 0 for a in mol.GetAtoms())
        or len(Chem.GetMolFrags(mol)) > 1
        or any(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds())
        or has_nonstandard_hydride_shape(mol)
    ):
        return
    if not name.rstrip(")]} ").endswith(_RADICAL_ENDINGS) or _PHENYL_ASSEMBLY_PARENT.search(name):
        raise UnsupportedStructure("the name does not cite the radical centre")


def _is_nucleoside_derivative(mol) -> bool:
    return has_nucleoside_name(mol) or oligonucleotide_name(mol) is not None or has_substituted_nucleoside_name(mol)


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


def _name_keeping_labels(mol, namer):
    """The name of `mol`, or None when `mol` has nuclide labels that this name leaves uncited."""
    name = namer(mol)
    return name if nuclides_cited(mol, name) else None


_NONCARBON_PEROXY_ESTER = Chem.MolFromSmarts("[#5,#14,#15,#16,#33,#34,#51,#52]OO[#6]")


def _name_isotope_label(mol) -> str | None:
    from ._isotope_labels import split_isotopes

    split = split_isotopes(mol)
    if split is not None and has_guanidine_shape(split[0]):
        try:
            return name_guanidine(mol)
        except UnsupportedStructure:
            pass
    if split is not None:
        name = name_labelled_cyanamide(mol)
        if name is not None:
            return name
    if split is not None and (
        has_amino_acid_shape(split[0])
        or has_inorganic_acid_derivative_shape(split[0])
        or has_noncarbon_oxoacid_shape(split[0])
        or has_functional_replacement_oxoacid_shape(split[0])
    ):
        try:
            name = _name_labelled_substituents(mol)
        except UnsupportedStructure:
            name = None
        if name is not None:
            return name
    if has_isotope_acid_hydrogen_shape(mol):
        try:
            return name_isotope_acid_hydrogen(mol)
        except UnsupportedStructure:
            pass
    try:
        name = name_polyfunctional(mol)
    except UnsupportedStructure:
        name = None
    if name is not None and nuclides_cited(mol, name):
        return name
    from ._nitrogen_parent_labels import (
        has_labelled_azoxy,
        has_labelled_diazene,
        has_labelled_hydrazone,
        name_labelled_azoxy,
        name_labelled_diazene,
        name_labelled_hydrazone,
    )

    for has_shape, namer in (
        (has_labelled_azoxy, name_labelled_azoxy),
        (has_labelled_diazene, name_labelled_diazene),
        (has_labelled_hydrazone, name_labelled_hydrazone),
        (has_isotope_diazonium_shape, name_isotope_diazonium),
        (has_isotope_acid_hydrogen_shape, name_isotope_acid_hydrogen),
        (has_isotope_hydride_shape, name_isotope_hydride),
        (has_isotope_substituted_hydride_shape, name_isotope_substituted_hydride),
        (has_isotope_one_carbon_shape, name_isotope_one_carbon),
        (has_isotope_acyl_halide_shape, name_isotope_acyl_halide),
        (has_radical_shape, hetero_radical_name),
    ):
        if has_shape(mol):
            error = None
            try:
                name = namer(mol)
            except UnsupportedStructure as caught:
                name, error = None, caught
            if name is not None:
                return name
            try:
                name = _name_labelled_substituents(mol)
            except UnsupportedStructure:
                if error is not None:
                    raise error from None
                raise
            if name is not None:
                return name
            if error is not None:
                raise error
    error = None
    try:
        name = _name_labelled_substituents(mol)
    except UnsupportedStructure as caught:
        name, error = None, caught
    if name is None and mol.HasSubstructMatch(_NONCARBON_PEROXY_ESTER):
        with contextlib.suppress(UnsupportedStructure):
            name = _labelled_namer(name_acid_derivative, oxygens=True)(mol)
    if name is None and error is not None:
        raise error
    return name


def _same_descriptors(mol, clean, index):
    """Whether every stereo element has the same CIP descriptor with and without the nuclides (P-92.1.4)."""
    return _isotopic_centres(mol, clean, index) is not None


def _isotopic_centres(mol, clean, index):
    """{atom of `clean`: descriptor} of the centres that only their nuclides make stereogenic (P-82.4.1) when every other
    stereo element has the same descriptor with and without the nuclides; None when a descriptor changes."""
    labelled, plain = Chem.Mol(mol), Chem.Mol(clean)
    assign_cip_labels(labelled)
    assign_cip_labels(plain)
    atoms = {index[a.GetIdx()]: a.GetProp("_CIPCode") for a in labelled.GetAtoms() if a.HasProp("_CIPCode") and a.GetIdx() in index}
    bonds = {
        frozenset((index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()])): b.GetProp("_CIPCode")
        for b in labelled.GetBonds()
        if b.HasProp("_CIPCode")
    }
    known = {a.GetIdx(): a.GetProp("_CIPCode") for a in plain.GetAtoms() if a.HasProp("_CIPCode")}
    if bonds != {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): b.GetProp("_CIPCode") for b in plain.GetBonds() if b.HasProp("_CIPCode")
    } or any(atoms.get(i) != code for i, code in known.items()):
        return None
    return {i: code for i, code in atoms.items() if i not in known}


def _name_labelled_substituents(mol):
    """Any parent named without its nuclides whose labelled atoms all sit in substituent groups that cite them."""
    from ._isotope_labels import split_isotopes
    from ._substituents import ISOTOPE_LABELS

    split = split_isotopes(mol)
    if split is None:
        raise UnsupportedStructure("this isotopically modified structure is not supported yet")
    clean, labels, index = split
    isotopic = {}
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms()):
        isotopic = _isotopic_centres(mol, clean, index)
        if isotopic is None:
            raise UnsupportedStructure("the nuclides change the stereodescriptors of this structure")
    context = {"labels": labels, "consumed": set(), "mol": clean}
    token = ISOTOPE_LABELS.set(context)
    cip_token = ISOTOPIC_CIP.set((clean.GetNumAtoms(), isotopic) if isotopic else None)
    try:
        name = name_polychalcogen_hydride(clean)
        if name is None or set(labels) - context["consumed"]:
            name, partial, failure = None, None, None
            for namer in (_oxoacid_name, name_chalcogen_imide, _name_mol, *_FALLBACK_NAMERS, name_carbonic_family, name_labelled_anhydride, name_polycarbonic_amide):
                context["consumed"].clear()
                try:
                    candidate = namer(clean)
                except UnsupportedStructure as error:
                    failure = error if namer is _name_mol else failure
                    continue
                if candidate is None:
                    continue
                if not set(labels) - context["consumed"]:
                    name, partial = candidate, None
                    break
                partial = partial or (candidate, set(context["consumed"]))
            if name is None and partial is None:
                raise failure or UnsupportedStructure("this isotopically modified structure is not supported yet")
            if name is None:
                name = partial[0]
                context["consumed"].clear()
                context["consumed"].update(partial[1])
    finally:
        ISOTOPIC_CIP.reset(cip_token)
        ISOTOPE_LABELS.reset(token)
    if set(labels) - context["consumed"]:
        name = _cite_parent_atom_nuclide(name, clean, labels, set(labels) - context["consumed"]) or name
    if isotopic and not _cites_every_stereo_element(mol, name):
        raise UnsupportedStructure("the stereocentre made by the nuclides is not cited by this name")
    return name if nuclides_cited(mol, name) else None


def _oxoacid_name(clean):
    """The name of a mononuclear oxoacid or its functional replacement analogue by the dedicated namers, which stand ahead of
    the phosphinic and phosphonic acid prefixes in the order of the unmodified names."""
    if _is_nucleoside_derivative(clean) or (carries_only_element_hydride_groups(clean) and has_carbonic_amide_parent(clean)):
        raise UnsupportedStructure("this structure is not named as a mononuclear oxoacid")
    for has_shape, namer in (
        (has_noncarbon_oxoacid_shape, name_noncarbon_oxoacid),
        (has_functional_replacement_oxoacid_shape, name_functional_replacement_oxoacid),
    ):
        if has_shape(clean):
            return namer(clean)
    raise UnsupportedStructure("this structure is not named as a mononuclear oxoacid")


_PARENT_ATOM_STEM = {
    5: "bor", 14: "sil", 15: "phosph", 16: "sulf|thi", 33: "ars", 34: "selen", 51: "stib", 52: "tellur",
}


_CYANATE_ANIONS = {"cyanate": 8, "thiocyanate": 16, "selenocyanate": 34, "tellurocyanate": 52}
_CYANATE_STEM = {16: "thi", 34: "selen", 52: "tellur"}


def _cite_cyanate_nuclides(name, clean, labels, open_atoms):
    """The nuclide of the chalcogen of a (thio)cyanate stands in front of its 'thio', 'seleno' or 'telluro' part, that of
    the carbon, nitrogen or the oxygen of cyanate in front of 'cyan' (P-82.2.1)."""
    from ._isotope_labels import descriptor

    chalcogen = _CYANATE_ANIONS[name]
    if any(labels[a]["H"] or not labels[a]["skeleton"] for a in open_atoms):
        return None
    at = {}
    for atom in open_atoms:
        is_chalcogen = clean.GetAtomWithIdx(atom).GetAtomicNum() == chalcogen
        offset = name.index(_CYANATE_STEM[chalcogen]) if is_chalcogen and chalcogen != 8 else name.index("cyan")
        at.setdefault(offset, {})[atom] = labels[atom]
    for offset in sorted(at, reverse=True):
        name = name[:offset] + descriptor(at[offset], {a: 1 for a in at[offset]}, True) + name[offset:]
    return name


def _cite_parent_atom_nuclide(name, clean, labels, open_atoms):
    """`name` with the nuclide of the one atom of its element, the atom the parent name is built on, in front of the stem
    that names that element (P-82.2.1, P-82.6.1.1); None when the labelled atoms are not the one atom of that element, the
    uncited nuclides then being refused by the final check; a name that does not delimit that stem exactly once is refused."""
    from ._isotope_labels import descriptor

    if name in _CYANATE_ANIONS:
        return _cite_cyanate_nuclides(name, clean, labels, open_atoms)
    if len(open_atoms) != 1:
        return None
    (atom,) = open_atoms
    z = clean.GetAtomWithIdx(atom).GetAtomicNum()
    entry = labels[atom]
    if (
        z not in _PARENT_ATOM_STEM
        or entry["H"]
        or not entry["skeleton"]
        or sum(1 for a in clean.GetAtoms() if a.GetAtomicNum() == z) != 1
    ):
        return None
    depth, outside = 0, []
    for char in name:
        depth += char in "([{"
        depth -= char in ")]}"
        outside.append(depth == 0)
    stems = [m for m in re.finditer(_PARENT_ATOM_STEM[z], name) if outside[m.start()]]
    if len(stems) != 1 or re.search(r"λ\d+-$", name[: stems[0].start()]):
        raise UnsupportedStructure("the parent stem of the labelled atom is not delimited in this name")
    at = stems[0].start()
    return name[:at] + descriptor({atom: entry}, {atom: 1}, True) + name[at:]


def _labelled_namer(namer, oxygens=False):
    """`namer` applied to the unlabelled structure, its organyl groups citing their nuclides and its acid parts named
    with theirs (P-82.2.4); None unless every nuclide appears in the name."""

    def name(mol):
        from ._isotope_labels import split_isotopes

        clean, labels, index = split_isotopes(mol)
        only_oxygen = sum(1 for atom in clean.GetAtoms() if atom.GetAtomicNum() == 8) == 1
        if any(
            clean.GetAtomWithIdx(a).GetAtomicNum() not in (6, 5, 14, 15, 16, 33, 34, 51, 52)
            and not (clean.GetAtomWithIdx(a).GetAtomicNum() == 8 and (only_oxygen or oxygens))
            and not (clean.GetAtomWithIdx(a).GetAtomicNum() in HALOGEN_PREFIXES and clean.GetAtomWithIdx(a).GetDegree() == 1)
            for a in labels
        ):
            return None
        isotopic = {}
        if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in clean.GetAtoms()):
            isotopic = _isotopic_centres(mol, clean, index)
            if isotopic is None:
                raise UnsupportedStructure("the nuclides change the stereodescriptors of this structure")
        token = ISOTOPE_LABELS.set({"labels": labels, "consumed": set(), "mol": clean})
        cip_token = ISOTOPIC_CIP.set((clean.GetNumAtoms(), isotopic) if isotopic else None)
        try:
            result = namer(clean)
        finally:
            ISOTOPIC_CIP.reset(cip_token)
            ISOTOPE_LABELS.reset(token)
        if isotopic and not _cites_every_stereo_element(mol, result):
            raise UnsupportedStructure("the stereocentre made by the nuclides is not cited by this name")
        return result if nuclides_cited(mol, result) else None

    return name


_MULTIPLIER_VALUE = {"di": 2, "bis": 2, "tri": 3, "tris": 3, "tetra": 4, "tetrakis": 4, "penta": 5, "hexa": 6}
_MULTIPLIED_GROUP = re.compile(r"(di|bis|tri|tris|tetra|tetrakis|penta|hexa)[\[{(]$")
_NUCLIDE_ITEM = re.compile(r"^(\d+)([A-Z][a-z]?)(\d*)$")


def _multiplied_factor(name, position):
    """How many times the group at `position` is repeated: the product of the multipliers of the enclosing brackets."""
    opened = []
    for index, char in enumerate(name[:position]):
        if char in "([{":
            opened.append(index)
        elif char in ")]}" and opened:
            opened.pop()
    factor = 1
    for index in opened:
        multiplier = _MULTIPLIED_GROUP.search(name[: index + 1])
        if multiplier:
            factor *= _MULTIPLIER_VALUE[multiplier.group(1)]
    return factor


_SYMMETRIC_ANHYDRIDE = re.compile(r"[^ ]+ anhydride")


def _cited_nuclides(name):
    table = Chem.GetPeriodicTable()
    cited = {}
    repeats = 2 if _SYMMETRIC_ANHYDRIDE.fullmatch(name) else 1
    for found in re.finditer(r"\(([^()]*)\)", name):
        group = found.group(1)
        if re.fullmatch(r"\d+H(,\d+H)*", group) and re.search(r"\d[a-z]?$", name[: found.start()]):
            continue
        factor = _multiplied_factor(name, found.start()) * repeats
        items = group.split(",")
        pending = []
        for item in items:
            locant, _, tail = item.rpartition("-")
            token = _NUCLIDE_ITEM.match(tail or locant)
            if token is None:
                pending.append(item)
                continue
            mass, symbol, count = int(token.group(1)), token.group(2), token.group(3)
            try:
                with rdBase.BlockLogs():
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


_ADDUCT_RATIO = re.compile(r" \((\d+(?:/\d+)+)\)$")


def _adduct_nuclides(name):
    """The nuclides cited by `name`, each component of an adduct counted as often as its ratio multiplies it (P-14.8.2)."""
    ratio = _ADDUCT_RATIO.search(name)
    if ratio is None:
        return _cited_nuclides(name)
    parts, counts = name[: ratio.start()].split("\u2014"), ratio.group(1).split("/")
    if len(parts) != len(counts):
        return _cited_nuclides(name)
    cited = {}
    for part, count in zip(parts, counts):
        for key, number in _cited_nuclides(part).items():
            cited[key] = cited.get(key, 0) + number * int(count)
    return cited


def nuclides_cited(mol, name):
    """Whether `name` cites every nuclide label of `mol`, as often as the structure carries it."""
    expected = {}
    for atom in mol.GetAtoms():
        if atom.GetIsotope():
            key = f"{atom.GetIsotope()}{atom.GetSymbol()}"
            expected[key] = expected.get(key, 0) + 1
    if not expected:
        return True
    cited = _adduct_nuclides(name)
    forced = FORCED_BRANCH_NAMES.get()
    given_names = [a.GetProp(PREFIX_PROP) for a in mol.GetAtoms() if a.HasProp(PREFIX_PROP)]
    given_names += [forced_name for forced_name, _ in forced[1].values()] if forced else []
    for given in given_names:
        if given in name:
            for key, count in _cited_nuclides(given).items():
                cited[key] = cited.get(key, 0) - count
    return all(key in name and cited.get(key, 0) == count for key, count in expected.items())


_NESTED_NAMES: dict = {}
_NESTED_NAMES_MAX = 4096


def _smiles_to_iupac_unabridged(smiles: str) -> str:
    # Acyl/substituent namers re-name the same fragment hundreds of times while ranking candidates.
    if not nested() or FORCED_BRANCH_NAMES.get():
        return _name_unabridged(smiles)
    key = (smiles, SYSTEMATIC_ACID_PROBE.get(), CITE_SUBSTITUENT_LOCANTS.get(), CATION_PARENT.get(), AMINIUM.get(), RING_CENTER.get())
    entry = _NESTED_NAMES.get(key)
    if entry is None:
        start = reason_count()
        try:
            outcome = (_name_unabridged(smiles), None)
        except Exception as error:
            outcome = (None, error)
        if len(_NESTED_NAMES) >= _NESTED_NAMES_MAX:
            _NESTED_NAMES.clear()
        entry = _NESTED_NAMES[key] = (outcome, tuple(reasons_since(start)))
    else:
        replay(entry[1])
    name, error = entry[0]
    if error is not None:
        raise error
    return name


def _name_unabridged(smiles: str) -> str:
    # The stereo, isotope and cationic-atom contexts of an enclosing name index the atoms of another molecule.
    isolated = BRANCH_STEREO.set({"atoms": {}, "bonds": {}, "used": set()}) if nested() and BRANCH_STEREO.get() else None
    unlabelled = ISOTOPE_LABELS.set(None) if nested() and ISOTOPE_LABELS.get() else None
    uncharged = CATION_RING_ATOM.set(None) if nested() and CATION_RING_ATOM.get() is not None else None
    enter()
    result = None
    try:
        result = _name_unabridged_body(smiles)
        return result
    finally:
        leave(result)
        if isolated is not None:
            BRANCH_STEREO.reset(isolated)
        if unlabelled is not None:
            ISOTOPE_LABELS.reset(unlabelled)
        if uncharged is not None:
            CATION_RING_ATOM.reset(uncharged)


def _name_unabridged_body(smiles: str) -> str:
    name = None
    lambda_token = None
    glycosyl_token = None
    try:
        parsed = _parse_smiles(smiles)
        if parsed is not None and parsed.HasSubstructMatch(_C_GLYCOSYL):
            glycosyl_token = SYSTEMATIC_FUSION.set(True)
        if parsed is not None and outermost() and Chem.MolToSmiles(parsed) in _METHYLBENZENES:
            return _METHYLBENZENES[Chem.MolToSmiles(parsed)]
        if parsed is not None and has_bare_inorganic_component(parsed):
            return mark(name_adduct(parsed, smiles_to_iupac), _NO_PIN_MIXED_ADDUCT)
        if parsed is not None and (azoxy_prefix := azoxy_prefix_name(parsed)) is not None:
            return azoxy_prefix
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
            try:
                return name_polycation(parsed)
            except UnsupportedStructure:
                pass
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
        if (
            parsed is not None
            and has_noncarbon_oxoacid_shape(parsed)
            and not _is_nucleoside_derivative(parsed)
            and not (carries_only_element_hydride_groups(parsed) and has_carbonic_amide_parent(parsed))
        ):
            name = _name_keeping_labels(parsed, name_noncarbon_oxoacid)
            if name is not None:
                return name
        if parsed is not None and has_functional_replacement_oxoacid_shape(parsed):
            name = _name_keeping_labels(parsed, name_functional_replacement_oxoacid)
            if name is not None:
                return name
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
        if parsed is not None and has_chalcogen_imide_shape(parsed):
            return name_chalcogen_imide(parsed)
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
        if parsed is not None and boron_acid_is_parent(parsed):
            with contextlib.suppress(UnsupportedStructure):
                return name_boron_acid(parsed)
            if _has_isotope_label(parsed):
                with contextlib.suppress(UnsupportedStructure):
                    name = _labelled_namer(name_boron_acid, oxygens=True)(parsed)
                    if name is not None:
                        return name
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
            carbonic_amide = (
                polycarbonic_amide_name(parsed)
                or cyanamide_name(parsed)
                or carbonohydrazonic_diamide_name(parsed)
                or hydrazine_dicarboxamide_name(parsed)
                or diazenecarbohydrazide_name(parsed)
            )
            if carbonic_amide is not None:
                return carbonic_amide
        if parsed is not None and has_dipolar_shape(parsed):
            return name_dipolar(parsed)
        if parsed is not None:
            sugar_carboxylate = sugar_carboxylate_name(parsed)
            if sugar_carboxylate is not None:
                return sugar_carboxylate
        if parsed is not None and has_sugar_lactone_shape(parsed):
            return name_sugar_lactone(parsed)
        if parsed is not None and has_sugar_lactam_shape(parsed):
            return name_sugar_lactam(parsed)
        if parsed is not None and has_sugar_ring_acid_shape(parsed):
            return name_sugar_ring_acid(parsed)
        if parsed is not None and has_ketoaldonic_chain_shape(parsed):
            return name_ketoaldonic_chain(parsed)
        if parsed is not None and has_glycoside_shape(parsed):
            return name_glycoside(parsed)
        if parsed is not None and has_sugar_alcohol_acid_shape(parsed):
            return name_sugar_alcohol_acid(parsed)
        if parsed is not None and has_substituted_chain_shape(parsed):
            return name_substituted_chain(parsed)
        if parsed is not None and sugar_acid_derivative_name(parsed) is not None:
            return sugar_acid_derivative_name(parsed)
        if parsed is not None and has_substituted_sugar_shape(parsed):
            return name_substituted_sugar(parsed)
        if parsed is not None and has_sugar_dianhydride_shape(parsed):
            return name_sugar_dianhydride(parsed)
        if parsed is not None and has_peptide_shape(parsed):
            return name_peptide(parsed)
        if parsed is not None and has_amino_acid_shape(parsed):
            return name_amino_acid(parsed)
        if parsed is not None and has_condensed_guanidine_shape(parsed):
            try:
                return name_condensed_guanidine(parsed)
            except UnsupportedStructure:
                pass
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
                name_inorganic_ester_parent,
                *((name_acid_derivative,) if in_acid_derivative_scope(parsed) else ()),
                name_boron_peroxy_ester,
                name_hetero_parent_acid,
                name_phosphorous_acid,
                *(
                    (_labelled_namer(name_polycarbonic, oxygens=True), _labelled_namer(name_carbonic_family), _labelled_namer(name_acid_derivative), _labelled_namer(name_sulfur_acid_ester, oxygens=True))
                    if _has_isotope_label(parsed)
                    else ()
                ),
            ):
                try:
                    candidate = namer(parsed)
                except UnsupportedStructure:
                    continue
                if candidate is None or not nuclides_cited(parsed, candidate):
                    continue
                # these namers do not cite stereodescriptors, so a flat name would silently drop the stereo
                if stereo_specified and namer is not name_acid_salt and not _cites_every_stereo_element(parsed, candidate):
                    continue
                return candidate
        if parsed is not None and parsed.GetRingInfo().NumRings() == 3 and sum(a.GetAtomicNum() == 8 for a in parsed.GetAtoms()) >= 8 and has_sugar_dianhydride_shape(parsed):
            return name_sugar_dianhydride(parsed)
        configured_spiro_union = None
        if parsed is not None and has_spiro_union_shape(parsed):
            try:
                name = name_spiro_union(parsed)
                if not _has_specified_stereo(parsed):
                    return name
                configured_spiro_union = name
            except UnsupportedStructure:
                pass
        if parsed is not None and has_inorganic_polyester_shape(parsed):
            return name_inorganic_polyester(parsed)
        beyond_preferred = None
        if parsed is not None:
            name = name_heteroacyclic(parsed)
            if name is not None:
                return name
            name = (
                name_chalcogen_chain_parent(parsed)
                or name_chalcogen_chain_multiplicative(parsed)
                or name_polychalcogen_hydride(parsed)
            )
            if name is not None:
                return name
            name = name_nitrogen_methylene_multiplicative(parsed)
            if name is not None:
                return name
            name = name_appendix3_skeleton(parsed) if glycosyl_token is None else None
            if name is not None:
                return name
            natural, operations = name_natural_product_ranked(parsed) if glycosyl_token is None else (None, 0)
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
            if has_isotope_steroid_shape(parsed):
                try:
                    return name_isotope_steroid(parsed)
                except UnsupportedStructure:
                    pass
            if configured_spiro_union is not None:
                return configured_spiro_union
            if not _has_ionic_centre(parsed) and has_chain_multiplicative_shape(parsed) and not has_phosphate_shape(parsed) and not has_borinic_acid_shape(parsed) and not has_boronic_acid_shape(parsed) and not has_inorganic_acid_derivative_shape(parsed):
                try:
                    name = name_polyfunctional(parsed)
                except UnsupportedStructure:
                    name = None
                if name is not None:
                    return name
            if _has_aromatic_oxo(parsed):
                try:
                    name = name_polyfunctional(parsed)
                except UnsupportedStructure:
                    name = None
                if name is not None:
                    return name
        hydro_fusion = parsed is not None and _has_hydro_fusion_system(parsed)
        try:
            name = _smiles_to_iupac_dispatch(smiles)
            if parsed is not None and _drops_anionic_charge(parsed, name):
                raise UnsupportedStructure("the negative charge of this structure is not cited by any supported name")
        except UnsupportedStructure as original:
            try:
                name = _run_fallbacks(smiles, original)
                if parsed is not None and _drops_anionic_charge(parsed, name, anywhere=True):
                    raise UnsupportedStructure("the negative charge of this structure is not cited by any supported name")
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
        if parsed is not None and _override_count(parsed) > len(_OVERRIDE_TOKENS.findall(name)):
            raise UnsupportedStructure("the allene stereochemistry of this structure is not cited by any supported name")
        visible = phantom_capped(parsed) if parsed is not None else None
        if visible is not None and _has_specified_stereo(visible):
            tokens = bool(_STEREO_TOKENS.search(name) or _STEREO_IN_RETAINED_NAME.search(name))
            if tokens and not _cites_every_double_bond(visible, name):
                raise UnsupportedStructure("a stereodefined double bond of this structure is not cited by any supported name")
            if tokens and not _cites_every_tetrahedral_center(visible, name):
                raise UnsupportedStructure("a stereocentre of this structure is not cited by any supported name")
            if tokens and name.startswith("("):
                cited = _engine_name(parsed)
                if cited is not None and cited.startswith("(") and cited != name:
                    strip = lambda text: re.sub(r"^\([^()]*\)-", "", text)
                    if strip(cited) == strip(name):
                        name = cited
            if not tokens or parsed.GetRingInfo().NumRings() == 0:
                cited = _engine_name(parsed)
                strip = (
                    (lambda text: _same_enclosing_marks(_DESCRIPTOR_GROUP.sub("", text)))
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
        if glycosyl_token is not None:
            SYSTEMATIC_FUSION.reset(glycosyl_token)


_SULFINYL_DESCRIPTOR = re.compile(r"\[\(([RS])\)-([a-z]+(?:sulfinyl|seleninyl|tellurinyl))\]([a-z]+)")


def _sulfinyl_descriptor_in_front(name):
    """P-93.3.4.1: with a single sulfinyl-type group on a simple parent the descriptor of its stereogenic atom leads the
    name, '(S)-(methanesulfinyl)ethane'."""
    match = _SULFINYL_DESCRIPTOR.fullmatch(name)
    return f"({match.group(1)})-({match.group(2)}){match.group(3)}" if match else name


_ANION_NAME_ENDING = re.compile(r"(?:ide|uide|ate|ite|ato|ido|elide)\b|(?:ide|uide|ate|ite)-|id(?:yl|ylidene|ylidyne)\b|-u?id-\d")


_ANION_WORD_ANYWHERE = re.compile(r"ato|ido|ide|ate|ite|uide|nidyl")


def _drops_anionic_charge(mol, name, anywhere=False) -> bool:
    if len(Chem.GetMolFrags(mol)) != 1 or sum(a.GetFormalCharge() for a in mol.GetAtoms()) >= 0:
        return False
    if (_ANION_WORD_ANYWHERE if anywhere else _ANION_NAME_ENDING).search(name) is None:
        return True
    if _has_carbanion(mol) and _CARBANION_WORD.search(name) is None:
        return True
    return bool(_ESTER_WORD_ENDING.search(name)) and _protonated_name(mol) == name


_CARBANION_WORD = re.compile(r"u?id(?:e(?![a-z])|-\d|e?diyl|yl|ylidene|ylidyne)")


def _has_carbanion(mol) -> bool:
    """A carbon or a nitrogen bonded twice carries its charge in an 'ide' word, never in an acid ending."""
    return any(
        a.GetFormalCharge() < 0
        and (a.GetAtomicNum() == 6 or (a.GetAtomicNum() == 7 and a.GetDegree() == 2))
        and not any(n.GetFormalCharge() > 0 for n in a.GetNeighbors())
        for a in mol.GetAtoms()
    )


_ESTER_WORD_ENDING = re.compile(r"^\S+ \S*ate$")


def _protonated_name(mol):
    """The name of `mol` with every anionic atom protonated: an ester name that equals it has lost the charge."""
    neutral = Chem.RWMol(mol)
    for atom in neutral.GetAtoms():
        if atom.GetFormalCharge() < 0:
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
            atom.SetNoImplicit(True)
    try:
        return smiles_to_iupac(Chem.MolToSmiles(neutral))
    except (UnsupportedStructure, ValueError):
        return None


_HYDRO_FUSION_RUNNING = set()
_STEREO_TOKENS = re.compile(
    r"(?<=[\d'′a-z⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])(?:[RSEZ]|seqCis|seqTrans)(?=[,)\]])|(?<=[\d'′a⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])[rs](?=[,)\]])|(?<=['′\d]\))[EZ](?=[,\]])|(?<=[RS],)[RS](?=[,)\]])|\((?:R|S|E|Z)\)|\((?:[RSEZ],)+[RSEZ]\)|\b(?:[DL]|alpha|beta)-|\((?:T|TS|SP|SS|TBPY|OC|SPY|TPR|PBPY|CU|SAPR|TPRS)-|cis-|trans-|rel-|rac-"
)


_DESCRIPTOR_GROUP = re.compile(
    r"[(\[](?:\d+[a-z]?[\u2032']*(?:\(\d+[\u2032']*\))?)?(?:[RSEZrsMPmp]|seqCis|seqTrans)"
    r"(?:,(?:\d+[a-z]?[\u2032']*(?:\(\d+[\u2032']*\))?)?(?:[RSEZrsMPmp]|seqCis|seqTrans))*[)\]]-"
)
def _same_enclosing_marks(text):
    """Dropping a descriptor group from a bracketed prefix lowers its enclosing marks (P-16.5.4.1.1)."""
    return text.translate(str.maketrans("[]{}", "()()"))


_FACE_DESCRIPTOR = re.compile(r"(?<=\d)[\u03b1\u03b2\u03be](?=[,\-)]|yl)")
_STEREO_IN_RETAINED_NAME = re.compile(r"inositol|(?:adenos|guanos|inos|xanthos|cytid|urid|thymid)in|" + _FACE_DESCRIPTOR.pattern)


_CENTER_CODE = re.compile(r"(?<=[\d'′a-z\u00b2\u00b3\u00b9\u2070-\u209f\u1d43-\u1d4a])[RSrs](?=[,)\]])|(?<=[RS],)[RS](?=[,)\]])|\((?:R|S)\)")
_CODE_WORDS = re.compile(r"bis[{\[(]|tris|tetrakis|\b[DL]-|[\u03b1\u03b2\u03be]|\brel-|\brac-|ose\b|inositol|ent-|\b(?:cis|trans)-|\((?:T|TS|SP|SS|TBPY|OC|SPY|TPR|PBPY|CU|SAPR|TPRS)-")


_OVERRIDE_TOKENS = re.compile(r"(?<=[\d'′a-z⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])[RSrsMPmp](?=[,)\]])|\([RSMP]\)")


def _override_count(mol) -> int:
    return sum(1 for a in mol.GetAtoms() if a.HasProp(STEREO_OVERRIDE))


_MULTIPLIED_GROUP_WORD = re.compile(r"\b(?:bis|tris|tetrakis)\b|(?:bis|tris|tetrakis)(?=[(\[{])|\b(?:di|tri|tetra)[(\[{]")


def _acid_repeats(name):
    """The acid of a symmetric anhydride is cited once for its two acyl groups."""
    return 2 if _SYMMETRIC_ANHYDRIDE.fullmatch(name) else 1


def _cites_every_tetrahedral_center(mol, name) -> bool:
    """A name without a D/L, alpha/beta or multiplying word needs one R or S for each specified tetrahedral centre."""
    if _CODE_WORDS.search(name):
        return True
    return len(_CENTER_CODE.findall(name)) * _acid_repeats(name) >= len(genuine_stereo_centres(mol))


def _cites_every_stereo_element(mol, name) -> bool:
    if _FACE_DESCRIPTOR.search(name):
        return True
    cited = len(_STEREO_TOKENS.findall(name))
    if _MULTIPLIED_GROUP_WORD.search(name):
        return cited > 0
    return cited * _acid_repeats(name) >= len(genuine_stereo_centres(mol)) + sum(
        1
        for s in Chem.FindPotentialStereo(mol)
        if s.specified == Chem.StereoSpecified.Specified and s.type != Chem.StereoType.Atom_Tetrahedral
    )


def _cites_every_double_bond(mol, name) -> bool:
    if _MULTIPLIED_GROUP_WORD.search(name):
        return True
    probe = Chem.Mol(mol)
    assign_cip_labels(probe)
    bonds = sum(
        1
        for s in Chem.FindPotentialStereo(probe)
        if s.type == Chem.StereoType.Bond_Double
        and s.specified == Chem.StereoSpecified.Specified
        and probe.GetBondWithIdx(s.centeredOn).HasProp("_CIPCode")
    )
    cited = len(re.findall(r"(?<=[\d'′a-z⁰¹²³⁴⁵⁶⁷⁸⁹ᵃᵇᶜᵈᵉᶠᵍʰ])(?:[EZ]|seqCis|seqTrans)(?=[,)\]])|(?<=['′\d]\))[EZ](?=[,\]])|\([EZ]\)", name))
    listed = sum(group.count(",") + 1 for group in re.findall(r"\(((?:[EZ],)+[EZ])\)", name))
    unlocated = sum(group.count(",") + 1 for group in re.findall(r"[(\[](?:[\d'\u2032]+[a-z]?[RSrs],)+((?:[EZ],)*[EZ])[)\]]", name))
    return (cited + listed + unlocated) * _acid_repeats(name) >= bonds


def _has_specified_stereo(mol) -> bool:
    return any(s.specified == Chem.StereoSpecified.Specified for s in Chem.FindPotentialStereo(mol))


def _engine_name(mol):
    return _hydro_fusion_name(mol, name_anion if any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()) else name_polyfunctional)


def _hydro_fusion_name(mol, namer=name_polyfunctional):
    key = Chem.MolToSmiles(mol)
    if key in _HYDRO_FUSION_RUNNING:
        return None
    _HYDRO_FUSION_RUNNING.add(key)
    try:
        return namer(mol)
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
            if is_hydro_fusion_system(mol, atoms) or (
                len(atoms) > len(ring) and system_type(mol, atoms) == "bridged_fused"
            ):
                return True
    return False


def _has_ionic_centre(mol) -> bool:
    """An anionic or cationic centre other than a nitro group, which a multiplicative name must keep as the parent."""
    from ._anion import anion_atoms, _is_nitro_nitrogen

    return bool(anion_atoms(mol)) or any(a.GetFormalCharge() > 0 and not _is_nitro_nitrogen(a) for a in mol.GetAtoms())


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


_FALLBACK_NAMERS = (
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
    name_hydrazine_with_heteroatom_groups,
    name_polyfunctional,
    name_guanidine,
    name_urea,
    name_thiourea,
    name_selenourea,
    name_tellurourea,
    name_carbamate_ester,
    _name_o_substituted_hydroxylamine,
    _name_hydride_onium,
    _name_chain_onium,
    name_anion,
    name_carbonic_ester_anion,
    name_ester_by_parts,
    name_chalcogen_chain_heterone,
)


_ACYCLIC_SULFONAMIDE = Chem.MolFromSmarts("[#6]S(=O)(=O)!@[NX3]")
_CARBON_CHALCOGENOL = Chem.MolFromSmarts("[#6][SX2H1,SeX2H1,TeX2H1]")


def _has_acid_senior_to_sulfonyl(mol):
    """A sulfonic, sulfinic, selenonic, seleninic, telluronic or tellurinic acid, or a carboxylic acid beside a sulfonamide."""
    return any(
        acid_group_at(mol, atom.GetIdx()) is not None for atom in mol.GetAtoms() if atom.GetAtomicNum() in (16, 34, 52)
    ) or (find_carboxylic_acid_carbon(mol) is not None and mol.HasSubstructMatch(_ACYCLIC_SULFONAMIDE))


def _name_standalone(mol):
    """An acid or amide is named as itself, whatever parent the molecule it was cut out of is being named around."""
    tokens = [(AMINIUM, AMINIUM.set(False)), (RING_CENTER, RING_CENTER.set(False)), (FORCED_PRINCIPAL, FORCED_PRINCIPAL.set(None))]
    try:
        name = name_polyfunctional(mol)
    finally:
        for var, token in reversed(tokens):
            var.reset(token)
    if name is None:
        raise UnsupportedStructure("this acid or amide has no parent hydride that the substitutive engine builds")
    return name


def _run_fallbacks(smiles, original):
    mol = _parse_smiles(smiles)
    key = Chem.MolToSmiles(mol)
    if key in _FALLBACKS_RUNNING:
        raise original
    _FALLBACKS_RUNNING.add(key)
    try:
        if has_ion_adduct_shape(mol):
            try:
                return mark(name_ion_adduct(mol, smiles_to_iupac), _NO_PIN_ION_ADDUCT)
            except UnsupportedStructure:
                pass
        skeletal_without_suffixes = None
        for skeletal in (name_skeletal_chain, name_hetero_macrocycle):
            try:
                name = skeletal(mol)
            except UnsupportedStructure:
                continue
            if name is not None and skeletal is name_skeletal_chain and amine_suffix_capacity(mol, adjacency(mol)):
                skeletal_without_suffixes = name
            elif name is not None:
                return name
        if skeletal_without_suffixes is None:
            try:
                skeletal_without_suffixes = name_heteroacyclic(mol, express_amine_suffixes=False)
            except UnsupportedStructure:
                pass
        for fallback in _FALLBACK_NAMERS:
            try:
                name = fallback(mol)
            except UnsupportedStructure:
                continue
            if name is not None and not _drops_anionic_charge(mol, name, anywhere=True):
                return name
        name = _name_via_fallbacks(mol)
        if name is not None:
            return name
        if skeletal_without_suffixes is not None:
            return mark(skeletal_without_suffixes, _NO_PIN_SKELETAL_AMINE)
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
        raise InvalidSmilesError(f"invalid SMILES: {smiles!r}")
    return _name_mol(mol)


def _name_via_fallbacks(mol):
    for contract in (contract_phosphanyl_groups, contract_hetero_groups_candidates, contract_nitro_groups):
        try:
            contracted = contract(mol)
            candidates = contracted if isinstance(contracted, list) else [contracted]
            names = []
            for candidate in candidates:
                if candidate is None:
                    continue
                try:
                    # a mapped placeholder survives SMILES, so the dedicated modules in front of the engine see it too
                    if contract is contract_nitro_groups:
                        names.append(smiles_to_iupac(Chem.MolToSmiles(candidate)))
                    else:
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
    aromatic_rings = [ring for ring in info.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    benzene = aromatic_rings[0] if len(aromatic_rings) == 1 else None
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
    modified_fullerene = name_fullerene_derivative(mol)
    if modified_fullerene is not None:
        return modified_fullerene
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
    if has_halogen_amide_shape(mol):
        try:
            return name_halogen_amide(mol)
        except UnsupportedStructure:
            pass

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
    if (
        has_metal_pair_shape(mol)
        and not has_chain_polyester_shape(mol)
        and not _has_senior_principal_group(mol)
        and not has_radical_shape(mol)
        and not has_carbonic_acid_shape(mol)
    ):
        # P-28.4.1: two identical rings joined by a bond are an assembly even when the joined atoms are Group 13-15.
        if mol.GetRingInfo().NumRings() == 2 and (ring_assembly_core := find_ring_assembly_core(mol)) is not None:
            try:
                return name_ring_assembly(mol, ring_assembly_core)
            except UnsupportedStructure:
                pass
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
        name = name_adduct(mol, smiles_to_iupac)
        return mark(name, _NO_PIN_MIXED_ADDUCT) if has_inorganic_component(mol) else name

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
    boron_acid = has_boronic_acid_shape(mol) or has_borinic_acid_shape(mol)
    acid_centre = (
        boron_acid
        or has_phosphorus_acid_shape(mol)
        or has_phosphanone_shape(mol)
        or has_phosphanimine_shape(mol)
    )
    if 3 <= mol.GetRingInfo().NumRings() <= 6 and not acid_centre:
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
    if mol.GetRingInfo().NumRings() == 2 and not acid_centre:
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
            name = namer(mol)
            if name is None:
                raise UnsupportedStructure("no supported name cites every nuclide of this structure")
            return name

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
        if any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
            for has_shape, namer in ((has_phosphate_shape, name_phosphate), (has_carbonic_acid_shape, name_carbonic_acid)):
                if has_shape(mol):
                    return namer(mol)
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

    # A urea or guanidine (carbonic acid amides, class 7b) outranks the amide of a phosphorus acid and any element hydride
    # group on its nitrogen (P-41); the phosphorus namers below take every molecule that contains phosphorus.
    if carries_only_element_hydride_groups(mol) and has_carbonic_amide_parent(mol):
        for has_shape, namer in ((has_urea_shape, name_urea), (has_guanidine_shape, name_guanidine)):
            if has_shape(mol):
                with contextlib.suppress(UnsupportedStructure):
                    return namer(mol)

    # An all-phosphorus, multi-atom skeleton (diphosphane, triphosphane,
    # ...) has no carbon at all.
    if has_phosphane_chain_shape(mol) and mol.GetNumAtoms() > 1:
        return name_phosphane_chain(mol)
    for has_shape, namer in (
        # the oxo acid groups carry hydroxyl oxygens `_phosphanone.py` does not expect, so they are routed first
        (has_phosphorus_acid_shape, name_phosphorus_acid),
        (has_phosphorus_acid_derivative_shape, name_phosphorus_acid_derivative),
        (has_azonic_acid_shape, name_azonic_acid),
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
        (has_phosphate_shape, name_phosphate),
        # Thio analogues of the phosphate, phosphonate and phosphinate esters carry sulfur on the phosphorus.
        (has_phosphorus_thioester_shape, name_phosphorus_thioester),
        (has_hydroxylamine_acid_shape, name_hydroxylamine_acid),
        (has_formazan_shape, name_formazan),
        # A phosphite ester (P-67.1.3.2's P(OR)3, no P=O) has three P-O-R
        # ester oxygens.
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
        # central carbon with two -O-R/-OH oxygens that `_ether.py` (which
        # rejects an oxygen bonded to more than one heavy atom outright)
        # does not expect -- routed here for the same early-ester reasoning
        # as sulfate/nitrate above.
        (has_carbonic_acid_shape, name_carbonic_acid),
        # A phosphine oxide (P-68.3.2.3.1's '-phosphanone' suffix, R-P(=O)<)
        # has its own phosphorus-bonded oxygen.
        (has_phosphanone_shape, name_phosphanone),
        (has_phosphanimine_shape, name_phosphanimine),
        # Thiophosphoric acid (P-67.1.2's own preselected infix-modified
        # oxoacid name) has a phosphorus with 4 substituents (=S plus three
        # -OH).
        (has_noncarbon_oxoacid_shape, name_noncarbon_oxoacid),
        (has_functional_replacement_oxoacid_shape, name_functional_replacement_oxoacid),
    ):
        if has_shape(mol):
            return namer(mol)

    # Phosphinine/phosphinoline/isophosphinoline (P-25's own P-ring
    # counterparts to pyridine/quinoline/isoquinoline) have P in a ring.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        if has_hetero_monocyclic_name(mol):
            return name_hetero_monocyclic(mol)
        if has_hetero_monocyclic_substituent_name(mol):
            return name_hetero_monocyclic_substituent(mol)
    for has_shape, namer in (
        (has_polyphosphane_shape, name_polyphosphane),
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
        # A boron atom (P-68's borane substitutive nomenclature) must be routed
        # here: none of the branches below recognize boron at all.
        (has_simple_borane_shape, name_simple_borane),
    ):
        if has_shape(mol):
            return namer(mol)

    # A Group 1/2 metal (Li/Na/K/Mg/Ca, P-69.3) is named additively; none of
    # the branches below recognize any of these elements at all.
    if has_group1_2_organometallic_shape(mol):
        return mark(name_group1_2_organometallic(mol), _NO_PIN_ORGANOMETALLIC)

    # None of the ring modules below understand a closed cage (P-27.2).
    if has_fullerene_name(mol):
        return name_fullerene(mol)
    if has_substituted_fullerene_cage(mol):
        require_defined_fullerene_numbering(mol, fullerene_cage_atoms(mol))
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
        (has_fused_ring_yne_shape, name_fused_ring_yne),
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
        # An oxo acid of sulfur, selenium or tellurium (P-65.3.1) has oxygens on its own atom that the ether/ester/
        # alcohol checks below do not understand, and a carboxylic acid outranks it and a sulfonamide (P-41,
        # Table 3.3), so the engine names these acids.
        (_has_acid_senior_to_sulfonyl, _name_standalone),
        # A sulfonyl group on a plain saturated ring nitrogen (e.g.
        # 1-methylsulfonylpiperidine) is a ring-as-parent construction
        # the sulfonyl checks below do not know -- must be routed here first (narrower than `has_ring_amine_shape` alone, so it
        # doesn't also preempt `_hidden_amide_ketone.py`'s unrelated
        # acyl-on-ring-nitrogen shape further down).
        (has_ring_amine_sulfonyl_shape, name_ring_amine),
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
        # oxygen-branch fallback and never even reach the nitrogen branch.
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
    # swallowed by that branch's unconditional fallback and
    # never reach the nitrogen branch further down. Gated on an oxygen
    # actually being present so this doesn't preempt a diazo group's own
    # C=N bond (checked later, in the nitrogen branch) -- diazo has no
    # oxygen at all, so it's unaffected either way, but this keeps the
    # check's intent explicit. A non-oxime imine (no oxygen) is still
    # correctly caught by has_simple_imine_shape's other call site further
    # down.
    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()) and has_simple_imine_shape(mol):
        return _name_standalone(mol)
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
        # before `has_amide_shape` further down, whose `name_amide`
        # unconditionally rejects any non-benzene ring.
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
        # An alkoxide anion (R-O(-), P-72.2.2.2.2) has a formal-charge -1
        # oxygen none of the neutral-oxygen checks below (or the
        # branch's own fallback) expect, so it must be routed first in this branch.
        if has_alkoxide_shape(mol):
            return name_alkoxide(mol)
        # Hydroxylamine (H2N-OH, P-68.3.1.1.1) and its O-substituted
        # derivatives (H2N-O-R) have a nitrogen the ether/carbonyl/alcohol
        # checks below don't expect at all, so it must be routed before all
        # of them. N-substituted forms (R-NH-OH) don't match this shape
        # (see _hydroxylamine.py's module docstring) and fall through to
        # the amine fallback below instead.
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
        # match any check above or the branch's own fallback below, so
        # it must be routed here too, right alongside its R-O-O-R' cousin.
        if has_hydroperoxide_shape(mol):
            return name_hydroperoxide(mol)
        # A carbamate's carbon (R-O-C(=O)-NH2) is simultaneously
        # ester-shaped (carbonyl + a second, carbon-bonded oxygen) and
        # amide-shaped (carbonyl + a primary -NH2) on the very same carbon,
        # so it must be routed before both has_ester_shape and
        # has_amide_shape below, either of which would otherwise misname it.
        if has_carbamate_shape(mol):
            return name_carbamate(mol)
        # A sulfonate anion's sulfur (R-SO3-, P-72.2.2.2.1.1) bears an
        # anionic oxygen (formal charge -1, no H) that every other module
        # here -- including `_sulfonic_acid.py`'s own neutral -SO3H check
        # further down -- would reject outright, so it must be routed
        # before those.
        if has_sulfonate_shape(mol):
            return name_sulfonate(mol)
        # A D-aldohexopyranose ring (a hemiacetal, no actual ring C=O
        # anywhere) is named as a sugar, not by the generic parent choice.
        if has_cyclic_aldopyranose_shape(mol):
            return name_cyclic_aldopyranose(mol)
        # Same reasoning for a furanose (5-membered) ring.
        if has_cyclic_aldofuranose_shape(mol):
            return name_cyclic_aldofuranose(mol)
        # A D-2-ketohexopyranose ring (e.g. fructopyranose) is a hemiketal,
        # not an actual ring ketone either.
        if has_cyclic_ketohexopyranose_shape(mol):
            return name_cyclic_ketohexopyranose(mol)
        if has_cyclic_ketohexofuranose_shape(mol):
            return name_cyclic_ketohexofuranose(mol)
        # Proline's pyrrolidine ring (P-103.1.2) is a plain secondary
        # cyclic amine with its own retained name.
        if has_proline_shape(mol):
            return name_proline(mol)
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
        if find_carboxylic_acid_carbon(mol) is not None:
            # P-41/Table 3.3: 'oic acid' outranks 'al', so a carboxylic acid
            # that also carries one or more aldehydes names the acid as the
            # suffix and demotes each aldehyde to an 'oxo' prefix.
            if has_aldehyde_carboxylic_acid_shape(mol):
                return name_aldehyde_carboxylic_acid(mol)
            # P-66.1.1.3.3: 'oic acid' also outranks 'amide', so a
            # carboxylic acid that also carries a coexisting primary amide
            # (an amic acid) names the acid as the suffix and demotes the
            # amide to 'amino'+'oxo' prefixes.
            if has_carboxylic_acid_amide_shape(mol):
                return name_carboxylic_acid_amide(mol)
            # P-103.1.1.1: histidine's imidazol-4-ylmethyl side chain gets its
            # own retained name + L/D descriptor.
            if has_histidine_shape(mol):
                return name_histidine(mol)
            return _name_standalone(mol)
        # A carbon bearing both a carbonyl oxygen and a primary-amide
        # nitrogen (-CONH2) is an amide (junior only to the acid/ester
        # suffixes above in Table 3.3) and must be routed before the
        # aldehyde/ketone checks below.
        # Urea (H2N-C(=O)-NH2) has a carbonyl carbon with two qualifying
        # nitrogens, the same shape `has_amide_shape` looks for -- routing
        # it first here lets `_urea.py` claim it before the amide branch.
        if has_urea_shape(mol):
            return name_urea(mol)
        if has_guanidine_shape(mol):
            return name_guanidine(mol)
        if has_amide_shape(mol):
            return _name_standalone(mol)
        # An imide (-C(=O)-NH-C(=O)-) is named as an N-acyl amide; it must be routed before the generic carbonyl
        # fallback below, which would otherwise misname each acyl carbon as a plain aldehyde.
        if has_imide_shape(mol):
            return _name_standalone(mol)
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
            # A plain open-chain aldose (P-102.5.1/P-102.5.2.2, #1039 M1
            # step 1) is itself aldehyde-shaped (a terminal -CHO), so it
            # must be routed before the generic aldehyde fallback below,
            # which would otherwise name it as a plain polyhydroxy-
            # aldehyde substitutive name instead of its carbohydrate name.
            if has_open_chain_aldose_shape(mol):
                return name_open_chain_aldose(mol)
            # A plain open-chain 2-ketose (P-102.5.2.1/P-102.5.2.2, #1039
            # M1 step 2) is itself ketone-shaped (a non-terminal
            # carbonyl), so it must be routed before the generic ketone
            # fallback below, which would otherwise name it as a plain
            # polyhydroxy-ketone substitutive name instead of its
            # carbohydrate name.
            if has_open_chain_2_ketose_shape(mol):
                return name_open_chain_2_ketose(mol)
            return _name_standalone(mol)
        # One of the 9 retained-name inositol stereoisomers (P-104.2.1) is
        # itself a plain cyclohexane-1,2,3,4,5,6-hexol shape, so it must be
        # routed before the generic systematic fallback below, which would
        # otherwise name it with CIP descriptors instead of its retained
        # name -- mirrors the open-chain-aldose/2-ketose carbohydrate-name-
        # before-generic-fallback pattern above.
        if has_inositol_shape(mol):
            return name_inositol(mol)
        return _name_standalone(mol)
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        # An aminide anion (-NH(-), P-72.2.2.2.3) has a formal-charge -1
        # nitrogen none of the neutral-nitrogen checks below (or
        # the amine fallback's own fallback) expect, so it must be routed first.
        if has_aminide_shape(mol):
            return name_aminide(mol)
        # An azide group (-N3, P-61.7) has no oxygen and three nitrogens
        # the amine fallback doesn't recognize at all, so it must be routed here
        # before the amine fallback for the same reason as the nitrile/imine checks
        # below.
        if has_azide_shape(mol):
            return name_azide(mol)
        # An isocyanide group (-NC, P-61.9) similarly has no oxygen and a
        # nitrogen shape (triple-bonded to a carbon) the amine fallback doesn't
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
        # routed before the amine fallback for the same reason as the checks above.
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
                return _name_standalone(mol)
            except UnsupportedStructure:
                pass
        # A hydrazine skeleton (H2N-NH2, P-68.3.1.2.1) has its own two
        # skeletal nitrogens (single-bonded, not double-bonded like
        # diazene above) with no carbon parent chain at all -- must be
        # routed before the amine fallback for the same reason as the checks
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
        # An amidine carbon (-C(=NH)NH2, P-66.4.1.1) has a C=N double bond that would otherwise look imine-shaped to
        # the check below; guanidine has the same imino+amino nitrogen shape and is claimed first.
        if has_guanidine_shape(mol):
            return name_guanidine(mol)
        if has_amidine_shape(mol):
            return _name_standalone(mol)
        # A plain (non-oxime) imine (C=N, P-62.3) has no oxygen, so it
        # reaches this branch alongside plain amines -- an oxime (which
        # does have an oxygen) is already caught by the has_simple_imine_shape
        # check much earlier, before the "any O atom" branch.
        if has_simple_imine_shape(mol):
            return _name_standalone(mol)
        # A plain saturated monocyclic amine whose sole ring nitrogen
        # carries one substituent (e.g. 1-methylpiperidine, or the N,N/
        # N,S piperazine/thiomorpholine siblings -- the N,O morpholine
        # sibling has an oxygen and is instead routed in the oxygen-gated
        # branch above) has the ring itself as the parent hydride, unlike
        # every other shape reaching the amine fallback below (which always
        # treats an acyclic chain as the parent) -- the amine fallback itself
        # explicitly defers this shape, so it must be routed here first.
        # These oxygen-free variants never collide with
        # `_hidden_amide_ketone.py`'s acyl-on-ring-nitrogen path, which
        # requires the substituent's own carbonyl oxygen and is routed
        # separately in the oxygen-gated branch above.
        if has_ring_amine_shape(mol):
            return name_ring_amine(mol)
        # The reverse shape: the ring nitrogen itself is plain (no
        # substituent), but a ring carbon bears an exocyclic primary
        # amine (e.g. piperidin-4-amine) -- the amine fallback also explicitly
        # defers any ring bearing a nitrogen heteroatom, so this must be
        # routed here first too, same reasoning as `has_ring_amine_shape`
        # just above.
        if has_hetero_ring_amine_shape(mol):
            return name_hetero_ring_amine(mol)
        return _name_standalone(mol)
    for has_shape, namer in (
        (has_chalcogen_peroxol_shape, name_hydroperoxide),
        (has_disulfide_shape, name_disulfide),  # A disulfide (R-S-S-R') has two sulfurs -- it would otherwise # look thiol-shaped to the check below (that check just looks for # the presence of any sulfur atom) -- must be routed here first.
        (has_sulfide_shape, name_sulfide),  # A plain -S- sulfide (P-63.2.1) has no suffix, so it must be routed # here before has_thiol_shape below: _thiol.py's validation rejects # a degree-2 sulfur outright (not a monovalent -SH), so a sulfide # would otherwise raise the wrong error there instead of being named.
        (has_diselenide_shape, name_diselenide),  # A diselenide (R-Se-Se-R') has two seleniums -- it would # otherwise look selenol-shaped to the check below (that check # just looks for the presence of any selenium atom) -- must be # routed here first.
        (has_selenide_shape, name_selenide),  # A plain -Se- selenide (P-63.2.1) has no suffix, so it must be # routed here before has_selenol_shape below for the same reason # as has_sulfide_shape above (that check doesn't look at degree # at all, so a degree-2 selenide would otherwise raise the wrong # error inside `name_selenol`'s degree-1 validation).
        (has_ditelluride_shape, name_ditelluride),  # A ditelluride (R-Te-Te-R') has two telluriums -- it would # otherwise look tellurol-shaped to the check below (that check # just looks for the presence of any tellurium atom) -- must be # routed here first.
        (has_telluride_shape, name_telluride),  # A plain -Te- telluride (P-63.2.1) has no suffix, so it must be # routed here before has_tellurol_shape below for the same reason # as has_selenide_shape above (that check doesn't look at degree # at all, so a degree-2 telluride would otherwise raise the wrong # error inside `name_tellurol`'s degree-1 validation).
    ):
        if has_shape(mol):
            return namer(mol)

    if mol.HasSubstructMatch(_CARBON_CHALCOGENOL):
        return _name_standalone(mol)

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
