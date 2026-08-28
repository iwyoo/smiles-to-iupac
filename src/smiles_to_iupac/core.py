from rdkit import Chem

from ._acyclic import name_acyclic_alkane
from ._acyl_halide import has_acyl_halide_shape, name_acyl_halide
from ._anhydride import has_anhydride_shape, name_anhydride
from ._carbamate import has_carbamate_shape, name_carbamate
from ._alcohol import name_alcohol
from ._aldehyde import name_aldehyde
from ._carboxylic_acid_amine import has_carboxylic_acid_amine_shape, name_carboxylic_acid_amine
from ._aldehyde_carboxylic_acid import (
    has_aldehyde_carboxylic_acid_shape,
    name_aldehyde_carboxylic_acid,
)
from ._aldehyde_ketone import has_aldehyde_ketone_shape, name_aldehyde_ketone
from ._acetal import has_acetal_shape, name_acetal
from ._amide import has_amide_shape, name_amide
from ._amidine import has_amidine_shape, name_amidine
from ._hydrazide import has_hydrazide_shape, name_hydrazide
from ._imide import has_imide_shape, name_imide
from ._amine import name_amine
from ._ammonium import has_ammonium_shape, name_ammonium
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
from ._carboxylic_acid import has_carboxylic_acid_shape, name_carboxylic_acid
from ._common import UnsupportedStructure, non_single_bonds
from ._cyclic import name_cycloalkane
from ._cyclic_unsaturated import find_cyclic_unsaturated_core, name_cyclic_unsaturated
from ._dihydro_aromatic import find_dihydronaphthalene_core, name_dihydronaphthalene
from ._ester import has_ester_shape, name_ester
from ._ether import has_ether_shape, name_ether
from ._fullerene import has_fullerene_name, name_fullerene
from ._androstane import has_androstane_name, name_androstane
from ._gonane import has_gonane_name, name_gonane
from ._heteroaromatic_fused import has_retained_heteroaromatic_fused_name, name_retained_heteroaromatic_fused
from ._hetero_monocyclic import has_hetero_monocyclic_name, name_hetero_monocyclic
from ._hydroxylamine import has_hydroxylamine_shape, name_hydroxylamine
from ._imine import has_simple_imine_shape, name_imine
from ._isotope import has_isotope_shape, name_isotope
from ._two_component_heterocycle_fusion import (
    has_two_component_heterocycle_fusion_name,
    name_two_component_heterocycle_fusion,
)
from ._ketone import name_ketone
from ._ketone_amide import has_ketone_amide_shape, name_ketone_amide
from ._ketone_ester import has_ketone_ester_shape, name_ketone_ester
from ._nitrile import has_nitrile_shape, name_nitrile
from ._azide import has_azide_shape, name_azide
from ._diazene import has_diazene_shape, name_diazene
from ._azine import has_azine_shape, name_azine
from ._hydrazine import has_hydrazine_shape, name_hydrazine
from ._hydrazone import has_hydrazone_shape, name_hydrazone
from ._diazo import has_diazo_shape, name_diazo
from ._isocyanate import has_isocyanate_shape, name_isocyanate
from ._isocyanide import has_isocyanide_shape, name_isocyanide
from ._nitro import has_nitro_shape, name_nitro
from ._nitroso import has_nitroso_shape, name_nitroso
from ._polycyclic import find_polycyclic_core, name_polycycloalkane
from ._cyclophane import has_cyclophane_name, name_cyclophane
from ._phosphane import has_simple_phosphane_shape, name_simple_phosphane
from ._phosphane_chain import has_phosphane_chain_shape, name_phosphane_chain
from ._peri_fused_aromatic import has_retained_peri_fused_name, name_retained_peri_fused
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
from ._sulfide import has_sulfide_shape, name_sulfide
from ._sulfinic_acid import has_sulfinic_acid_shape, name_sulfinic_acid
from ._sulfonic_acid import has_sulfonic_acid_shape, name_sulfonic_acid
from ._sulfone import has_sulfone_shape, name_sulfone
from ._sulfoxide import has_sulfoxide_shape, name_sulfoxide
from ._selenol import has_selenol_shape, name_selenol
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

    # A sulfonic acid (-SO3H, P-65.3.1) has three oxygens on its own sulfur,
    # so it must be routed here before the plain "any O atom" branch below --
    # none of the ether/ester/carboxylic-acid/aldehyde/ketone/alcohol checks
    # in that branch understand a sulfur-centered oxygen cluster at all.
    if has_sulfonic_acid_shape(mol):
        return name_sulfonic_acid(mol)
    # A sulfinic acid (-SO2H, P-65.3.1) has two oxygens on its own sulfur --
    # the same reasoning as sulfonic acid above -- so it too must be routed
    # before the plain "any O atom" branch.
    if has_sulfinic_acid_shape(mol):
        return name_sulfinic_acid(mol)
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

    if any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms()):
        # Hydroxylamine (H2N-OH, P-68.3.1.1.1) and its O-substituted
        # derivatives (H2N-O-R) have a nitrogen the ether/carbonyl/alcohol
        # checks below don't expect at all, so it must be routed before all
        # of them. N-substituted forms (R-NH-OH) don't match this shape
        # (see _hydroxylamine.py's module docstring) and fall through to
        # name_amine below instead.
        if has_hydroxylamine_shape(mol):
            return name_hydroxylamine(mol)
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
        # A carbon bearing both a carbonyl oxygen and a second, carbon-bonded
        # oxygen is an ester (-COO-), which must be routed before the
        # carboxylic-acid/aldehyde/ketone checks below: its carbonyl half
        # would otherwise look aldehyde/ketone-shaped, and (for a rejected,
        # out-of-scope case) its non-carbonyl oxygen would never satisfy the
        # carboxylic acid module's hydroxyl (O-H) requirement anyway.
        if has_ester_shape(mol):
            # P-41/Table 3.3: 'oate' outranks 'one', so an ester whose acyl
            # chain also carries one or more ketones names the ester as the
            # suffix and demotes each ketone to an 'oxo' prefix instead of
            # `_ester.py`'s own "coexisting oxygen" rejection.
            if has_ketone_ester_shape(mol):
                return name_ketone_ester(mol)
            return name_ester(mol)
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
        # A nitrile nitrogen (-C#N, P-66.5) has no oxygen, so it reaches this
        # branch alongside plain amines; it must be routed here before
        # name_amine, which doesn't recognize a triple-bonded nitrogen at all.
        if has_nitrile_shape(mol):
            return name_nitrile(mol)
        # An amidine carbon (-C(=NH)NH2, P-66.4.1.1) has a C=N double bond
        # that would otherwise look imine-shaped to the check below (a
        # different, unrelated interpretation of the same C=N bond) -- must
        # be routed here first.
        if has_amidine_shape(mol):
            return name_amidine(mol)
        # A plain (non-oxime) imine (C=N, P-62.3) has no oxygen, so it
        # reaches this branch alongside plain amines -- an oxime (which
        # does have an oxygen) is already caught by the has_simple_imine_shape
        # check much earlier, before the "any O atom" branch.
        if has_simple_imine_shape(mol):
            return name_imine(mol)
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
    if has_selenol_shape(mol):
        # A selenol (-SeH, P-63.1.1) is the next chalcogen analogue after
        # a thiol -- has neither O, N, nor S, so it only reaches this
        # branch once all three are ruled out above.
        return name_selenol(mol)

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
