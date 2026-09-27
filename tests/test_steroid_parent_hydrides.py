from smiles_to_iupac import smiles_to_iupac


def test_all_seven_parent_hydrides():
    # Each canonical SMILES cross-checked against its PubChem CID's
    # ConnectivitySMILES (see module docstring for the CIDs and the
    # side-chain/methyl removal chain confirming each skeleton reduces to
    # the previous one).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("C1CCC2C(C1)CCC3C2CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"
    assert smiles_to_iupac("CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "ergostane"


def test_natural_configuration_parent_hydrides_resolve_to_retained_name():
    # Each skeleton's own PubChem CID's current IsomericSMILES (CID 6857523
    # gonane, 6857536 androstane, 5460658 estrane, 439513 pregnane, 6857459
    # cholane, 6857534 cholestane, 6857535 ergostane) -- the natural
    # ring-fusion (and, where present, side-chain) configuration that
    # essentially every real-world instance of these compounds actually
    # has -- now resolves to the bare retained name, no descriptor needed
    # (Rule 3S-2.2/2.3/2.4 defines each retained name as this one fixed
    # configuration).
    assert (
        smiles_to_iupac("C1CCC2CC[C@H]3[C@@H]4CCC[C@H]4CC[C@@H]3[C@H]2C1") == "gonane"
    )
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@]4([C@H]3CC2)C")
        == "androstane"
    )
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@H]4[C@H]3CC2") == "estrane"
    )
    assert (
        smiles_to_iupac("CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C")
        == "pregnane"
    )
    assert (
        smiles_to_iupac("CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C")
        == "cholane"
    )
    assert (
        smiles_to_iupac("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "cholestane"
    )
    assert (
        smiles_to_iupac("C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "ergostane"
    )


def test_c24_epimer_skeletons_resolve_to_retained_name_without_collision():
    # campestane (CID 6857532) and stigmastane (CID 6857438) are the C24
    # epimers of ergostane and poriferastane (CID 6857528) respectively --
    # same constitution, opposite configuration at that one carbon. Each
    # resolves to its own distinct retained name, not its epimer's or
    # ergostane's, confirming the epic's collision concern doesn't
    # materialize (see module docstring).
    assert (
        smiles_to_iupac("C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "campestane"
    )
    assert (
        smiles_to_iupac("CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C")
        == "poriferastane"
    )
    assert (
        smiles_to_iupac("CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C")
        == "stigmastane"
    )


def test_non_natural_stereo_specified_parent_hydrides_still_fall_through_to_von_baeyer():
    # A partially-specified stereo input that doesn't match either the
    # stereo-free entries or any natural-configuration entry above still
    # falls through to the general von Baeyer engine, since it produces a
    # different canonical SMILES from every lookup entry. (This function
    # used to also assert this for a second, differently-partially-
    # specified androstane SMILES -- that one turned out to actually BE
    # androstane's own fully-C5-specified natural form once fetched fresh
    # from its dedicated "5alpha-androstane" CID, see
    # test_fully_c5_specified_parent_hydrides_resolve_to_retained_name;
    # the assertion was simply wrong, not a regression.)
    assert (
        smiles_to_iupac("C1CCC2[C@H](C1)CCC3C2CCC4C3CCC4")
        == "tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_extra_methyl_beyond_each_parent_hydride_falls_through_to_von_baeyer():
    # A methyl (or, for gonane, both angular methyls) beyond what each
    # skeleton's own definition allows changes the canonical SMILES, so it
    # correctly falls through to the general von Baeyer engine as a plain
    # substituted hydrocarbon instead of being misrecognized as the next
    # skeleton up.
    assert (
        smiles_to_iupac("CC12CCCC1(C)CCC1C2CCC2C1CCCC2")
        == "11,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C3CC(C)C4CCCCC4(C3CC2)C")
        == "2,8,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C1CC(C)C3CCCCC3C1CC2")
        == "8,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCC1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-ethyl-2,12,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(pentan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)C(C)CC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-(4,5-dimethylhexan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_fully_c5_specified_parent_hydrides_resolve_to_retained_name():
    # 8 of the 10 registered natural-configuration entries left the C5
    # ring-fusion stereocenter unspecified (a PubChem data-depiction quirk
    # on those particular CIDs -- see module docstring); each skeleton's
    # own dedicated "5alpha-<name>" PubChem CID fully specifies it instead
    # (94144, 6857525, 6857465, 2723895, 164641, 6857533, 188005, 5491914).
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "androstane"
    )
    assert (
        smiles_to_iupac("C1CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H]2CC[C@H]4[C@H]3CCC4") == "gonane"
    )
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@H]4[C@H]3CC2") == "estrane"
    )
    assert (
        smiles_to_iupac("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C")
        == "cholestane"
    )
    assert (
        smiles_to_iupac("C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C")
        == "ergostane"
    )
    assert (
        smiles_to_iupac("C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C")
        == "campestane"
    )
    assert (
        smiles_to_iupac("CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C")
        == "poriferastane"
    )
    assert (
        smiles_to_iupac("CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C")
        == "stigmastane"
    )


def test_5beta_androstane_resolves_via_c5_cip_computation():
    # etiocholane (5-beta-androstane), PubChem CID 6857462 -- differs from
    # the natural (5-alpha) entry only at the C5 ring-fusion stereocenter,
    # recognized via the general per-locant alpha/beta mechanism (not a
    # hardcoded row -- see module docstring/`_alpha_beta_stereo_prefix`).
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC[C@@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "5beta-androstane"
    )


def test_androstane_other_diastereomer_gets_its_own_alpha_beta_citation():
    # A ring-fusion stereocenter *other* than C5 differing from the
    # natural configuration is a genuinely different diastereomer, not
    # simply "the 5-beta epimer" -- but the general per-locant mechanism
    # (unlike the old single-case hardcoded one) still names it correctly
    # rather than falling through to the general von Baeyer engine.
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "5alpha,8alpha-androstane"
    )


def test_blue_book_worked_example_pregnane_multi_locant_inversion():
    # The Blue Book's own P-101.2.6.1.1 worked example: pregnane with C9
    # and C10 inverted from the natural configuration, plus C5 (always
    # cited once specified) -- "5beta,9beta,10alpha-pregnane"
    # (`tmp/bluebook/P10.txt` ~245), reproduced exactly here by inverting
    # locants 5/9/10 of the natural-configuration entry.
    assert (
        smiles_to_iupac("CC[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@@]4(C)[C@@H]3CC[C@]12C")
        == "5beta,9beta,10alpha-pregnane"
    )
