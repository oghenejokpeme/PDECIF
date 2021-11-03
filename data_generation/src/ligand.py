import os
from rdkit import Chem
from rdkit.Chem import GraphDescriptors
from rdkit.Chem import EState
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import MolSurf
from rdkit.Chem import Crippen
from rdkit.Chem import QED
from rdkit.Chem import Fragments
import math

HEADER = ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1',
              'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 
              'Chi4n', 'Chi4v', 'Kappa1', 'Kappa2', 'Kappa3', 'HallKierAlpha',
              'EstateVSA1', 'EstateVSA10', 'EstateVSA11', 'EstateVSA2', 
              'EstateVSA3', 'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 
              'EstateVSA7', 'EstateVSA8', 'EstateVSA9', 'Estate10', 'Estate4',
              'Estate5', 'Estate8', 'Estate9', 'ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 
              'FpDensityMorgan3', 'HeavyAtomMolWt', 'MaxAbsPartialCharge',
              'MaxPartialCharge', 'MinAbsPartialCharge', 'MinPartialCharge', 
              'MolWt', 'NumRadicalElectrons', 'NumValenceElectrons',
              'FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount',
              'NumAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds',
              'NumAromaticRings', 'NumSaturatedRings', 'NumAliphaticRings',
              'NumAromaticHeteroCycles', 'NumAromaticCarbocycles', 
              'NumSaturatedHeteroCycles', 'NumSaturatedCarbocycles',
              'NumAliphaticHeteroCycles', 'NumAliphaticCarbocycles',
              'RingCount', 'MaxAbsEstateIndex',
              'MaxEstateIndex', 'MinAbsEstateIndex', 'MinEstateIndex',
              'LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 
              'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 
              'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8',
              'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 
              'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9',
              'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 
              'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 
              'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'Crippen_MolLogP', 'Crippen_MolMR', 'QED',
              'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 
              'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 
              'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 
              'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 
              'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 
              'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine',
              'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether',
              'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide',
              'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone',
              'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 
              'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 
              'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 
              'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 
              'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 
              'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea'
              ]

def graph_properties(m):
    gbala  = GraphDescriptors.BalabanJ(m)
    gbert  = GraphDescriptors.BertzCT(m)
    gchi0  = GraphDescriptors.Chi0(m)
    gchi0n = GraphDescriptors.Chi0n(m)
    gchi0v = GraphDescriptors.Chi0v(m)
    gchi1  = GraphDescriptors.Chi1(m)
    gchi1n = GraphDescriptors.Chi1n(m)
    gchi1v = GraphDescriptors.Chi1v(m)
    gchi2n = GraphDescriptors.Chi2n(m)
    gchi2v = GraphDescriptors.Chi2v(m)
    gchi3n = GraphDescriptors.Chi3n(m)
    gchi3v = GraphDescriptors.Chi3v(m)
    gchi4n = GraphDescriptors.Chi4n(m)
    gchi4v = GraphDescriptors.Chi4v(m)
    gkap1  = GraphDescriptors.Kappa1(m)
    gkap2  = GraphDescriptors.Kappa2(m)
    gkap3  = GraphDescriptors.Kappa3(m)
    ghalk  = GraphDescriptors.HallKierAlpha(m)

    header = ['BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1',
              'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 
              'Chi4n', 'Chi4v', 'Kappa1', 'Kappa2', 'Kappa3', 'HallKierAlpha']

    features = (gbala, gbert, gchi0, gchi0n, gchi0v, gchi1, gchi1n, 
                gchi1v, gchi2n, gchi2v, gchi3n, gchi3v, gchi4n, 
                gchi4v, gkap1, gkap2, gkap3, ghalk)
    
    return header, features

def estatevsa_properties(m):
    evas1  = EState.EState_VSA.EState_VSA1(m)
    evas10 = EState.EState_VSA.EState_VSA10(m)
    evas11 = EState.EState_VSA.EState_VSA11(m)
    evas2  = EState.EState_VSA.EState_VSA2(m)
    evas3  = EState.EState_VSA.EState_VSA3(m)
    evas4  = EState.EState_VSA.EState_VSA4(m)
    evas5  = EState.EState_VSA.EState_VSA5(m)
    evas6  = EState.EState_VSA.EState_VSA6(m)
    evas7  = EState.EState_VSA.EState_VSA7(m)
    evas8  = EState.EState_VSA.EState_VSA8(m)
    evas9  = EState.EState_VSA.EState_VSA9(m)
    est10  = EState.EState_VSA.VSA_EState10(m)
    est4   = EState.EState_VSA.VSA_EState4(m)
    est5   = EState.EState_VSA.VSA_EState5(m)
    est8   = EState.EState_VSA.VSA_EState8(m)
    est9   = EState.EState_VSA.VSA_EState9(m)

    header = ['EstateVSA1', 'EstateVSA10', 'EstateVSA11', 'EstateVSA2', 
              'EstateVSA3', 'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 
              'EstateVSA7', 'EstateVSA8', 'EstateVSA9', 'Estate10', 'Estate4',
              'Estate5', 'Estate8', 'Estate9']
    
    features = (evas1, evas10, evas11, evas2, evas3, evas4, evas5, evas6, 
                evas7, evas8, evas9, est10, est4, est5, est8, est9)
    
    return header, features

def general_properties(m):
    emolwt  = Descriptors.ExactMolWt(m)
    morgan1 = Descriptors.FpDensityMorgan1(m)
    morgan2 = Descriptors.FpDensityMorgan2(m)
    morgan3 = Descriptors.FpDensityMorgan3(m)
    atomwt  = Descriptors.HeavyAtomMolWt(m)
    apchrge = Descriptors.MaxAbsPartialCharge(m)
    mpchrge = Descriptors.MaxPartialCharge(m)
    miapchrge = Descriptors.MinAbsPartialCharge(m)
    mipchrge = Descriptors.MinPartialCharge(m)
    molwt   = Descriptors.MolWt(m)
    radelec = Descriptors.NumRadicalElectrons(m)
    valelec = Descriptors.NumValenceElectrons(m)

    header = ['ExactMolWt', 'FpDensityMorgan1', 'FpDensityMorgan2', 
              'FpDensityMorgan3', 'HeavyAtomMolWt', 'MaxAbsPartialCharge',
              'MaxPartialCharge', 'MinAbsPartialCharge', 'MinPartialCharge', 
              'MolWt', 'NumRadicalElectrons', 'NumValenceElectrons']
    
    features = (emolwt, morgan1, morgan2, morgan3, atomwt, apchrge,
                mpchrge, miapchrge, mipchrge, molwt, radelec, valelec)
    
    return header, features 

def lipinksi_properties(m):
    frac   = Lipinski.FractionCSP3(m)
    heavac = Lipinski.HeavyAtomCount(m)
    nhoc   = Lipinski.NHOHCount(m)
    noc    = Lipinski.NOCount(m)
    hacc   = Lipinski.NumHAcceptors(m)
    hdon   = Lipinski.NumHDonors(m)
    hetatm = Lipinski.NumHeteroatoms(m)
    rot    = Lipinski.NumRotatableBonds(m)
    ar     = Lipinski.NumAromaticRings(m)
    sat    = Lipinski.NumSaturatedRings(m)
    aliph  = Lipinski.NumAliphaticRings(m)
    arhet  = Lipinski.NumAromaticHeterocycles(m)
    arcar  = Lipinski.NumAromaticCarbocycles(m)
    sathet = Lipinski.NumSaturatedHeterocycles(m)
    satcar =Lipinski.NumSaturatedCarbocycles(m)
    alihet = Lipinski.NumAliphaticHeterocycles(m)
    alicar = Lipinski.NumAliphaticCarbocycles(m)           
    rings  = Lipinski.RingCount(m)

    header = ['FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount',
              'NumAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds',
              'NumAromaticRings', 'NumSaturatedRings', 'NumAliphaticRings',
              'NumAromaticHeteroCycles', 'NumAromaticCarbocycles', 
              'NumSaturatedHeteroCycles', 'NumSaturatedCarbocycles',
              'NumAliphaticHeteroCycles', 'NumAliphaticCarbocycles',
              'RingCount']

    features = (frac, heavac, nhoc, noc, hacc, hdon, hetatm, rot, ar, sat,
                aliph, arhet, arcar, sathet, satcar, alihet, alicar, rings)
    
    return header, features

def estate_properties(m):
    mxaei  = EState.MaxAbsEStateIndex(m)
    mai    = EState.MaxEStateIndex(m)
    mnaei  = EState.MinAbsEStateIndex(m)
    mnei   = EState.MinEStateIndex(m)

    header = ['MaxAbsEstateIndex', 'MaxEstateIndex', 
              'MinAbsEstateIndex', 'MinEstateIndex']
    
    features = (mxaei, mai, mnaei, mnei)

    return header, features

def labute_properties(m):
    labuteasa  = MolSurf.LabuteASA(m)
    peoe_vsa1  = MolSurf.PEOE_VSA1(m)
    peoe_vsa10 = MolSurf.PEOE_VSA10(m)
    peoe_vsa11 = MolSurf.PEOE_VSA11(m)
    peoe_vsa12 = MolSurf.PEOE_VSA12(m)
    peoe_vsa13 = MolSurf.PEOE_VSA13(m)
    peoe_vsa14 = MolSurf.PEOE_VSA14(m)
    peoe_vsa2 = MolSurf.PEOE_VSA2(m)
    peoe_vsa3 = MolSurf.PEOE_VSA3(m)
    peoe_vsa4 = MolSurf.PEOE_VSA4(m)
    peoe_vsa5 = MolSurf.PEOE_VSA5(m)
    peoe_vsa6 = MolSurf.PEOE_VSA6(m)
    peoe_vsa7 = MolSurf.PEOE_VSA7(m)
    peoe_vsa8 = MolSurf.PEOE_VSA8(m)
    peoe_vsa9 = MolSurf.PEOE_VSA9(m)
    smr_vsa1  = MolSurf.SMR_VSA1(m)
    smr_vsa10 = MolSurf.SMR_VSA10(m)
    smr_vsa2  = MolSurf.SMR_VSA2(m)
    smr_vsa3  = MolSurf.SMR_VSA3(m)
    smr_vsa4  = MolSurf.SMR_VSA4(m)
    smr_vsa5  = MolSurf.SMR_VSA5(m)
    smr_vsa6  = MolSurf.SMR_VSA6(m)
    smr_vsa7  = MolSurf.SMR_VSA7(m)
    smr_vsa8  = MolSurf.SMR_VSA8(m)
    smr_vsa9  = MolSurf.SMR_VSA9(m)
    slogp_vsa1  = MolSurf.SlogP_VSA1(m)
    slogp_vsa10 = MolSurf.SlogP_VSA10(m)
    slogp_vsa11 = MolSurf.SlogP_VSA11(m)
    slogp_vsa12 = MolSurf.SlogP_VSA12(m)
    slogp_vsa2  = MolSurf.SlogP_VSA2(m)
    slogp_vsa3  = MolSurf.SlogP_VSA3(m)
    slogp_vsa4  = MolSurf.SlogP_VSA4(m)
    slogp_vsa5  = MolSurf.SlogP_VSA5(m)
    slogp_vsa6  = MolSurf.SlogP_VSA6(m)
    slogp_vsa7  = MolSurf.SlogP_VSA7(m)
    slogp_vsa8  = MolSurf.SlogP_VSA8(m)
    slogp_vsa9  = MolSurf.SlogP_VSA9(m)
    tpsa =  MolSurf.TPSA(m)

    header = ['LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 
              'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 
              'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8',
              'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 
              'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9',
              'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 
              'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 
              'SlogP_VSA8', 'SlogP_VSA9', 'TPSA']
    
    features = (labuteasa, peoe_vsa1, peoe_vsa10, peoe_vsa11, peoe_vsa12, peoe_vsa13,
                peoe_vsa14, peoe_vsa2, peoe_vsa3, peoe_vsa4, peoe_vsa5, peoe_vsa6,
                peoe_vsa7, peoe_vsa8, peoe_vsa9, smr_vsa1, smr_vsa10, smr_vsa2,
                smr_vsa3, smr_vsa4, smr_vsa5, smr_vsa6, smr_vsa7, smr_vsa8, smr_vsa9,
                slogp_vsa1, slogp_vsa10, slogp_vsa11, slogp_vsa12, slogp_vsa2, slogp_vsa3,
                slogp_vsa4, slogp_vsa5, slogp_vsa6, slogp_vsa7, slogp_vsa8, slogp_vsa9,
                tpsa)
    
    return header, features

def crippen_properties(m):
    mollogp = Crippen.MolLogP(m)
    molmr   = Crippen.MolMR(m)

    header = ['Crippen_MolLogP', 'Crippen_MolMR']
    features = (mollogp, molmr)

    return header, features

def qed_propeties(m):
    qed = QED.qed(m)

    return ['QED'], (qed,)

def fragment_properties(m):
    fr_Al_COO          = Fragments.fr_Al_COO(m)
    fr_Al_OH           = Fragments.fr_Al_OH(m)
    fr_Al_OH_noTert    = Fragments.fr_Al_OH_noTert(m)
    fr_ArN             = Fragments.fr_ArN(m)
    fr_Ar_COO          = Fragments.fr_Ar_COO(m)
    fr_Ar_N            = Fragments.fr_Ar_N(m)
    fr_Ar_NH           = Fragments.fr_Ar_NH(m)
    fr_Ar_OH           = Fragments.fr_Ar_OH(m)
    fr_COO             = Fragments.fr_COO(m)
    fr_COO2            = Fragments.fr_COO2(m)
    fr_C_O             = Fragments.fr_C_O(m)
    fr_C_O_noCOO       = Fragments.fr_C_O_noCOO(m)
    fr_C_S             = Fragments.fr_C_S(m)
    fr_HOCCN           = Fragments.fr_HOCCN(m)
    fr_Imine           = Fragments.fr_Imine(m)
    fr_NH0             = Fragments.fr_NH0(m)
    fr_NH1             = Fragments.fr_NH1(m)
    fr_NH2             = Fragments.fr_NH2(m)
    fr_N_O             = Fragments.fr_N_O(m)
    fr_Ndealkylation1  = Fragments.fr_Ndealkylation1(m)
    fr_Ndealkylation2  = Fragments.fr_Ndealkylation2(m)
    fr_Nhpyrrole       = Fragments.fr_Nhpyrrole(m)
    fr_SH              = Fragments.fr_SH(m)
    fr_aldehyde        = Fragments.fr_aldehyde(m)
    fr_alkyl_carbamate = Fragments.fr_alkyl_carbamate(m)
    fr_alkyl_halide    = Fragments.fr_alkyl_halide(m)
    fr_allylic_oxid    = Fragments.fr_allylic_oxid(m)
    fr_amide           = Fragments.fr_amide(m)
    fr_amidine         = Fragments.fr_amidine(m)
    fr_aniline         = Fragments.fr_aniline(m)
    fr_aryl_methyl     = Fragments.fr_aryl_methyl(m)
    fr_azide           = Fragments.fr_azide(m)
    fr_azo             = Fragments.fr_azo(m)
    fr_barbitur        = Fragments.fr_barbitur(m)
    fr_benzene         = Fragments.fr_benzene(m)
    fr_benzodiazepine  = Fragments.fr_benzodiazepine(m)
    fr_bicyclic        = Fragments.fr_bicyclic(m)
    fr_diazo           = Fragments.fr_diazo(m)
    fr_dihydropyridine = Fragments.fr_dihydropyridine(m)
    fr_epoxide         = Fragments.fr_epoxide(m)
    fr_ester           = Fragments.fr_ester(m)
    fr_ether           = Fragments.fr_ether(m)
    fr_furan           = Fragments.fr_furan(m)
    fr_guanido         = Fragments.fr_guanido(m)
    fr_halogen         = Fragments.fr_halogen(m)
    fr_hdrzine         = Fragments.fr_hdrzine(m)
    fr_hdrzone         = Fragments.fr_hdrzone(m)
    fr_imidazole       = Fragments.fr_imidazole(m)
    fr_imide           = Fragments.fr_imide(m)
    fr_isocyan         = Fragments.fr_isocyan(m)
    fr_isothiocyan     = Fragments.fr_isothiocyan(m)
    fr_ketone          = Fragments.fr_ketone(m)
    fr_ketone_Topliss  = Fragments.fr_ketone_Topliss(m)
    fr_lactam          = Fragments.fr_lactam(m)
    fr_lactone         = Fragments.fr_lactone(m)
    fr_methoxy         = Fragments.fr_methoxy(m)
    fr_morpholine      = Fragments.fr_morpholine(m)
    fr_nitrile         = Fragments.fr_nitrile(m)
    fr_nitro           = Fragments.fr_nitro(m)
    fr_nitro_arom      = Fragments.fr_nitro_arom(m)
    fr_nitro_arom_nonortho = Fragments.fr_nitro_arom_nonortho(m)
    fr_nitroso         = Fragments.fr_nitroso(m)
    fr_oxazole         = Fragments.fr_oxazole(m)
    fr_oxime           = Fragments.fr_oxime(m)
    fr_para_hydroxylation = Fragments.fr_para_hydroxylation(m)
    fr_phenol          = Fragments.fr_phenol(m)
    fr_phenol_noOrthoHbond = Fragments.fr_phenol_noOrthoHbond(m)
    fr_phos_acid       = Fragments.fr_phos_acid(m)
    fr_phos_ester      = Fragments.fr_phos_ester(m)
    fr_piperdine       = Fragments.fr_piperdine(m)
    fr_piperzine       = Fragments.fr_piperzine(m)
    fr_priamide        = Fragments.fr_priamide(m)
    fr_prisulfonamd    = Fragments.fr_prisulfonamd(m)
    fr_pyridine        = Fragments.fr_pyridine(m)
    fr_quatN           = Fragments.fr_quatN(m)
    fr_sulfide         = Fragments.fr_sulfide(m)
    fr_sulfonamd       = Fragments.fr_sulfonamd(m)
    fr_sulfone         = Fragments.fr_sulfone(m)
    fr_term_acetylene  = Fragments.fr_term_acetylene(m)
    fr_tetrazole       = Fragments.fr_tetrazole(m)
    fr_thiazole        = Fragments.fr_thiazole(m)
    fr_thiocyan        = Fragments.fr_thiocyan(m)
    fr_thiophene       = Fragments.fr_thiophene(m)
    fr_unbrch_alkane   = Fragments.fr_unbrch_alkane(m)
    fr_urea            = Fragments.fr_urea(m)

    header = ['fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 
              'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 
              'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 
              'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 
              'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 
              'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine',
              'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether',
              'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide',
              'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone',
              'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 
              'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 
              'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 
              'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 
              'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 
              'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea']

    features = (fr_Al_COO, fr_Al_OH, fr_Al_OH_noTert, fr_ArN, fr_Ar_COO, fr_Ar_N, fr_Ar_NH, 
                fr_Ar_OH, fr_COO, fr_COO2, fr_C_O, fr_C_O_noCOO, fr_C_S, fr_HOCCN, fr_Imine, 
                fr_NH0, fr_NH1, fr_NH2, fr_N_O, fr_Ndealkylation1, fr_Ndealkylation2, 
                fr_Nhpyrrole, fr_SH, fr_aldehyde, fr_alkyl_carbamate, fr_alkyl_halide, 
                fr_allylic_oxid, fr_amide, fr_amidine, fr_aniline, fr_aryl_methyl, fr_azide, 
                fr_azo, fr_barbitur, fr_benzene, fr_benzodiazepine, fr_bicyclic, fr_diazo, 
                fr_dihydropyridine, fr_epoxide, fr_ester, fr_ether, fr_furan, fr_guanido, 
                fr_halogen, fr_hdrzine, fr_hdrzone, fr_imidazole, fr_imide, fr_isocyan, 
                fr_isothiocyan, fr_ketone, fr_ketone_Topliss, fr_lactam, fr_lactone, fr_methoxy, 
                fr_morpholine, fr_nitrile, fr_nitro, fr_nitro_arom, fr_nitro_arom_nonortho, 
                fr_nitroso, fr_oxazole, fr_oxime, fr_para_hydroxylation, fr_phenol, 
                fr_phenol_noOrthoHbond, fr_phos_acid, fr_phos_ester, fr_piperdine, fr_piperzine, 
                fr_priamide, fr_prisulfonamd, fr_pyridine, fr_quatN, fr_sulfide, fr_sulfonamd, 
                fr_sulfone, fr_term_acetylene, fr_tetrazole, fr_thiazole, fr_thiocyan, fr_thiophene,
                fr_unbrch_alkane, fr_urea)
    
    return header, features

def read_molecule_ids(subset):
    filepath = '../input/ligands/' + subset + '/'
    filenames = [name.split('.')[0] for name in os.listdir(filepath)]
    filenames.sort()

    return filenames

def read_ligand(subset, molecule_id):
    filepath = '../input/ligands/' + subset + '/' + molecule_id + '.mol'
    return Chem.MolFromMolFile(filepath)

def add_features(func, feature_map, ligand):
    features, values = func(ligand)
    for feature, value in zip(features, values):
        feature_map[feature] = value

def get_ligand_features(ligand):
    features = {}
    add_features(graph_properties, features, ligand)
    add_features(estatevsa_properties, features, ligand)
    add_features(general_properties, features, ligand)
    add_features(lipinksi_properties, features, ligand)
    add_features(estate_properties, features, ligand)
    add_features(labute_properties, features, ligand)
    add_features(crippen_properties, features, ligand)
    add_features(qed_propeties, features, ligand)
    add_features(fragment_properties, features, ligand)

    return features

def build_ligand_dataset(subset):
    molecule_ids = read_molecule_ids(subset)

    with open('../output/ligand_' + subset.lower() + '.csv', 'w') as f:
        header = 'Complex,' + ','.join(HEADER) + '\n'
        f.write(header)
        for molecule_id in molecule_ids:
            ligand = read_ligand(subset, molecule_id)
            ligand_features = get_ligand_features(ligand)
            values = [str(round(ligand_features[feature], 3)) for feature in HEADER]
            entry = molecule_id + ',' + ','.join(values)
            f.write(entry + '\n')

def find_missing(subset):
    molecule_ids = read_molecule_ids(subset)

    for molecule_id in molecule_ids:
        ligand = read_ligand(subset, molecule_id)
        ligand_features = get_ligand_features(ligand)

        print(molecule_id)
        for feature in HEADER:
            value = ligand_features[feature]
            if math.isnan(value):
                print(feature, value, type(value))
        print()

def main():
    subsets = ['CASF-07', 'CASF-13', 'CASF-16', 'CASF-19']

    for subset in subsets:
        print(subset)
        build_ligand_dataset(subset)

if __name__ == '__main__':
    main()