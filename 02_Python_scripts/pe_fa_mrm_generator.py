#PE-FA MRM Generator 2.0
import pandas as pd


fatty_acids = ['14:0', '14:1', '16:0', '16:1', '18:0', '18:1', '18:2', '18:3', '20:0', '20:1', '20:2', '20:3', '20:4', '20:5', '22:4', '22:5', '22:6']      # FAs pulled from Lipidyzer


def mrm_output(fa1):

    carbon = 12.000000
    hydrogen = 1.007825
    oxygen = 15.994915
    phosphorus = 30.973762
    nitrogen = 14.003074

    fa_carbon_no = []
    fa_desat_no = []
    fa_desat_no_mass = []

    # converting string to int for calculations below
    for i in fa1:
        fa_carbon_no.append(int(i[:2]))
        fa_desat_no.append(int(i[-1]))
        fa_desat_no_mass.append(int(i[-1]))
        
    

    pe_carbon = []
    pe_desat = []
    pe_desat_mass = []
    

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            pe_carbon += [fa_carbon_no[0] + i]

        fa_carbon_no.remove(fa_carbon_no[0])
            

    # to get the total fa desat number for label
    while fa_desat_no != []:
        for i in fa_desat_no:
            pe_desat += [fa_desat_no[0] + i]

        fa_desat_no.remove(fa_desat_no[0])
            

    # to get the total fa desat number for mass calculation
    while fa_desat_no_mass != []:
        for i in fa_desat_no_mass:
            pe_desat_mass += [fa_desat_no_mass[0] + i]

        fa_desat_no_mass.remove(fa_desat_no_mass[0])
            


    # to output the pe id (i.e. PE-(fa C no.):(fa desat no.))
    pe_id = []
    pe_M1_mass = []

    id_counter = len(pe_carbon)
    id_number = len(pe_carbon)

    
    while id_counter != 0:
        pe_id += ["PE-" + str(pe_carbon[id_number - id_counter]) + ":" + str(pe_desat[id_number - id_counter])]
        
        pe_M1_mass += [(carbon * (pe_carbon[id_number - id_counter] + 5)) + (hydrogen * ((pe_carbon[id_number - id_counter] * 2) - (pe_desat_mass[id_number - id_counter] * 2) + 10)) + (oxygen * 8) + nitrogen + phosphorus - hydrogen]

        id_counter -= 1
    


    # to output the FA combinations for each PE (this is the same as the original TAG/FA combinator)
    pe = []
        
    while fa1 != []:
        for i in fa1:
            pe += [fa1[0] + '/' + i]

        fa1.remove(fa1[0])




    # calculating each possible MRM transition (Ac adducts)
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []

    for i in pe:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_combo.clear()
    


    pe_M1_mass_counter = len(pe_M1_mass)
    pe_M1_mass_number = len(pe_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pe_M1_mass_counter -= 1
    
    pe_M1_mass_counter = len(pe_M1_mass)
    pe_M1_mass_number = len(pe_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pe_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2
    all_fatty_acids = fatty_acid_1 + fatty_acid_2
    all_pe_id = pe_id + pe_id
    all_pe_M1_mass = pe_M1_mass + pe_M1_mass
    all_pe_FA_label = []

    pe_FA_label_counter = len(all_pe_id)
    pe_FA_label_number = len(all_pe_id)

    for i in all_pe_id:
        all_pe_FA_label += [i + '-FA' + all_fatty_acids[pe_FA_label_number - pe_FA_label_counter]]

        pe_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'PE ID': pe_id, 'FA composition': pe, 'M1 mass (-H)': pe_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'PE ID': all_pe_id, 'FA loss': all_fatty_acids, 'PE-FA label': all_pe_FA_label, 'M1 mass (-H)': all_pe_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)