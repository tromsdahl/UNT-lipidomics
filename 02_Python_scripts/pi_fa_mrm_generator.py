#PI-FA MRM Generator 2.0
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
        
    

    pi_carbon = []
    pi_desat = []
    pi_desat_mass = []
    

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            pi_carbon += [fa_carbon_no[0] + i]

        fa_carbon_no.remove(fa_carbon_no[0])
            

    # to get the total fa desat number for label
    while fa_desat_no != []:
        for i in fa_desat_no:
            pi_desat += [fa_desat_no[0] + i]

        fa_desat_no.remove(fa_desat_no[0])
            

    # to get the total fa desat number for mass calculation
    while fa_desat_no_mass != []:
        for i in fa_desat_no_mass:
            pi_desat_mass += [fa_desat_no_mass[0] + i]

        fa_desat_no_mass.remove(fa_desat_no_mass[0])
            


    # to output the pi id (i.e. PI-(fa C no.):(fa desat no.))
    pi_id = []
    pi_M1_mass = []

    id_counter = len(pi_carbon)
    id_number = len(pi_carbon)

    
    while id_counter != 0:
        pi_id += ["PI-" + str(pi_carbon[id_number - id_counter]) + ":" + str(pi_desat[id_number - id_counter])]
        
        pi_M1_mass += [(carbon * (pi_carbon[id_number - id_counter] + 9)) + (hydrogen * ((pi_carbon[id_number - id_counter] * 2) - (pi_desat_mass[id_number - id_counter] * 2) + 15)) + (oxygen * 13) + phosphorus - hydrogen]

        id_counter -= 1
    


    # to output the FA combinations for each PI (this is the same as the original TAG/FA combinator)
    pi = []
        
    while fa1 != []:
        for i in fa1:
            pi += [fa1[0] + '/' + i]

        fa1.remove(fa1[0])




    # calculating each possible MRM transition (Ac adducts)
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []

    for i in pi:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_combo.clear()
    


    pi_M1_mass_counter = len(pi_M1_mass)
    pi_M1_mass_number = len(pi_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pi_M1_mass_counter -= 1
    
    pi_M1_mass_counter = len(pi_M1_mass)
    pi_M1_mass_number = len(pi_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pi_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2
    all_fatty_acids = fatty_acid_1 + fatty_acid_2
    all_pi_id = pi_id + pi_id
    all_pi_M1_mass = pi_M1_mass + pi_M1_mass
    all_pi_FA_label = []

    pi_FA_label_counter = len(all_pi_id)
    pi_FA_label_number = len(all_pi_id)

    for i in all_pi_id:
        all_pi_FA_label += [i + '-FA' + all_fatty_acids[pi_FA_label_number - pi_FA_label_counter]]

        pi_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'PI ID': pi_id, 'FA composition': pi, 'M1 mass (-H)': pi_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'PI ID': all_pi_id, 'FA loss': all_fatty_acids, 'PI-FA label': all_pi_FA_label, 'M1 mass (-H)': all_pi_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)