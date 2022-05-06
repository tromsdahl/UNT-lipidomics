#PC-FA MRM Generator 2.0
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
        
    

    pc_carbon = []
    pc_desat = []
    pc_desat_mass = []
    

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            pc_carbon += [fa_carbon_no[0] + i]

        fa_carbon_no.remove(fa_carbon_no[0])
            

    # to get the total fa desat number for label
    while fa_desat_no != []:
        for i in fa_desat_no:
            pc_desat += [fa_desat_no[0] + i]

        fa_desat_no.remove(fa_desat_no[0])
            

    # to get the total fa desat number for mass calculation
    while fa_desat_no_mass != []:
        for i in fa_desat_no_mass:
            pc_desat_mass += [fa_desat_no_mass[0] + i]

        fa_desat_no_mass.remove(fa_desat_no_mass[0])
            


    # to output the pc id (i.e. PC-(fa C no.):(fa desat no.))
    pc_id = []
    pc_M1_mass = []

    id_counter = len(pc_carbon)
    id_number = len(pc_carbon)

    
    while id_counter != 0:
        pc_id += ["PC-" + str(pc_carbon[id_number - id_counter]) + ":" + str(pc_desat[id_number - id_counter])]
        
        pc_M1_mass += [(carbon * (pc_carbon[id_number - id_counter] + 8)) + (hydrogen * ((pc_carbon[id_number - id_counter] * 2) - (pc_desat_mass[id_number - id_counter] * 2) + 15)) + (oxygen * 8) + nitrogen + phosphorus + ((oxygen * 2) + (carbon * 2) + (hydrogen * 4))]

        id_counter -= 1
    


    # to output the FA combinations for each PC (this is the same as the original TAG/FA combinator)
    pc = []
        
    while fa1 != []:
        for i in fa1:
            pc += [fa1[0] + '/' + i]

        fa1.remove(fa1[0])




    # calculating each possible MRM transition (Ac adducts)
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []

    for i in pc:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_combo.clear()
    


    pc_M1_mass_counter = len(pc_M1_mass)
    pc_M1_mass_number = len(pc_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pc_M1_mass_counter -= 1
    
    pc_M1_mass_counter = len(pc_M1_mass)
    pc_M1_mass_number = len(pc_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pc_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2
    all_fatty_acids = fatty_acid_1 + fatty_acid_2
    all_pc_id = pc_id + pc_id
    all_pc_M1_mass = pc_M1_mass + pc_M1_mass
    all_pc_FA_label = []

    pc_FA_label_counter = len(all_pc_id)
    pc_FA_label_number = len(all_pc_id)

    for i in all_pc_id:
        all_pc_FA_label += [i + '-FA' + all_fatty_acids[pc_FA_label_number - pc_FA_label_counter]]

        pc_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'PC ID': pc_id, 'FA composition': pc, 'M1 mass (+AcO-)': pc_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'PC ID': all_pc_id, 'FA loss': all_fatty_acids, 'PC-FA label': all_pc_FA_label, 'M1 mass (+AcO-)': all_pc_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)