# TAG-FA MRM generator 2.0
import pandas as pd

fatty_acids = ['14:0', '14:1', '16:0', '16:1', '18:0', '18:1', '18:2', '18:3', '20:0', '20:1', '20:2', '20:3', '20:4', '20:5', '22:4', '22:5', '22:6']      # FAs from lipidyzer


def mrm_output(fa1):

    carbon = 12.000000
    hydrogen = 1.007825
    oxygen = 15.994915
    ammonium = 18.03382555

    fa_carbon_no = []
    fa_desat_no = []

    # converting string to int for calculations below
    for i in fa1:
        fa_carbon_no.append(int(i[:2]))
        fa_desat_no.append(int(i[-1]))
        
    

    tag_carbon_mass = []
    tag_desat_mass = []
    final_tag_carbon = []
    final_tag_desat = []

    fa_carbon_no_2 = []
    fa_desat_no_2 = []
    
    m = len(fa_carbon_no)
    n = len(fa_desat_no)

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            tag_carbon_mass += [fa_carbon_no[0] + i]

        fa_carbon_no_2.append(fa_carbon_no[0])
        fa_carbon_no.remove(fa_carbon_no[0])
            
        
    for i in fa_carbon_no_2:
        for x in tag_carbon_mass:
            final_tag_carbon += [i + x]
    
        del tag_carbon_mass[:m]

        m -= 1

    # to get the total fa desat number
    while fa_desat_no != []:
        for i in fa_desat_no:
            tag_desat_mass += [fa_desat_no[0] + i]

        fa_desat_no_2.append(fa_desat_no[0])
        fa_desat_no.remove(fa_desat_no[0])
            
        
    for i in fa_desat_no_2:
        for x in tag_desat_mass:
            final_tag_desat += [i + x]
    
        del tag_desat_mass[:n]

        n -= 1



    # to output the TAG id (i.e. TAG-(fa C no.):(fa desat no.))
    TAG_id = []
    TAG_M1_mass = []

    id_counter = len(final_tag_carbon)
    id_number = len(final_tag_carbon)

    
    while id_counter != 0:
        TAG_id += ["TAG-" + str(final_tag_carbon[id_number - id_counter]) + ":" + str(final_tag_desat[id_number - id_counter])]
        
        TAG_M1_mass += [(carbon * (final_tag_carbon[id_number - id_counter] + 3)) + (hydrogen * ((final_tag_carbon[id_number - id_counter] * 2) - (final_tag_desat[id_number - id_counter] * 2) + 2)) + (oxygen * 6) + ammonium]

        id_counter -= 1
    

    
    # to output the FA combinations for each TAG (this is the same as the original TAG/FA combinator) 
    fa1_2 = []
    fa2 = []
    tag = []

    t = len(fa1)

       
    while fa1 != []:
        for i in fa1:
            fa2 += [fa1[0] + '/' + i]

        fa1_2.append(fa1[0])
        fa1.remove(fa1[0])
            
        
    for i in fa1_2:
        for x in fa2:
            tag += [i + "/" + x]
    
        del fa2[:t]

        t -= 1



    # calculating each possible MRM transition
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_3 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []
    M3_fa3 = []

    for i in tag:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_3 += [fatty_acid_combo[2]]
        fatty_acid_combo.clear()


    TAG_M1_mass_counter = len(TAG_M1_mass)
    TAG_M1_mass_number = len(TAG_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [(TAG_M1_mass[TAG_M1_mass_number - TAG_M1_mass_counter]) - ((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen + ammonium)]

        TAG_M1_mass_counter -= 1
    
    TAG_M1_mass_counter = len(TAG_M1_mass)
    TAG_M1_mass_number = len(TAG_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [(TAG_M1_mass[TAG_M1_mass_number - TAG_M1_mass_counter]) - ((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen + ammonium)]

        TAG_M1_mass_counter -= 1

    TAG_M1_mass_counter = len(TAG_M1_mass)
    TAG_M1_mass_number = len(TAG_M1_mass)
    for i in fatty_acid_3:
        M3_fa3 += [(TAG_M1_mass[TAG_M1_mass_number - TAG_M1_mass_counter]) - ((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen + ammonium)]

        TAG_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2 + M3_fa3
    all_fatty_acids = fatty_acid_1 + fatty_acid_2 + fatty_acid_3
    all_TAG_id = TAG_id + TAG_id + TAG_id
    all_TAG_M1_mass = TAG_M1_mass + TAG_M1_mass + TAG_M1_mass
    all_TAG_FA_label = []

    TAG_FA_label_counter = len(all_TAG_id)
    TAG_FA_label_number = len(all_TAG_id)

    for i in all_TAG_id:
        all_TAG_FA_label += [i + '-FA' + all_fatty_acids[TAG_FA_label_number - TAG_FA_label_counter]]

        TAG_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'TAG ID': TAG_id, 'FA composition': tag, 'M1 mass (+NH4)': TAG_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2, 'M3 -FA3 mass': M3_fa3})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'TAG ID': all_TAG_id, 'FA loss': all_fatty_acids, 'TAG-FA label': all_TAG_FA_label, 'M1 mass (+NH4)': all_TAG_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)