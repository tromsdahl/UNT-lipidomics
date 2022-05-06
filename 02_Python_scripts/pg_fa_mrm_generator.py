#PG-FA MRM Generator 2.0
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
        
    

    pg_carbon = []
    pg_desat = []
    pg_desat_mass = []
    

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            pg_carbon += [fa_carbon_no[0] + i]

        fa_carbon_no.remove(fa_carbon_no[0])
            

    # to get the total fa desat number for label
    while fa_desat_no != []:
        for i in fa_desat_no:
            pg_desat += [fa_desat_no[0] + i]

        fa_desat_no.remove(fa_desat_no[0])
            

    # to get the total fa desat number for mass calculation
    while fa_desat_no_mass != []:
        for i in fa_desat_no_mass:
            pg_desat_mass += [fa_desat_no_mass[0] + i]

        fa_desat_no_mass.remove(fa_desat_no_mass[0])
            


    # to output the pg id (i.e. PG-(fa C no.):(fa desat no.))
    pg_id = []
    pg_M1_mass = []

    id_counter = len(pg_carbon)
    id_number = len(pg_carbon)

    
    while id_counter != 0:
        pg_id += ["PG-" + str(pg_carbon[id_number - id_counter]) + ":" + str(pg_desat[id_number - id_counter])]
        
        pg_M1_mass += [(carbon * (pg_carbon[id_number - id_counter] + 6)) + (hydrogen * ((pg_carbon[id_number - id_counter] * 2) - (pg_desat_mass[id_number - id_counter] * 2) + 11)) + (oxygen * 10) + phosphorus - hydrogen]

        id_counter -= 1
    


    # to output the FA combinations for each PG (this is the same as the original TAG/FA combinator)
    pg = []
        
    while fa1 != []:
        for i in fa1:
            pg += [fa1[0] + '/' + i]

        fa1.remove(fa1[0])




    # calculating each possible MRM transition (Ac adducts)
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []

    for i in pg:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_combo.clear()
    


    pg_M1_mass_counter = len(pg_M1_mass)
    pg_M1_mass_number = len(pg_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pg_M1_mass_counter -= 1
    
    pg_M1_mass_counter = len(pg_M1_mass)
    pg_M1_mass_number = len(pg_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen)]

        pg_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2
    all_fatty_acids = fatty_acid_1 + fatty_acid_2
    all_pg_id = pg_id + pg_id
    all_pg_M1_mass = pg_M1_mass + pg_M1_mass
    all_pg_FA_label = []

    pg_FA_label_counter = len(all_pg_id)
    pg_FA_label_number = len(all_pg_id)

    for i in all_pg_id:
        all_pg_FA_label += [i + '-FA' + all_fatty_acids[pg_FA_label_number - pg_FA_label_counter]]

        pg_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'PG ID': pg_id, 'FA composition': pg, 'M1 mass (-H)': pg_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'PG ID': all_pg_id, 'FA loss': all_fatty_acids, 'PG-FA label': all_pg_FA_label, 'M1 mass (-H)': all_pg_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)