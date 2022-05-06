#DAG-FA MRM Generator 2.0
import pandas as pd


fatty_acids = ['16:0', '18:3', '18:2', '18:1', '18:0', '20:2', '20:1', '20:0', '22:1', '22:0']      # Arabidopsis https://plantfadb.org/datasets?plant_id=17354


def mrm_output(fa1):

    carbon = 12.000000
    hydrogen = 1.007825
    oxygen = 15.994915
    ammonium = 18.03382555

    fa_carbon_no = []
    fa_desat_no = []
    fa_desat_no_mass = []

    # converting string to int for calculations below
    for i in fa1:
        fa_carbon_no.append(int(i[:2]))
        fa_desat_no.append(int(i[-1]))
        fa_desat_no_mass.append(int(i[-1]))
        


    dag_carbon = []
    dag_desat = []
    dag_desat_mass = []
    

    # to get the total fa carbon number
    while fa_carbon_no != []:
        for i in fa_carbon_no:
            dag_carbon += [fa_carbon_no[0] + i]

        fa_carbon_no.remove(fa_carbon_no[0])
            

    # to get the total fa desat number for label
    while fa_desat_no != []:
        for i in fa_desat_no:
            dag_desat += [fa_desat_no[0] + i]

        fa_desat_no.remove(fa_desat_no[0])
            

    # to get the total fa desat number for mass calculation
    while fa_desat_no_mass != []:
        for i in fa_desat_no_mass:
            dag_desat_mass += [fa_desat_no_mass[0] + i]

        fa_desat_no_mass.remove(fa_desat_no_mass[0])
            

    
    # to output the dag id (i.e. DAG-(fa C no.):(fa desat no.))
    dag_id = []
    dag_M1_mass = []

    id_counter = len(dag_carbon)
    id_number = len(dag_carbon)


    while id_counter != 0:
        dag_id += ["DAG-" + str(dag_carbon[id_number - id_counter]) + ":" + str(dag_desat[id_number - id_counter])]
        
        dag_M1_mass += [(carbon * (dag_carbon[id_number - id_counter] + 3)) + (hydrogen * ((dag_carbon[id_number - id_counter] * 2) - (dag_desat_mass[id_number - id_counter] * 2) + 4)) + (oxygen * 5) + ammonium]

        id_counter -= 1
    
    

    # to output the FA combinations for each DAG (this is the same as the original TAG/FA combinator)    
    dag = []

    while fa1 != []:
        for i in fa1:
            dag += [fa1[0] + '/' + i]

        fa1.remove(fa1[0])



    # calculating each possible MRM transition
    fatty_acid_1 = []
    fatty_acid_2 = []
    fatty_acid_combo = []
    
    M3_fa1 = []
    M3_fa2 = []

    for i in dag:
        fatty_acid_combo = i.split('/')
        fatty_acid_1 += [fatty_acid_combo[0]]
        fatty_acid_2 += [fatty_acid_combo[1]]
        fatty_acid_combo.clear()
    

    dag_M1_mass_counter = len(dag_M1_mass)
    dag_M1_mass_number = len(dag_M1_mass)
    for i in fatty_acid_1:
        M3_fa1 += [(dag_M1_mass[dag_M1_mass_number - dag_M1_mass_counter]) - ((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen + ammonium)]

        dag_M1_mass_counter -= 1
    
    dag_M1_mass_counter = len(dag_M1_mass)
    dag_M1_mass_number = len(dag_M1_mass)
    for i in fatty_acid_2:
        M3_fa2 += [(dag_M1_mass[dag_M1_mass_number - dag_M1_mass_counter]) - ((int(i[0:2]) * carbon) + (hydrogen * ((int(i[0:2]) * 2) - (int(i[-1]) * 2))) + (oxygen * 2) - hydrogen + ammonium)]

        dag_M1_mass_counter -= 1



    # Final MRM table output
    all_M3_transitions = M3_fa1 + M3_fa2
    all_fatty_acids = fatty_acid_1 + fatty_acid_2
    all_dag_id = dag_id + dag_id
    all_dag_M1_mass = dag_M1_mass + dag_M1_mass
    all_dag_FA_label = []

    dag_FA_label_counter = len(all_dag_id)
    dag_FA_label_number = len(all_dag_id)

    for i in all_dag_id:
        all_dag_FA_label += [i + '-FA' + all_fatty_acids[dag_FA_label_number - dag_FA_label_counter]]

        dag_FA_label_counter -= 1



    MRM_table = pd.DataFrame({'DAG ID': dag_id, 'FA composition': dag, 'M1 mass (+NH4)': dag_M1_mass, 'M3 -FA1 mass': M3_fa1, 'M3 -FA2 mass': M3_fa2})
    print(MRM_table)
    MRM_table.to_csv('MRM_table.csv')
 

    MRM_condensed_table = pd.DataFrame({'DAG ID': all_dag_id, 'FA loss': all_fatty_acids, 'DAG-FA label': all_dag_FA_label, 'M1 mass (+NH4)': all_dag_M1_mass, 'M3 -FA mass': all_M3_transitions})
    MRM_condensed_table_unique = MRM_condensed_table.drop_duplicates()
    
    print(MRM_condensed_table)
    print(MRM_condensed_table_unique)
    
    MRM_condensed_table.to_csv('MRM_condensed_table.csv')
    MRM_condensed_table_unique.to_csv('MRM_condensed_table_unique.csv')


mrm_output(fatty_acids)