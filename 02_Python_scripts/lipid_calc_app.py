#!/usr/bin/env python3

import re

C = 12.000000
C13 = 13.003355
H = 1.007825
O = 15.994915
P = 30.973762
N = 14.003074
NH4 = 18.03382555
Na = 22.989769
K = 38.963707
S = 31.972072



def TAG_MW(molec_fa_c, molec_fa_d):
    
    TAG = (C * (molec_fa_c + 3)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 2)) + (O * 6)
    TAG_H = TAG + H
    TAG_Na = TAG + Na
    TAG_K = TAG + K
    TAG_NH4 = TAG + NH4

    print('Monoisotopic mass: ' + str(round(TAG, 5)))
    print('[M + H]^+: ' + str(round(TAG_H, 5)))
    print('[M + Na]^+: ' + str(round(TAG_Na, 5)))
    print('[M + K]^+: ' + str(round(TAG_K, 5)))
    print('[M + NH4]^+: ' + str(round(TAG_NH4, 5)))


def DAG_MW(molec_fa_c, molec_fa_d):
    DAG = (C * (molec_fa_c + 3)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 4)) + (O * 5)
    DAG_H = DAG + H
    DAG_Na = DAG + Na
    DAG_K = DAG + K
    DAG_NH4 = DAG + NH4

    print('Monoisotopic mass: ' + str(round(DAG, 5)))
    print('[M + H]^+: ' + str(round(DAG_H, 5)))
    print('[M + Na]^+: ' + str(round(DAG_Na, 5)))
    print('[M + K]^+: ' + str(round(DAG_K, 5)))
    print('[M + NH4]^+: ' + str(round(DAG_NH4, 5)))


def MAG_MW(molec_fa_c, molec_fa_d):
    MAG = (C * (molec_fa_c + 3)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 6)) + (O * 4)
    MAG_H = MAG + H
    MAG_Na = MAG + Na
    MAG_K = MAG + K
    MAG_NH4 = MAG + NH4

    print('Monoisotopic mass: ' + str(round(MAG, 5)))
    print('[M + H]^+: ' + str(round(MAG_H, 5)))
    print('[M + Na]^+: ' + str(round(MAG_Na, 5)))
    print('[M + K]^+: ' + str(round(MAG_K, 5)))
    print('[M + NH4]^+: ' + str(round(MAG_NH4, 5)))


def FA_MW(molec_fa_c, molec_fa_d):
    FA = (C * molec_fa_c) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2))) + (O * 2)
    FA_H = FA + H
    FA_neg_H = FA - H
    FA_Na = FA + Na
    FA_K = FA + K
    FA_NH4 = FA + NH4

    print('Monoisotopic mass: ' + str(round(FA, 5)))
    print('[M + H]^+: ' + str(round(FA_H, 5)))
    print('[M - H]^-: ' + str(round(FA_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(FA_Na, 5)))
    print('[M + K]^+: ' + str(round(FA_K, 5)))
    print('[M + NH4]^+: ' + str(round(FA_NH4, 5)))


def PC_MW(molec_fa_c, molec_fa_d):

    PC = (C * (molec_fa_c + 8)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 16)) + (O * 8) + N + P
    PC_H = PC + H
    PC_Na = PC + Na
    PC_K = PC + K
    PC_NH4 = PC + NH4

    print('Monoisotopic mass: ' + str(round(PC, 5)))
    print('[M + H]^+: ' + str(round(PC_H, 5)))
    print('[M + Na]^+: ' + str(round(PC_Na, 5)))
    print('[M + K]^+: ' + str(round(PC_K, 5)))
    print('[M + NH4]^+: ' + str(round(PC_NH4, 5)))


def PE_MW(molec_fa_c, molec_fa_d):
    
    PE = (C * (molec_fa_c + 5)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 10)) + (O * 8) + N + P
    PE_H = PE + H
    PE_neg_H = PE - H
    PE_Na = PE + Na
    PE_K = PE + K
    PE_NH4 = PE + NH4

    print('Monoisotopic mass: ' + str(round(PE, 5)))
    print('[M + H]^+: ' + str(round(PE_H, 5)))
    print('[M - H]^-: ' + str(round(PE_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(PE_Na, 5)))
    print('[M + K]^+: ' + str(round(PE_K, 5)))
    print('[M + NH4]^+: ' + str(round(PE_NH4, 5)))


def PI_MW(molec_fa_c, molec_fa_d):

    PI = (C * (molec_fa_c + 9)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 15)) + (O * 13) + P
    PI_H = PI + H
    PI_neg_H = PI - H
    PI_Na = PI + Na
    PI_K = PI + K
    PI_NH4 = PI + NH4

    print('Monoisotopic mass: ' + str(round(PI, 5)))
    print('[M + H]^+: ' + str(round(PI_H, 5)))
    print('[M - H]^-: ' + str(round(PI_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(PI_Na, 5)))
    print('[M + K]^+: ' + str(round(PI_K, 5)))
    print('[M + NH4]^+: ' + str(round(PI_NH4, 5)))


def PG_MW(molec_fa_c, molec_fa_d):

    PG = (C * (molec_fa_c + 6)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 11)) + (O * 10) + P
    PG_H = PG + H
    PG_Na = PG + Na
    PG_K = PG + K
    PG_NH4 = PG + NH4

    print('Monoisotopic mass: ' + str(round(PG, 5)))
    print('[M + H]^+: ' + str(round(PG_H, 5)))
    print('[M + Na]^+: ' + str(round(PG_Na, 5)))
    print('[M + K]^+: ' + str(round(PG_K, 5)))
    print('[M + NH4]^+: ' + str(round(PG_NH4, 5)))


def PS_MW(molec_fa_c, molec_fa_d):
    PS = (C * (molec_fa_c + 6)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 10)) + (O * 10) + P + N
    PS_H = PS + H
    PS_neg_H = PS - H
    PS_Na = PS + Na
    PS_K = PS + K
    PS_NH4 = PS + NH4

    print('Monoisotopic mass: ' + str(round(PS, 5)))
    print('[M + H]^+: ' + str(round(PS_H, 5)))
    print('[M - H]^-: ' + str(round(PS_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(PS_Na, 5)))
    print('[M + K]^+: ' + str(round(PS_K, 5)))
    print('[M + NH4]^+: ' + str(round(PS_NH4, 5)))


def PA_MW(molec_fa_c, molec_fa_d):
    PA = (C * (molec_fa_c + 3)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 5)) + (O * 8) + P
    PA_H = PA + H
    PA_neg_H = PA - H
    PA_Na = PA + Na
    PA_K = PA + K
    PA_NH4 = PA + NH4

    print('Monoisotopic mass: ' + str(round(PA, 5)))
    print('[M + H]^+: ' + str(round(PA_H, 5)))
    print('[M - H]^-: ' + str(round(PA_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(PA_Na, 5)))
    print('[M + K]^+: ' + str(round(PA_K, 5)))
    print('[M + NH4]^+: ' + str(round(PA_NH4, 5)))


def MGDG_MW(molec_fa_c, molec_fa_d):
    MGDG = (C * (molec_fa_c + 9)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 14)) + (O * 10)
    MGDG_H = MGDG + H
    MGDG_Na = MGDG + Na
    MGDG_K = MGDG + K
    MGDG_NH4 = MGDG + NH4

    print('Monoisotopic mass: ' + str(round(MGDG, 5)))
    print('[M + H]^+: ' + str(round(MGDG_H, 5)))
    print('[M + Na]^+: ' + str(round(MGDG_Na, 5)))
    print('[M + K]^+: ' + str(round(MGDG_K, 5)))
    print('[M + NH4]^+: ' + str(round(MGDG_NH4, 5)))


def DGDG_MW(molec_fa_c, molec_fa_d):
    DGDG = (C * (molec_fa_c + 15)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 24)) + (O * 15)
    DGDG_H = DGDG + H
    DGDG_Na = DGDG + Na
    DGDG_K = DGDG + K
    DGDG_NH4 = DGDG + NH4

    print('Monoisotopic mass: ' + str(round(DGDG, 5)))
    print('[M + H]^+: ' + str(round(DGDG_H, 5)))
    print('[M + Na]^+: ' + str(round(DGDG_Na, 5)))
    print('[M + K]^+: ' + str(round(DGDG_K, 5)))
    print('[M + NH4]^+: ' + str(round(DGDG_NH4, 5)))


def SQDG_MW(molec_fa_c, molec_fa_d):
    SQDG = (C * (molec_fa_c + 9)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 14)) + (O * 12) + S
    SQDG_H = SQDG + H
    SQDG_neg_H = SQDG - H
    SQDG_Na = SQDG + Na
    SQDG_K = SQDG + K
    SQDG_NH4 = SQDG + NH4

    print('Monoisotopic mass: ' + str(round(SQDG, 5)))
    print('[M + H]^+: ' + str(round(SQDG_H, 5)))
    print('[M - H]^-: ' + str(round(SQDG_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(SQDG_Na, 5)))
    print('[M + K]^+: ' + str(round(SQDG_K, 5)))
    print('[M + NH4]^+: ' + str(round(SQDG_NH4, 5)))


def LPC_MW(molec_fa_c, molec_fa_d):
    LPC = (C * (molec_fa_c + 8)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 18)) + (O * 7) + N + P
    LPC_H = LPC + H
    LPC_Na = LPC + Na
    LPC_K = LPC + K
    LPC_NH4 = LPC + NH4

    print('Monoisotopic mass: ' + str(round(LPC, 5)))
    print('[M + H]^+: ' + str(round(LPC_H, 5)))
    print('[M + Na]^+: ' + str(round(LPC_Na, 5)))
    print('[M + K]^+: ' + str(round(LPC_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPC_NH4, 5)))


def LPE_MW(molec_fa_c, molec_fa_d):
    LPE = (C * (molec_fa_c + 5)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 12)) + (O * 7) + N + P
    LPE_H = LPE + H
    LPE_Na = LPE + Na
    LPE_K = LPE + K
    LPE_NH4 = LPE + NH4

    print('Monoisotopic mass: ' + str(round(LPE, 5)))
    print('[M + H]^+: ' + str(round(LPE_H, 5)))
    print('[M + Na]^+: ' + str(round(LPE_Na, 5)))
    print('[M + K]^+: ' + str(round(LPE_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPE_NH4, 5)))


def LPI_MW(molec_fa_c, molec_fa_d):
    LPI = (C * (molec_fa_c + 9)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 17)) + (O * 12) + P
    LPI_H = LPI + H
    LPI_neg_H = LPI - H
    LPI_Na = LPI + Na
    LPI_K = LPI + K
    LPI_NH4 = LPI + NH4

    print('Monoisotopic mass: ' + str(round(LPI, 5)))
    print('[M + H]^+: ' + str(round(LPI_H, 5)))
    print('[M - H]^-: ' + str(round(LPI_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(LPI_Na, 5)))
    print('[M + K]^+: ' + str(round(LPI_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPI_NH4, 5)))


def LPG_MW(molec_fa_c, molec_fa_d):

    LPG = (C * (molec_fa_c + 6)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 13)) + (O * 9) + P
    LPG_H = LPG + H
    LPG_Na = LPG + Na
    LPG_K = LPG + K
    LPG_NH4 = LPG + NH4

    print('Monoisotopic mass: ' + str(round(LPG, 5)))
    print('[M + H]^+: ' + str(round(LPG_H, 5)))
    print('[M + Na]^+: ' + str(round(LPG_Na, 5)))
    print('[M + K]^+: ' + str(round(LPG_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPG_NH4, 5)))


def LPS_MW(molec_fa_c, molec_fa_d):

    LPS = (C * (molec_fa_c + 6)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 12)) + (O * 9) + P + N
    LPS_H = LPS + H
    LPS_neg_H = LPS - H
    LPS_Na = LPS + Na
    LPS_K = LPS + K
    LPS_NH4 = LPS + NH4

    print('Monoisotopic mass: ' + str(round(LPS, 5)))
    print('[M + H]^+: ' + str(round(LPS_H, 5)))
    print('[M - H]^-: ' + str(round(LPS_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(LPS_Na, 5)))
    print('[M + K]^+: ' + str(round(LPS_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPS_NH4, 5)))


def LPA_MW(molec_fa_c, molec_fa_d):

    LPA = (C * (molec_fa_c + 3)) + (H * ((molec_fa_c * 2) - (molec_fa_d * 2) + 7)) + (O * 7) + P
    LPA_H = LPA + H
    LPA_neg_H = LPA - H
    LPA_Na = LPA + Na
    LPA_K = LPA + K
    LPA_NH4 = LPA + NH4

    print('Monoisotopic mass: ' + str(round(LPA, 5)))
    print('[M + H]^+: ' + str(round(LPA_H, 5)))
    print('[M - H]^-: ' + str(round(LPA_neg_H, 5)))
    print('[M + Na]^+: ' + str(round(LPA_Na, 5)))
    print('[M + K]^+: ' + str(round(LPA_K, 5)))
    print('[M + NH4]^+: ' + str(round(LPA_NH4, 5)))


#def PC-13C_MW(molec_fa_c, molec_fa_d):
#def TAG-13C_MW(molec_fa_c, molec_fa_d):


while True:

    try:

        molec_input = input("Enter molecular name (e.g. TAG-54:3); enter 'H' for help; enter 'Q' to exit: ")
        
        lipid_id = re.split('-|:', molec_input)
        lipid_class = lipid_id[0]
        lipid_fa_c = int(lipid_id[1])
        lipid_fa_d = int(lipid_id[2])
    
    except IndexError:
        
        if molec_input.upper() == 'Q':
            
            print('All done!')
            break
        
        elif molec_input.upper() == 'H':
            print('Available lipid classes include: TAG, DAG, MAG, FA, PC, PE, PI, PG, PS, PA, LPC, LPE, LPI, LPG, LPS, LPA, MGDG, DGDG, SQDG\nUse the following format to enter (without parantheses): (lipid class)-(FA carbons):(FA unsaturations)')
            continue
        
        else:
        
            print('Invalid.')
            continue
    
    except ValueError:
        
        if molec_input.upper() == 'Q':
            
            print('All done!')
            break

        elif molec_input.upper() == 'H':
            print('Available lipid classes include: TAG, DAG, MAG, FA, PC, PE, PI, PG, PS, PA, LPC, LPE, LPI, LPG, LPS, LPA, MGDG, DGDG, SQDG\nUse the following format to enter (without parantheses): (lipid class)-(FA carbons):(FA unsaturations)')
            continue
        
        else:
        
            print('Invalid.')
            continue

    if lipid_class.upper() == 'TAG':

        TAG_MW(lipid_fa_c,lipid_fa_d)
    
    elif lipid_class.upper() == 'DAG':

        DAG_MW(lipid_fa_c,lipid_fa_d)

    elif lipid_class.upper() == 'MAG':

        MAG_MW(lipid_fa_c,lipid_fa_d)

    elif lipid_class.upper() == 'FA':

        FA_MW(lipid_fa_c,lipid_fa_d)

    elif lipid_class.upper() == 'PC':

        PC_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'PE':

        PE_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'PI':

        PI_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'PG':

        PG_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'PS':

        PS_MW(lipid_fa_c, lipid_fa_d)
    
    elif lipid_class.upper() == 'PA':

        PA_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'MGDG':

        MGDG_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'DGDG':

        DGDG_MW(lipid_fa_c, lipid_fa_d)
    
    elif lipid_class.upper() == 'SQDG':

        SQDG_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'LPC':

        LPC_MW(lipid_fa_c, lipid_fa_d)
    
    elif lipid_class.upper() == 'LPE':

        LPE_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'LPI':

        LPI_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'LPG':

        LPG_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'LPS':

        LPS_MW(lipid_fa_c, lipid_fa_d)

    elif lipid_class.upper() == 'LPA':

        LPA_MW(lipid_fa_c, lipid_fa_d)

    else:

        print('Not available.')