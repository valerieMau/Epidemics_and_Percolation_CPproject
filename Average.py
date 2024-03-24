import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

# import from other files
from Lattice_additions import Bethe, total_number_for_nShells, random_connections, super_spreader, combine_lattice
# from Plot_Layout import hierarchy_pos, update_node_color
from Epidemic_algorithm import SIR_only_bond, SIR_bond_site_vaccination_with_time 


def average(network, p_trans, p_reinf_vac, 
            p_reinf_reco, Num_vac_per_Tstep, 
            T, Time_till_recovery, death_rate, Time_till_reinfec,
            vac_thresh, t_vac_develop, pP_notVac, pP_antiVac, patient_zero,n_average):
    I_all = []
    S_all = []
    R_all = []
    V_all = []
    D_all = []
    NI_all = []
    NS_all = []
    NR_all = []
    NV_all = []
    ND_all = []
    for n in tqdm (range(n_average),  desc="Loading"):
        # S, I, R, t, _ = SIR_bond_Scontrolled(network, p_trans, T, Time_till_recovery)
        S,I,R,V,D, t, color, NS, NI, NR, ND = SIR_bond_site_vaccination_with_time(network, p_trans, p_reinf_vac, 
                                                            p_reinf_reco, Num_vac_per_Tstep, 
                                                            T, Time_till_recovery, death_rate, Time_till_reinfec,
                                                            vac_thresh, t_vac_develop, pP_notVac, pP_antiVac, patient_zero)
        I_all.append(I)
        S_all.append(S)
        R_all.append(R)
        V_all.append(V)
        D_all.append(D)
        NI_all.append(NI)
        NS_all.append(NS)
        NR_all.append(NR)
        ND_all.append(ND)

    mean_I=np.average(np.array(I_all),axis=0)
    mean_S=np.average(np.array(S_all),axis=0)
    mean_R=np.average(np.array(R_all),axis=0)
    mean_V=np.average(np.array(V_all),axis=0)
    mean_D=np.average(np.array(D_all),axis=0)
    
    mean_NI=np.average(np.array(NI_all),axis=0)
    mean_NS=np.average(np.array(NS_all),axis=0)
    mean_NR=np.average(np.array(NR_all),axis=0)
    mean_ND=np.average(np.array(ND_all),axis=0)
    
    std_I=np.std(np.array(I_all),axis=0)
    std_S=np.std(np.array(S_all),axis=0)
    std_R=np.std(np.array(R_all),axis=0)
    std_V=np.std(np.array(V_all),axis=0)
    std_D=np.std(np.array(D_all),axis=0)
    
    std_NI=np.std(np.array(NI_all),axis=0)
    std_NS=np.std(np.array(NS_all),axis=0)
    std_NR=np.std(np.array(NR_all),axis=0)
    std_ND=np.std(np.array(ND_all),axis=0)
    
    Result_Dict = {'t': t,
                   'S': mean_S,
                   'I': mean_I,
                   'R': mean_R,
                   'D': mean_D,
                   'V': mean_V,
                   'NS': mean_NS,
                   'NI': mean_NI,
                   'NR': mean_NR,
                   'ND': mean_ND,
                   'std_S': std_S,
                   'std_I': std_I,
                   'std_R': std_R,
                   'std_D': std_D,
                   'std_V': std_V,
                   'std_NS': std_NS,
                   'std_NI': std_NI,
                   'std_NR': std_NR,
                   'std_ND': std_ND}

    return Result_Dict


def average_bond(network, p_trans, T, Time_till_recovery,patient_zero, n_average):
    I_all = []
    S_all = []
    R_all = []
    for n in tqdm (range(n_average),  desc="Loading"):
        S, I, R, t, _, _ = SIR_only_bond(network, p_trans, T, Time_till_recovery, patient_zero)

        I_all.append(I)
        S_all.append(S)
        R_all.append(R)

    mean_I=np.average(np.array(I_all),axis=0)
    mean_S=np.average(np.array(S_all),axis=0)
    mean_R=np.average(np.array(R_all),axis=0)

    std_I=np.std(np.array(I_all),axis=0)
    std_S=np.std(np.array(S_all),axis=0)
    std_R=np.std(np.array(R_all),axis=0)
    
    
    Result_Dict = {'t': t,
                   'S': mean_S,
                   'I': mean_I,
                   'R': mean_R,
                   'std_S': std_S,
                   'std_I': std_I,
                   'std_R': std_R}

    return Result_Dict


