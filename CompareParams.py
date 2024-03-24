import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import pandas as pd
import csv

# import from other files
from lattice_additions import Bethe, total_number_for_nShells, random_connections, super_spreader, combine_lattice
# from Plot_Layout import hierarchy_pos, update_node_color
from Epidemic_algorithm import SIR_only_bond, SIR_bond_site_vaccination_with_time 
from Average import average

# genrate a 2D Graph 
grid_size= 151
population = grid_size*grid_size
print(population)
grid_2D = nx.grid_2d_graph(grid_size,grid_size, periodic=True)
node_grid_2D = list(grid_2D.nodes)

#we need to convert the labels of the nodes to integers from 0 to N such that the SIR simulation works for this lattice
grid_2D_renamed = nx.convert_node_labels_to_integers(grid_2D, first_label=0) 

# get a list of all nodes in the grid
grid_2D_renamed_nodes = list(grid_2D_renamed.nodes)

# extract the position of the nodes from the grid without renamed nodes
pos_tuple = [(y,-x) for x,y in grid_2D.nodes()]
pos = { grid_2D_renamed_nodes[i]: pos_tuple[i] for i in range(len(grid_2D_renamed_nodes)) }

grid_2D_renamed_copy=grid_2D_renamed.copy()
grid_2D_combined=combine_lattice(grid_2D_renamed_copy,n_shortcuts=int(population*0.1), n_spreader=int(population*0.01), n_spreader_edges=8)


# Average for different Vacc_per_day

save_name = 'InstaVacc_n50'

Params_Dict = {'p_trans' : 0.10, 
               'p_reinf_vac': 0.03, 
                'p_reinf_reco':0.01,  
                'Num_vac_per_Tstep': 1/80,
                'T':300, 
                'Time_till_recovery':7, 
                'death_rate': 0.05,
                'Time_till_reinfec': 15,
                'vac_thresh': 0, 
                't_vac_develop':0, 
                'pP_notVac': 0.1,
                'pP_antiVac': 0}

# save params in csv
f = open('Data/InstaVacc/Params'+save_name+ '.csv', "w")
w = csv.writer(f)

# loop over dictionary keys and values
for key, val in Params_Dict.items():
    # write every key and value to file
    w.writerow([key, val])
    
f.close()

# Vacc_per_day = np.array([0, 1/320,1/160, 1/80, 3/80])
# Vacc_file_name= {str(Vacc_per_day[0]): '0',
#                 str(Vacc_per_day[1]): '1_320',
#                 str(Vacc_per_day[2]): '1_160',
#                 str(Vacc_per_day[3]): '1_80',
#                 str(Vacc_per_day[4]): '3_80'}



# for Vrate in Vacc_per_day:
#     Dict_res = average(network=grid_2D_combined, 
#                     p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
#                     p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
#                     T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
#                     Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=Params_Dict['vac_thresh'], 
#                     t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'], pP_antiVac=Params_Dict['pP_antiVac'], 
#                     patient_zero = int((grid_size**2)*1/2),n_average=1)

#     df = pd.DataFrame(Dict_res)

#     df.to_csv('Data/VaccPerDay/'+ save_name + '_' + str(Vacc_file_name[str(Vrate)]) + '.csv', index=True)


############### anti Vaxxers #########################################

# percentage_anti_vaxxers = np.array([0, 0.05, 0.1, 0.3, 0.5, 0.7, 1])
# AntiVacc_file_name= {str(percentage_anti_vaxxers[0]): '0',
#                 str(percentage_anti_vaxxers[1]): '005',
#                 str(percentage_anti_vaxxers[2]): '01',
#                 str(percentage_anti_vaxxers[3]): '03',
#                 str(percentage_anti_vaxxers[4]): '05',
#                 str(percentage_anti_vaxxers[5]): '07',
#                 str(percentage_anti_vaxxers[6]): '1'}

# percentage_anti_vaxxers = np.array([0.1])
# AntiVacc_file_name= {str(percentage_anti_vaxxers[0]): '01'}



# for pav in percentage_anti_vaxxers:
#     Dict_res = average(network=grid_2D_combined, 
#                     p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
#                     p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
#                     T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
#                     Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=Params_Dict['vac_thresh'], 
#                     t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'], pP_antiVac=pav, 
#                     patient_zero = int((grid_size**2)*1/2),n_average=1000)

#     df = pd.DataFrame(Dict_res)

#     df.to_csv('Data/AntiVacc/'+ save_name + '_' + str(AntiVacc_file_name[str(pav)]) + '.csv', index=True)


##################### Vacc Devlopement time ###########################################

# Vacc_start = np.array([0, 5, 10, 15, 20])
# Vacc_file_name= {str(Vacc_start[0]): '0',
#                 str(Vacc_start[1]): '5',
#                 str(Vacc_start[2]): '10',
#                 str(Vacc_start[3]): '15',
#                 str(Vacc_start[4]): '20'}


# for t_s in Vacc_start:
#     Dict_res = average(network=grid_2D_combined, 
#                     p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
#                     p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
#                     T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
#                     Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=Params_Dict['vac_thresh'], 
#                     t_vac_develop=t_s, pP_notVac=Params_Dict['pP_notVac'],pP_antiVac=Params_Dict['pP_antiVac'], 
#                     patient_zero = int((grid_size**2)*1/2),n_average=50)

#     df = pd.DataFrame(Dict_res)

#     df.to_csv('Data/VaccStart/'+ save_name + '_' + str(Vacc_file_name[str(t_s)]) + '.csv', index=True)


############################# Vacc Threshold #######################################################

# Vacc_thresh = np.array([0.0001, 0.001, 0.01 , 0.05, 0.1, 1])
# Vacc_file_name= {str(Vacc_thresh[0]): '00001',
#                 str(Vacc_thresh[1]): '0001',
#                 str(Vacc_thresh[2]): '001',
#                 str(Vacc_thresh[3]): '005',
#                 str(Vacc_thresh[4]): '01',
#                 str(Vacc_thresh[5]): '1'}

# for thresh in Vacc_thresh:
#     Dict_res = average(network=grid_2D_combined, 
#                     p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
#                     p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
#                     T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
#                     Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=thresh, 
#                     t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'],pP_antiVac=Params_Dict['pP_antiVac'], 
#                     patient_zero = int((grid_size**2)*1/2), n_average=50)


#     df = pd.DataFrame(Dict_res)

#     df.to_csv('Data/VaccThresh/'+ save_name + '_' + str(Vacc_file_name[str(thresh)]) + '.csv', index=True)
# so
    

######################################## vacc against death ##############################################################
    

# # p_inf_vac = np.array([0, 0.01, 0.03, 0.05, 0.07, 0.1])
# # InfVacc_file_name= {str(p_inf_vac[0]): '0',
# #                 str(p_inf_vac[1]): '001',
# #                 str(p_inf_vac[2]): '003',
# #                 str(p_inf_vac[3]): '005',
# #                 str(p_inf_vac[4]): '007',
# #                 str(p_inf_vac[5]): '1'}

# p_inf_vac = np.array([1])
# InfVacc_file_name= {str(p_inf_vac[0]): 'NoVacc'}




# for p in p_inf_vac:
#     Dict_res = average(network=grid_2D_combined, 
#                     p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
#                     p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
#                     T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
#                     Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=p, 
#                     t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'], pP_antiVac=Params_Dict['pP_antiVac'], 
#                     patient_zero = int((grid_size**2)*1/2),n_average=50)

#     df = pd.DataFrame(Dict_res)

#     df.to_csv('Data/VaccReinf/'+ save_name + '_' + str(InfVacc_file_name[str(p)]) + '.csv', index=True)


################################## insta Impfung ###########################################################

# Vacc_per_day = np.array([0, 0.2, 0.35, 0.4, 0.45, 0.5])
# Vacc_file_name= {str(Vacc_per_day[0]): '0',
#                 str(Vacc_per_day[1]): '02',
#                 str(Vacc_per_day[2]): '035',
#                 str(Vacc_per_day[3]): '04',
#                 str(Vacc_per_day[4]): '045',
#                 str(Vacc_per_day[5]): '05'}

# Vacc_per_day = np.array([0.35, 0.4, 0.45])
# Vacc_file_name= {str(Vacc_per_day[0]): '035',
#                 str(Vacc_per_day[1]): '04',
#                 str(Vacc_per_day[2]): '045'}

Vacc_per_day = np.array([0.55, 0.6])
Vacc_file_name= {str(Vacc_per_day[0]): '055',
                str(Vacc_per_day[1]): '06'}



for Vrate in Vacc_per_day:
    Dict_res = average(network=grid_2D_combined, 
                    p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
                    p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Vrate, 
                    T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], 
                    Time_till_reinfec=Params_Dict['Time_till_reinfec'],vac_thresh=Params_Dict['vac_thresh'], 
                    t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'], pP_antiVac=Params_Dict['pP_antiVac'], 
                    patient_zero = int((grid_size**2)*1/2),n_average=50)

    df = pd.DataFrame(Dict_res)

    df.to_csv('Data/InstaVacc/'+ save_name + '_' + str(Vacc_file_name[str(Vrate)]) + '.csv', index=True)
    