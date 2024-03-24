import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import csv
import pandas as pd

# import from other files
from Lattice_additions import Bethe, total_number_for_nShells, random_connections, combine_lattice
from Epidemic_algorithm import SIR_bond_site_vaccination_with_time, SIR_only_bond

# Animation Plot, if 1 show animation
Ani_plot = 1

# genrate a 2D Graph 
grid_size= 151
population = grid_size*grid_size
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

print(population)




######### define paramter ###########

save_name = '2d_grid_2'

Params_Dict = {'p_trans' : 0.10, 
               'p_reinf_vac': 0.03, 
                'p_reinf_reco':0.01,  
                'Num_vac_per_Tstep': 1/80,
                'T':150, 
                'Time_till_recovery':7, 
                'death_rate': 0.01,
                'Time_till_reinfec': 15,
                'vac_thresh': 0.01, 
                't_vac_develop':15, 
                'pP_notVac': 0.1,
                'pP_antiVac': 0.1}

##########################

# save params in csv
# f = open('Data/2DGrid_Simulation/'+save_name+ '.csv', "w")
# w = csv.writer(f)

# # loop over dictionary keys and values
# for key, val in Params_Dict.items():
#     # write every key and value to file
#     w.writerow([key, val])
    
# f.close()

# df_params = pd.DataFrame(Params_Dict, index=[0])

# df_params.to_csv('Data/2DGrid_Simulation/'+save_name+ '.csv', index=True)

########################################

S,I,R,V,D, t, color, NS, NI, NR, ND = SIR_bond_site_vaccination_with_time(network=grid_2D_combined, p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
                                                            p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
                                                            T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], Time_till_reinfec=Params_Dict['Time_till_reinfec'],
                                                            vac_thresh=Params_Dict['vac_thresh'], t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'],pP_antiVac=Params_Dict['pP_antiVac'], patient_zero = int((grid_size**2)*1/2))






# plot evolution of S,I,R,V,D with time

# plt.figure(1)
# plt.errorbar(t, S, label='S', fmt='.', color='green')   
# plt.errorbar(t, I, label='I', fmt='.', color='red')  
# plt.errorbar(t, R, label='R', fmt='.', color='blue') 
# # plt.errorbar(t, V, label='V', fmt='.', color='orange')   
# # plt.errorbar(t, D, label='D', fmt='.', color='black')  
# plt.xlabel('time')
# plt.ylabel('population')
# plt.grid()
# plt.legend()
# plt.show() 

# plt.figure(2)
# plt.errorbar(t, NS, label='S', fmt='.', color='magenta')   
# plt.errorbar(t, NI, label='I', fmt='.', color='red')  
# plt.errorbar(t, NR, label='R', fmt='.', color='blue')    
# plt.errorbar(t, ND, label='D', fmt='.', color='black')  
# plt.xlabel('time')
# plt.ylabel('population')
# plt.grid()
# plt.legend()
# plt.title('nodes without vaccination')
# plt.show() 

#show animation 


if Ani_plot == 1:
    # Initialize the figure
    fig, ax = plt.subplots()
    # get the nodes colors
    SIR_bond_animation_colors = color


    def update_node_color(frame, network, SIR_bond_animation_colors, pos, edge_colors=None):
        ax.clear()
        node_color_frame = SIR_bond_animation_colors[frame]
        
        # you can either draw the network with edges
        # nx.draw(network, pos, with_labels=False, node_color=node_color_frame, edge_color=edge_colors, node_size=10)
        # or only draw the nodes
        nx.draw_networkx_nodes(network, pos, node_color=node_color_frame,node_size=10)
        ax.set_title(f'Time step {frame}')

    # start the animation
    ani = FuncAnimation(fig, update_node_color, frames=100, fargs=[grid_2D_renamed,SIR_bond_animation_colors, pos], interval=100)
    # ani.save('Data/2DGrid_Simulation/Animation'+save_name+'.gif',  writer='Pillow')
    plt.show()

