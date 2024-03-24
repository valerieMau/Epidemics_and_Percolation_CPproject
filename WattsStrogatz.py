import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
import csv
import pandas as pd

# import from other files
from Lattice_additions import Bethe, total_number_for_nShells, random_connections, combine_lattice, super_spreader
from Epidemic_algorithm import SIR_bond_site_vaccination_with_time, SIR_only_bond
from Average import average

Ani_plot = 1

population =20
WattsStrog = nx.newman_watts_strogatz_graph(population, 4 , 0.1)
WattsStrog_combined = super_spreader(WattsStrog, n_spreader=int(population*0.01), n_edges=8)

pos = nx.circular_layout(WattsStrog)

nx.draw(WattsStrog_combined, pos, node_color='grey',node_size=60)
plt.show()

Params_Dict = {'p_trans' : 0.10, 
               'p_reinf_vac': 0.03, 
                'p_reinf_reco':0.01,  
                'Num_vac_per_Tstep': 1/80,
                'T':160, 
                'Time_till_recovery':7, 
                'death_rate': 0.05,
                'Time_till_reinfec': 15,
                'vac_thresh': 0.001, 
                't_vac_develop':15, 
                'pP_notVac': 0.1,
                'pP_antiVac': 0.1,
                'k': 4,
                'population': population,
                'p_rand_con': 0.1}


S,I,R,V,D, t, color, NS, NI, NR, ND = SIR_bond_site_vaccination_with_time(network=WattsStrog_combined, p_trans=Params_Dict['p_trans'], p_reinf_vac=Params_Dict['p_reinf_vac'], 
                                                            p_reinf_reco=Params_Dict['p_reinf_reco'], Num_vac_per_Tstep=Params_Dict['Num_vac_per_Tstep'], 
                                                            T=Params_Dict['T'], Time_till_recovery=Params_Dict['Time_till_recovery'], death_rate=Params_Dict['death_rate'], Time_till_reinfec=Params_Dict['Time_till_reinfec'],
                                                            vac_thresh=Params_Dict['vac_thresh'], t_vac_develop=Params_Dict['t_vac_develop'], pP_notVac=Params_Dict['pP_notVac'],pP_antiVac=Params_Dict['pP_antiVac'], patient_zero = 0)



if Ani_plot == 1:
    # Initialize the figure
    fig, ax = plt.subplots()
    # get the nodes colors
    SIR_bond_animation_colors = color


    def update_node_color(frame, network, SIR_bond_animation_colors, pos, edge_colors=None):
        ax.clear()
        node_color_frame = SIR_bond_animation_colors[frame]
        
        # you can either draw the network with edges
        nx.draw(network, pos, with_labels=False, node_color=node_color_frame, node_size=10)
        # or only draw the nodes
        # nx.draw_networkx_nodes(network, pos, node_color=node_color_frame,node_size=10)
        ax.set_title(f'Time step {frame}')

    # start the animation
    ani = FuncAnimation(fig, update_node_color, frames=100, fargs=[WattsStrog_combined,SIR_bond_animation_colors, pos], interval=100)
    # ani.save('Data/2DGrid_Simulation/Animation'+save_name+'.gif',  writer='Pillow')
    plt.show()

