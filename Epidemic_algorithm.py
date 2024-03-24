import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

# import from other files
from Lattice_additions import Bethe, total_number_for_nShells, random_connections, combine_lattice


def SIR_bond_site_vaccination_with_time(network, p_trans, p_reinf_vac, p_reinf_reco, 
                                        Num_vac_per_Tstep, T, Time_till_recovery, death_rate, 
                                        Time_till_reinfec, vac_thresh, t_vac_develop, pP_notVac , 
                                        pP_antiVac, patient_zero = 0):
    """
    We combine here bond and site percolation. The site percolation is used for describing the vaccination rate.
    Here we define the vaccinated people in the beginning and dont say that the vaccination starts at a certain time step
    
    Args:
        network : the lattice
        p_trans : transmissibility
        p_reinf_vac: probability to become reinfected after vaccination
        p_reinf_reco: probability to become reinfected after recovering
        vac_rate: this is the maximal vaccination rate
        vac_thresh: threshold of infected people for starting vaccination
        Num_vac_per_Tstep: Number of vaccinations per time step --> aktuell ist es implementiert, dass das ein Anteil der Gesamtbevölkerung ist
        T: number of time steps
        Time_till_recovery: recovery time
        vac_start: after how many timsteps vaccination starts
        t_vac_develop: time it takes to develope the vaccination 
        pP_antiVac: nodes that refuse vaccination
        

    """
    #remove random seed used in the grid to make sure that simulations are different
    np.random.seed(None)
    
    # death rate for vacc nodes
    vacc_death_rate = death_rate * 0.3
    
    # we assign a color to every state: I: 'firebrick', S: 'silver', R: 'steelblue' for the animation
    nodes_color = []
    nodes_color_timestep0 = [] #since we start our for loop at t=2 we have to keep track of the colors for the first time step putside of the loop
    
    # we want to keep track of the number of infected I, susceptible S and recovefirebrick R people
    i = 0
    s = 0
    r = 0 
    v = 0
    d = 0
    
    ns = 0
    nd = 0
    ni = 0
    nr = 0
    
    ############### time step 0 #############
    
    # assign every node a state S or not S depending on the susceptibility
    notVac_nodes = np.random.choice(len(network.nodes), size=int(pP_notVac*len(network.nodes)), replace=False)
    antiVac_nodes = np.random.choice(len(network.nodes), size=int(pP_antiVac*len(network.nodes)), replace=False)
    for node in network.nodes:
        # decide which nodes can not be vaccinated   
        if node in notVac_nodes:
            network.nodes[node]['canVacc'] = False
            nodes_color_timestep0.append('snow')  
            ns += 1        
        else:
            network.nodes[node]['canVacc'] = True
            nodes_color_timestep0.append('silver') 
            
        if node in antiVac_nodes:
            network.nodes[node]['antiVacc'] = True
        else:
            network.nodes[node]['antiVacc'] = False
            
        
        s += 1   
        network.nodes[node]['state'] = 'S' #susceptible   
        network.nodes[node]['wasVacc'] = False         
        network.nodes[node]['neighbors'] = [n for n in network.neighbors(node)]
        network.nodes[node]['bonds'] = [b for b in network.edges(node)]
     

    
    
    nodes_color.append(nodes_color_timestep0)
    # the first time step:
    I = [i]
    S = [s]
    R = [r]
    V = [v]
    D = [d]
    
    NS = [ns]
    ND = [nd]
    NI = [ni]
    NR = [nr]
    

    ############### time step 1 #############

    # node 0 is patient 0
    network.nodes[patient_zero]['state'] = 'I'
    network.nodes[patient_zero]['TimeOfInfection'] = 0
    s -= 1
    i += 1
    if network.nodes[patient_zero]['canVacc'] == False:
        ni += 1
        ns -= 1
    # for timestep 1 the colors are the same as for timestep 0, only difference is that patient 0 is now infected -> gets color firebrick
    nodes_color_timestep1 = nodes_color_timestep0.copy()
    nodes_color_timestep1[patient_zero] = 'firebrick' 
    
    nodes_color.append(nodes_color_timestep1)

    I.append(i)
    R.append(r)
    S.append(s)
    V.append(v)
    D.append(d)
    
    NS.append(ns)
    ND.append(nd)
    NI.append(ni)
    NR.append(nr)

    
    # time steps for simulation
    time_steps = np.arange(0,T)
    # variable for starting vaccination
    start_vac = 0
    t_till_vac = 0 
     
    ############ time step 2 to T ################
    
    # for t in time_steps[2:]:
    # for t in tqdm (np.arange(2,T), desc="Loading"):
    for t in np.arange(2,T):
    
        # need to save the states for the t+1 time step in an different network, such that we update the network only after the end of the current timestep
        network_tplus1 = network.copy(network)

        # this list stores the node colors in an extra list for the current time steps
        node_colors_t = []
        
        # loop over every node and check if an infected person changes its state to recovefirebrick
        # and check if a susceptible person has infected neighbors
        for node in network.nodes:         
               
            if network.nodes[node]['state'] == 'S':
                if network.nodes[node]['canVacc'] == False:
                    node_colors_t.append('snow') 
                else:
                    node_colors_t.append('silver') 
                # check if a susceptible person has infected neighbors
                # count the number of infected neighbors
                Number_infected_neighbors = 0

                for nb_neighbor in range(len(network.nodes[node]['neighbors'])):
                    neighbor = network.nodes[node]['neighbors'][nb_neighbor]
                    neighbor_node_bond = network.nodes[node]['bonds'][nb_neighbor]

                    if network.nodes[neighbor]['state'] == 'I':
                        Number_infected_neighbors += 1
                    
                if np.random.uniform(0,1) < 1-(1-p_trans)**Number_infected_neighbors:
                    # then the node will become infected
                    network_tplus1.nodes[node]['state'] = 'I'
                    network_tplus1.nodes[node]['TimeOfInfection'] = 1
                    # node_colors_t.append('firebrick')
                    i += 1
                    s -= 1
                    
                    if network.nodes[node]['canVacc'] == False:
                        ni += 1
                        ns -= 1
                    
                
            elif network.nodes[node]['state'] == 'R':
                node_colors_t.append('steelblue')
                #check if the 'cooldown' for recovefirebrick nodes still applies
                if network.nodes[node]['TimeTillReinfection'] < Time_till_reinfec:
                    network_tplus1.nodes[node]['TimeTillReinfection'] += 1
                    # node_colors_t.append('steelblue')
                else:
                    # check if there are infected neighbors
                    Number_infected_neighbors = 0

                    for nb_neighbor in range(len(network.nodes[node]['neighbors'])):
                        neighbor = network.nodes[node]['neighbors'][nb_neighbor]
                        neighbor_node_bond = network.nodes[node]['bonds'][nb_neighbor]

                        if network.nodes[neighbor]['state'] == 'I':
                            Number_infected_neighbors += 1
                    
                    # there is a small probability that recovefirebrick nodes can become reinfected        
                    if np.random.uniform(0,1) < 1-(1-p_reinf_reco)**Number_infected_neighbors:
                        network_tplus1.nodes[node]['state'] = 'I'
                        network_tplus1.nodes[node]['TimeOfInfection'] = 1
                        # node_colors_t.append('firebrick')
                        i += 1
                        r -= 1
                        
                        if network.nodes[node]['canVacc'] == False:
                            ni += 1
                            nr -= 1
                        
                
            elif network.nodes[node]['state'] == 'V':
                node_colors_t.append('darkorange')
                # check if there are infected neighbors
                Number_infected_neighbors = 0

                for nb_neighbor in range(len(network.nodes[node]['neighbors'])):
                    neighbor = network.nodes[node]['neighbors'][nb_neighbor]
                    neighbor_node_bond = network.nodes[node]['bonds'][nb_neighbor]

                    if network.nodes[neighbor]['state'] == 'I':
                        Number_infected_neighbors += 1
                
                # there is a small probability that vaccinated nodes can become infected        
                if np.random.uniform(0,1) < 1-(1-p_reinf_vac)**Number_infected_neighbors:
                    network_tplus1.nodes[node]['state'] = 'I'
                    network_tplus1.nodes[node]['TimeOfInfection'] = 1
                    # node_colors_t.append('firebrick')
                    i += 1
                    v -= 1 

            elif network.nodes[node]['state'] == 'D':
                node_colors_t.append('black')
                
                
            # check if nodes have recovefirebrick 
            elif network.nodes[node]['state'] == 'I':
                node_colors_t.append('firebrick')
                if network.nodes[node]['TimeOfInfection'] < Time_till_recovery:
                    network_tplus1.nodes[node]['TimeOfInfection'] += 1
                    # node_colors_t.append('firebrick')
                else:
                    #check if infected node dies
                    if network.nodes[node]['wasVacc'] == True:
                        death_rate_ = vacc_death_rate
                    else:
                        death_rate_ = death_rate
                        
                    if np.random.uniform(0,1) < death_rate_:
                        network_tplus1.nodes[node]['state'] = 'D'
                        d += 1
                        i -= 1
                        # node_colors_t.append('black')
                        
                        if network.nodes[node]['canVacc'] == False:
                            nd += 1
                            ni -= 1
                    else:
                        #if it doesn't die, it recovers
                        network_tplus1.nodes[node]['state'] = 'R'
                        #introduce a 'cooldown' for recovefirebrick nodes before they can become reinfected
                        network_tplus1.nodes[node]['TimeTillReinfection'] = 1
                        r += 1
                        i -= 1
                        # node_colors_t.append('steelblue')
                        
                        if network.nodes[node]['canVacc'] == False:
                            nr += 1
                            ni -= 1
                    
        if (i+r) > vac_thresh*len(network.nodes):
            start_vac = 1 
            
        if start_vac == 1:
            t_till_vac += 1
        
        # if start_vac == 1 and v+r < vac_rate * len(network.nodes) and t_till_vac > t_vac_develop:
        if start_vac == 1 and t_till_vac >= t_vac_develop:
        # if t <=2:
            # print(v)
            # we only add new vaccinated people until the maximal vaccination rate is achieved 
        
            susc_nodes = [ node for node in network_tplus1.nodes if network_tplus1.nodes[node]['state'] == 'S' 
                          and network_tplus1.nodes[node]['canVacc'] == True and network.nodes[node]['antiVacc'] == False] # only susceptible person will be vaccinated 
            
            # need this if loop because otherwise random.choice() gives problems if the sample size is smaller than the size of susc_nodes
            if Num_vac_per_Tstep*len(network.nodes) < len(susc_nodes):
                size = round(Num_vac_per_Tstep*len(network.nodes))
                # print(size)
            else:
                size = len(susc_nodes)
            
            vac_node = np.random.choice(susc_nodes, size=size, replace=False)
            
            for node in vac_node:
                network_tplus1.nodes[node]['state'] = 'V'
                network_tplus1.nodes[node]['wasVacc'] = True
                s -= 1
                v += 1
                    
        
        # save the colors of the current time step         
        nodes_color.append(node_colors_t)   
        
                  

        # update the network
        network = network_tplus1       
        
        # save number of suscptible, infected and recovefirebrick persons in the current time step
        S.append(s)
        I.append(i)
        R.append(r)
        V.append(v)
        D.append(d)
        
        NS.append(ns)
        ND.append(nd)
        NI.append(ni)
        NR.append(nr)
        
        
    return S,I,R,V,D, time_steps, (nodes_color), NS, NI, NR, ND




def SIR_only_bond(network, p_trans, T, Time_till_recovery, patient_zero=0, n_shells=None, n_con=None):
    
    np.random.seed(None)
    # network = Network_bond_perco(network, p_trans)
    
    # we assign a color to every state: I: 'firebrick', S: 'silver', R: 'steelblue' for the animation
    nodes_color = []
    nodes_color_timestep0 = [] #since we start our for loop at t=2 we have to keep track of the colors for the first time step putside of the loop

    # at the beginning every node is susceptible
    for node in network.nodes:
        network.nodes[node]['state'] = 'S'
        network.nodes[node]['neighbors'] = [n for n in network.neighbors(node)]
        network.nodes[node]['bonds'] = [b for b in network.edges(node)]
        nodes_color_timestep0.append('silver')
    
    nodes_color.append(nodes_color_timestep0)

    # node 0 is patient 0
    network.nodes[patient_zero]['state'] = 'I'
    network.nodes[patient_zero]['TimeOfInfection'] = 0
    # for timestep 1 the colors are the same as for timestep 0, only difference is that patient 0 is nor infected -> gets color firebrick
    nodes_color_timestep1 = nodes_color_timestep0.copy()
    nodes_color_timestep1[patient_zero] = 'firebrick' 
    
    nodes_color.append(nodes_color_timestep1)

    # we want to keep track of the number of infected I, susceptible S and recovefirebrick R people
    i = 1
    s = len(network.nodes)-1
    r = 0 
    
    # these are the first two time steps
    I = [0, i]
    S = [s+1, s]
    R = [0, 0]
    
    # time steps for simulation
    time_steps = np.arange(0,T)
     
    # bar = Bar('Processing', max=T) 
    for t in time_steps[2:]:
    # for t in tqdm (np.arange(2,T), desc="Loading", ncols=75):
        # need to save the states for the t+1 time step in an different network, such that we update the network only after the end of the current timestep
        network_tplus1 = network.copy(network)

        # this list stores the node colors in an extra list for the current time steps
        node_colors_t = []
        
        # loop over every node and check if an infected person changes its state to recovefirebrick
        # and check if a susceptible person has infected neighbors
        for node in network.nodes:
               
            # check if a susceptible person has infected neighbors
            if network.nodes[node]['state'] == 'S':
                node_colors_t.append('silver')
                # count the number of infected neighbors
                Number_infected_neighbors = 0

                for nb_neighbor in range(len(network.nodes[node]['neighbors'])):
                    neighbor = network.nodes[node]['neighbors'][nb_neighbor]
                    neighbor_node_bond = network.nodes[node]['bonds'][nb_neighbor]

                    if network.nodes[neighbor]['state'] == 'I' : #and network.edges[neighbor_node_bond]['state'] == 'not_blocked':
                        Number_infected_neighbors += 1
                    
                if np.random.uniform(0,1) <= 1-(1-p_trans)**Number_infected_neighbors:
                # if Number_infected_neighbors>0:
                    # then the node will become infected
                    network_tplus1.nodes[node]['state'] = 'I'
                    network_tplus1.nodes[node]['TimeOfInfection'] = 1
                    # node_colors_t.append('firebrick')
                    i += 1
                    s -= 1
                # else:
                    # the node stays susceptible
                    # node_colors_t.append('silver')
                
            elif network.nodes[node]['state'] == 'R':
                node_colors_t.append('steelblue')
                
            # check if nodes have recovefirebrick 
            elif network.nodes[node]['state'] == 'I':
                node_colors_t.append('firebrick')
                
                if network.nodes[node]['TimeOfInfection'] < Time_till_recovery:
                    network_tplus1.nodes[node]['TimeOfInfection'] += 1
                    # node_colors_t.append('firebrick')
                else:
                    network_tplus1.nodes[node]['state'] = 'R'
                    # network.nodes[node]['state'] = 'R'
                    r += 1
                    i -= 1
                    # node_colors_t.append('steelblue')
        
        # save the colors of the current time step         
        nodes_color.append(node_colors_t)                

        # update the network
        network = network_tplus1       
        
        # save number of suscptible, infected and recovefirebrick persons in the current time step
        S.append(s)
        I.append(i)
        R.append(r)
        
        perco = 0 
        if n_shells != None:
            for node in list(np.arange(total_number_for_nShells(n_shells-1,n_con),total_number_for_nShells(n_shells,n_con))):
                if network.nodes[node]['state'] == 'I' or network.nodes[node]['state'] == 'R':
                    perco += 1
        
    return S,I,R, time_steps, (nodes_color), perco
