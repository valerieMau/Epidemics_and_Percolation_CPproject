import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm



def Bethe(nShells, z):
    """_summary_

    Args:
        nShells (_type_): number of shells in Bethe lattice
        z (_type_): number of connections per node  
    """
    G = nx.Graph() #generate a Graph
    G.add_nodes_from(np.arange(0,z+1)) #add the nodes of the zeroth and first shell 
    edge_list = [(0,i) for i in range(1,z+1)]
    G.add_edges_from(edge_list,color='black') #add edges from 0 node to every node of first shell

    depth = nShells -2 #because shell 1 and 2 were already generated -> so depth is the depth for the lattice without the first and second shell
    parent_beginning = 1 #we define parent nodes for each shell
    parent_end = len(G.nodes)

    while parent_end <= total_number_for_nShells(nShells-1,z):
        #the total number of nodes is given by z* sum[(from k=0 to k=depth)(2**k)] + z+1  (for the first and second shell) =z* 2*(2**depth -1 ) +z+1
        for i in range(parent_beginning,parent_end):
            # for each parent node in the current shell we add z nodes and corresponding edges to the parent
            parent = i
            Glen_before_adding_nodes = len(G)
            G.add_nodes_from(np.arange(Glen_before_adding_nodes,Glen_before_adding_nodes+z-1))
            edge_list_parent = [(parent,i) for i in range(Glen_before_adding_nodes,Glen_before_adding_nodes+z-1)]
            G.add_edges_from(edge_list_parent,color='black')
        
        #start a new shell -> new parent nodes
        parent_beginning = parent_end 
        parent_end = len(G.nodes)
    
    return G  


def total_number_for_nShells(nshells, number_con):
    """
    Calculate the total number of nodes in a Bethe lattice with nshells and number_con connections between nodes

    Args:
        nshells (_type_): number of shells of bethe lattice 
        number_con (_type_): number of connections from one node
    """
    total_num = 0
    for n in range(1,nshells+1):
        n_= number_nodes_in_given_shell_BetheLattice(n, number_con)
        total_num = total_num + n_

    
    return total_num

def number_nodes_in_given_shell_BetheLattice(nshell, number_con):
    """
    Calculate the  number of nodes in the one shell
    We count the shells from 1 to n

    Args:
        nshells (_type_): number of shells of bethe lattice 
        number_con (_type_): number of connections from one node
    """
    if nshell < 1:
        print("Seid ihr dumm?")
        n_in_this_shell=None
    elif nshell == 1:
        n_in_this_shell=1
    else:
        n_in_this_shell=number_con*(number_con-1)**(nshell-2)

    return n_in_this_shell


def random_connections(network, p):
    #adds edges between random nodes of a network with probability p
    
    for m in tqdm (np.arange(0,len(network.nodes)), desc="Loading", ncols=75):
        for n in range(m+1, len(network.nodes)):
            r=np.random.uniform()
            if r<=p:
                network.add_edge(m,n,color='pink')#grey
    
    return network

def faster_random_connections(network, n_shortcuts):
    np.random.seed(42)
    for i in range(n_shortcuts):
        n12=np.random.choice(len(network.nodes),size=2, replace=False)
        network.add_edge(n12[0],n12[1],color='pink')
    return network


def super_spreader(network, n_spreader, n_edges):
    #adds n_spreader superspreaders, who get n_edges extra connections to random nodes across the grid
    np.random.seed(420)
    spreader=np.random.choice(len(network.nodes),size=n_spreader, replace=False)
    # print("Spreader nodes:",spreader)

    for j in tqdm(range(n_spreader), desc="Loading", ncols=75):
        s_nodes=(np.random.choice(len(network.nodes),size=n_edges, replace=False))
        # print(s_nodes)
        for k in s_nodes:
            if spreader[j]==k:
                continue
            else:
                network.add_edge(spreader[j],k,color='lightgrey')

    return network

# def Bond
        
def combine_lattice(network,n_shortcuts,n_spreader,n_spreader_edges, p_random_connections=0):
    # this gives our normal Bethe Lattice random connections and we can implement super spreader
    
    network_copy=network.copy()
    # network_rand=random_connections(network_copy,p_random_connections)
    network_rand=faster_random_connections(network_copy,n_shortcuts)
    network_copy_2=network_rand.copy()
    network_rand_super=super_spreader(network_copy_2,n_spreader,n_spreader_edges)

    return network_rand_super
     

def critical_threshold(network,n_con,n_shells,pure_Bethe=False):
    #k: degree of the vertex (number of edges connected to a node)
    k=[]
    for node in network.nodes:
        #print(type(node))
        if node in list(np.arange(total_number_for_nShells(n_shells-1,n_con),total_number_for_nShells(n_shells,n_con))):
        # if node in list(np.arange(10,20)):
            k.append(len(network.edges(node))+(n_con-1))
        else:
            k.append(len(network.edges(node)))

    k=np.array(k)

    if pure_Bethe==True:
        k=k[:total_number_for_nShells(n_shells-1,n_con)]
        print("Länge k:",len(k))
        print("total number for nShells-1:",total_number_for_nShells(n_shells-1,n_con))

    k_mean=np.mean(k)
    k_squared_mean=np.mean(k**2)
    print(k,k_mean,k_squared_mean)
    p_c=k_mean/(k_squared_mean-k_mean)
    return p_c
    
    
          
        
        
    
        
        


    
    