# import wntr
import networkx as nx
# import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import math
# from collections import Counter

def link_density(uG):
    link_density=nx.density(uG)
    return (link_density)

def central_point_dominance(uG):
    bet_cen = nx.betweenness_centrality(nx.Graph(uG)) # not implemented for multigraph 
    bet_cen = list(bet_cen.values())
    cpd = sum(max(bet_cen) - np.array(bet_cen))/(len(bet_cen)-1)
    return (cpd)

def average_path_length(G):
   apl=nx.average_shortest_path_length(G, weight=None)
   return (apl)

def meshedness_coefficient(G):
    mc = float(G.number_of_edges() - G.number_of_nodes() + 1)/(2*G.number_of_nodes()-5)
    return(mc)

def algebraic_connectivity(uG):
    eig = nx.laplacian_spectrum(uG)
    eig = np.sort(eig)
    ac = eig[1]
    return (ac)

def spectral_gap(uG):
    eig = nx.adjacency_spectrum(uG)
    sg = abs(eig[0] - eig[1])
    return (sg)


