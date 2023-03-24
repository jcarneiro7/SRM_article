# import wntr
import networkx as nx
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
import math
from collections import Counter



def entropy_awumah(G, sources=None, sinks=None):
    if G.is_directed() == False:
        return
    if sources is None:
        sources = [key for key,value in nx.get_node_attributes(G,'type').items() if value == 'Reservoir']
    if sinks is None:
        sinks = G.nodes()
    S = {}
    Q = {}
    for nodej in sinks:
        # print('Nó: ',nodej)
        if nodej in sources:
            S[nodej] = 0 # nodej is the source
            continue

        sp = [] # simple path
        # print('type: ', G.nodes[nodej]['type'])
        if G.nodes[nodej]['type']  == 'Junction':
            for source in sources:
                if nx.has_path(G, source, nodej):
                    simple_paths = nx.all_simple_paths(G,source,target=nodej)
                    sp = sp + ([p for p in simple_paths])
            # print('Simple paths: ', sp)
        if len(sp) == 0:
            S[nodej] = np.nan # nodej is not connected to any sources
            continue
        # print('length simple paths',len(sp))
        # "dtype=object" is needed to create an array from a list of lists with differnet lengths
        sp = np.array(sp, dtype=object)

        # Uj = set of nodes on the upstream ends of links incident on node j
        Uj = G.predecessors(nodej)
        # qij = flow in link from node i to node j
        qij = []
        # aij = number of equivalent independent paths through the link from node i to node j
        aij = []
        for nodei in Uj:
            # print('upstream node: ',nodei)
            mask = np.array([nodei in path for path in sp])
            # NDij = number of paths through the link from node i to node j
            NDij = sum(mask)
            if NDij == 0:
                continue
            temp = sp[mask]
            # MDij = links in the NDij path
            MDij = [(t[idx],t[idx+1]) for t in temp for idx in range(len(t)-1)]

            flow = 0
            for link in G[nodei][nodej].keys():
                flow = flow + G[nodei][nodej][link]['weight']
            qij.append(flow)
            # print('qij: ', qij)
            # dk = degree of link k in MDij
            dk = Counter()
            for elem in MDij:
                # divide by the numnber of links between two nodes
                dk[elem] += 1/len(G[elem[0]][elem[1]].keys())
            # print(dk)
            V = np.array(list(dk.values()))
            aij.append(NDij*(1-float(sum(V - 1))/sum(V)))
            # print('aij: ', aij)

        Q[nodej] = sum(qij) # Total flow into node j

        # Equation 7
        S[nodej] = 0
        for idx in range(len(qij)):
            if Q[nodej] != 0 and qij[idx]/Q[nodej] > 0:
                S[nodej] = S[nodej] - \
                    qij[idx]/Q[nodej]*math.log(qij[idx]/Q[nodej]) + \
                    qij[idx]/Q[nodej]*math.log(aij[idx])
    # print(S)
    Q0 = sum(nx.get_edge_attributes(G, 'weight').values())
    # print(Q0)
    # print(Q)
    # Equation 3
    S_ave = 0
    for nodej in sinks:
        if not np.isnan(S[nodej]):
            if nodej not in sources:
                if Q[nodej]/Q0 > 0:
                    S_ave = S_ave + \
                        (Q[nodej]*S[nodej])/Q0 - \
                        Q[nodej]/Q0*math.log(Q[nodej]/Q0)
                        
    S = pd.Series(S) # convert S to a series
    # print(S_ave)
    return [S, S_ave]

def entropy_tanyimboh(G, demand,flowrate, sources=None, sinks=None):
    if G.is_directed() == False:
        return
    if sources is None:
        sources = [key for key,value in nx.get_node_attributes(G,'type').items() if value == 'Reservoir']
    if sinks is None:
        sinks = G.nodes()
    # print(sinks.values())

    sources_demand=abs(demand.loc[sources])
    T=sum(sources_demand) # total demand
    S_0=0
    # Primeira parte equação: Sources => S_0
    for node_s in sources:
        Q_s = sources_demand.loc[node_s]
        # Equation 2 and 3
        S_0 = S_0 - Q_s/T*math.log(Q_s/T)

    # Segunda parte da equação -> Nodes S_j
    S = {}
    S_j={}
    S_j_1={}
    S_j_2={}

    Q_j={}
    T_j={}

    for nodej in sinks:
        # print(nodej)
        if nodej in sources:
            S[nodej] = 0 # nodej is the source
            continue
        sp = [] # simple path
        if G.nodes[nodej]['type']  == 'Junction':
            for source in sources:
                if nx.has_path(G, source, nodej):
                    simple_paths = nx.all_simple_paths(G,source,target=nodej)
                    sp = sp + ([p for p in simple_paths])
        if len(sp) == 0:
            S_j[nodej] = 0 # nodej is not connected to any sources
            continue

        # Demand of node j
        Q_j[nodej] = demand.loc[nodej]

        # Flow into node j (chega)
        Uj = G.predecessors(nodej)
        T_j_aux = []  
        for nodei in Uj:
            mask = np.array([nodei in path for path in sp])
            # NDij = number of paths through the link from node i to node j
            NDij = sum(mask)
            if NDij == 0:
                continue
            flow = 0
            for link in G[nodei][nodej].keys():
                flow = flow + flowrate.loc[link]
            T_j_aux.append(flow)
        T_j[nodej]=sum(T_j_aux)
        if T_j[nodej]<0:
            S_j[nodej]=0
            continue


        # Primeira parte da equação S_j
        S_j_1[nodej]=0
        if Q_j[nodej]>0 and Q_j[nodej]/T_j[nodej] > 0:
            S_j_1[nodej]=Q_j[nodej]/T_j[nodej]*math.log(Q_j[nodej]/T_j[nodej])

        # Flow emanated from node j
        Dj=G.successors(nodej)
        q_jk=[]
        for nodek in Dj:
            mask = np.array([nodej in path for path in sp])
            # NDij = number of paths through the link from node i to node j
            NDij = sum(mask)
            if NDij == 0:
                continue            
            flow=0
            for link in G[nodej][nodek].keys():
                flow = flowrate.loc[link] ### REVER!!!!!!!!!!!!!!!! acho que não deve ser a somar!!! -> antigo flow = flow + flowrate.loc[link]+*            q_jk.append(flow)            
            q_jk.append(flow)
        # Segunda parte da equação S_j
        S_j_2[nodej] = 0
        for idx in range(len(q_jk)):
            if q_jk[idx] > 0 and q_jk[idx]/T_j[nodej] > 0:
                S_j_2[nodej] = S_j_2[nodej] + \
                    q_jk[idx]/T_j[nodej]*math.log(q_jk[idx]/T_j[nodej])

        S_j[nodej]=T_j[nodej] * (S_j_1[nodej] + S_j_2[nodej])
        if  np.isnan(S_j[nodej]):
            S_j[nodej]=0
    S = S_0 - 1/T * sum(S_j.values())
    return [S]    

def diameter_sensitive_flow_entropy (G, demand,flowrate, velocity, sources=None, sinks=None):
    c=0.1 # velocity constant 

    if G.is_directed() == False:
        return
    if sources is None:
        sources = [key for key,value in nx.get_node_attributes(G,'type').items() if value == 'Reservoir']
    if sinks is None:
        sinks = G.nodes()
    # print(sinks.values())

    sources_demand=abs(demand.loc[sources])
    T=sum(sources_demand)
    S_0=0
    # Primeira parte equação -> Sources = S_0
    for nodei in sources:
        Q_i = sources_demand.loc[nodei]
        # Equation 2 and 3
        S_0 = S_0 - Q_i/T*math.log(Q_i/T)

    # Segunda parte da equação -> Nodes S_j
    S = {}
    S_j={}
    S_j_1={}
    S_j_2={}

    Q_j={}
    T_j={}

    for nodej in sinks:
        # print(nodej)
        if nodej in sources:
            S[nodej] = 0 # nodej is the source
            continue
        sp = [] # simple path
        if G.nodes[nodej]['type']  == 'Junction':
            for source in sources:
                if nx.has_path(G, source, nodej):
                    simple_paths = nx.all_simple_paths(G,source,target=nodej)
                    sp = sp + ([p for p in simple_paths])
        if len(sp) == 0:
            S_j[nodej] = 0 # nodej is not connected to any sources
            continue

        # Demand of node j
        Q_j[nodej] = demand.loc[nodej]
        # Flow into node j (chega)
        Uj = G.predecessors(nodej)
        T_j_aux = []  
        for nodei in Uj:
            mask = np.array([nodei in path for path in sp])
            # NDij = number of paths through the link from node i to node j
            NDij = sum(mask)
            if NDij == 0:
                continue
            flow = 0
            for link in G[nodei][nodej].keys():
                flow = flow + flowrate.loc[link]
            T_j_aux.append(flow)
        T_j[nodej]=sum(T_j_aux)
        if T_j[nodej]<0:
            S_j[nodej]=0
            continue


        # Primeira parte da equação S_j
        S_j_1[nodej]=0
        if Q_j[nodej]>0 and Q_j[nodej]/T_j[nodej] > 0:
            S_j_1[nodej]=Q_j[nodej]/T_j[nodej]*math.log(Q_j[nodej]/T_j[nodej])

        # Flow emanated from node j
        Dj=G.successors(nodej)
        q_jk=[]
        v_jk=[]
        for nodek in Dj:
            mask = np.array([nodej in path for path in sp])
            # NDij = number of paths through the link from node i to node j
            NDij = sum(mask)
            if NDij == 0:
                continue            
            flow=0
            for link in G[nodej][nodek].keys():
                flow = flowrate.loc[link] 
                vel = velocity.loc[link]
            q_jk.append(flow)  
            v_jk.append(vel)            

        # Segunda parte da equação S_j
        S_j_2[nodej] = 0
        for idx in range(len(q_jk)):
            if q_jk[idx] > 0.00001 and q_jk[idx]/T_j[nodej] > 0: # q_jk[idx] > 0.00001 para remover python errors
                S_j_2[nodej] = S_j_2[nodej] + \
                    (c/v_jk[idx]) * \
                    (q_jk[idx]/T_j[nodej]*math.log(q_jk[idx]/T_j[nodej]))
            
                # print(c/v_jk[idx])
                # print(q_jk[idx]/T_j[nodej]*math.log(q_jk[idx]/T_j[nodej]))

        # print('S_j_1: ', S_j_1[nodej])
        # print(nodej)
        # print('vel',v_jk)
        # print('S_j_2: ',S_j_2[nodej])
        S_j[nodej]=T_j[nodej] * (S_j_1[nodej] + S_j_2[nodej])

        if  np.isnan(S_j[nodej]):
            S_j[nodej]=0

    dsfe = S_0 - 1/T * sum(S_j.values())
    # print('S_0: ', S_0)
    # print('S_j: ',1/T * sum(S_j.values()))
    # print('DSFE: ',dsfe)
    return [dsfe]    