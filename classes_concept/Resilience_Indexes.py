import wntr
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import Counter


def resilience_index(wn,results,p_req,teste):
    node_elevation = wn.query_node_attribute('elevation')
    gamma = 9810
    list_of_reservoirs = wn.reservoir_name_list
    list_of_pumps = wn.pump_name_list
    nodes_with_demand = wn.query_node_attribute('base_demand', np.greater,0)
    list_of_nodes = list(nodes_with_demand.index)

    ri_1 = 0
    ri_res = 0
    ri_pumps = 0
    ri_req = 0

    for r in list_of_reservoirs:
        # ri_res += results.node['head'].loc[:,r] * (abs(results.node['demand'].loc[:,r]))
        ri_res += results.node['head'].loc[:,r] * (-(results.node['demand'].loc[:,r]))

    for p in list_of_pumps:
        ri_pumps += abs(results.link['flowrate'].loc[:,p]) * abs(results.link['headloss'].loc[:,p]) / gamma

    for n in list_of_nodes:
        ri_1 += results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        ri_req += results.node['demand'].loc[:,n] * (p_req+node_elevation.loc[n])
        
        # Grafico
        ri_1_node_aux = results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        ri_1_node=ri_1_node_aux.iloc[0]
        no = wn.get_node(n)
        no.initial_quality=float(ri_1_node)
    
    # # plot distribuicao espacial
    # ax = wntr.graphics.plot_network(wn, node_attribute='initial_quality', node_colorbar_label='Node resilience importance',node_size=40,node_cmap='turbo',node_range=[0, 0.3])
    # ax.set_box_aspect(4)
    # l, b, w, h = ax.get_position().bounds
    # ax.set_position([l-0.2, b, w, h])
    # plt.savefig('Results/RI_'+teste+'.png',dpi=300)
    # plt.show()

    # Metrica Resiliencia
    ri=ri_1 / (ri_res + ri_pumps - ri_req)
    return (ri.loc[0])

def network_resilience_index(wn,results,p_req,uG,teste):
    node_elevation = wn.query_node_attribute('elevation')
    gamma = 9810
    list_of_reservoirs = wn.reservoir_name_list
    list_of_pumps = wn.pump_name_list
    nodes_with_demand = wn.query_node_attribute('base_demand', np.greater,0)
    list_of_nodes = list(nodes_with_demand.index)
    diameters = wn.query_link_attribute('diameter')

    ri_1 = 0
    ri_res = 0
    ri_pumps = 0
    ri_req = 0

    for r in list_of_reservoirs:
        ri_res += results.node['head'].loc[:,r] * (abs(results.node['demand'].loc[:,r]))

    for p in list_of_pumps:
        ri_pumps += abs(results.link['flowrate'].loc[:,p]) * abs(results.link['headloss'].loc[:,p]) / gamma

    for n in list_of_nodes:
        links = list(uG.edges(keys=True))
        links_pd = pd.DataFrame(links)
        node_links = links_pd.iloc[:,0:2]
        idx = []
        idx_aux = node_links.index[node_links.iloc[:,0] == n].tolist()
        idx_aux2 = node_links.index[node_links.iloc[:,1] == n].tolist()
        idx.extend(idx_aux)
        idx.extend(idx_aux2)
        pipes=links_pd.iloc[idx,2]
        diams=[]
        for p in pipes:
            # d = diameters[p]
            # print(d)
            diams.append(diameters[p])
            # print(diams)

        Ui = sum(diams) / (len(diams) * max(diams))

        ri_1 += Ui * results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        ri_req += results.node['demand'].loc[:,n] * (p_req+node_elevation.loc[n])

        # Grafico
        nri_1_node_aux = Ui * results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        nri_1_node=nri_1_node_aux.iloc[0]
        no = wn.get_node(n)
        no.initial_quality=float(nri_1_node)
    
    # # plot distribuicao espacial
    # ax = wntr.graphics.plot_network(wn, node_attribute='initial_quality', node_colorbar_label='Node resilience importance',node_size=40,node_cmap='turbo',node_range=[0, 0.3])
    # ax.set_box_aspect(4)
    # l, b, w, h = ax.get_position().bounds
    # ax.set_position([l-0.2, b, w, h])
    # plt.savefig('Results/NRI_'+teste+'.png',dpi=300)
    # plt.show()

    # Métrica resiliencia
    nri = ri_1 / (ri_res + ri_pumps - ri_req)
    return (nri.loc[0])

def modified_resilience_index(wn,results,p_req):
    node_elevation = wn.query_node_attribute('elevation')
    nodes_with_demand = wn.query_node_attribute('base_demand', np.greater,0)
    list_of_nodes = list(nodes_with_demand.index)

    ri_1 = 0
    ri_req = 0
    for n in list_of_nodes:
        ri_1 += results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        ri_req += results.node['demand'].loc[:,n] * (p_req+node_elevation.loc[n])
        # print(results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n])))
        # print(results.node['demand'].loc[:,n] * (p_req+node_elevation.loc[n]))
        # # print(ri_1)
        # # print(ri_req)
        # print('aqui')

    mri = (ri_1 / ri_req) * 100
    return (mri.loc[0])

def target_hydraulic_resilience_index(wn,results,p_req,p_min): # p_req é o p_target
    nodes_with_demand = wn.query_node_attribute('base_demand', np.greater,0)
    list_of_nodes = list(nodes_with_demand.index)
    thri_aux_1 = 0
    thri_aux_2 = 0
    thri_aux_3 = 0
    for n in list_of_nodes:
        thri_aux_1 += results.node['demand'].loc[:,n] * results.node['pressure'].loc[:,n]
        thri_aux_2 += results.node['demand'].loc[:,n] * p_min
        thri_aux_3 += results.node['demand'].loc[:,n] * p_req 
    # print('QP',thri_aux_1)
    # print('QPtarget',thri_aux_3)
    # print('QPmin',thri_aux_2)
    thri = (thri_aux_1 - thri_aux_2) / (thri_aux_3 - thri_aux_2)
    return(thri.loc[0])

def weighted_resilience_index (wn,results,p_req,teste):
    node_elevation = wn.query_node_attribute('elevation')
    gamma = 9810
    list_of_reservoirs = wn.reservoir_name_list
    list_of_pumps = wn.pump_name_list
    nodes_with_demand = wn.query_node_attribute('base_demand', np.greater,0)
    list_of_nodes = list(nodes_with_demand.index)
    diameters = wn.query_link_attribute('diameter')

    start_nodes = wn.query_link_attribute('start_node_name')
    end_nodes = wn.query_link_attribute('end_node_name')
    
    wri_1 = 0
    wri_res = 0
    wri_pumps = 0
    wri_req = 0
    k_t = 0
    k_i = 0
    k_u = 0
    for r in list_of_reservoirs:
        wri_res += results.node['head'].loc[:,r] * (abs(results.node['demand'].loc[:,r]))

    for p in list_of_pumps:
        wri_pumps += abs(results.link['flowrate'].loc[:,p]) * abs(results.link['headloss'].loc[:,p]) / gamma

    # Auxiliar de calculo do coeficiente de importancia
    flows_tot=[]
    # flowrate=(results.link['flowrate'])
    for n in list_of_nodes:
        pipes_start=list(start_nodes[start_nodes==(n)].index.values)
        pipes_end=list(end_nodes[end_nodes==(n)].index.values)
        flows=[]
        for p in pipes_start: # nó n é start node, só entra água no nó se o flow for negativo
            # print(results.link['flowrate'].loc[0,p])
            if results.link['flowrate'].loc[0,p] < 0:
                flow=abs(results.link['flowrate'].loc[:,p])
                flows.extend(flow)
        for p in pipes_end: # nó n é end node, só entra água no nó se o flow for positivo
            if results.link['flowrate'].loc[0,p] > 0:
                flow=results.link['flowrate'].loc[:,p]
                flows.extend(flow)
        # print('aqui')
        # flows_aux=sum(flows)
        flows_tot.append(sum(flows))

    for n in list_of_nodes:
        # print(n)

        pipes_start=list(start_nodes[start_nodes==(n)].index.values)
        pipes_end=list(end_nodes[end_nodes==(n)].index.values)
        diams=[]
        flows=[]   
        for p in pipes_start: # nó n é start node, só entra água no nó se o flow for negativo

            if results.link['flowrate'].loc[0,p] < 0:
                diams.append(diameters[p])
                flow=abs(results.link['flowrate'].loc[:,p])
                flows.extend(flow)
        for p in pipes_end: # nó n é end node, só entra água no nó se o flow for positivo
            if results.link['flowrate'].loc[0,p] > 0:
                diams.append(diameters[p])
                flow=results.link['flowrate'].loc[:,p]
                flows.extend(flow)

        np_i = len(diams)
        k_t = 0.5 + ((np_i-1)/np_i)
        k_i = sum(flows)/(max(flows_tot))
        squared_diams = [number ** 2 for number in diams]
        k_u = sum(squared_diams) / (min([np_i,2]) * (max(diams)**2))

        wri_1 += k_i * k_t * k_u * results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        wri_req += results.node['demand'].loc[:,n] * (p_req+node_elevation.loc[n])

        # print("k_i: ", k_i)
        # print(n)
        # print("k_t: ", k_t)       
        # print("k_u: ", k_u)
        # print("flow: ", results.node['demand'].loc[:,n])
        # print("Current head: ",results.node['head'].loc[:,n])
        # print("Required head: ",(p_req+node_elevation.loc[n]))
        # print("wri_1 no: ", k_i * k_t * k_u * results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n])))
        # print(wri_req)

        # Para obter a distribuicao espacial da importancia de cada no para a resiliencia
        wri_1_node_aux = k_i * k_t * k_u * results.node['demand'].loc[:,n] * (results.node['head'].loc[:,n]-(p_req+node_elevation.loc[n]))
        wri_1_node=wri_1_node_aux.iloc[0]
        no = wn.get_node(n)
        no.initial_quality=float(wri_1_node)

    wri = (wri_1) / (wri_res + wri_pumps - wri_req)    

    # plot distribuicao espacial
    ax = wntr.graphics.plot_network(wn, node_attribute='initial_quality', node_colorbar_label='Node resilience importance',node_size=40,node_cmap='turbo',node_range=[0, 0.1])
    ax.set_box_aspect(4)
    l, b, w, h = ax.get_position().bounds
    ax.set_position([l-0.2, b, w, h])
    plt.savefig('Results/WRI_'+teste+'.png',dpi=300)

    return (wri.loc[0])