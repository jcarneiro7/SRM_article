import wntr
import networkx as nx
import numpy as np
import pandas as pd

from classes_concept import Graph_Theory_Metrics, Entropy_metrics, Resilience_Indexes


teste='2_1'
if teste=='1_0':
    filename='Networks/concept_1/Concept_1'
elif teste=='1_1':
    filename='Networks/concept_1/Concept_1_1'
elif teste=='1_2':
    filename='Networks/concept_1/Concept_1_2'
elif teste=='1_3':
    filename='Networks/concept_1/Concept_1_3'
elif teste == '1_4':
    filename='Networks/concept_1/Concept_1_4' 
elif teste == '1_5':
    filename='Networks/concept_1/Concept_1_5' 

elif teste=='2_0':
    filename='Networks/concept_2/Concept_2'
elif teste=='2_1':
    filename='Networks/concept_2/Concept_2_1'
elif teste=='2_2':
    filename='Networks/concept_2/Concept_2_2'
elif teste=='2_3':
    filename='Networks/concept_2/Concept_2_3'
elif teste == '2_4':
    filename='Networks/concept_2/Concept_2_4' 
elif teste == '2_5':
    filename='Networks/concept_2/Concept_2_5' 
elif teste=='3_0':
    filename='Networks/concept_3/Concept_3'
elif teste=='3_1':
    filename='Networks/concept_3/Concept_3_1'
elif teste=='3_2':
    filename='Networks/concept_3/Concept_3_2'
elif teste=='3_3':
    filename='Networks/concept_3/Concept_3_3' 
elif teste == '3_4':
    filename='Networks/concept_3/Concept_3_4' 
elif teste == '3_5':
    filename='Networks/concept_3/Concept_3_5' 

inp_file=filename+'.inp'
wn=wntr.network.WaterNetworkModel(inp_file)
sim = wntr.sim.EpanetSimulator(wn)
results = sim.run_sim()
G_top = wn.get_graph()
uG_top = G_top.to_undirected() # undirected multigraph

## Resilience Index Metrics
p_req = 20
p_min = 10

ri = Resilience_Indexes.resilience_index(wn,results,p_req,teste)
nri = Resilience_Indexes.network_resilience_index(wn,results,p_req,uG_top,teste)
mri = Resilience_Indexes.modified_resilience_index(wn,results,p_req)
thri = Resilience_Indexes.target_hydraulic_resilience_index(wn,results,p_req,p_min)
wri = Resilience_Indexes.weighted_resilience_index (wn,results,p_req,teste)
print('finish resilience indexes metrics')

list_of_pipes = wn.pipe_name_list
list_of_nodes = wn.junction_name_list
pressao_original = results.node['pressure'].loc[:,list_of_nodes]

# ## Entropy metrics
flowrate_orig = results.link['flowrate']
list_of_times = list(flowrate_orig.index)

entropia_awumah=[]
entropia_tanyimboh=[]
entropia_dsfe=[]

for t in list_of_times: # loop time
    wn_graph =  wntr.network.WaterNetworkModel(inp_file)
    for p in list_of_pipes:     # loop pipes 
        if flowrate_orig.loc[t,p]<0:
            pipe_orig=wn_graph.get_link(p)
            new_pipe_name=pipe_orig.name
            new_pipe_no_final=pipe_orig.start_node_name
            new_pipe_no_inicial=pipe_orig.end_node_name
            new_pipe_diameter=pipe_orig.diameter
            new_pipe_length=pipe_orig.length
            new_pipe_roughness=pipe_orig.roughness
            wn_graph.remove_link(pipe_orig.name)
            wn_graph.add_pipe(new_pipe_name, start_node_name=new_pipe_no_inicial, end_node_name=new_pipe_no_final,
            length=new_pipe_length, diameter=new_pipe_diameter, roughness=new_pipe_roughness, minor_loss=0)
    wn_graph.write_inpfile('entropia_new_direction_new.inp', version=2.2)
    sim = wntr.sim.EpanetSimulator(wn_graph)
    results = sim.run_sim() 
    flowrate = results.link['flowrate'].loc[t,:]
    demand = results.node['demand'].loc[t,:]
    velocity = results.link['velocity'].loc[t,:]
    G = wn_graph.get_graph(link_weight=flowrate)
    entropy, awumah = Entropy_metrics.entropy_awumah(G)
    tanyimboh = Entropy_metrics.entropy_tanyimboh(G,demand,flowrate)    
    dsfe = Entropy_metrics.diameter_sensitive_flow_entropy (G, demand,flowrate, velocity)

    entropia_awumah.append(awumah)
    entropia_tanyimboh.extend(tanyimboh)
    entropia_dsfe.extend(dsfe)

entropia_awumah=np.array(entropia_awumah)
entropia_tanyimboh=np.array(entropia_tanyimboh)
entropia_dsfe=np.array(entropia_dsfe)
print('finish entropy metrics')

## Graph Theory Metrics
ld = Graph_Theory_Metrics.link_density(uG_top)
cpd = Graph_Theory_Metrics.central_point_dominance(uG_top)
apl = Graph_Theory_Metrics.average_path_length(uG_top)
mc = Graph_Theory_Metrics.meshedness_coefficient(G_top)
ac = Graph_Theory_Metrics.algebraic_connectivity(uG_top)
sg = Graph_Theory_Metrics.spectral_gap(uG_top)
print('finish graph theory')

## save results csv
results = [ri, nri, mri, thri, wri, entropia_awumah[0], entropia_tanyimboh[0], entropia_dsfe[0], ld, cpd, apl, mc, ac, sg]
resultados=pd.DataFrame(results)
resultados.to_csv('Results/res_'+teste+'.csv')
print('finished')

