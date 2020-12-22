from network_models import sir_network, seir_network

import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import pandas as pd
import random

# TODO: must change models to jax + numpyro

def run_simulation(graph, model, num_sims=10, t_sim=300):
    df_sims = []
    for iter_i in tqdm(range(num_sims)):
        (df_graph, df_states) = model(graph, len(graph.nodes), t_sim, num_infected_init=None)
        df_states['sim_id'] = iter_i
        df_sims.append(df_states)
    df_all = pd.concat(df_sims)

    return (df_all, df_graph)



graph = nx.scale_free_graph(500).to_undirected()

(df_dynamics_sir, df_graph_sir) = run_simulation(graph, sir_network, num_sims=10, t_sim=100)
(df_dynamics_seir, df_graph_seir) = run_simulation(graph, seir_network, num_sims=10, t_sim=100)

df_all_m = df_all_SIR.copy().reset_index()
df_all_m = df_all_m.groupby(['index']).mean()


ax = df_all_m.plot(y='I')
df_all_m.plot(ax=ax, y='R')




graph = nx.scale_free_graph(500).to_undirected()
pos = nx.spring_layout(graph, iterations=100) #,prog='twopi',args='')
nx.draw( G=graph, pos=pos,
        node_size=5,
        node_color= 'black',#'#e35340',
        edge_color='gray',
        width=.05,
        edge_cmap=plt.cm.Blues, with_labels=False)
plt.show()



plt.show()
