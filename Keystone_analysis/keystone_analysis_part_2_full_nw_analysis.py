import networkx as nx
import pandas as pd
from networkx.algorithms import community

otu_table = pd.read_csv('Dataframes/all_otus.tsv', sep='\t')
all_otu_list = otu_table.columns.to_list()
all_otu_list = all_otu_list[1:]

file_path = 'Networks/wl_nw_all_otus.gml'
network = nx.read_gml(file_path)

print('Network is read')

for source, target, attributes in network.edges(data=True):
    weight = attributes.get('weight', 0)  
    if weight > 0:
        attributes['InteractionType'] = 'Positive'
    else:
        attributes['InteractionType'] = 'Negative'

nodes_to_remove = [node for node, data in network.nodes(data=True) if data.get('mv') == 1]
network.remove_nodes_from(nodes_to_remove)

densities = [nx.density(network)]
num_nodes = [nx.number_of_nodes(network)]
num_edges = [nx.number_of_edges(network)]
subnetwork_counts = [nx.number_connected_components(network)]

nw_properties = {node: [degree] for node, degree in network.degree()}

df_nw_properties = pd.DataFrame.from_dict(nw_properties, orient='index', columns=['Degree'])
df_nw_properties.loc['Density'] = densities
df_nw_properties.loc['num_nodes'] = num_nodes
df_nw_properties.loc['num_edges'] = num_edges
df_nw_properties.loc['subnetwork_counts'] = subnetwork_counts

df_nw_properties.to_csv('Results/Full_nw_properties_1_degree_density.tsv', sep='\t')

betweenness = nx.betweenness_centrality(network)
closeness = nx.closeness_centrality(network)

df_nw_centralities = pd.DataFrame({
    'Betweenness': pd.Series(betweenness),
    'Closeness': pd.Series(closeness)
})


df_nw_centralities.to_csv('Results/Full_nw_properties_2_centrality.tsv', sep='\t')

node_communities = community.label_propagation_communities(network)
community_data = []

for i, module in enumerate(node_communities):
    for node in module:
        community_data.append({'Network_Index': 1, 'Module': i+1, 'Node': node})

community_df = pd.DataFrame(community_data)

community_df.to_csv('Results/Full_nw_properties_3_module_numbers.tsv', sep='\t')
