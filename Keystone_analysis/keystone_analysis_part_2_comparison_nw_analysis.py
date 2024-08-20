import networkx as nx
import pandas as pd
from networkx.algorithms import community

otu_table = pd.read_csv('Dataframes/all_otus.tsv', sep='\t')
all_otu_list = otu_table.columns.to_list()
all_otu_list = all_otu_list[1:]

# Read in the dropout networks
networks_flashweave = {}

for name in all_otu_list:
    networks_flashweave[name] = []
    file_path = f'Networks/wl_nw_{name}.gml'
    network = nx.read_gml(file_path)
    networks_flashweave[name].append(network)

print('Dropout networks are read')

for dict_network_name, dict_network_list in networks_flashweave.items():
    for dict_network in dict_network_list:
        for source, target, attributes in dict_network.edges(data=True):
            weight = attributes.get('weight', 0)  
            if weight > 0:
                attributes['InteractionType'] = 'Positive'
            else:
                attributes['InteractionType'] = 'Negative'
        nodes_to_remove = [node for node, data in dict_network.nodes(data=True) if data.get('mv') == 1]
        dict_network.remove_nodes_from(nodes_to_remove)

print('Dict Interaction Types are changed')

# Read in the full networks
file_path = f'Networks/wl_nw_all_otus.gml'
full_flashweave = nx.read_gml(file_path)
 
print('Full networks are read')


for source, target, attributes in full_flashweave.edges(data=True):
    weight = attributes.get('weight', 0)  
    if weight > 0:
        attributes['InteractionType'] = 'Positive'
    else:
        attributes['InteractionType'] = 'Negative'
nodes_to_remove = [node for node, data in dict_network.nodes(data=True) if data.get('mv') == 1]
full_flashweave.remove_nodes_from(nodes_to_remove)

print('Full NW Interaction Types are changed')#

# Compare the full networks
comparison_results = []

for dict_network_name, dict_network_list in networks_flashweave.items():
    for dict_network in dict_network_list:
    # Compare Edge numbers
        all_edges = full_flashweave.number_of_edges()
        degree_dropout = full_flashweave.degree(dict_network_name)
        theoretical_edges = all_edges - degree_dropout
        edge_num_dropout = dict_network.number_of_edges()

        # Identical Edges without Edge Attributes
        shared_edges = len(set(dict_network.edges()).intersection(set(full_flashweave.edges())))

        # Identical Edges with Edge attributes
        edges1 = set((u, v, data['InteractionType']) for u, v, data in full_flashweave.edges(data=True))
        edges2 = set((u, v, data['InteractionType']) for u, v, data in dict_network.edges(data=True))
        
        identical_edges = 0
        for edge1 in edges1:
            if edge1 in edges2:
                identical_edges += 1

        comparison_results.append({
            'Dict_Network': dict_network_name,
            'List_Network': 'full_network',
            'Theoretical_edge_num': theoretical_edges,
            'Real_edge_num': edge_num_dropout,
            'Edge_num_difference': theoretical_edges - edge_num_dropout,
            'Shared_Edges': shared_edges,
            'Shared_Edges_with_attributes' : identical_edges
        })

        print(f'{dict_network_name} appended')

comparison_df = pd.DataFrame(comparison_results)

comparison_df.to_csv('Results/Comparison_full_and_dropout_nws.tsv', sep='\t')