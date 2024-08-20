import networkx as nx
import pandas as pd
from networkx.algorithms import community
from collections import Counter

otu_table = pd.read_csv('Dataframes/all_otus.tsv', sep='\t')
all_otu_list = otu_table.columns.to_list()
all_otu_list = all_otu_list[1:]

# Read in the dropout networks
networks_flashweave = {}

otu_names = []
densities = []
num_nodes_list = []
num_edges_list = []
num_subnetworks_list = []
modularity_data = []

for name in all_otu_list:
    networks_flashweave[name] = []
    file_path = f'Networks/wl_nw_{name}.gml'
    network = nx.read_gml(file_path)
    networks_flashweave[name].append(network)

    # Calculate network properties
    density = nx.density(network)
    num_nodes = nx.number_of_nodes(network)
    num_edges = nx.number_of_edges(network)
    num_subnetworks = nx.number_connected_components(network)

    # Store properties in lists
    otu_names.append(name)
    densities.append(density)
    num_nodes_list.append(num_nodes)
    num_edges_list.append(num_edges)
    num_subnetworks_list.append(num_subnetworks)

    communities = community.greedy_modularity_communities(network)
    modularity_coef = community.modularity(network, communities)
    modularity_data.append({'OTU_Name': name, 'Modularity_Coefficient': modularity_coef})

modularity_df = pd.DataFrame(modularity_data)
modularity_df.to_csv('Results/Dropout_nw_properties_2_modularity.tsv', sep='\t', index=False)

# Create DataFrame from properties
network_properties_df = pd.DataFrame({
    'OTU_Name': otu_names,
    'Density': densities,
    'Num_Nodes': num_nodes_list,
    'Num_Edges': num_edges_list,
    'Num_Subnetworks': num_subnetworks_list
})

# Save the DataFrame to CSV
network_properties_df.to_csv('Results/Dropout_nw_properties_1_degree_density.csv', index=False)

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

all_edges_with_types = []

for dict_network_list in networks_flashweave.values():
    for dict_network in dict_network_list:
        for source, target, attributes in dict_network.edges(data=True):
            weight = attributes.get('weight', 0)  
            if weight > 0:
                interaction_type = 'Positive'
            else:
                interaction_type = 'Negative'
            all_edges_with_types.append(((source, target), interaction_type)) 

unique_edges_with_types = set(all_edges_with_types)
edge_counts = Counter(all_edges_with_types)

df_edge_counts = pd.DataFrame.from_dict(edge_counts, orient='index').reset_index()
df_edge_counts.columns = ['Edge', 'Count']

df_edge_counts[['Source', 'Target']] = pd.DataFrame(df_edge_counts['Edge'].tolist(), index=df_edge_counts.index)
df_edge_counts.to_csv('Results/Dropout_nw_edge_counts.tsv', sep='\t')
