import networkx as nx
import pandas as pd
import random
from multiprocessing import Pool

# Function to generate a random network
def generate_random_network(all_otu_list, pos_percentage, num_edges):
    random_network = nx.Graph()
    random_network.add_nodes_from(all_otu_list)

    edge_set = set()
    while len(edge_set) < num_edges:
        node1, node2 = random.sample(all_otu_list, 2)
        if (node1, node2) not in edge_set and (node2, node1) not in edge_set:
            interaction_type = 'Positive' if random.random() < pos_percentage else 'Negative'
            edge_set.add((node1, node2, (('InteractionType', interaction_type),)))

    random_network.add_edges_from((u, v, dict(attr)) for u, v, attr in edge_set)
    return random_network

# Function to generate random networks with dropouts
def generate_dropout_networks(params):
    node_to_remove, all_otu_list, pos_percentage, density, num_networks = params
    node_list_with_dropout = all_otu_list.copy()
    node_list_with_dropout.remove(node_to_remove)

    num_nodes_without = len(node_list_with_dropout)
    total_possible_edges = (num_nodes_without * (num_nodes_without - 1)) / 2
    num_edges = int(total_possible_edges * density)

    dropout_networks = [generate_random_network(node_list_with_dropout, pos_percentage, num_edges) for _ in range(num_networks)]
    return node_to_remove, dropout_networks

def compare_networks(params):
    node_to_remove, dropout_networks, random_full_networks = params
    comparison_results = []
    for dropout_network in dropout_networks:
        for full_network in random_full_networks:
            all_edges = full_network.number_of_edges()
            degree_dropout = full_network.degree[node_to_remove]
            theoretical_edges = all_edges - degree_dropout
            edge_num_dropout = dropout_network.number_of_edges()

            betweenness_centrality = nx.betweenness_centrality(full_network).get(node_to_remove, 0)
            closeness_centrality = nx.closeness_centrality(full_network).get(node_to_remove, 0)

            shared_edges = len(set(dropout_network.edges()).intersection(set(full_network.edges())))
            edges_full = set((u, v, tuple(data.items())) for u, v, data in full_network.edges(data=True))
            edges_dropout = set((u, v, tuple(data.items())) for u, v, data in dropout_network.edges(data=True))
            identical_edges_with_attributes = len(edges_full.intersection(edges_dropout))

            subnetwork_count = nx.number_connected_components(dropout_network)

            comparison_results.append({
                'Dict_Network': node_to_remove,
                'Theoretical_edge_num': theoretical_edges,
                'Degree': degree_dropout,
                'Betweenness_C': betweenness_centrality,
                'Closeness_C': closeness_centrality,
                'Real_edge_num': edge_num_dropout,
                'Edge_num_difference': theoretical_edges - edge_num_dropout,
                'Shared_Edges': shared_edges,
                'Shared_Edges_with_attributes': identical_edges_with_attributes,
                'Num_Subnetwork': subnetwork_count
            })
    return comparison_results

# Read OTU table
otu_table = pd.read_csv('Dataframes/all_otus.tsv', sep='\t')
all_otu_list = otu_table.columns.to_list()[1:]

# Read in the full network
network = nx.read_gml('Networks/wl_nw_all_otus.gml')
for source, target, attributes in network.edges(data=True):
    weight = attributes.get('weight', 0)
    attributes['InteractionType'] = 'Positive' if weight > 0 else 'Negative'

# Remove nodes with 'mv' attribute set to 1
nodes_to_remove = [node for node, data in network.nodes(data=True) if data.get('mv') == 1]
network.remove_nodes_from(nodes_to_remove)

# Set Parameters
density = nx.density(network)
positive_edges = sum(1 for _, _, data in network.edges(data=True) if data.get('InteractionType') == 'Positive')
total_edges = network.number_of_edges()
pos_percentage = positive_edges / total_edges
num_nodes = len(all_otu_list)
total_possible_edges = (num_nodes * (num_nodes - 1)) / 2
num_edges = int(total_possible_edges * density)
num_networks = 100

# Generate random full networks
random_full_networks = [generate_random_network(all_otu_list, pos_percentage, num_edges) for _ in range(num_networks)]

# Generate dropout networks using multiprocessing
with Pool() as pool:
    dropout_params = [(node, all_otu_list, pos_percentage, density, num_networks) for node in all_otu_list]
    results = pool.map(generate_dropout_networks, dropout_params)
    random_nw_with_dropouts = {node: networks for node, networks in results}

# Compare networks using multiprocessing
with Pool() as pool:
    compare_params = [(node, networks, random_full_networks) for node, networks in random_nw_with_dropouts.items()]
    comparison_results = pool.map(compare_networks, compare_params)

comparison_df = pd.DataFrame([item for sublist in comparison_results for item in sublist])
comparison_df.to_csv('Results/Result_comparison_random_networks.tsv', sep='\t', index=False)

