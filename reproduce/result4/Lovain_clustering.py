import networkx as nx
import community as community_louvain
dir_list = {
    'Exp2_1_2_t12_eeisp':'/home/liu/43_Gene2Role/Gene2role/codes/data/singleCell/glioblastoma/splitMatrix/Exp2_1_2_t12_eeisp.edgelist',
    'GMP_0_network':'/home/liu/43_Gene2Role/Gene2role/codes/data/mult-omics/GMP_0_network.edgelist',
}
for project, dir_ in dir_list.items():
    df = pd.read_csv(dir_, sep='\t', header=None)
    df = df[df[2] > 0]
    G = nx.from_pandas_edgelist(df,0,1,edge_attr=None, create_using=nx.Graph)
    connected_components = nx.connected_components(G)
    # Find the largest connected component
    largest_component = max(connected_components, key=len)
    # Create a new graph containing only the largest component
    largest_component_graph = G.subgraph(largest_component)
    connected_components = nx.connected_components(G)
    for x in connected_components:
        print (x)
    # Perform community detection using the Louvain method
    partition = community_louvain.best_partition(largest_component_graph, resolution=0.8)
    print (community_louvain.modularity(partition, largest_component_graph))
    df_out = pd.DataFrame(partition, index=['community']).T
    df_out.insert(0, "gene", df_out.index)
    print (df_out['community'].unique())
    df_out.to_csv(f"Louvain_{project}.csv", index=None, header=True, sep='\t')
    # Print the detected communities
    # print("Communities:")
    # for node, community in partition.items():
    #     print(f"Node {node} belongs to community {community}")