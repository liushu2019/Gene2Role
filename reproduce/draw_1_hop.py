import networkx as nx
import pandas as pd
from collections import deque
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
from matplotlib.colors import ListedColormap
colors = ["#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"]
# positions = [0.0, 0.25, 0.5, 1.0]  # Positions for each color in the range [0, 1]
plt.rcParams['font.family'] = 'Ubuntu Mono'
custom_cmap = ListedColormap(colors)

def get_k_hop_structure(edgelist, id, path, k=1):
    df = pd.read_csv(edgelist, sep='\t', header=None)
    edgelist_out = pd.DataFrame()
    queue = deque()
    queue.append(id)
    local_queue = deque()
    left_nodes = 1
    processed_nodes = set()
    bk_k = k
    layers = {x:[] for x in range(0,k+1)}
    layers_p = {x:[] for x in range(0,k+1)}
    layers_n = {x:[] for x in range(0,k+1)}

    while queue:
        left_nodes -= 1
        node = queue.popleft()
        if node in processed_nodes:
            continue
        processed_nodes.add(node)
        layers[k].append(node)

        df_tmp = df[(df[0] == node) | (df[1] == node)]
        df_tmp_p = df_tmp[df_tmp[2] > 0]
        df_tmp_n = df_tmp[df_tmp[2] < 0]
        if k < bk_k:
            layers_p[k] += list((set(df_tmp_p[0]) | set(df_tmp_p[1]))-processed_nodes - set(layers_n[k+1]))
            layers_n[k] += list((set(df_tmp_n[0]) | set(df_tmp_n[1]))-processed_nodes - set(layers_p[k+1]))
        else:
            layers_p[k] += list((set(df_tmp_p[0]) | set(df_tmp_p[1]))-processed_nodes)
            layers_n[k] += list((set(df_tmp_n[0]) | set(df_tmp_n[1]))-processed_nodes)
            
        local_queue.extend((set(df_tmp[0]) | set(df_tmp[1])))
        # edgelist_out = pd.concat([edgelist_out, df_tmp])
        if left_nodes == 0:
            k -= 1
            if k < 0:
                break
            queue.extend(set(local_queue) - processed_nodes)
            left_nodes = len(queue)
            #print (left_nodes, queue)
    g = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=True)
    g = nx.subgraph(g, processed_nodes)
    #print (len(processed_nodes))
    #print (g)
    pos = nx.spring_layout(g, center=(0,0))
    pos[id]=(0,0)
    df_plot = nx.to_pandas_edgelist(g, source=0, target=1 )
    df_plot_other = df_plot[(df_plot[0]!=id) & (df_plot[1]!=id)]
    df_plot_n_neg = df_plot[((df_plot[0]==id) | (df_plot[1]==id)) & (df_plot[2]<0) ]
    df_plot_n_pos = df_plot[((df_plot[0]==id) | (df_plot[1]==id)) & (df_plot[2]>0) ]
    nx.draw_networkx(g, pos, with_labels=False, nodelist=[],#nodelist=list(set(g.nodes) - set([id])),
                node_size=15,
        # font_size=22,
        node_color='skyblue',
        # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
        edgelist=tuple(zip(df_plot_other[0],df_plot_other[1])),
        width=0.2,
        edge_color=[colors[2] if x[2] > 0 else 
                    colors[3] for x in df_plot_other.values],
    )
    nx.draw_networkx(g, pos, with_labels=False, nodelist=[],#nodelist=list(set(g.nodes) - set([id])),
                node_size=15,
        # font_size=22,
        node_color='skyblue',
        # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
        edgelist=tuple(zip(df_plot_n_pos[0],df_plot_n_pos[1])),
        width=2,
        edge_color=[colors[2]],
    )
    nx.draw_networkx(g, pos, with_labels=False, nodelist=[],#nodelist=list(set(g.nodes) - set([id])),
                node_size=15,
        # font_size=22,
        node_color='skyblue',
        # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
        edgelist=tuple(zip(df_plot_n_neg[0],df_plot_n_neg[1])),
        width=2,
        edge_color=[colors[3]],
    )
    # nx.draw_networkx(g, pos, with_labels=False, nodelist=[],#nodelist=list(set(g.nodes) - set([id])),
    #             node_size=15,
    #     # font_size=22,
    #     node_color='skyblue',
    #     # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
    #     width=[1.5 if(x[0] == id or x[1] == id) else 
    #                 0.2 for x in nx.to_pandas_edgelist(g, source=0, target=1 ).values],
    #     edge_color=[colors[2] if x[2] > 0 else 
    #                 colors[3] for x in nx.to_pandas_edgelist(g, source=0, target=1 ).values],
    # )
    # nx.draw_networkx(g, pos, with_labels=False, nodelist=[],#nodelist=list(set(g.nodes) - set([id])),
    #             node_size=15,
    #     # font_size=22,
    #     node_color='skyblue',
    #     # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
    #     width=[1.5 if(x[0] == id or x[1] == id) else 
    #                 0.2 for x in nx.to_pandas_edgelist(g, source=0, target=1 ).values],
    #     edge_color=[colors[2] if x[2] > 0 else 
    #                 colors[3] for x in nx.to_pandas_edgelist(g, source=0, target=1 ).values],
    # )
    # print (layers)

    node_color_list = dict(zip(list(set (df_plot[0] ) | set(df_plot[1])), [4]*(df_plot.shape[0])))
    node_color_list[id]=0
    node_color_list.update(dict(zip(list((set (df_plot_n_pos[0] ) | set(df_plot_n_pos[1])) - set([id])), [2]*(df_plot_n_pos.shape[0]))))
    node_color_list.update(dict(zip(list((set (df_plot_n_neg[0] ) | set(df_plot_n_neg[1])) - set([id])), [3]*(df_plot_n_neg.shape[0]))))

    for nodes in layers.items():
        # print (nodes)
        # print ([1+ (max(layers.keys())-nodes[0])*2 if n in layers_n[max(layers.keys())-nodes[0]] else 0+ (max(layers.keys())-nodes[0])*2 for n in nodes[1]])
        nx.draw_networkx(g, pos, with_labels=False, nodelist=nodes[1], font_size=5,
                    node_size=15+50*nodes[0],
            # font_size=22,
            node_color=[colors[node_color_list[n]] for n in nodes[1]],
            # width=[weight_map[x] if x in weight_map else 1 for x in g.edges ],
            edgelist=[],
            # labels=f'{max(layers.keys())-nodes[0]}-hop'
        )
    plt.gca().set_axis_off()

    legend_labels = {0:{'label':'0-hop', 'color':colors[0]},
                    1:{'label':'1-hop (+)', 'color':colors[2]},
                    2:{'label':'1-hop (-)', 'color':colors[3]},}
    if bk_k >=2:
        legend_labels[3] = {'label':'2-hop (+/-)', 'color':colors[4]}
    for group, info in legend_labels.items():
        plt.scatter([], [], c=info['color'], label=info['label'])
    plt.legend()
    legend_labels = {}
    plt.savefig(path + str(id) + '.pdf',dpi=300, bbox_inches='tight')
    plt.close()
    #plt.savefig(f'sub_{os.path.basename(edgelist).split(".")[0]}_{id}_{bk_k}.pdf',dpi=500, bbox_inches='tight')

def batch_plot_1_hop(edgelist, index_list, path, gene_name, edgelist_list):
    file_path = path + gene_name + '/'
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    for id in index_list:
        if id in edgelist_list:
            get_k_hop_structure(edgelist, id, file_path)
        else:
            print(str(id))


#edgelist = str(sys.argv[1])
#id = int(sys.argv[2])
#k = int(sys.argv[3])
#get_k_hop_structure(edgelist, id, '/Users/xinzeng/Desktop/research/role_singlecell/codes/jupyter-notebook/', k)
