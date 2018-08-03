import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import pickle
import matplotlib.colors as colors
import matplotlib.cm as cmx
import json


if __name__ == "__main__":

    num_args = len(sys.argv)

    threshold = 90
    plot_heatmaps = False
    num_connected_components_displayed = 5
    # num_connected_components_displayed = 22
    # num_connected_components_displayed = 40

    if num_args > 1:
        threshold = int(sys.argv[1])
    if num_args > 2:
        num_connected_components_displayed = int(sys.argv[2])
    if num_args > 3:
        plot_heatmaps = bool(sys.argv[3])

    data_dir = "data_top10_motifs"
    # data_dir = "mini_data_top10_motifs"

    print("opening reactions_bag.p...")
    reactions_bag = pickle.load(open(data_dir + "/reactions_bag.p", "rb"))
    print("...done")
    print("opening heatmaps.p...")
    heatmaps = pickle.load(open(data_dir + "/heatmaps.p", "rb"))
    print("...done")
    print("opening column_names.p...")
    column_names = pickle.load(open(data_dir + "/column_names.p", "rb"))
    print("...done")
    print("opening motif_ids_and_labels.json...")
    with open('motif_ids_and_labels.json') as f:
        motif_ids_and_labels = json.load(f)
    print("...done")

    # motif_ids_and_labels_dict for use in displaying legend
    motif_ids_and_labels_dict = {}
    for result in motif_ids_and_labels["results"]["bindings"]:
        motif_primary_id = result["MotifPrimaryId"]["value"]
        motif_label = result["MotifLabel"]["value"]
        motif_ids_and_labels_dict[motif_primary_id] = motif_label

    reactions_heatmap_items = heatmaps[0]
    motifs_heatmap_items = heatmaps[1]
    glycan_distance_cols = column_names[0]
    motif_distance_cols = column_names[1]

    reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                   orient='index', columns=glycan_distance_cols)
    motif_distance_df = pd.DataFrame.from_dict(motifs_heatmap_items,
                                                   orient='index', columns=motif_distance_cols)

    # Plot heatmaps
    if plot_heatmaps:
        print("plotting heatmaps...")
        sns.set()
        plt.subplot(1,2,1)

        ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True)
        plt.title('pairwise glycan similarity by reactions')

        plt.subplot(1,2,2)
        ax_motifs = sns.heatmap(motif_distance_df, annot=False, fmt=".2f", square=True)
        plt.title('pairwise glycan similarity by motifs')

        plt.show()
        print("...done")

    # Apply threshold
    # reaction_distance_array = reaction_distance_df.values
    # # thresholded_reactions_bool = reaction_distance_df.values > threshold
    #
    # reaction_distance_df.values = np.argwhere(reaction_distance_array > threshold)

    # print("")
    # print("reaction_distance_df:")
    # print(reaction_distance_df)
    # get indices of values in distance dataframe above a threshold
    thresholded_reactions = []
    print("applying threshold to glycan distances...")
    for col in reaction_distance_df.columns:
        # print("col:", col)
        # print("reaction_distance_df[col]:", reaction_distance_df[col])
        col_array = np.zeros((1, len(reaction_distance_df[col])))
        i = 0
        for item in reaction_distance_df[col]:
            # print("item:", item)
            if item > threshold:
                col_array[0, i] = 1
            i += 1

        thresholded_reactions.append((col, col_array[0][:]))
        # threshold_indices.append((col, reaction_distance_df[reaction_distance_df[col] > threshold].index.tolist()))
        # print("start")
        # print([reaction_distance_df[col]]==reaction_distance_df[reaction_distance_df[col] > threshold].index.tolist())
        # print("end")
    print("...done")

    # print("thresholded_reactions:")
    # print(thresholded_reactions)

    cols = [x[0] for x in thresholded_reactions]
    # print("cols:", cols)
    thresholded_reactions_df = pd.DataFrame.from_items(thresholded_reactions,
                                                       orient='index', columns=cols)

    # print("reaction_distance_df.values:")
    # print(reaction_distance_df.values)
    #
    # thresholded_reactions_bool = reaction_distance_df.values > threshold
    #
    # print("thresholded_reactions_bool:")
    # print(thresholded_reactions_bool)
    # print("reaction_distance_df.keys():")
    # print(reaction_distance_df.keys())
    #
    # print("reactions_bag:")
    # print(reactions_bag)


    print("creating network graph from thresholded reactions...")
    # G = nx.from_pandas_adjacency(reaction_distance_df)
    G = nx.from_pandas_adjacency(thresholded_reactions_df)
    print("...done")

    print("getting connected_components from graph...")
    # len_c = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    connected_components = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    print("connected_components:", connected_components)
    print("...done")

    # print("removing nodes with zero edges...")
    # zero_edge_nodes = []
    # for n, nbrs in G.adjacency():
    #     # print(n)
    #     # print(nbrs)
    #     # print(len(nbrs))
    #     if len(nbrs) == 1:
    #         # print("zero edges")
    #         zero_edge_nodes.append(n)
    #     # for nbr, eattr in nbrs.items():
    #     #     print(nbr)
    #     #     print(eattr)
    #         # data = eattr['weight']
    #         # if data < 0.5: print('(%d, %d, %.3f)' % (n, nbr, data))
    #
    # # node_degrees = nx.degree(G)
    # # for node in node_degrees:
    # #     if node[1] == 0:
    # #         G.remove_node(node[0])
    #
    # G.remove_edges_from(zero_edge_nodes)
    # print("...done")

    motifs = []
    for reaction in reactions_bag:
        sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
        if sorted_motifs not in motifs:
            motifs.append(sorted_motifs)

    # print("motifs:", motifs)
    print("len(motifs):", len(motifs))

    node_colors = range(len(motifs))
    print("creating graph layout...")
    # pos = nx.spring_layout(G)
    # pos = nx.circular_layout(G)
    # pos = nx.shell_layout(G)
    pos = nx.fruchterman_reingold_layout(G) # similar to spring
    # pos = nx.kamada_kawai_layout(G) # might work
    # pos = nx.spectral_layout(G) # everything appears on a dense horizontal line, not much use
    print("...done")

    print("applying colours to nodes in network...")

    tab20b = plt.cm.tab20b
    cNorm = colors.Normalize(vmin=node_colors[0], vmax=node_colors[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=tab20b)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(1, 1, 1)

    for label in node_colors:
        motif_labels = []
        for motif_id in motifs[label]:
            motif_labels.append(motif_ids_and_labels_dict[motif_id])
        ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=',||| '.join(motif_labels))
        # ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label='abc')

    # node_label_dict = {}
    for node_color in node_colors:
        nodelist = []
        for reaction in reactions_bag:
            sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
            if sorted_motifs == motifs[node_color]:
                nodelist.append(reaction)
                # G.nodes[reaction]['motifs'] = sorted_motifs
                # node_label_dict[reaction] = '-'.join(sorted_motifs)
                # print(G.nodes[reaction])

        nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                               node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                               nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)
        # nx.draw_networkx_labels(G, pos, nodelist=nodelist, label="test "+str(ii), font_size=10)

    # nx.draw_networkx_labels(G, pos, nodelist=nodelist, font_size=10, labels=node_label_dict, ax=ax)
    print("...done")

    print("drawing network graph...")
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
    plt.legend(loc='upper right', fontsize="xx-small")
    # plt.axis('off')
    # plt.show()
    print("...done")

    print("drawing connected components...")
    fig = plt.figure(2)
    ax = fig.add_subplot(1, 1, 1)
    motifs = []

    print("Motifs for top", num_connected_components_displayed, "connected components:")
    for i in range(0, num_connected_components_displayed):
        motifs_for_conn_comp = []
        for reaction in connected_components[i]:
            sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
            if sorted_motifs not in motifs:
                motifs.append(sorted_motifs)
            if sorted_motifs not in motifs_for_conn_comp:
                motifs_for_conn_comp.append(sorted_motifs)
            # motifs_for_conn_comp = motifs_for_conn_comp + sorted_motifs
        # print(i+1,"^", ', '.join([motif_ids_and_labels_dict[x] for x in motifs_for_conn_comp]))
        for sorted_mts in motifs_for_conn_comp:
            print(i+1,"^", ', '.join([motif_ids_and_labels_dict[x] for x in sorted_mts]))

    print("(connected components) len(motifs):", len(motifs))
    node_colors = range(len(motifs))

    tab20b = plt.cm.tab20b
    cNorm = colors.Normalize(vmin=node_colors[0], vmax=node_colors[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=tab20b)

    for label in node_colors:
        motif_labels = []
        for motif_id in motifs[label]:
            motif_labels.append(motif_ids_and_labels_dict[motif_id])
        ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=',||| '.join(motif_labels))

    motif_set_counts = dict.fromkeys(node_colors, 0)

    for i in range(0, num_connected_components_displayed):

        # fig = plt.figure(i+2)
        # ax = fig.add_subplot(6, 4, i+1)

        for node_color in node_colors:
            nodelist = []
            for reaction in connected_components[i]:
                sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
                if sorted_motifs == motifs[node_color]:
                    nodelist.append(reaction)

            if len(nodelist) > 0:
                motif_set_counts[node_color] += 1

            # nx.draw_networkx_nodes(G, pos, with_labels=False, font_weight='bold',
            #                    nodelist=nodelist, ax=ax)
            nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                                   node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                                   nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)

            nx.draw_networkx_edges(G, pos, edgelist=G.edges(nbunch=nodelist), width=1.0, alpha=0.5, ax=ax)

        plt.legend(loc='upper right', fontsize="xx-small")
        # plt.axis('off')
        # plt.show()

    print("...done")
    motif_set_counts_sorted_by_value = dict(sorted(motif_set_counts.items(),
                                                   key=lambda kv: kv[1],
                                                   reverse=True))
    print("motif_set_counts:", motif_set_counts_sorted_by_value)
    for motif_set_idx in motif_set_counts_sorted_by_value:
        print(', '.join([motif_ids_and_labels_dict[x] for x in motifs[motif_set_idx]]), ':', motif_set_counts_sorted_by_value[motif_set_idx])
    plt.show()
