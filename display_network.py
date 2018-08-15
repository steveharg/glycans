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
    num_connected_components_displayed = 5
    plot_heatmaps = False
    # num_connected_components_displayed = 22
    # num_connected_components_displayed = 40

    # display_motifs = [['G00055MO'],
    #                     ['G00055MO', 'G00057MO'],
    #                     ['G00055MO', 'G00068MO'],
    #                     ['G00046MO', 'G00055MO'],
    #                     ['G00056MO'],
    #                     ['G00026MO', 'G00055MO'],
    #                     ['G00026MO'],
    #                     ['G00046MO', 'G00055MO', 'G00068MO'],
    #                     ['G00055MO', 'G00056MO'],
    #                     ['G00057MO'],
    #                     ['G00027MO', 'G00055MO'],
    #                     ['G00056MO', 'G00068MO'],
    #                     ['G00026MO', 'G00046MO', 'G00055MO'],
    #                     ['G00026MO', 'G00046MO', 'G00055MO', 'G00057MO'],
    #                     ['G00026MO', 'G00055MO', 'G00068MO'],
    #                     ['G00026MO', 'G00056MO'],
    #                     ['G00026MO', 'G00057MO'],
    #                     ['G00026MO', 'G00055MO', 'G00056MO'],
    #                     ['G00026MO', 'G00055MO', 'G00057MO'],
    #                     ['G00027MO', 'G00046MO', 'G00055MO'],
    #                     ['G00056MO', 'G00057MO']]  # display all motif sets

    # similar?
    # display_motifs = [['G00056MO'], ['G00057MO']]  # 2 & 4 - Neo Lactosamine & LacDiNAc
    # display_motifs = [['G00056MO'], ['G00055MO']]  # 2 & 3 - Neo Lactosamine & Lactosamine motif
    display_motifs = [['G00046MO', 'G00055MO', 'G00068MO'], ['G00055MO', 'G00068MO']]  # 16 & 10 - Galalpha1-3Gal epitope, Lactosamine motif, Blood group H & Lactosamine motif, Blood group H

    # re-check
    # display_motifs = [['G00056MO'], ['G00055MO', 'G00056MO']]  # 2 & 12 - Neo Lactosamine & Lactosamine motif, Neo Lactosamine
    # display_motifs = [['G00026MO', 'G00055MO', 'G00056MO'], ['G00026MO', 'G00055MO', 'G00057MO']]  # 18 & 19 - N-Glycan core basic, Lactosamine motif, Neo Lactosamine & N-Glycan core basic, Lactosamine motif, LacDiNAc
    # display_motifs = [['G00026MO', 'G00057MO'], ['G00026MO', 'G00055MO']]  # 5 & 7 - N-Glycan core basic, LacDiNAc & N-Glycan core basic, Lactosamine motif
    # display_motifs = [['G00026MO', 'G00055MO'], ['G00027MO', 'G00055MO']]  # 7 & 8 - N-Glycan core basic, Lactosamine motif & N-Glycan truncated motif. First GlcpNAC cut off, Lactosamine motif
    # display_motifs = [['G00046MO', 'G00055MO'], ['G00055MO']]  # 9 & 3 - Galalpha1-3Gal epitope, Lactosamine motif & Lactosamine motif
    # display_motifs = [['G00057MO'], ['G00055MO', 'G00057MO']]  # 4 & 13 LacDiNAc & Lactosamine motif, LacDiNAc

    # dissimilar?
    # display_motifs = [['G00026MO'], ['G00055MO']]  # 1 & 3 N-Glycan core basic & Lactosamine motif
    # display_motifs = [['G00057MO'], ['G00026MO', 'G00057MO']]  # 4 & 5 – LacDiNAc & N-Glycan core basic, LacDiNAc

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
    print("opening G.p...")
    G = pickle.load(open(data_dir + "/G.p", "rb"))
    print("...done")
    print("opening pos" + str(threshold) + ".p...")
    pos = pickle.load(open(data_dir + "/pos_" + str(threshold) + ".p", "rb"))
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
        ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=' & '.join(motif_labels))
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
    leg = plt.legend(loc='upper right', fontsize="x-small")
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    plt.axis('off')
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
        if motifs[label] in display_motifs:
            ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=' & '.join(motif_labels))
        # else:
        #     ax.plot([0], [0], color='Gray', label=' & '.join(motif_labels))
    ax.plot([0], [0], color='Gray', label='All other motif sets')


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

            # color this set of nodes in gray first, then re-color if we've specified them in the
            # display_motifs list
            nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                                   node_color='gray',
                                   nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)
            print(motifs[node_color])
            if motifs[node_color] in display_motifs:
                nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                                       node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                                       nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)
            # else:
            #     nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
            #                            node_color='gray',
            #                            nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)


            nx.draw_networkx_edges(G, pos, edgelist=G.edges(nbunch=nodelist), width=1.0, alpha=0.5, ax=ax)

        # leg = plt.legend(loc='upper right', fontsize="x-small")
        leg = plt.legend(fontsize="x-small")
        for line in leg.get_lines():
            line.set_linewidth(4.0)

        # Construct figure title
        plot_title = ''
        for motif_id_list in display_motifs:
            motif_names = []
            for motif_id in motif_id_list:
                motif_names.append(motif_ids_and_labels_dict[motif_id])
            plot_title += ', '.join(motif_names)
            plot_title += ' & '
        plot_title = plot_title[0:-3]
        plt.title(plot_title)
        plt.axis('off')
        # plt.show()

    print("...done")
    motif_set_counts_sorted_by_value = dict(sorted(motif_set_counts.items(),
                                                   key=lambda kv: kv[1],
                                                   reverse=True))
    print("motif_set_counts:", motif_set_counts_sorted_by_value)
    for motif_set_idx in motif_set_counts_sorted_by_value:
        print(', '.join([motif_ids_and_labels_dict[x] for x in motifs[motif_set_idx]]), ':', motif_set_counts_sorted_by_value[motif_set_idx])
    plt.show()
