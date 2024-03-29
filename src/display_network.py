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
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    num_args = len(sys.argv)

    plot_heatmaps = False
    if num_args > 1:
        plot_heatmaps = bool(sys.argv[1])

    threshold = config['DEFAULT'].getint('Threshold', fallback=90)
    num_connected_components_displayed = config['DEFAULT'].getint('NumConnectedComponentsDisplayed', fallback=5)
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')

    if include_zero_motif_glycans:
        if use_reaction_quantities:
            data_dir = config['DEFAULT']['QuantifiedReactionsZeroMotifsDataDir']
        else:
            data_dir = config['DEFAULT']['SetOfReactionsZeroMotifsDataDir']
    else:
        if use_reaction_quantities:
            data_dir = config['DEFAULT']['QuantifiedReactionsDataDir']
        else:
            data_dir = config['DEFAULT']['SetOfReactionsDataDir']

    print("opening reactions_bag.p...")
    reactions_collection = pickle.load(open(data_dir + "/reactions_bag.p", "rb"))
    print("...done")
    print("opening G_" + str(threshold) + ".p...")
    G = pickle.load(open(data_dir + "/G_" + str(threshold) + ".p", "rb"))
    print("...done")
    print("opening pos_" + str(threshold) + ".p...")
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

    # Plot heatmaps
    if plot_heatmaps:

        print("opening heatmaps.p...")
        heatmaps = pickle.load(open(data_dir + "/heatmaps.p", "rb"))
        print("...done")
        print("opening column_names.p...")
        column_names = pickle.load(open(data_dir + "/column_names.p", "rb"))
        print("...done")

        reactions_heatmap_items = heatmaps[0]
        motifs_heatmap_items = heatmaps[1]
        glycan_distance_cols = column_names[0]
        motif_distance_cols = column_names[1]

        reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                      orient='index', columns=glycan_distance_cols)
        motif_distance_df = pd.DataFrame.from_dict(motifs_heatmap_items,
                                                   orient='index', columns=motif_distance_cols)

        print("plotting heatmaps...")
        sns.set()
        # plt.subplot(1,2,1)

        ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True, xticklabels=False,
                                   yticklabels=False)
        plt.title('Glycan Similarity by Reactions')

        # plt.subplot(1,2,2)
        # ax_motifs = sns.heatmap(motif_distance_df, annot=False, fmt=".2f", square=True)
        # plt.title('glycan similarity by motifs')

        plt.show()
        print("...done")


    print("getting connected_components from graph...")
    connected_components = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    print("connected_components:", connected_components)
    print("...done")

    motif_sets_from_all_glycans = []
    for reaction in reactions_collection:
        sorted_motifs = sorted(reactions_collection[reaction]['motifs'])
        if sorted_motifs not in motif_sets_from_all_glycans:
            motif_sets_from_all_glycans.append(sorted_motifs)

    print("len(motif_sets_from_all_glycans):", len(motif_sets_from_all_glycans))

    node_colors = range(len(motif_sets_from_all_glycans))

    print("applying colours to nodes in network...")

    tab20b = plt.cm.tab20b
    cNorm = colors.Normalize(vmin=node_colors[0], vmax=node_colors[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=tab20b)

    # draw the full network
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(1, 1, 1)

    # for label in node_colors:
    #     motif_labels = []
    #     for motif_id in motif_sets_from_all_glycans[label]:
    #         motif_labels.append(motif_ids_and_labels_dict[motif_id])
    #     ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=' & '.join(motif_labels))

    for node_color in node_colors:
        nodelist = []
        for reaction in reactions_collection:
            sorted_motifs = sorted(reactions_collection[reaction]['motifs'])
            if sorted_motifs == motif_sets_from_all_glycans[node_color]:
                nodelist.append(reaction)

        nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                               node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                               nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)
        node_motif_dict = dict.fromkeys(nodelist, ', '.join([motif_ids_and_labels_dict[x] for x in motif_sets_from_all_glycans[node_color]]))
        nx.set_node_attributes(G, node_motif_dict, 'motifs')

    print("attribute_assortativity_coefficient:", nx.attribute_assortativity_coefficient(G, 'motifs'))

    print("...done")

    print("drawing network graph...")
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
    leg = plt.legend(loc='upper right', fontsize="x-small")
    for line in leg.get_lines():
        line.set_linewidth(4.0)
    plt.axis('off')
    print("...done")

    # draw only the top (num_connected_components_displayed) connected components
    print("drawing connected components...")
    fig = plt.figure(2)
    ax = fig.add_subplot(1, 1, 1)
    plt.axis('off')
    motif_sets_from_most_conn_comps = []

    print("Motifs for top", num_connected_components_displayed, "connected components:")
    for i in range(0, num_connected_components_displayed):
        motifs_for_conn_comp = []
        for reaction in connected_components[i]:
            sorted_motifs = sorted(reactions_collection[reaction]['motifs'])
            if sorted_motifs not in motif_sets_from_most_conn_comps:
                motif_sets_from_most_conn_comps.append(sorted_motifs)
            if sorted_motifs not in motifs_for_conn_comp:
                motifs_for_conn_comp.append(sorted_motifs)
            # motifs_for_conn_comp = motifs_for_conn_comp + sorted_motifs
        # print(i+1,"^", ', '.join([motif_ids_and_labels_dict[x] for x in motifs_for_conn_comp]))
        for sorted_mts in motifs_for_conn_comp:
            print(i+1,"^", ', '.join([motif_ids_and_labels_dict[x] + ' ' + x for x in sorted_mts]))

    display_motifs = motif_sets_from_most_conn_comps
    # display_motifs = [['G00071MO']]
    # display_motifs = [['G00026MO'], ['G00055MO'], ['G00026MO', 'G00055MO']]
    # display_motifs = [['G00026MO'], ['G00026MO', 'G00055MO']]

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
    # display_motifs = [['G00056MO'], ['G00055MO']]  # 2 & 3 - Neo Lactosamine & Lactosamine motif

    # re-check
    # display_motifs = [['G00026MO', 'G00055MO', 'G00056MO'], ['G00026MO', 'G00055MO', 'G00057MO']]  # 18 & 19 - N-Glycan core basic, Lactosamine motif, Neo Lactosamine & N-Glycan core basic, Lactosamine motif, LacDiNAc
    # display_motifs = [['G00026MO', 'G00057MO'], ['G00026MO', 'G00055MO']]  # 5 & 7 - N-Glycan core basic, LacDiNAc & N-Glycan core basic, Lactosamine motif
    # display_motifs = [['G00026MO', 'G00055MO'], ['G00027MO', 'G00055MO']]  # 7 & 8 - N-Glycan core basic, Lactosamine motif & N-Glycan truncated motif. First GlcpNAC cut off, Lactosamine motif

    # display_motifs = [['G00046MO', 'G00055MO', 'G00068MO'], ['G00055MO', 'G00068MO']]  # 16 & 10 - Galalpha1-3Gal epitope, Lactosamine motif, Blood group H & Lactosamine motif, Blood group H
    # display_motifs = [['G00056MO'], ['G00055MO', 'G00056MO']]  # 2 & 12 - Neo Lactosamine & Lactosamine motif, Neo Lactosamine
    # display_motifs = [['G00046MO', 'G00055MO'], ['G00055MO']]  # 9 & 3 - Galalpha1-3Gal epitope, Lactosamine motif & Lactosamine motif
    # display_motifs = [['G00057MO'], ['G00055MO', 'G00057MO']]  # 4 & 13 LacDiNAc & Lactosamine motif, LacDiNAc

    # dissimilar?
    display_motifs = [['G00026MO'], ['G00055MO']]  # 1 & 3 N-Glycan core basic & Lactosamine motif
    # display_motifs = [['G00057MO'], ['G00026MO', 'G00057MO']]  # 4 & 5 – LacDiNAc & N-Glycan core basic, LacDiNAc
    # display_motifs = [['G00056MO'], ['G00057MO']]  # 2 & 4 - Neo Lactosamine & LacDiNAc

    print("(connected components) len(motifs):", len(motif_sets_from_most_conn_comps))
    node_colors = range(len(motif_sets_from_most_conn_comps))

    tab20b = plt.cm.tab20b
    cNorm = colors.Normalize(vmin=node_colors[0], vmax=node_colors[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=tab20b)

    # for label in node_colors:
    #     motif_labels = []
    #     for motif_id in motif_sets_from_most_conn_comps[label]:
    #         motif_labels.append(motif_ids_and_labels_dict[motif_id])
    #     if motif_sets_from_most_conn_comps[label] in display_motifs:
    #         ax.plot([0], [0], color=scalarMap.to_rgba(node_colors[label]), label=' & '.join(motif_labels))
    # ax.plot([0], [0], color='Gray', label='All other motif sets')

    motif_set_counts = dict.fromkeys(node_colors, 0)

    for i in range(0, num_connected_components_displayed):

        for node_color in node_colors:
            nodelist = []
            for reaction in connected_components[i]:
                sorted_motifs = sorted(reactions_collection[reaction]['motifs'])
                if sorted_motifs == motif_sets_from_most_conn_comps[node_color]:
                    nodelist.append(reaction)

            if len(nodelist) > 0:
                motif_set_counts[node_color] += 1

            # color this set of nodes in gray first, then re-color if we've specified them in the
            # display_motifs list
            nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                                   node_color='gray',
                                   nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)
            if motif_sets_from_most_conn_comps[node_color] in display_motifs:
                nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                                       node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                                       nodelist=nodelist, cmap=plt.cm.tab20b, vmin=node_colors[0], vmax=node_colors[-1], ax=ax)

            nx.draw_networkx_edges(G, pos, edgelist=G.edges(nbunch=nodelist), width=1.0, alpha=0.5, ax=ax)

        leg = plt.legend(fontsize="medium")
        # leg = plt.legend(fontsize="small")
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
        # plt.title(plot_title)
        plt.axis('off')

    print("...done")
    motif_set_counts_sorted_by_value = dict(sorted(motif_set_counts.items(),
                                                   key=lambda kv: kv[1],
                                                   reverse=True))
    print("motif_set_counts:", motif_set_counts_sorted_by_value)

    for motif_set_idx in motif_set_counts_sorted_by_value:
        print(', '.join([motif_ids_and_labels_dict[x] for x in motif_sets_from_most_conn_comps[motif_set_idx]]), ':', motif_set_counts_sorted_by_value[motif_set_idx])
    plt.show()
