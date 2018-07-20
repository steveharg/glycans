import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys
from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances

if __name__ == "__main__":

    num_args = len(sys.argv)

    threshold = 90
    num_datafile_lines = None
    plot_heatmaps = False

    if num_args > 1:
        threshold = int(sys.argv[1])
    if num_args > 2:
        num_datafile_lines = int(sys.argv[2])
    if num_args > 3:
        plot_heatmaps = bool(sys.argv[3])

    reactions_bag = parse_glytoucan_results('glycans.json', num_datafile_lines)
    heatmaps, column_names = calc_jaccard_distances(reactions_bag)

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
        sns.set()
        plt.subplot(1,2,1)

        ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True)
        plt.title('pairwise glycan similarity by reactions')

        plt.subplot(1,2,2)
        ax_motifs = sns.heatmap(motif_distance_df, annot=False, fmt=".2f", square=True)
        plt.title('pairwise glycan similarity by motifs')

        plt.show()

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



    # G = nx.from_pandas_adjacency(reaction_distance_df)
    G = nx.from_pandas_adjacency(thresholded_reactions_df)

    motifs = []
    for reaction in reactions_bag:
        sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
        if sorted_motifs not in motifs:
            motifs.append(sorted_motifs)

    # print("motifs:", motifs)

    node_colors = range(len(motifs))
    pos = nx.spring_layout(G)

    for node_color in node_colors:
        nodelist = []
        for reaction in reactions_bag:
            sorted_motifs = sorted(reactions_bag[reaction]['motifs'])
            if sorted_motifs == motifs[node_color]:
                nodelist.append(reaction)

        nx.draw_networkx_nodes(G, pos, with_labels=True, font_weight='bold',
                               node_color=np.zeros((len(nodelist)))+node_colors[node_color],
                               nodelist=nodelist, cmap=plt.cm.Reds, vmin=node_colors[0], vmax=node_colors[-1])

    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
    plt.show()