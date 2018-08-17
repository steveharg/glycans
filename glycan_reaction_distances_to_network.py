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
    if num_args > 1:
        threshold = int(sys.argv[1])

    # data_dir = "data_top10_motifs"
    # data_dir = "mini_data_top10_motifs"
    data_dir = "data_top10_with_reac_quants_motifs"

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

    # Apply threshold
    thresholded_reactions = []
    print("applying threshold to glycan distances...")
    for col in reaction_distance_df.columns:
        col_array = np.zeros((1, len(reaction_distance_df[col])))
        i = 0
        for item in reaction_distance_df[col]:
            if item > threshold:
                col_array[0, i] = 1
            i += 1

        thresholded_reactions.append((col, col_array[0][:]))
    print("...done")

    cols = [x[0] for x in thresholded_reactions]
    # print("cols:", cols)
    thresholded_reactions_df = pd.DataFrame.from_items(thresholded_reactions,
                                                       orient='index', columns=cols)

    print("creating network graph from thresholded reactions...")
    G = nx.from_pandas_adjacency(thresholded_reactions_df)
    print("...done")

    print("creating graph layout...")
    # pos = nx.spring_layout(G)
    # pos = nx.circular_layout(G)
    # pos = nx.shell_layout(G)
    pos = nx.fruchterman_reingold_layout(G) # similar to spring
    # pos = nx.kamada_kawai_layout(G) # might work
    # pos = nx.spectral_layout(G) # everything appears on a dense horizontal line, not much use
    print("...done")

    print("Pickling G and pos...")
    pickle.dump(G, open(data_dir + "/G.p", "wb"))
    pickle.dump(pos, open(data_dir + "/pos_" + str(threshold) + ".p", "wb"))
    print("...done")
