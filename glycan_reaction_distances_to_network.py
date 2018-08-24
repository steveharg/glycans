import networkx as nx
import numpy as np
import pandas as pd
import sys
import pickle
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    threshold = config['DEFAULT'].getint('Threshold', fallback=90)
    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    heatmaps_pickle_file = config['DEFAULT']['HeatmapsPickleFile']
    col_names_pickle_file = config['DEFAULT']['ColNamesPickleFile']

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

    print("opening " + heatmaps_pickle_file + "...")
    heatmaps = pickle.load(open(data_dir + "/" + heatmaps_pickle_file, "rb"))
    print("...done")
    print("opening " + col_names_pickle_file + "...")
    column_names = pickle.load(open(data_dir + "/" + col_names_pickle_file, "rb"))
    print("...done")

    reactions_heatmap_items = heatmaps[0]
    motifs_heatmap_items = heatmaps[1]
    glycan_distance_cols = column_names[0]
    motif_distance_cols = column_names[1]

    reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                   orient='index', columns=glycan_distance_cols)
    motif_distance_df = pd.DataFrame.from_dict(motifs_heatmap_items,
                                                   orient='index', columns=motif_distance_cols)

    if use_reaction_quantities:
        print("opening G_" + str(threshold) + ".p from " + config['DEFAULT']['SetOfReactionsDataDir'] + "/...")
        G_set = pickle.load(open(config['DEFAULT']['SetOfReactionsDataDir'] + "/G_" + str(threshold) + ".p", "rb"))
        print("...done")
        num_edges_in_g_set = len(G_set.edges)
        print("num_edges_in_g_set:", num_edges_in_g_set)
        # G_set_ordered_edges = sorted(G_set.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)
        # print("G_set_ordered_edges:", G_set_ordered_edges[0:200])
        # sys.exit(0)

        glycan_name_and_idx = dict(zip(range(0, len(glycan_distance_cols)), glycan_distance_cols))
        G_list = nx.from_numpy_array(reaction_distance_df.values)
        G_list2 = nx.relabel_nodes(G_list, glycan_name_and_idx)
        g_list_edges = G_list2.edges(data=True)
        # print("unordered g_list_edges:", g_list_edges)
        # ordered_edges = sorted(g_list_edges, key=lambda edge: edge['weight'])
        ordered_edges = sorted(g_list_edges, key=lambda x: x[2]['weight'], reverse=True)
        print("ordered_edges:", ordered_edges)
        truncated_edges_list = ordered_edges[0:num_edges_in_g_set]
        print("truncated_edges_list:", truncated_edges_list)
        sys.exit(0)
    else:
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
        print("Pickling G_" + str(threshold) + "...")
        pickle.dump(G, open(data_dir + "/G_" + str(threshold) + ".p", "wb"))

    print("creating graph layout...")
    # pos = nx.spring_layout(G)
    # pos = nx.circular_layout(G)
    # pos = nx.shell_layout(G)
    pos = nx.fruchterman_reingold_layout(G) # similar to spring
    # pos = nx.kamada_kawai_layout(G) # might work
    # pos = nx.spectral_layout(G) # everything appears on a dense horizontal line, not much use
    print("...done")

    print("Pickling pos_" + str(threshold) + "...")
    pickle.dump(pos, open(data_dir + "/pos_" + str(threshold) + ".p", "wb"))
    print("...done")
