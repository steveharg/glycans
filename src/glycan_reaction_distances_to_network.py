import networkx as nx
import numpy as np
import pandas as pd
import sys
import pickle
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')


def apply_threshold_to_reaction_distance_df(reaction_dis_df, threshld):
    threshld_reactions = []
    print("applying threshold to glycan distances...")
    for col in reaction_dis_df.columns:
        col_array = np.zeros((1, len(reaction_dis_df[col])))
        i = 0
        for item in reaction_dis_df[col]:
            if item > threshld:
                col_array[0, i] = 1
            i += 1

            threshld_reactions.append((col, col_array[0][:]))
    print("...done")

    cols = [x[0] for x in threshld_reactions]
    thresholded_reactions_df = pd.DataFrame.from_items(threshld_reactions,
                                                       orient='index', columns=cols)
    return thresholded_reactions_df


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
        # open the graph G previously created using the reaction set method
        print("opening G_" + str(threshold) + ".p from " + config['DEFAULT']['SetOfReactionsDataDir'] + "/...")
        G_set = pickle.load(open(config['DEFAULT']['SetOfReactionsDataDir'] + "/G_" + str(threshold) + ".p", "rb"))
        print("...done")
        # count its number of edges
        num_edges_in_g_set = len(G_set.edges)
        print("num_edges_in_g_set:", num_edges_in_g_set)

        # create a new graph from the full glycan similarity matrix created using the reaction quantities method,
        # applying the similarity measure as edge weights
        G_list = nx.from_numpy_array(reaction_distance_df.values)

        # the nodes of this new graph only have the indices of the similarity matrix as labels, so relabel using glycan
        # IDs
        glycan_name_and_idx = dict(zip(range(0, len(glycan_distance_cols)), glycan_distance_cols))
        G_list_relabelled = nx.relabel_nodes(G_list, glycan_name_and_idx)

        # sort the edges by weight
        G_list_edges = G_list_relabelled.edges(data=True)
        ordered_edges = sorted(G_list_edges, key=lambda x: x[2]['weight'], reverse=True)

        # create a new graph using only the same number of (highest weight) edges present in the original 'reaction set
        # method' graph
        G = nx.Graph()
        G.add_edges_from(ordered_edges[0:num_edges_in_g_set])

        # save the new graph. the threshold value in the filename refers to the threshold used to produce
        # the corresponding graph created using the reaction set method, from which we read the number of edges
        # to be used in this graph
        print("Pickling G_" + str(threshold) + ".p...")
        pickle.dump(G, open(data_dir + "/G_" + str(threshold) + ".p", "wb"))
        print("...done")

    else:
        # Apply threshold
        thresholded_reactions_df = apply_threshold_to_reaction_distance_df(reaction_distance_df, threshold)

        print("creating network graph from thresholded reactions...")
        G = nx.from_pandas_adjacency(thresholded_reactions_df)
        print("...done")
        print("Pickling G_" + str(threshold) + "...")
        pickle.dump(G, open(data_dir + "/G_" + str(threshold) + ".p", "wb"))
        print("...done")

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
