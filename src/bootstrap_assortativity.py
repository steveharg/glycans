import pickle
import configparser
import networkx as nx
import json
import numpy as np
import pandas as pd
from glycan_reaction_distances_to_network import apply_threshold_to_reaction_distance_df
from calc_jaccard_distances import calc_jaccard_distances

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')


if __name__ == "__main__":

    threshold = config['DEFAULT'].getint('Threshold', fallback=90)
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')
    reactions_collection_pickle_file = config['DEFAULT']['ReactionsCollectionPickleFile']

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

    print("opening " + reactions_collection_pickle_file + "...")
    reactions_collection = pickle.load(open(data_dir + "/" + reactions_collection_pickle_file, "rb"))
    print("...done")

    print("opening motif_ids_and_labels.json...")
    with open('motif_ids_and_labels.json') as f:
        motif_ids_and_labels = json.load(f)
    print("...done")

    num_bootstraps = 10
    use_bootstrap = True
    assortativity_coeffs = []

    motif_ids_and_labels_dict = {}
    for result in motif_ids_and_labels["results"]["bindings"]:
        motif_primary_id = result["MotifPrimaryId"]["value"]
        motif_label = result["MotifLabel"]["value"]
        motif_ids_and_labels_dict[motif_primary_id] = motif_label

    for n in range(0, num_bootstraps):
        print("bootstrap", n + 1, "of", num_bootstraps)
        heatmaps, column_names, reactions_collection_btsrp = calc_jaccard_distances(reactions_collection, use_reaction_quantities, use_bootstrap)
        reaction_distance_df = pd.DataFrame.from_dict(heatmaps[0],
                                                      orient='index', columns=column_names[0])
        thresholded_reactions_df = apply_threshold_to_reaction_distance_df(reactions_collection_btsrp, threshold)

        print("creating network graph from thresholded reactions...")
        G = nx.from_pandas_adjacency(thresholded_reactions_df)
        print("...done")

        motif_sets_from_all_glycans = []
        for reaction in reactions_collection_btsrp:
            sorted_motifs = sorted(reactions_collection_btsrp[reaction]['motifs'])
            if sorted_motifs not in motif_sets_from_all_glycans:
                motif_sets_from_all_glycans.append(sorted_motifs)

        node_colors = range(len(motif_sets_from_all_glycans))
        for node_color in node_colors:
            nodelist = []
            for reaction in reactions_collection_btsrp:
                sorted_motifs = sorted(reactions_collection_btsrp[reaction]['motifs'])
                if sorted_motifs == motif_sets_from_all_glycans[node_color]:
                    nodelist.append(reaction)

            node_motif_dict = dict.fromkeys(nodelist, ', '.join(
                [motif_ids_and_labels_dict[x] for x in motif_sets_from_all_glycans[node_color]]))
            nx.set_node_attributes(G, node_motif_dict, 'motifs')

        assortativity_coeff = nx.attribute_assortativity_coefficient(G, 'motifs')
        print("attribute_assortativity_coefficient:", assortativity_coeff)
        assortativity_coeffs.append(assortativity_coeff)

assortativity_coeff_mean = np.mean(assortativity_coeffs)
assortativity_coeff_std = np.std(assortativity_coeffs)

print("assortativity_coeff_mean:", assortativity_coeff_mean)
print("assortativity_coeff_std:", assortativity_coeff_std)