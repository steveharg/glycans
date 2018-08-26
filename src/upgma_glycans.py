import pickle
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor, DistanceCalculator
from Bio import Phylo
import json
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    heatmaps_pickle_file = config['DEFAULT']['HeatmapsPickleFile']
    col_names_pickle_file = config['DEFAULT']['ColNamesPickleFile']
    motif_ids_and_labels_json_file = config['DEFAULT']['MotifIdsAndLabelsJsonFile']
    glycans_json_file = config['DEFAULT']['GlycansJsonFile']
    tree_with_motif_names_pickle_file = config['DEFAULT']['TreeWithMotifNamesPickleFile']

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

    # create a glycan id to motif ids dictionary
    glycan_motifs = {}

    print("opening " + motif_ids_and_labels_json_file + "...")
    with open(motif_ids_and_labels_json_file) as f:
        motif_ids_and_labels = json.load(f)
    print("...done")

    # motif_ids_and_labels_dict
    motif_ids_and_labels_dict = {}
    for result in motif_ids_and_labels["results"]["bindings"]:
        motif_primary_id = result["MotifPrimaryId"]["value"]
        motif_label = result["MotifLabel"]["value"]
        motif_ids_and_labels_dict[motif_primary_id] = motif_label

    with open(glycans_json_file) as f:
        results = json.load(f)

    for result in results["results"]["bindings"]:
        glycan_id = result["PrimaryId"]["value"]
        if glycan_id not in glycan_motifs:
            glycan_motifs[glycan_id] = []
        if "MotifPrimaryId" in result:
            motif_id = result["MotifPrimaryId"]["value"]
            if motif_id not in glycan_motifs[glycan_id]:
                glycan_motifs[glycan_id].append(motif_id)

    # replace motif ids with motif names
    for glycan_id in glycan_motifs:
        motif_ids = sorted(glycan_motifs[glycan_id])
        motif_names = []
        for motif_id in motif_ids:
            motif_names.append(motif_ids_and_labels_dict[motif_id])
        glycan_motifs[glycan_id] = ', '.join(motif_names)

    reactions_heatmap_items = heatmaps[0]
    glycan_distance_cols = column_names[0]

    # replace glycan_ids in glycan_distance_cols with motif name lists
    glycan_distance_cols_temp = []
    for glycan_id in glycan_distance_cols:
        print("appending:", glycan_motifs[glycan_id])
        glycan_distance_cols_temp.append(glycan_motifs[glycan_id] + '(' + glycan_id + ')')
        if glycan_id == "G00001NT":
            print('no motifs')

    glycan_distance_cols = glycan_distance_cols_temp

    reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                  orient='index', columns=glycan_distance_cols)

    reaction_distance_df_vals = reaction_distance_df.get_values()
    reaction_distance_lower_tri = np.multiply(reaction_distance_df_vals, np.tril(np.ones(reaction_distance_df_vals.shape)))
    reaction_distance_lower_tri = (100 - reaction_distance_lower_tri)/100

    i = 1
    reaction_distance_lower_tri_list = []
    for row in reaction_distance_lower_tri:
        reaction_distance_lower_tri_list.append(row[0:i].tolist())
        i += 1

    print('calculating DistanceMatrix...')
    dmm = DistanceMatrix(names=glycan_distance_cols, matrix=reaction_distance_lower_tri_list)
    print('...done')

    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dmm)

    print('pickling tree...')
    pickle.dump(tree, open(data_dir + "/" + tree_with_motif_names_pickle_file, "wb"))
    print('...done')

    # Phylo.draw(tree)