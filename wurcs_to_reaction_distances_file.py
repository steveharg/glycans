import sys
import pickle
import json
from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances
from filter_by_common_motifs import filter_by_common_motifs

if __name__ == "__main__":

    num_args = len(sys.argv)

    include_zero_motif_glycans = False
    num_datafile_lines = None # None means use all datafile lines
    use_reaction_quantities = True

    if num_args > 1:
        include_zero_motif_glycans = 'True' == sys.argv[1]
    if num_args > 2:
        num_datafile_lines = int(sys.argv[2])

    reactions_bag = parse_glytoucan_results('glycans.json', num_datafile_lines, use_reaction_quantities=use_reaction_quantities)
    reactions_bag = filter_by_common_motifs(reactions_bag, include_zero_motif_glycans)
    heatmaps, column_names = calc_jaccard_distances(reactions_bag, use_reaction_quantities=use_reaction_quantities)

    reaction_distances = [reactions_bag, heatmaps, column_names]

    if include_zero_motif_glycans:
        if use_reaction_quantities:
            data_dir = "data_top10_with_reac_quants_and_zero_motifs"
        else:
            data_dir = "data_top10_and_zero_motifs"
    else:
        if use_reaction_quantities:
            data_dir = "data_top10_with_reac_quants_motifs"
        else:
            data_dir = "data_top10_motifs"

    # data_dir = "mini_data_top10_motifs"

    pickle.dump(reactions_bag, open(data_dir + "/reactions_bag.p", "wb"))
    pickle.dump(heatmaps, open(data_dir + "/heatmaps.p", "wb"))
    pickle.dump(column_names, open(data_dir + "/column_names.p", "wb"))
