import sys
import pickle
import json
from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances
from filter_by_common_motifs import filter_by_common_motifs

if __name__ == "__main__":

    num_args = len(sys.argv)

    num_datafile_lines = None

    if num_args > 1:
        num_datafile_lines = int(sys.argv[1])

    reactions_bag = parse_glytoucan_results('glycans.json', num_datafile_lines)
    reactions_bag = filter_by_common_motifs(reactions_bag)
    heatmaps, column_names = calc_jaccard_distances(reactions_bag)

    reaction_distances = [reactions_bag, heatmaps, column_names]

    data_dir = "data_top10_motifs"

    pickle.dump(reactions_bag, open(data_dir + "/reactions_bag.p", "wb"))
    pickle.dump(heatmaps, open(data_dir + "/heatmaps.p", "wb"))
    pickle.dump(column_names, open(data_dir + "/column_names.p", "wb"))
