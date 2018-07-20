import sys
import pickle
from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances

if __name__ == "__main__":

    num_args = len(sys.argv)

    num_datafile_lines = None

    if num_args > 1:
        num_datafile_lines = int(sys.argv[1])

    reactions_bag = parse_glytoucan_results('glycans.json', num_datafile_lines)
    heatmaps, column_names = calc_jaccard_distances(reactions_bag)

    reaction_distances = [reactions_bag, heatmaps, column_names]
    pickle.dump(reaction_distances, open("reaction_distances.p", "wb"))
