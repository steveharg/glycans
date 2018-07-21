import sys
import pickle
import json
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

    pickle.dump(reactions_bag, open("reactions_bag.p", "wb"))
    pickle.dump(heatmaps, open("heatmaps.p", "wb"))
    pickle.dump(column_names, open("column_names.p", "wb"))

    # with open('reactions_bag.json', 'w') as outfile:
    #     json.dump(reactions_bag, outfile)
    #
    # with open('heatmaps.json', 'w') as outfile:
    #     json.dump(heatmaps, outfile)
    #
    # with open('column_names.json', 'w') as outfile:
    #     json.dump(column_names, outfile)