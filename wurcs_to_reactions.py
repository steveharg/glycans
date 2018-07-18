from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances
from plot_glycan_heatmaps import plot_glycan_heatmaps

if __name__ == "__main__":
    reactions_bag = parse_glytoucan_results('glycans.json')
    heatmaps, column_names = calc_jaccard_distances(reactions_bag)
    plot_glycan_heatmaps(heatmaps, column_names)
