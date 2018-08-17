import sys
import numpy as np
import time
from collections import Counter

np.random.seed(0)


def calc_jaccard_distances(reactions_bag, use_reaction_quantities=True):

    start = time.time()

    # Calculate Jaccard distances between all structures (both, and separately, using reactions and motifs)
    glycan_distances=np.zeros((len(reactions_bag),len(reactions_bag)))
    reactions_heatmap_items = {}

    motif_distances=np.zeros((len(reactions_bag),len(reactions_bag)))
    motifs_heatmap_items = {}

    glycan_idx1 = 0

    if use_reaction_quantities:
        for glycan1 in reactions_bag:
            glycan_idx2 = 0
            for glycan2 in reactions_bag:
                reactions1 = reactions_bag[glycan1]['reactions']
                reactions2 = reactions_bag[glycan2]['reactions']
                motifs1 = reactions_bag[glycan1]['motifs']
                motifs2 = reactions_bag[glycan2]['motifs']

                reactions1_set = set(reactions1)
                reactions2_set = set(reactions2)
                counter_reactions1_set = Counter(reactions1_set)
                counter_reactions2_set = Counter(reactions2_set)
                unique_reactions_in_either_set = reactions1_set | reactions2_set
                unique_reactions_in_both_sets = reactions1_set & reactions2_set
                jaccard_entries = {}
                for x in unique_reactions_in_both_sets:
                    jaccard_entries[x] = min(counter_reactions1_set[x], counter_reactions2_set[x]) / max(counter_reactions1_set[x], counter_reactions2_set[x])

                motifs1_set = set(motifs1)
                motifs2_set = set(motifs2)
                motifs1_exclude_none_set = set(motifs1) - set(['None'])
                motifs2_exclude_none_set = set(motifs2) - set(['None'])
                unique_motifs_in_either_set = motifs1_set | motifs2_set
                unique_motifs_in_both_sets = motifs1_exclude_none_set & motifs2_exclude_none_set

                glycan_distances[glycan_idx1][glycan_idx2] = 100 * sum(jaccard_entries.values())/len(reactions1_set | reactions2_set)
                motif_distances[glycan_idx1][glycan_idx2] = 100 * len(unique_motifs_in_both_sets) / len(
                    unique_motifs_in_either_set)
                glycan_idx2 += 1

            reactions_heatmap_items[glycan1] = glycan_distances[glycan_idx1][:]
            motifs_heatmap_items[glycan1] = motif_distances[glycan_idx1][:]
            glycan_idx1 += 1

            if not np.mod(glycan_idx1, 100):
                sys.stdout.write("\rglycan_idx1: %i of %i" % (glycan_idx1, len(reactions_bag)))
                sys.stdout.flush()

    else:
        for glycan1 in reactions_bag:
            glycan_idx2 = 0
            for glycan2 in reactions_bag:
                reactions1 = reactions_bag[glycan1]['reactions']
                reactions2 = reactions_bag[glycan2]['reactions']
                motifs1 = reactions_bag[glycan1]['motifs']
                motifs2 = reactions_bag[glycan2]['motifs']

                reactions1_set = set(reactions1)
                reactions2_set = set(reactions2)
                unique_reactions_in_either_set = reactions1_set | reactions2_set
                unique_reactions_in_both_sets = reactions1_set & reactions2_set

                motifs1_set = set(motifs1)
                motifs2_set = set(motifs2)
                motifs1_exclude_none_set = set(motifs1) - set(['None'])
                motifs2_exclude_none_set = set(motifs2) - set(['None'])
                unique_motifs_in_either_set = motifs1_set | motifs2_set
                unique_motifs_in_both_sets = motifs1_exclude_none_set & motifs2_exclude_none_set

                glycan_distances[glycan_idx1][glycan_idx2] = 100 * len(unique_reactions_in_both_sets) /  len(unique_reactions_in_either_set)
                motif_distances[glycan_idx1][glycan_idx2] = 100 * len(unique_motifs_in_both_sets) / len(unique_motifs_in_either_set)
                glycan_idx2 += 1

            reactions_heatmap_items[glycan1] = glycan_distances[glycan_idx1][:]
            motifs_heatmap_items[glycan1] = motif_distances[glycan_idx1][:]
            glycan_idx1 += 1

            if not np.mod(glycan_idx1, 100):
                sys.stdout.write("\rglycan_idx1: %i of %i" % (glycan_idx1, len(reactions_bag)))
                sys.stdout.flush()

    heatmaps = [reactions_heatmap_items, motifs_heatmap_items]
    column_names = [list(reactions_bag.keys()), list(reactions_bag.keys())]

    end = time.time()
    print()
    print("elapsed jaccard time:", end-start)

    return heatmaps, column_names
