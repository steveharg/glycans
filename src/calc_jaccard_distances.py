import sys
import numpy as np
import time
from collections import Counter

np.random.seed(0)

def bootstrap_reactions(reactions1, reactions2):

    reactions1_set = set(reactions1)
    reactions2_set = set(reactions2)

    union_reactions_set_list = list(reactions1_set | reactions2_set)
    num_reactions = len(union_reactions_set_list)

    btstrp_reactions1 = []
    btstrp_reactions2 = []

    for i in range(0, num_reactions):
        rnd_idx = np.random.randint(0, num_reactions)
        if union_reactions_set_list[rnd_idx] in reactions1_set:
            btstrp_reactions1.append(union_reactions_set_list[rnd_idx])
        if union_reactions_set_list[rnd_idx] in reactions2_set:
            btstrp_reactions2.append(union_reactions_set_list[rnd_idx])

    return btstrp_reactions1, btstrp_reactions2

def calc_jaccard_distances(reactions_collection, use_reaction_quantities, use_bootstrap = False):

    start = time.time()

    # Calculate Jaccard distances between all structures (both, and separately, using reactions and motifs)
    glycan_distances=np.zeros((len(reactions_collection),len(reactions_collection)))
    reactions_heatmap_items = {}

    motif_distances=np.zeros((len(reactions_collection),len(reactions_collection)))
    motifs_heatmap_items = {}
    if use_bootstrap:
        reactions_collection_btsrp = dict.fromkeys(reactions_collection.keys())
    else:
        reactions_collection_btsrp = None

    glycan_idx1 = 0

    if use_reaction_quantities:
        for glycan1 in reactions_collection:
            glycan_idx2 = 0
            for glycan2 in reactions_collection:
                reactions1 = reactions_collection[glycan1]['reactions']
                reactions2 = reactions_collection[glycan2]['reactions']
                motifs1 = reactions_collection[glycan1]['motifs']
                motifs2 = reactions_collection[glycan2]['motifs']

                reactions1_set = set(reactions1)
                reactions2_set = set(reactions2)
                counter_reactions1 = Counter(reactions1)
                counter_reactions2 = Counter(reactions2)
                unique_reactions_in_either_set = reactions1_set | reactions2_set
                unique_reactions_in_both_sets = reactions1_set & reactions2_set
                jaccard_entries = {}
                for x in unique_reactions_in_both_sets:
                    jaccard_entries[x] = min(counter_reactions1[x], counter_reactions2[x]) / max(counter_reactions1[x], counter_reactions2[x])

                motifs1_set = set(motifs1)
                motifs2_set = set(motifs2)
                motifs1_exclude_none_set = set(motifs1) - set(['None'])
                motifs2_exclude_none_set = set(motifs2) - set(['None'])
                unique_motifs_in_either_set = motifs1_set | motifs2_set
                unique_motifs_in_both_sets = motifs1_exclude_none_set & motifs2_exclude_none_set

                glycan_distances[glycan_idx1][glycan_idx2] = 100 * sum(jaccard_entries.values())/len(unique_reactions_in_either_set)
                motif_distances[glycan_idx1][glycan_idx2] = 100 * len(unique_motifs_in_both_sets) / len(
                    unique_motifs_in_either_set)
                glycan_idx2 += 1

            reactions_heatmap_items[glycan1] = glycan_distances[glycan_idx1][:]
            motifs_heatmap_items[glycan1] = motif_distances[glycan_idx1][:]
            glycan_idx1 += 1

            if not np.mod(glycan_idx1, 100):
                sys.stdout.write("\rglycan_idx1: %i of %i" % (glycan_idx1, len(reactions_collection)))
                sys.stdout.flush()

    else:
        for glycan1 in reactions_collection:
            glycan_idx2 = 0
            for glycan2 in reactions_collection:
                reactions1 = reactions_collection[glycan1]['reactions']
                reactions2 = reactions_collection[glycan2]['reactions']
                motifs1 = reactions_collection[glycan1]['motifs']
                motifs2 = reactions_collection[glycan2]['motifs']

                if use_bootstrap:
                    reactions1, reactions2 = bootstrap_reactions(reactions1, reactions2)

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
                sys.stdout.write("\rglycan_idx1: %i of %i" % (glycan_idx1, len(reactions_collection)))
                sys.stdout.flush()

    heatmaps = [reactions_heatmap_items, motifs_heatmap_items]
    column_names = [list(reactions_collection.keys()), list(reactions_collection.keys())]

    end = time.time()
    print()
    print("elapsed jaccard time:", end-start)

    return heatmaps, column_names, reactions_collection_btsrp
