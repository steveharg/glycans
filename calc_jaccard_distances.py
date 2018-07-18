import sys
import numpy as np

np.random.seed(0)


def calc_jaccard_distances(reactions_bag):

    # Calculate Jaccard distances between all structures (both, and separately, using reactions and motifs)
    glycan_distances=np.zeros((len(reactions_bag),len(reactions_bag)))
    glycan_distance_cols = []
    reactions_heatmap_items = []

    motif_distances=np.zeros((len(reactions_bag),len(reactions_bag)))
    motif_distance_cols = []
    motifs_heatmap_items = []

    glycan_idx1 = 0
    for glycan1 in reactions_bag:
        glycan_distance_cols.append(glycan1)
        motif_distance_cols.append(glycan1)
        # glycan_distance_cols.append('-'.join(reactions_bag[glycan1]['motifs']))
        glycan_idx2 = 0
        for glycan2 in reactions_bag:
            reactions1 = reactions_bag[glycan1]['reactions']
            reactions2 = reactions_bag[glycan2]['reactions']
            motifs1 = reactions_bag[glycan1]['motifs']
            motifs2 = reactions_bag[glycan2]['motifs']
            reactions_match_count = 0
            motifs_match_count = 0
            unique_reactions_in_both = []
            unique_motifs_in_both = []
            for reaction1 in reactions1:
                for reaction2 in reactions2:
                    if reaction1 not in unique_reactions_in_both:
                        unique_reactions_in_both.append(reaction1)
                    if reaction2 not in unique_reactions_in_both:
                        unique_reactions_in_both.append((reaction2))

                    if reaction1 == reaction2:
                        reactions_match_count += 1

            for motif1 in motifs1:
                for motif2 in motifs2:
                    if motif1 not in unique_motifs_in_both:
                        unique_motifs_in_both.append(motif1)
                    if motif2 not in unique_motifs_in_both:
                        unique_motifs_in_both.append((motif2))

                    if motif1 == motif2 and motif1 != 'None':
                        motifs_match_count += 1

            glycan_distances[glycan_idx1][glycan_idx2] = 100 * reactions_match_count / len(unique_reactions_in_both)
            motif_distances[glycan_idx1][glycan_idx2] = 100 * motifs_match_count / len(unique_motifs_in_both)
            glycan_idx2 += 1

        reactions_heatmap_items.append((glycan1, glycan_distances[glycan_idx1][:]))
        motifs_heatmap_items.append((glycan1, motif_distances[glycan_idx1][:]))
        glycan_idx1 += 1
        sys.stdout.write("\rglycan_idx1: %i of %i" % (glycan_idx1, len(reactions_bag)))
        sys.stdout.flush()

        heatmaps = [reactions_heatmap_items, motifs_heatmap_items]
        column_names = [glycan_distance_cols, motif_distance_cols]

    return heatmaps, column_names
