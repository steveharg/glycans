import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def plot_glycan_heatmaps(heatmaps, column_names):
    reactions_heatmap_items = heatmaps[0]
    motifs_heatmap_items = heatmaps[1]
    glycan_distance_cols = column_names[0]
    motif_distance_cols = column_names[1]

    sns.set()

    # Plot heatmaps
    reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                   orient='index', columns=glycan_distance_cols)
    motif_distance_df = pd.DataFrame.from_dict(motifs_heatmap_items,
                                                   orient='index', columns=motif_distance_cols)

    plt.subplot(1,2,1)

    ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True)
    plt.title('pairwise glycan similarity by reactions')

    plt.subplot(1,2,2)
    ax_motifs = sns.heatmap(motif_distance_df, annot=False, fmt=".2f", square=True)
    plt.title('pairwise glycan similarity by motifs')

    plt.show()

    return reaction_distance_df, motif_distance_df