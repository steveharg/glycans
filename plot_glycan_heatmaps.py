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
    reaction_distance_df = pd.DataFrame.from_items(reactions_heatmap_items[::-1],
                                                   orient='index', columns=glycan_distance_cols)
    motif_distance_df = pd.DataFrame.from_items(motifs_heatmap_items[::-1],
                                                   orient='index', columns=motif_distance_cols)

    plt.subplot(1,2,1)

    # ax = sns.heatmap(ani, annot=True, fmt=".2f", square=True)
    ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True)
    plt.title('pairwise glycan similarity by reactions')

    plt.subplot(1,2,2)
    ax_motifs = sns.heatmap(motif_distance_df, annot=False, fmt=".2f", square=True)
    plt.title('pairwise glycan similarity by motifs')

    plt.show()

    # fig = ax.get_figure()
    # fig.savefig('strain_ANIs2.png')