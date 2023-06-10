import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data_dir = "data_top10_motifs"
# data_dir = "data_top10_with_reac_quants_motifs"

heatmaps = pickle.load(open(data_dir + "/heatmaps.p", "rb"))
column_names = pickle.load(open(data_dir + "/column_names.p", "rb"))
reactions_heatmap_items = heatmaps[0]
glycan_distance_cols = column_names[0]
reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                   orient='index', columns=glycan_distance_cols)

data_dir_q = "data_top10_with_reac_quants_motifs"
# data_dir = "data_top10_with_reac_quants_motifs_random"

heatmaps_q = pickle.load(open(data_dir_q + "/heatmaps.p", "rb"))
column_names_q = pickle.load(open(data_dir_q + "/column_names.p", "rb"))
reactions_heatmap_items_q = heatmaps_q[0]
glycan_distance_cols_q = column_names_q[0]
reaction_distance_df_q = pd.DataFrame.from_dict(reactions_heatmap_items_q,
                                                   orient='index', columns=glycan_distance_cols_q)

df_diffs = reaction_distance_df - reaction_distance_df_q

stop_ere = 902

print("plotting heatmaps...")
sns.set()
plt.subplot(1, 3, 1)

ax_reactions = sns.heatmap(reaction_distance_df, annot=False, fmt=".2f", square=True)
plt.title('pairwise glycan similarity by reactions')

plt.subplot(1, 3, 2)
ax_reactions = sns.heatmap(reaction_distance_df_q, annot=False, fmt=".2f", square=True)
plt.title('pairwise glycan similarity by reactions with quantities')

plt.subplot(1, 3, 3)
ax_reactions = sns.heatmap(df_diffs, annot=False, fmt=".2f", square=True)
plt.title('df_diffs')

plt.show()
print("...done")

stop_ere = 902