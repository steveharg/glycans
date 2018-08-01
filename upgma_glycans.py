import pickle
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor, DistanceCalculator
from Bio import Phylo

if __name__ == "__main__":

    data_dir = "data_top10_motifs"
    # data_dir = "mini_data_top10_motifs"

    print("opening heatmaps.p...")
    heatmaps = pickle.load(open(data_dir + "/heatmaps.p", "rb"))
    print("...done")
    print("opening column_names.p...")
    column_names = pickle.load(open(data_dir + "/column_names.p", "rb"))
    print("...done")

    reactions_heatmap_items = heatmaps[0]
    glycan_distance_cols = column_names[0]

    reaction_distance_df = pd.DataFrame.from_dict(reactions_heatmap_items,
                                                  orient='index', columns=glycan_distance_cols)

    reaction_distance_df_vals = reaction_distance_df.get_values()

    reaction_distance_lower_tri = np.multiply(reaction_distance_df_vals, np.tril(np.ones(reaction_distance_df_vals.shape)))

    reaction_distance_lower_tri = (100 - reaction_distance_lower_tri)/100

    print(reaction_distance_lower_tri[:][0])
    print(reaction_distance_lower_tri[:][1])

    i = 1
    reaction_distance_lower_tri_list = []
    for row in reaction_distance_lower_tri:
        reaction_distance_lower_tri_list.append(row[0:i].tolist())
        i += 1

    # reaction_distance_lower_tri_list = reaction_distance_lower_tri.tolist()

    print(reaction_distance_lower_tri_list[0])
    print(reaction_distance_lower_tri_list[1])

    # print(reaction_distance_lower_tri.tolist())

    dmm = DistanceMatrix(names=glycan_distance_cols, matrix=reaction_distance_lower_tri_list)

    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dmm)

    Phylo.draw(tree)

    print("dfs")