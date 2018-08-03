import matplotlib.pyplot as plt

def filter_by_common_motifs(reactions_bag, include_zero_motif_glycans):

    print("len(reactions_bag) before removing any glycans with zero or more than one motif:", len(reactions_bag))

    print("counting motifs from multi motif glycans...")
    motif_counts_from_multi_motif_glycans = {}
    for primary_id in reactions_bag:
        for motif in reactions_bag[primary_id]["motifs"]:
            if motif in motif_counts_from_multi_motif_glycans:
                motif_counts_from_multi_motif_glycans[motif] += 1
            else:
                motif_counts_from_multi_motif_glycans[motif] = 1

    multi_glycan_motif_counts_sorted_by_value = dict(sorted(motif_counts_from_multi_motif_glycans.items(), key=lambda kv: kv[1], reverse=True))

    plt.figure(1)
    plt.bar(list(multi_glycan_motif_counts_sorted_by_value.keys())[1:22], list(multi_glycan_motif_counts_sorted_by_value.values())[1:22])
    plt.xticks(rotation=45)
    plt.title("top 20 motifs from multi-motif glycans")

    print("counting motifs from single motif glycans...")
    motif_counts_from_single_motif_glycans = {}
    for primary_id in reactions_bag:
        if len(reactions_bag[primary_id]["motifs"]) == 1:
            motif = reactions_bag[primary_id]["motifs"][0]
            if motif in motif_counts_from_single_motif_glycans:
                motif_counts_from_single_motif_glycans[motif] += 1
            else:
                motif_counts_from_single_motif_glycans[motif] = 1

    single_glycan_motif_counts_sorted_by_value = dict(sorted(motif_counts_from_single_motif_glycans.items(), key=lambda kv: kv[1], reverse=True))

    plt.figure(2)
    plt.bar(list(single_glycan_motif_counts_sorted_by_value.keys())[1:22], list(single_glycan_motif_counts_sorted_by_value.values())[1:22])
    plt.xticks(rotation=45)
    plt.title("top 20 motifs from single-motif glycans")

    plt.figure(3)
    ax = plt.subplot(211)
    ax.bar(list(multi_glycan_motif_counts_sorted_by_value.keys())[1:22], list(multi_glycan_motif_counts_sorted_by_value.values())[1:22],
            color='b', align="center")
    plt.xticks(rotation=45)
    plt.title("top 20 motifs from multi- and single-motif glycans")

    single_motif_counts_from_corresponding_multi_motif_glycans = []
    for x in list(multi_glycan_motif_counts_sorted_by_value.keys())[1:22]:
        if x in motif_counts_from_single_motif_glycans:
            single_motif_counts_from_corresponding_multi_motif_glycans.append(motif_counts_from_single_motif_glycans[x])
        else:
            single_motif_counts_from_corresponding_multi_motif_glycans.append(0)

    ax = plt.subplot(212)
    ax.bar(list(multi_glycan_motif_counts_sorted_by_value.keys())[1:22], single_motif_counts_from_corresponding_multi_motif_glycans,
            color='g', align="edge")
    plt.xticks(rotation=45)


    plt.figure(4)
    ax = plt.subplot(211)
    ax.bar(list(single_glycan_motif_counts_sorted_by_value.keys())[1:22], list(single_glycan_motif_counts_sorted_by_value.values())[1:22],
            color='b', align="center")
    plt.xticks(rotation=45)
    plt.title("top 20 motifs from multi- and single-motif glycans")
    multi_motif_counts_from_corresponding_single_motif_glycans = []
    for x in list(single_glycan_motif_counts_sorted_by_value.keys())[1:22]:
        if x in motif_counts_from_multi_motif_glycans:
            multi_motif_counts_from_corresponding_single_motif_glycans.append(motif_counts_from_multi_motif_glycans[x])
        else:
            multi_motif_counts_from_corresponding_single_motif_glycans.append(0)

    ax = plt.subplot(212)
    ax.bar(list(single_glycan_motif_counts_sorted_by_value.keys())[1:22], multi_motif_counts_from_corresponding_single_motif_glycans,
            color='g', align="edge")
    plt.xticks(rotation=45)
    plt.show()

    # print("removing glycans with zero or more than one motif from reactions_bag...")
    #
    # trimmed_reactions_bag = {}
    # for primary_id in reactions_bag:
    #     if reactions_bag[primary_id]["motifs"] != ['None'] and len(reactions_bag[primary_id]["motifs"]) == 1:
    #         trimmed_reactions_bag[primary_id] = reactions_bag[primary_id]
    #
    # print("removing", len(reactions_bag) - len(trimmed_reactions_bag),
    #       "glycans which had either zero or more than one motif from reactions_bag")
    # reactions_bag = trimmed_reactions_bag
    #
    # print("counting motifs...")
    # motif_counts = {}
    # for primary_id in reactions_bag:
    #     if reactions_bag[primary_id]["motifs"][0] in motif_counts:
    #         motif_counts[reactions_bag[primary_id]["motifs"][0]] += 1
    #     else:
    #         motif_counts[reactions_bag[primary_id]["motifs"][0]] = 1
    #
    # motif_counts_sorted_by_value = dict(sorted(motif_counts.items(), key=lambda kv: kv[1], reverse=True))
    #
    # plt.bar(motif_counts_sorted_by_value.keys(), motif_counts_sorted_by_value.values())
    # plt.xticks(rotation=90)
    # plt.show()

    # Filter the reactions_bag so that it only includes glycans containing motifs from the top 10 motifs (as found from single motif glycans)
    trimmed_reactions_bag = {}

    if include_zero_motif_glycans:
        top_10_motifs = list(single_glycan_motif_counts_sorted_by_value.keys())[0:11]
    else:
        top_10_motifs = list(single_glycan_motif_counts_sorted_by_value.keys())[1:11]

    for primary_id in reactions_bag:
        all_motifs_in_top10 = True
        for motif in reactions_bag[primary_id]["motifs"]:
            if motif not in top_10_motifs:
                all_motifs_in_top10 = False
                break

        if all_motifs_in_top10:
            trimmed_reactions_bag[primary_id] = reactions_bag[primary_id]

    return trimmed_reactions_bag
