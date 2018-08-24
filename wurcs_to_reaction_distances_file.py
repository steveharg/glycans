import pickle
from parse_glytoucan_results import parse_glytoucan_results
from calc_jaccard_distances import calc_jaccard_distances
from filter_by_common_motifs import filter_by_common_motifs
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')
    num_datafile_lines = config['DEFAULT'].getint('NumDatafileLines')  # None means use all datafile lines
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    glycans_json_file = config['DEFAULT']['GlycansJsonFile']
    reactions_collection_pickle_file = config['DEFAULT']['ReactionsCollectionPickleFile']
    heatmaps_pickle_file = config['DEFAULT']['HeatmapsPickleFile']
    col_names_pickle_file = config['DEFAULT']['ColNamesPickleFile']

    if include_zero_motif_glycans:
        if use_reaction_quantities:
            data_dir = config['DEFAULT']['QuantifiedReactionsZeroMotifsDataDir']
        else:
            data_dir = config['DEFAULT']['SetOfReactionsZeroMotifsDataDir']
    else:
        if use_reaction_quantities:
            data_dir = config['DEFAULT']['QuantifiedReactionsDataDir']
        else:
            data_dir = config['DEFAULT']['SetOfReactionsDataDir']

    reactions_collection = parse_glytoucan_results(glycans_json_file, num_datafile_lines, use_reaction_quantities)
    reactions_collection = filter_by_common_motifs(reactions_collection, include_zero_motif_glycans)
    heatmaps, column_names = calc_jaccard_distances(reactions_collection, use_reaction_quantities)

    reaction_distances = [reactions_collection, heatmaps, column_names]

    pickle.dump(reactions_collection, open(data_dir + "/" + reactions_collection_pickle_file, "wb"))
    pickle.dump(heatmaps, open(data_dir + "/" + heatmaps_pickle_file, "wb"))
    pickle.dump(column_names, open(data_dir + "/" + col_names_pickle_file, "wb"))
