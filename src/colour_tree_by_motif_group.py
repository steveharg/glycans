from Bio import Phylo
import pickle
import seaborn as sns
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')


def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents


def colour_internal_node(parents, node):
    child1 = node.clades[0]
    child2 = node.clades[1]
    child1_colour = child1.color
    child2_colour = child2.color
    if child1_colour is None or child2_colour is None:
        return
    elif child1_colour.to_hex() == child2_colour.to_hex():
        node.color = child1_colour
    else:
        return

    node_parent = parents[node]

    if node_parent is not None:
        colour_internal_node(parents, node_parent)

if __name__ == "__main__":

    include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')
    use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
    tree_with_motif_names_pickle_file = config['DEFAULT']['TreeWithMotifNamesPickleFile']
    tree_coloured_by_motif_names_phyloxml_file = config['DEFAULT']['TreeColouredByMotifNamesPhyloXmlFile']

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

    tree = pickle.load(open(data_dir + "/" + tree_with_motif_names_pickle_file, "rb"))

    # get distinct motif groups, and correct for lack of space before opening bracket
    distinct_motif_groups = {}
    for term in tree.get_terminals():
        name_split = term.name.split('(')
        print(name_split[0])
        if name_split[0] not in distinct_motif_groups:
            distinct_motif_groups[name_split[0]] = ''
        term.name = name_split[0] + ' (' + name_split[1]

    print('')
    print('-----------------')
    print('')

    for motif_group in sorted(distinct_motif_groups):
        print(motif_group)

    print('')
    print('-----------------')
    print('')

    for term in tree.get_terminals():
        print(term)

    pal = sns.color_palette('colorblind', len(distinct_motif_groups))

    hex_pal = pal.as_hex()

    i = 0
    for distinct_motif_group in distinct_motif_groups:
        distinct_motif_groups[distinct_motif_group] = hex_pal[i]
        i += 1

    # Colour the terminal nodes by motif group
    for term in tree.get_terminals():
        motif_group = term.name.split(' (')[0]
        term.color = distinct_motif_groups[motif_group]

    # Colour internal nodes where they have same colour children
    parents = all_parents(tree)
    for term in tree.get_terminals():
        term_parent = parents[term]
        colour_internal_node(parents, term_parent)

    # save the tree
    Phylo.write(tree, data_dir + "/" + tree_coloured_by_motif_names_phyloxml_file, 'phyloxml')