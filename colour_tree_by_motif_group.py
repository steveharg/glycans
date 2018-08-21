from Bio import Phylo
from Bio.Phylo.PhyloXML import Property
import pickle
import seaborn as sns

data_dir = 'data_top10_motifs'

tree = pickle.load(open(data_dir + "/tree_with_motif_names.p", "rb"))
# trees = Phylo.parse('/run/media/sah217/DATA/dev/python/glycans/data_top10_motifs/tree_with_motif_names_colour.xml', 'phyloxml')
# tree = next(trees)

# get distinct motif groups, and correct for lack of space before opening bracket
distinct_motif_groups = {}
for term in tree.get_terminals():
    name_split = term.name.split('(')
    print(name_split[0])
    if name_split[0] not in distinct_motif_groups:
        distinct_motif_groups[name_split[0]] = ''
    term.name = name_split[0] + ' (' + name_split[1]

# print(term_names)

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
# pal = sns.color_palette('RdBu', len(distinct_motif_groups))

hex_pal = pal.as_hex()  # do we need to convert to hex???

i = 0
for distinct_motif_group in distinct_motif_groups:
    distinct_motif_groups[distinct_motif_group] = hex_pal[i]
    i += 1

# Colour the terminal nodes by motif group
for term in tree.get_terminals():
    motif_group = term.name.split(' (')[0]
    # term.properties = [Property(value=distinct_motif_groups[motif_group],
    #                            ref='style:font_color', applies_to='node', datatype="xsd:token")]
    term.color = distinct_motif_groups[motif_group]

# save the tree
Phylo.write(tree, data_dir + "/tree_coloured_by_motif_names.xml", 'phyloxml')