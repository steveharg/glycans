import pickle
import networkx as nx
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

use_reaction_quantities = True if config['DEFAULT']['ReactionCountMethod'] == 'List' else False
include_zero_motif_glycans = config['DEFAULT'].getboolean('IncludeZeroMotifGlycans')

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

threshold = config['DEFAULT'].getint('Threshold', fallback=90)
G = pickle.load(open(data_dir + "/G_" + str(threshold) + ".p", "rb"))
reactions_collection = pickle.load(open(data_dir + "/reactions_bag.p", "rb"))

connected_components = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
conn_comp_freqs = {}

for conn_comp in connected_components:
    if len(conn_comp) in conn_comp_freqs:
        conn_comp_freqs[len(conn_comp)] += 1
    else:
        conn_comp_freqs[len(conn_comp)] = 1

print(conn_comp_freqs)

nodes = G.nodes

zero_conn_nodes = 0
degree_ones = False
node_conn_freqs = {}
for node in nodes:
    if G.degree[node] == 2:
        zero_conn_nodes += 1

    # print("node:", node)
    # print(G.adj[node])
    # print(G.degree[node])
    if G.degree[node] == 1:
        degree_ones = True


print("zero_conn_nodes:", zero_conn_nodes, " out of ", len(nodes), "nodes")

# zero_motif_glycans_exist = False
# for reaction in reactions_collection:
#     # print(reactions_collection[reaction]['motifs'])
#     print(reaction)
#     print(len(reactions_collection[reaction]['motifs']))
#     if len(reactions_collection[reaction]['motifs']) == 0:
#         zero_motif_glycans_exist = True
#
# print(zero_motif_glycans_exist)