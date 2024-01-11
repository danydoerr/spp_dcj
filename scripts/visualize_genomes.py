#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit

# import from third-party packages
import os, locale
locale.setlocale(locale.LC_ALL, 'en_US')
if not os.environ.get('DISPLAY', None):
    import matplotlib; matplotlib.use('Agg')
from matplotlib import pylab as plt
from itertools import chain
from collections import defaultdict
import networkx as nx
import csv

# import from own packages
import data_utils as du

COLOR_MAP = {du.ETYPE_ADJ: 'black', du.ETYPE_ID: 'gray' }
COLORS_COMP = ['orange', 'yellowgreen', 'lightseagreen', 'cornflowerblue', 'orchid']
LINESTYLE_MAP = {du.ETYPE_ADJ: 'solid', du.ETYPE_ID: 'dotted' }


def constructGenomeGraph(adjacencies, genes):
    ''' constructs genome graphs from input adjacencies'''

    G = nx.MultiGraph()

    for (g1, ext1), (g2, ext2) in adjacencies:
        if ext1 != du.EXTR_CAP or ext2 != du.EXTR_CAP:
            G.add_edge((g1, ext1), (g2, ext2), type=du.ETYPE_ADJ)

    extremities = list(chain(*zip(*adjacencies)))
    for gene in set(genes).difference(map(lambda x: x[0], filter(lambda x: x[1] == du.EXTR_CAP, extremities))):
        G.add_edge((gene, du.EXTR_HEAD), (gene, du.EXTR_TAIL),
                type=du.ETYPE_ID)
    return G


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('-i', '--highlight', type=open,
            help='highlight adjacencies in visualization')
    parser.add_argument('-o', '--omit_components', type=int, default=0,
            help='Omit components of with x genes or less')

    args = parser.parse_args()

    # load data
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies)

    highlights = None
    if args.highlight:
        highlights = du.parseAdjacencies(args.highlight)['adjacencies']

    # construct genomes and draw figures
    n_species = len(candidateAdjacencies['species'])
    plt.figure(figsize=(10, 5 * (n_species+1)//2))
    for i, (species, adjacencies) in enumerate(candidateAdjacencies[ \
            'adjacencies'].items()):
        # construct genome graph
        genes = candidateAdjacencies['genes'][species]
        G = constructGenomeGraph(adjacencies, genes)
        k = 0
        node_color = list()
        nodes = list()
        for c in sorted(nx.connected_components(G), key=len, reverse=True):
            if len(c) <= args.omit_components * 2:
                print(f'omitting component of size {len(c)//2} in plotting {species}', file=stderr)
                G.remove_nodes_from(c)
            else:
                nodes.extend(c)
                node_color.extend([COLORS_COMP[k % len(COLORS_COMP)]]* len(c))
                k += 1


        # plot figure
        plt.subplot((n_species + 1) // 2, 2, i+1)
        layout = nx.spring_layout(G, iterations=200)

        edges = list(G.edges(data = True))

        edge_color, edges = zip(*sorted(zip(list(map(lambda x: \
                COLOR_MAP[x[2]['type']], edges)), edges)))

#        line_style, edges = zip(*sorted(zip(list(map(lambda x: \
#                LINESTYLE_MAP[x[2]['type']], edges)), edges)))

        if highlights and species in highlights:
            highlight_edges  = list(x for x in highlights[species] if
                    G.has_edge(x[0], x[1]))
            nx.draw_networkx(G, pos=layout, edgelist=highlight_edges, \
                    edge_color = 'red', with_labels=False, node_size=0, \
                    width=4)

        nx.draw_networkx(G, pos=layout, nodelist=nodes, node_color=node_color, \
                         edgelist=edges, edge_color = edge_color, \
                         with_labels=False, font_size=1, node_size=2, \
                         width=list(map(lambda x: x == 'black' and 2 or 1, edge_color)))
#        nx.draw_networkx(G, pos=layout, edgelist=edges, edge_color='black', \
#                font_size=10, node_color='white', node_size=5, \
#                width=list(map(lambda x: x == 'dotted' and 0.4 or 0.5, \
#                line_style)), style=line_style, font_color='black', \
#                font_family='serif', labels=dict(((v, '.') for v in G.nodes)),
#                with_labels=True, verticalalignment='baseline',)
        plt.title(species)
        ax = plt.gca()
        ax.axis('off')
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

    plt.tight_layout(h_pad=4)
    # plot to output
    with os.fdopen(stdout.fileno(), 'wb', closefd=False) as out:
        plt.savefig(out, format='pdf')
        out.flush()


