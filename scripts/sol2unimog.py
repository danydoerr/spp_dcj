#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from itertools import chain
from sys import stdout, stderr, exit
import csv

# import from own packages
import data_utils as du

import networkx as nx


def calculateDCJindelDistance(vars_, n_star):

    z = 0
    s = 0
    t = 0
    o = 0

    for var_, val in vars_.items():
        val = int(round(val))
        if var_.startswith('z') and val > .9:
            z += 1
        if var_.startswith('s') and val > .9:
            s += 1
        if var_.startswith('t') and val > .9:
            t += 1
    return n_star - z + t/2 +s


def _add_edge(G, id_, **attr):

    try:
        u, v, _ = id_.split('_')
    except ValueError:
        u, v = id_.split('_')
    exists = False
    if G.has_edge(u, v):
        for data in G[u][v].values():
            if data['id'] == id_:
                exists = True
                for k, v in attr.items():
                    data[k] = v
                break
    if not exists:
        G.add_edge(u, v, id=id_, **attr)

def checkConsistency(G):
    for v, vdata in G.nodes(data=True):
        id_ = vdata['id']
        if G.degree(v) != 2:
            raise Exception('{0}:{1}.{2} ({4}) has degree {3}'.format(
                *vdata['id'], G.degree(v), v))
        if vdata.get('o', 0) == 1 and vdata['type'] != du.VTYPE_CAP:
            raise Exception('{0}:{1}.{2} ({3}) is not a telomere but ' + \
                    'has o-label'.format(*vdata['id'], v))
        hasAdj = False
        hasExtId = False

        for u in G.neighbors(v):
            for data in G[u][v].values():
                if data['type'] == du.ETYPE_ADJ:
                    if hasAdj:
                        raise Exception(('{0}:{1}.{2} ({3}) is incident ' + \
                                'to two adjacency edges').format(*vdata['id'], \
                                v))
                    hasAdj = True
                elif hasExtId:
                    raise Exception(('{0}:{1}.{2} ({3}) is incident to two' + \
                            'ext/id edges').format(*vdata['id'], v))
                else:
                    hasExtId = True


def annotateGraph(G, id2ext):

    for v, vdata in G.nodes(data = True):
        ext = id2ext[v]
        vdata['id'] = ext
        if ext[1].startswith('t'):
            vdata['type'] = du.VTYPE_CAP
        else:
            vdata['type'] = du.VTYPE_EXTR

    for u, v, data in G.edges(data=True):
        gName1, g1, ext1 = G.nodes[u]['id']
        gName2, g2, ext2 = G.nodes[v]['id']
        if gName1 != gName2:
            data['type'] = du.ETYPE_EXTR
        elif data['id'].count('_') == 2:
            data['type'] = du.ETYPE_ID
        else:
            data['type'] = du.ETYPE_ADJ


def constructGraph(vars_):

    G = nx.MultiGraph()

    for var_, val in vars_.items():
        val = int(round(val))
        if val:
            if var_.startswith('x'):
                _add_edge(G, var_[1:])
            if var_.startswith('z') and val > .9:
                G.add_node(var_[1:var_.find('_')], z = 1)
            if var_.startswith('t') and val > .9:
                _add_edge(G, var_[1:var_.rfind('_')], t = 1)
            if var_.startswith('o') and val > .9:
                G.add_node(var_[1:], o = 1)

    return G

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sol_file', type=open,
            help='solution file of GUROBI optimizer')
    parser.add_argument('id_to_extremity_map', type=open,
            help='mapping between node IDs and extremities')

    args = parser.parse_args()

    out = stdout
    #
    # load data
    #

    # id-to-extremity mapping
    id2ext = dict()
    for line in csv.reader(args.id_to_extremity_map, delimiter = '\t'):
        if line:
            id2ext[line[0]] = tuple(line[1:])

    adjacenciesList, indelList, matchingList, _, vars_ = \
            du.parseSOL(args.sol_file, id2ext)

    G = constructGraph(vars_)
    annotateGraph(G, id2ext)
    # checkConsistency(G)


#    from matplotlib import pylab as plt
#    G = nx.Graph()
#    for gName, adjs in chain(adjacenciesList.items(), indelList.items()):
#        for (g1, ext1), (g2, ext2) in adjs:
#            G.add_edge((gName, g1, ext1), (gName, g2, ext2))
#    G.add_edges_from(map(lambda x: (x[0] + (du.EXTR_HEAD, ), x[1] +
#        (du.EXTR_HEAD, )), matchingList))
#    G.add_edges_from(map(lambda x: (x[0] + (du.EXTR_TAIL, ), x[1] +
#        (du.EXTR_TAIL, )), matchingList))
#    nx.draw(G)
#    plt.show()
#    import pdb; pdb.set_trace() 

    telomeres = {x: set() for x in adjacenciesList.keys()}
    for (gName, g, ext) in id2ext.values():
        if ext == du.EXTR_CAP:
            telomeres[gName].add(g)
    genomes = du.adjacencies2unimog(adjacenciesList, matchingList, telomeres)
    genomes.sort()
    n_star = len([x for x in matchingList if not x[1][1].startswith('t')]) + len([x for x in matchingList if x[1][1].startswith('t')])/2
    #p_star = len(matchingList) - n_star
    ddcj = calculateDCJindelDistance(vars_, n_star)
    print(f'DCJ INDEL distance is {ddcj}', file=stderr)

    # write genomes in UniMoG format 
    gene2str = lambda x: x[0] == du.ORIENT_NEGATIVE and f'-{x[1]}' or x[1]
    for gName, chrs in genomes:
        print(f'>{gName}', file = out)
        for ctype, chr_ in chrs:
            print(' '.join(map(gene2str, chr_)) + f' {ctype}', file = out)

