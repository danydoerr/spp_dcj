#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from sys import stdout, stderr, exit
from itertools import product, combinations, chain, repeat
from functools import reduce
from collections import defaultdict
from math import ceil
import logging
import csv

# import from third-party packages
import networkx as nx

# import from own packages
import data_utils as du
import spp_dcj as spp

#
# global variables
#

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def writeConstraints(G, out):

    vars_ = set()
    for v, vdata in G.nodes(data = True):
        const = ''
        for u in G.neighbors(v):
            for edata in G[u][v].values():
                if edata['type'] in {du.ETYPE_ADJ, du.ETYPE_ID}:
                    if const:
                        const += ' + '
                    var = 'x{}'.format(edata['id'])
                    vars_.add(var)
                    const += var
        if vdata['type'] == du.VTYPE_EXTR:
            const += ' = 1'
        else:
            const += f' - o{v} = 0'
            vars_.add(f'o{v}')
        print(const, file=out)
    return vars_


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=open,
            help='phylogenetic tree as parent-child relation table')
    parser.add_argument('candidateAdjacencies', type=open,
            help='candidate adjacencies of the genomes in the phylogeny')
    parser.add_argument('-m', '--output_id_mapping', type=FileType('w'),
            help='writs a table with ID-to-gene extremity mapping')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    # load & process input data
    LOG.info('loading species tree from {}'.format(args.tree.name))
    speciesTree = du.parseTree(args.tree)

    LOG.info('loading candidate adjacencies from {}'.format(
        args.candidateAdjacencies.name))
    candidateAdjacencies = du.parseAdjacencies(args.candidateAdjacencies)

    # add telomeres
    telomeres = spp.identifyCandidateTelomeres(candidateAdjacencies,
            spp.ADJ_TRUST_THRESHOLD)

    # construct adjacency graphs
    genes = candidateAdjacencies['genes']
    adjacencies = candidateAdjacencies['adjacencies']
    weights = candidateAdjacencies['weights']

    ext2id = du.IdManager()
    LOG.info(('constructing relational diagrams for all {} branches of ' + \
            'the tree').format(len(speciesTree)))
    relationalDiagrams = du.constructRelationalDiagrams(speciesTree,
            adjacencies, telomeres, weights, defaultdict(dict), genes, ext2id)

    graphs = relationalDiagrams['graphs']

    # construct & output ILP
    out = stdout

    print('Maximize\nobj:\nSubject To', file=out)
    vars_ = set()
    for graph in graphs.values():
        vars_.update(writeConstraints(graph, out))

    print('Binaries')
    for var in vars_:
        print(var, file=out)


    if args.output_id_mapping:
        LOG.info('writing ID-to-gene extremity mapping to {}'.format(
            args.output_id_mapping.name))
        idMap = ext2id.getMap()
        out_table = list()
        for k, v in idMap.items():
            out_table.append((str(v), k[0], k[1][0], k[1][1]))
        out_table.sort(key = lambda x: int(x[0]))
        print('\n'.join(map(lambda x: '\t'.join(x), out_table)),
                file=args.output_id_mapping)

    LOG.info('DONE')
    out.write('end\n')

