#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from sys import stdout, stderr, exit
import logging
import csv
import re

#
# global variables
#

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

PAT_T = re.compile('^t(\d+)_(\d+)_(\d+)$')

def check_c10(active_vars, species_tree, id2ext):

    for var in active_vars:
        m = PAT_T.match(var)
        if m:
            id1, id2, gid = map(int, m.groups())
            # check that the edge is not an active adjacency edge:
            if 'x{}_{}'.format(id1, id2) in active_vars and \
                    (id2ext[id1][:2] != id2ext[id2][:2] or \
                    id2ext[id1][0] != species_tree[i][1]):
                LOG.error(('variable {} (edge {}-{}) violates constraint ' + \
                        'c10').format(var, id2ext[id1], id2ext[id2]))

#    for u, v, data in G.edges(data = True):
#        if data['type'] != du.ETYPE_ADJ or G.nodes[u]['genome'] != 0:
#            out.write('t{}_{} = 0\n'.format(data['id'], i))



if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sol_file', type=open,
            help='Solution of SPP_DCJ in Gurobi SOL format')
    parser.add_argument('id_to_extremity_map', type=open,
            help='mapping between node IDs and extremities')
    parser.add_argument('species_tree', type=open,
            help='Sepcies tree in tabular child-parent format')
    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    # load & process input data
    LOG.info('loading id-to-extremity mapping from {}'.format(
        args.id_to_extremity_map.name))

    # id-to-extremity mapping
    id2ext = dict()
    for line in csv.reader(args.id_to_extremity_map, delimiter = '\t'):
        if line:
            id2ext[int(line[0])] = tuple(line[1:])

    LOG.info('loading species tree from {}'.format(args.species_tree.name))

    species_tree = list()
    for line in csv.reader(args.species_tree, delimiter = '\t'):
        if not line or line[0].startswith('#'):
            continue
        species_tree.append(tuple(line[:2]))
    species_tree.sort()

    LOG.info('loading ILP solution from {}'.format(args.sol_file.name))
    active_vars = set()

    for line in csv.reader(args.sol_file, delimiter = ' '):
        if not line or line[0].startswith('#'):
            continue
        if int(line[1]) == 1:
            active_vars.add(line[0])

    check_c10(active_vars, species_tree, id2ext)

    LOG.info('DONE')
