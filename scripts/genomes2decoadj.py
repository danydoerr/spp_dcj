#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit

# import from own packages
import data_utils as du
import zombi_utils as zu


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF, description =
            'Output adjacencies of leaf genomes of a given phylogeny in ' + \
                    'DeCo format')
    parser.add_argument('species_tree', type=open,
                        help='species tree of input genomes')
    parser.add_argument('genome', nargs = '+', type=open,
                        help='gene order in ZOMBI genome format')
    parser.add_argument('-c', '--circular', action = 'store_true',
                        help='make genome ciruclar')
    args = parser.parse_args()


    tree = du.parseTree(args.species_tree)
    leaves = [v for v, isLeaf in du.getLeaves(tree).items() if isLeaf]

    gene_orders = zu.parseGeneOrderFiles(args.genome, leaves)

    out = stdout
    for leaf, go in gene_orders.items():
        for i in range(len(go)-1):
            print('{0}_{1} {0}_{2} {3} {4}'.format(leaf, go[i][0],
                go[i+1][0], go[i][2], go[i+1][2]), file = out)

        if args.circular:
            print('{0}_{1} {0}_{2} {3} {4}'.format(leaf, go[-1][0],
                go[0][0], go[-1][2], go[0][2]), file = out)

    if not len(gene_orders):
        print('WARNING: no matching leaf genome found in phylogeny!', file =
                stderr)
