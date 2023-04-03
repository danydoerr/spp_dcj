#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit, maxsize
from itertools import combinations, chain
from collections import defaultdict

from matplotlib import pylab as plt

# import from own packages
import data_utils as du

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('unimog_file', type=open,
            help='genomes in UNIMOG format')

    args = parser.parse_args()

    #
    # load data
    #
    genomes = du.parseUniMoG(args.unimog_file)

    for (A, chrsA), (B, chrsB) in combinations(genomes, 2):
        f, ax = plt.subplots(figsize=(10, 10))


        # construct positional map for genome B
        posB = defaultdict(list)
        for p, g in enumerate(map(lambda x: x[1].split('_')[0], chain(*map(lambda x: x[1], chrsB)))):
            posB[g].append(p)

        X = list()
        Y = list()
        for x, (_, g) in enumerate(chain(*map(lambda x: x[1], chrsA))):
            for y in posB.get(g.split('_')[0]):
                X.append(x)
                Y.append(y)

        ax.scatter(X, Y, marker='.')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(visible=False)
        ax.set_ylabel(B)
        ax.set_xlabel(A)
        ax.set_title(f'dotplot {A} vs {B}')

        plt.savefig(f'{A}_{B}.pdf', format='pdf')
        plt.close()
