#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
from itertools import combinations
from os.path import basename, join
from glob import glob
import logging
import csv

from ReconciledTreeIO import recPhyloXML_parser

#
# global variables
#

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)


def parseReconciledGeneTrees(zombi_folder):
    parser = recPhyloXML_parser()
    trees = dict()
    for f in glob(join(zombi_folder, 'G', 'Gene_trees', '*_rec.xml')):
        LOG.info(f'parsing {f}')
        trees[basename(f).split('_')[0]] = parser.parse(f).recTrees[0]
    return trees


def getTrueOrthologs(tree, familyID):
    res = set()

    for u, v in combinations(tree.get_leaves(), 2):
        a = tree.get_common_ancestor(u, v)
        e = a.getEvents()
        if len(e) > 1:
            LOG.fatal(f'expected only one event, got more than two for ' + \
                    f'node {a.name} in tree {tree.name}')
            exit(1)
        if e[0].eventCode == 'S':
            ugenome, uid = u.name.split('_')
            vgenome, vid = v.name.split('_')
            id1 = f'{ugenome}@{familyID}_{uid}'
            id2 = f'{vgenome}@{familyID}_{vid}'
            res.add(id1 < id2 and (id1, id2) or (id2, id1))
    return res


def readLoc2GeneIDMap(zombi_folder):
    res = dict()
    for f in glob(join(zombi_folder, 'G', 'Genomes', '*_GENOME.tsv')):
        genome = basename(f).split('_', 1)[0]
        isHeader = True
        for line in csv.reader(open(f), delimiter='\t'):
            if isHeader:
                isHeader = False
                continue
            geneID = f'{genome}@{line[1]}_{line[3]}'
            locus = f'{genome}.{line[0]}'
            res[locus] = geneID
    return res


def parseGraph(data, loc2geneid, true_orthologs):

    fp = 0
    fn = 0
    tp = 0

    for line in csv.reader(data, delimiter='\t'):
        if line and line[0].startswith('#'):
            continue
        l1 = loc2geneid[line[0]]
        l2 = loc2genei[line[1]]
        orth = l1 < l2 and (l1, l2) or (l2, l1)
        if orth in true_orthologs:
            tp += 1
        else:
            fp += 1
    fn = len(true_orthologs) - tp - fp
    return tp, fp, fn


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('zombi_folder', type=str,
            help='ZOMBI ouput folder')
    parser.add_argument('graph', type=open,
            help='ProteinOrtho/PoFF graph output file')

    args = parser.parse_args()
    out = stdout

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    LOG.info(f'parse reconciled trees from {args.zombi_folder}/G/Gene_trees')
    trees = parseReconciledGeneTrees(args.zombi_folder)

    LOG.info(f'extracting ortholog assignments of {len(trees)} families')
    true_orthologs = set()
    for tid, t in trees.items():
        true_orthologs.update(getTrueOrthologs(t, tid))

    LOG.info(f'reading genome files from {args.zombi_folder}/G/Genomes')
    loc2geneid = readLoc2GeneIDMap(args.zombi_folder)

    LOG.info(f'parsing orthology assignments from {args.graph.name}')
    tp, fp, fn = parseGraph(args.graph, loc2geneid, true_orthologs)
    print(f'TP: {tp}', file=out)
    print(f'FP: {fp}', file=out)
    print(f'FN: {fn}', file=out)

