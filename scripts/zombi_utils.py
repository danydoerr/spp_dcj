#!/usr/bin/env python3

from sys import stdout, stderr, exit
from os.path import basename
import csv
import re


import data_utils as du


PAT_GENOME_NAME = re.compile('^(\w+)_GENOME.tsv$')


def readGeneOrder(data):
    ''' reads ZOMBI gene order table assuming the followign columns:
    POSITION    GENE_FAMILY     ORIENTATION     GENE_ID'''

    res = list()

    isHeader = True
    for line in csv.reader(data, delimiter = '\t'):
        if isHeader:
            isHeader = False
            continue
        res.append((int(line[0]), int(line[1]), line[2], int(line[3])))

    return res


def geneorder2adjacencies(gene_order, isCircular = False):
    ''' constructs list of adjacencies from ZOMBI gene order list '''

    res = list()

    for i in range(len(gene_order)-1):
        gene1 = gene_order[i]
        gene2 = gene_order[i+1]
        ext1 = ('_'.join((str(gene1[1]), str(gene1[3]))),
                du.SIGN2EXT_1[gene1[2]])
        ext2 = ('_'.join((str(gene2[1]), str(gene2[3]))),
                du.SIGN2EXT_2[gene2[2]])
        if ext1 > ext2:
            ext1, ext2 = ext2, ext1
        res.append((ext1, ext2))

    if isCircular:
        gene1 = gene_order[-1]
        gene2 = gene_order[0]
        ext1 = ('_'.join((str(gene1[1]), str(gene1[3]))),
                du.SIGN2EXT_1[gene1[2]])
        ext2 = ('_'.join((str(gene2[1]), str(gene2[3]))),
                du.SIGN2EXT_2[gene2[2]])
        if ext1 > ext2:
            ext1, ext2 = ext2, ext1
        res.append((ext1, ext2))

    return res


def parseGeneOrderFiles(genome_files, filterGenomes = None):
    ''' reads ZOMBI gene orders from a list of files, but only from those that
    correspond to leaf genomes'''

    res = dict()
    for g in genome_files:
        name = basename(g.name)
        m = PAT_GENOME_NAME.match(name)
        if m:
            name = m.group(1)
        else:
            name = name.rsplit('.', 1)[0]

        if filterGenomes == None or name in filterGenomes:
            noOutput = False
            go = readGeneOrder(g)
            res[name] = go
    return res


def mapLeafsetsToNodes(tree):
    ''' given a tree, the function returns the unique mapping between leaf sets
    and vertices such that each leaf set L maps to the vertex whose induced
    subtree contains exactly L'''

    # result is a dictionary that maps leaf sets to vertices of the tree
    res = dict()

    parent = dict(tree)
    leaves = du.getLeaves(tree)

    # leaves are mapped onto themselves
    node2leafset = dict()

    # initialization:
    # iterate over all vertices (yes, leaves contains *all* vertices of the
    # tree) and intialize table 
    for v in leaves:
        node2leafset[v] = list()


    # tree bottom-up traversal
    for leaf, isLeaf in leaves.items():
        if isLeaf:
            # bottom-up
            queue = [leaf]
            while queue:
                v = queue.pop()
                node2leafset[v].append(leaf)
                # proceed to next higher level in tree unless we reached the
                # root node
                if v in parent:
                    queue.append(parent[v])

    for v, leafset in node2leafset.items():
        res[tuple(sorted(leafset))] = v

    return res


def readGeneContent(data):

    res = dict()

    for line in csv.reader(data, delimiter = '\t'):
        species = line[0]
        res[species] = line[1:]


if __name__ == '__main__':
    print('this is a module', file = stdout)

