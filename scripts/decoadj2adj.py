#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
from itertools import chain, product, combinations
from collections import deque
import csv
import re

# import from own packages
import data_utils as du
import zombi_utils as zu
from newick_parser import parse_tree, Leaf


PAT_FAM = re.compile('^>family (\d+)$')
PAT_ANC_GENE = re.compile('^(\d+)\|(\d+)$')

def readDecoAncestralGenes(data):
    res = dict()

    for line in csv.reader(data, delimiter = ' '):
        if len(line) > 2:
            m = PAT_ANC_GENE.match(line[1])
            if m:
                fam_id = int(m.group(1))
                if fam_id not in res:
                    res[fam_id] = list()
                res[fam_id].append((line[1], tuple(line[2:])))
    mx = max(res.keys())
    return [x in res and res[x] or [] for x in range(mx+1)]


def readDecoRecGeneTrees(data):

    res = dict()

    line = data.readline()
    while line:
        m = PAT_FAM.match(line)
        if not m:
            raise SyntaxError(('expected line match expression "%s", but ' + \
                    'received "%s"') %(m.pattern, line))
        fam_id = int(m.group(1))
        tree = parse_tree(data.readline())
        res[fam_id] = tree

        try:
            line = data.readline()
        except StopIteration:
            line = None

    mx = max(res.keys())
    return [x in res and res[x] or [] for x in range(mx+1)]


def getAdjacenciesFromDeCo(data, deco2zombi):

    adjacencies = dict()
    weights = dict()
    for line in csv.reader(data, delimiter = ' '):
        if line[1].find('|') >= 0:
            fam_id1 = int(PAT_ANC_GENE.match(line[1]).group(1))
            fam_id2 = int(PAT_ANC_GENE.match(line[2]).group(1))
            line[1] = deco2zombi[fam_id1][line[1]]
            line[2] = deco2zombi[fam_id2][line[2]]

        species, pos1 = line[1].split('@', 1)
        sp2,     pos2 = line[2].split('@', 1)
        if species != sp2:
            print(('FATAL: excpected extremities of adjacency to be part ' + \
                    'of the same genome, but got {}-{}').format(line[1], \
                    line[2]))
            exit(1)
        sign1 = line[3]
        sign2 = line[4]

        ext1 = (pos1, du.SIGN2EXT_1[sign1])
        ext2 = (pos2, du.SIGN2EXT_2[sign2])

        if ext1 > ext2:
            ext1, ext2 = ext2, ext1
        if species not in adjacencies:
            adjacencies[species] = set()
        adjacencies[species].add((ext1, ext2))
        weights[(ext1, ext2)] = float(line[-1])

    return adjacencies, weights


def mapDecoIDsToZombi(ancestral_genes, tree):

    res = dict((x[0], None) for x in ancestral_genes)

    # intialize map with leafset of gene tree
    for v in tree.getLeaves():
        label = v.label[:v.label.find('|')]
        res[label] = label

    # value is a triple of: 1) node name, 2) node object
    desc2anc = dict((v.label[:v.label.find('|')], \
            (v.getAncestor().label[:v.getAncestor().label.find('|')], \
            v.getAncestor())) for v in tree.getNodes() if v.getAncestor())

    # bottom-up processing of gene tree
    queue = deque(ancestral_genes)

    # failsafe: loose upper bound on number of iterations
    N = (len(queue) * (len(queue)+1)) // 2
    c = 0
    while queue and c < N:
        v, descendants = queue.popleft()
        # ancestral gene deletions lead to genes not being part of the result
        # set. They will be mapped through their association to surviving
        # siblings
        if all(map(lambda x: res[x], descendants)):
            # identify ancestor
            # first, identify ancestor with surviving descendants
            ancestors = set()
            for desc in descendants:
                anc = desc2anc[res[desc]]
                while '|Dup|' in anc[1].label:
                    anc = desc2anc[anc[0]]
                ancestors.add(anc[0])
            if len(ancestors) != 1:
                print('FATAL: assumed 1 ancestor, found {}'.format(ancestors),
                        file = stderr)
                exit(1)
            # - map deco id to zombi id
            res[v] = ancestors.pop()
        else:
            queue.append((v, descendants))
        c += 1

    if queue:
        print('FATAL: number of iteration exceeded in method' + \
                'mapDecoIDsToZombi', file = stderr)
        exit(1)

    return res


def complementCandidateAdjacencySet(adjacencies, weights):

    # XXX following assumptions are made:
    # - genome is circular
    # - single-gene chromosomes are not permitted

    # process each genome separately
    for species, adjs in adjacencies.items():
        # count number of incident candidate adjacencies of each extremity
        extremities = dict()
        for (gene1, ext1), (gene2, ext2) in adjs:
            # initialize data structure for unvisited genes
            for e in product((gene1, gene2), (du.EXTR_HEAD, du.EXTR_TAIL)):
                if e not in extremities:
                    extremities[e] = 0
            # count
            extremities[(gene1, ext1)] += 1
            extremities[(gene2, ext2)] += 1

        unsaturated = [e for e, c in extremities.items() if not c]

        for ext1, ext2 in combinations(unsaturated, 2):
            if ext1[0] != ext2[0]:
                if ext1 > ext2:
                    ext1, ext2 = ext2, ext1
                adjacencies[species].add((ext1, ext2))
                weights[(ext1, ext2)] = 0.0


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF, description =
            'Output adjacencies of leaf genomes of a given phylogeny in ' + \
                    'DeCo format')
    parser.add_argument('deco_adjacencies', type=open,
                        help='inferred adjacencies from DeCo')
    parser.add_argument('deco_genes', type=open,
                        help='mapping of inferred genes to descendants')
    parser.add_argument('deco_reconciliations', type=open,
                        help='file containing the reconciliated gene ' + \
                                'trees of DeCo in newick format' )
    parser.add_argument('-c', '--complement', action='store_true',
                        help='complement adjacency set to saturate open ' + \
                                'extremities')
    args = parser.parse_args()

    out = stdout

    genes = readDecoAncestralGenes(args.deco_genes)
    gene_trees = readDecoRecGeneTrees(args.deco_reconciliations)

    deco2zombi = list()
    for gene_set, gene_tree in zip(genes, gene_trees):
        deco2zombi.append(mapDecoIDsToZombi(gene_set, gene_tree))

    #
    # process adjacencies
    #
    # get candidate adjacencies from DeCo to output
    adjacencies, weights = getAdjacenciesFromDeCo(args.deco_adjacencies,
            deco2zombi)

    # complement adjacency set
    if args.complement:
        complementCandidateAdjacencySet(adjacencies, weights)

    du.writeAdjacencies(adjacencies, weights, stdout)
