#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
import data_utils as du
import csv

#
# COMPARING TWO SOLUTIONS
#

def SDDistance(adjacencies1,adjacencies2,speciesList,noiseInLeaves=False,leavesDict={}):
    ''' Compute the number of TP, FP, FN per species assuming adjacencies 2 is the ground truth'''
    adjacencies1Aux = {}
    adjacencies2Aux = {}
    TP,FP,FN        = {species:0 for species in speciesList},{species:0 for species in speciesList},{species:0 for species in speciesList}
    for species in speciesList:
        if species not in leavesDict:
            continue
        if ((not leavesDict[species]) or (noiseInLeaves and leavesDict[species])):
            adjacencies1Aux[species] = []
            adjacencies2Aux[species] = []
            for [ext1,ext2] in adjacencies1[species]:
                if ext1<ext2:
                    adjacencies1Aux[species].append([ext1,ext2])
                else:
                    adjacencies1Aux[species].append([ext2,ext1])
            for [ext1,ext2] in adjacencies2[species]:
                if ext1<ext2:
                    adjacencies2Aux[species].append([ext1,ext2])
                else:
                    adjacencies2Aux[species].append([ext2,ext1])
            commonAdjacencies = [[ext1,ext2] for [ext1,ext2] in adjacencies1Aux[species] if [ext1,ext2]  in adjacencies2Aux[species]]
            specific1Adjacencies = [[ext1,ext2] for [ext1,ext2] in adjacencies1Aux[species] if [ext1,ext2] not in adjacencies2Aux[species]]
            specific2Adjacencies = [[ext1,ext2] for [ext1,ext2] in adjacencies2Aux[species] if [ext1,ext2] not in adjacencies1Aux[species]]
            TP[species] = len(commonAdjacencies)
            FP[species] = len(specific1Adjacencies)
            FN[species] = len(specific2Adjacencies)

    precision = sum([float(TP[species]) for species in speciesList]) / sum([float(TP[species]+FP[species]) for species in speciesList])
    recall    = sum([float(TP[species]) for species in speciesList]) / sum([float(TP[species]+FN[species]) for species in speciesList])

    return({'TP':TP,'FP':FP,'FN':FN,'precision':precision,'recall':recall})

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=open,
                        help='phylogenetic tree as parent-child relation table')
    parser.add_argument('predictedAdjacencies', type=open,
                        help='predicted of the genomes in the phylogeny')
    parser.add_argument('trueAdjacencies', type=open,
                        help='true adjacencies of the genomes in the phylogeny')
    parser.add_argument('-l', '--noiseInLeaves', action = 'store_true',
                        help='adds noise to the leaves')
    args = parser.parse_args()

    speciesTree          = du.parseTree(args.tree)
    leavesDict           = du.getLeaves(speciesTree)
    predictedAdjacencies = du.parseAdjacencies(args.predictedAdjacencies)
    trueAdjacencies      = du.parseAdjacencies(args.trueAdjacencies)
    speciesList          = trueAdjacencies['species']

    RES = SDDistance(predictedAdjacencies['adjacencies'],trueAdjacencies['adjacencies'],speciesList,
                     noiseInLeaves=args.noiseInLeaves,leavesDict=leavesDict)

    print('#Precision\tRecall\n',RES['precision'],'\t',RES['recall'])
    print('#Species\tTP\tFP\tFN')
    for species in speciesList:
        if species in leavesDict and ((not leavesDict[species]) or \
                (args.noiseInLeaves and leavesDict[species])):
            print(species,'\t',RES['TP'][species],'\t',RES['FP'][species],'\t',RES['FN'][species])

