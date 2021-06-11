#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
import csv
import re

# import from own packages
import data_utils as du


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sol_file', type=open,
            help='solution file of GUROBI optimizer')
    parser.add_argument('id_to_extremity_map', type=open,
            help='mapping between node IDs and extremities')

    args = parser.parse_args()

    #
    # load data
    #

    # id-to-extremity mapping
    id2ext = dict()
    for line in csv.reader(args.id_to_extremity_map, delimiter = '\t'):
        if line:
            id2ext[line[0]] = tuple(line[1:])

    adjacenciesList, _, weightsDict, _, _, _ = du.parseSOL(args.sol_file, id2ext)

    # write adjacencies
    du.writeAdjacencies(adjacenciesList, weightsDict, stdout)
