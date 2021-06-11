#!/usr/bin/env python3

from sys import stdout,stderr,exit
from optparse import OptionParser
from newick_parser import parse_tree_iterator, Branch
from tree_span import calculateSpan
from copy import deepcopy

def rescale_absolute(tree, max_length):
    span = calculateSpan(tree)
    s = max_length/float(span)
    return rescale(tree, s)

def rescale(tree, scale_factor):

    res = deepcopy(tree)

    stack = list()
    stack.append(res.subtree)

    while stack:
        x = stack.pop()
        if x.length:
            x.length *= scale_factor
        if type(x) == Branch:
            stack.extend(x.subtrees)

    return res


if __name__ == '__main__':
    usage = 'usage: %prog [options] <NEWICK FILE>'
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--scale_factor', dest='scale_factor',
            help='Scale factor of distances in tree',
            type=float, default=0, metavar='FLOAT')
    parser.add_option('-a', '--absolute_length', dest='absolute',
            help='Absolute length of maximal distance in tree',
            type=float, default=0, metavar='FLOAT')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit(1)
    if not ((options.absolute > 0) ^ (options.scale_factor > 0)):
        print('!! Specify either scale factor or absolute length with ' + \
                'strictly positive number', file = stderr)
        exit(1)

    for tree in parse_tree_iterator(open(args[0])):
        if options.absolute > 0:
            print(rescale_absolute(tree, options.absolute), file = stdout)
        else:
            print(rescale(tree, options.scale_factor), file = stdout)


