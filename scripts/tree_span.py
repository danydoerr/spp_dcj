#!/usr/bin/env python3

from sys import argv, exit, stdout
from newick_parser import parse_tree_iterator

def constructDistTable(tree):
    res = {}
    leaves = [(l, [(0, l.label)]) for l in tree.getLeaves()]
    while len(leaves) > 2:
        hasMatch = False
        for i in range(len(leaves)-1):
            for j in range(i+1, len(leaves)):
                if leaves[i][0].getAncestor() == leaves[j][0].getAncestor():
                    joinItems(leaves, i, j, res)
                    hasMatch = True
                    break;
            if hasMatch:
                break
    joinItems(leaves, 0, 1, res)
    return res

def joinItems(leaves, i, j, res):
    l1 = leaves.pop(j)
    l2 = leaves.pop(i)
    anc = l1[0].getAncestor()
    for d1, k in l1[1]:
        for d2, l in l2[1]:
            key = [k, l]
            key.sort()
            res[tuple(key)] = d1 + d2 + l1[0].length + l2[0].length
    k = [(d+l1[0].length, n) for d, n in l1[1]]
    l = [(d+l2[0].length, n) for d, n in l2[1]]
    leaves.append((anc, k+l))


def calculateSpan(tree):
    dists = constructDistTable(tree.subtree)
    return max(dists.values())

if __name__ == '__main__':
    if len(argv) != 2:
        print('\tusage: {} <NEWICK FILE>'.format(argv[0]), file = stdout)
        exit(1)
    for tree in parse_tree_iterator(open(argv[1])):
        print(calculateSpan(tree), file = stdout)

