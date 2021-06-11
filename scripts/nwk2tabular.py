#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit

# import from own packages
from newick_parser import parse_tree, Branch


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree', type=open, help='tree in newick format')
    args = parser.parse_args()

    tree = parse_tree(args.tree)

    out = stdout

    print('#Node\tParent', file = out)
    for v in tree.getNodes():
        u = v.getAncestor()
        if u != None:
            print('\t'.join((v.label, u.label)), file = out)


