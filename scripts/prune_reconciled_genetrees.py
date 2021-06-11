#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
from ReconciledTreeIO import recPhyloXML_parser as recPhyloParser



def pruneTree(tree):

    flags = dict((l.name, l.getEvent(0).eventCode != 'L') for l in tree.get_leaves())

    # bottom-up traversal: identify sub-trees that are entirely non-surviving
    for node in tree.traverse('postorder'):
        anc = node.up
        if anc:
            if anc.name not in flags:
                flags[anc.name] = False
            flags[anc.name] |= flags[node.name]

    # top-down traversal: remove sub-trees that are entirely non-surviving
    for node in tree.traverse('preorder'):
        if node.getEvent(0).eventCode in ('S', 'D') and not flags[node.name]:
            node.up.remove_child(node)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF, description =
            'Remove subtrees that did not survive in extant species.')
    parser.add_argument('reconciled_gene_tree', type=open,
                        help='reconciled gene tree in recPhyloXML format')
    args = parser.parse_args()

    out = stdout

    parser = recPhyloParser()
    forest =  parser.parse(args.reconciled_gene_tree.name)

    for tree in forest:
        pruneTree(tree)
        print(tree.getTreeRecPhyloXML(), file = out)

