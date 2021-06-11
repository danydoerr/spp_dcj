#!/usr/bin/env python3

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
from os.path import basename
import csv
import re


PAT_GENOME_FILE = re.compile('^(\w+)_GENOME\.tsv$')

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('genome', nargs='+', type=open,
            help='Converter for ZOMBI genome to OUR gene order format')
    args = parser.parse_args()

    out = stdout

    print('#species\tscaffold\tgene\torientation', file = out)

    for f in args.genome:
        m = PAT_GENOME_FILE.match(basename(f.name))
        if not m:
            print(('File name should match regular expression {}, but ' + \
                    'instead hads form {}, exiting').format( \
                    PAT_GENOME_FILE.pattern, basename(f.name)), file = stderr)
            exit(1)
        species = m.group(1)
        isHeader = True

        #
        # THIS CODE ASSUMES THAT INPUT GENOMES ARE SORTED BY GENE POSITION
        #
        prev_pos = None
        for line in csv.reader(f, delimiter = '\t'):
            if isHeader:
                isHeader = False
                continue
            if len(line) != 4:
                print('Unknown format, expected 4 columns, exiting', file =
                        stderr)
                exit(1)
            pos, fam, orient, geneID = line
            if prev_pos != None and prev_pos + 1 != int(pos):
                print(('Assumption that input genome {} is sorted by gene' + \
                        'position has been violated, exiting').format(f.name),
                        file = stderr)
                exit(1)
            prev_pos = int(pos)
            print('\t'.join((species, '1', '@'.join((species, '_'.join((fam,
                geneID)))), orient)), file = out)



