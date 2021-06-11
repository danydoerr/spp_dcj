#!/usr/bin/env python3
from sys import argv, stdout
import csv
import os

import matplotlib.pylab as plt

X = [int(row[0]) for row in csv.reader(open(argv[1]), delimiter='\t')]
spp = [float(row[1]) for row in csv.reader(open(argv[1]), delimiter='\t')]
ding = [float(row[1]) for row in csv.reader(open(argv[2]), delimiter='\t')]

plt.figure()
plt.plot(X, spp, 'x')
plt.plot(X, ding, 'o')
plt.yscale('log')
plt.legend(('SPP DCJ', 'DING'))
plt.ylabel('seconds')
plt.xlabel('# rearrangements + indels')
plt.title('Runtime of SPP DCJ vs DING (n = 1000, 100 linear genomes)')
with os.fdopen(stdout.fileno(), 'wb', closefd=False) as out:
    plt.savefig(out, format='pdf')
    out.flush()
