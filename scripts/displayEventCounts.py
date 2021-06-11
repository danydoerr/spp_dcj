#!/usr/bin/env python3
#
# copy and paste from Aniket Mane' code
#

# import from built-in packages
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, FileType
from sys import stdout, stderr, exit
from os.path import join
from os import listdir
import csv

import pandas as pd


def readTreeEvents(data):
    ''' Obtaining ancestral and extant species '''
    ancestral = []
    extant = []
    rev_tree_dict = {}

    for line in csv.reader(data, delimiter = '\t'):
        if line[0].startswith('#'):
            continue
        if line[1] == 'S':
            species = line[2].split(';')[0]
            ancestral.append(species)

            ch1, ch2 = line[2].split(';')[1], line[2].split(';')[2]
            rev_tree_dict[ch1] = species
            rev_tree_dict[ch2] = species
        elif line[1] == 'F':
            species = line[2]
            extant.append(species)
    return ancestral, extant, rev_tree_dict


def readEventCounts(events_folder, rev_tree_dict):

    event_counts = pd.DataFrame(columns=['Branch','Dup_events',
        'Dup_genes','Loss_events', 'Lost_genes', 'Inversions', 'Transpositions'
        ])

    #Key = branches, Value = (D,L,list of inv timestamps, list of transl
    #timestamps)
    branch_dict = {}
    total = [0,0,0,0,0,0]
    for filename in listdir(events_folder):
        ch = filename.split('_')[0]
        if ch != 'Root':
            par = rev_tree_dict[ch]
            branch = (par,ch)
            branch_dict[branch] = [0,0,0,0,0,0]

        dup_events, dups = 0, 0
        loss_events, losses = 0, 0
        invs, transps = 0, 0

        prev_time = 0
        prev_event = None
        with open(join(events_folder, filename)) as data:
            for line in csv.reader(data, delimiter = '\t'):
                if line[0].startswith('#'):
                    continue

                curr_time = line[0]
                curr_event = line[1]
                if curr_time == prev_time:
                    if curr_event == 'D':
                        dups += 1
                    elif curr_event == 'L':
                        losses += 1
                else:
                    if curr_event == 'D':
                        dups += 1
                        dup_events += 1
                    elif curr_event == 'L':
                        losses += 1
                        loss_events += 1
                    elif curr_event == 'I':
                        invs += 1
                    elif curr_event == 'P':
                        transps += 1

                prev_event = curr_event
                prev_time = curr_time

        branch_dict[branch] = [dup_events,dups,loss_events,losses,invs,transps]
        for k in range(len(total)):
            total[k] += branch_dict[branch][k]

    branch_dict['Overall'] = total
    for branch in branch_dict:
        event_counts = event_counts.append({ \
                'Branch':           branch, \
                'Dup_events':       branch_dict[branch][0], \
                'Dup_genes':        branch_dict[branch][1], \
                'Loss_events':      branch_dict[branch][2], \
                'Lost_genes':       branch_dict[branch][3], \
                'Inversions':       branch_dict[branch][4], \
                'Transpositions':   branch_dict[branch][5] \
                },ignore_index=True)

    return event_counts


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('tree_events', type=open,
            help='Events.tsv of ZOMBI experiment')
    parser.add_argument('events_folder', type=str,
            help='events folder of ZOMBI simulation')
    args = parser.parse_args()

    ancestral, extant, rev_tree_dict = readTreeEvents(args.tree_events)
    events_table = readEventCounts(args.events_folder, rev_tree_dict)

    pd.set_option('expand_frame_repr', True)
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_colwidth', 50)
    pd.set_option('display.max_rows', 1000)
    pd.set_option('display.max_columns', 100)
    print(events_table, file = stdout)

