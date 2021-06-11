#!/usr/bin/env python

#
# TODO:
#
#   - sample DCJ cut positions by adjcacencies, not by #chrs+adj
#   - allow duplications in circular chromosomes to go beyond the "last" position
#   - allow deletions in circular chromosomes to go beyond "last" position

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF, \
        FileType
from os.path import basename, dirname, abspath, isfile, join
from sys import stdout, stderr, exit, setrecursionlimit, maxint
from cStringIO import StringIO
from copy import deepcopy, copy
from functools import partial
from operator import add
import logging

import numpy as np
import networkx as nx



from trees import parse_tree, Leaf
from trees import reroot
from trees import rescale_absolute

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

DIRECTION_CRICK_STRAND = '+'
DIRECTION_WATSON_STRAND = '-'
DIRECTION_BOTH_STRAND = '|'

DEFAULT_ZIPF_S = 3.5

setrecursionlimit(10000000)

# construct join functions: calling the function joins the given gene with
# the loose extremity
joinHead = lambda x, end: setattr(end, 'head', x)
joinTail = lambda x, end: setattr(end, 'tail', x)


class Gene:
    def __init__(self, name, head=None, tail=None):
        self.name = name
        self.head = head
        self.tail = tail

    def __str__(self):
        return str(self.name)

    def next(self, incoming=None):
        if incoming == self.head:
            return self.tail
        elif incoming == self.tail:
            return self.head
        else:
            raise LookupError(("%s is neither connected to head (%s) nor tail (%s) of " \
                    + "gene %s") %(incoming, self.head, self.tail, self.name))

    def orient(self, incoming):

        # if self.tail == self.head (in circular chromosomes of size 2, then
        # return DIRECTION_BOTH_STRAND
        if self.tail == self.head:
            return DIRECTION_BOTH_STRAND
        if incoming == self.tail:
            return DIRECTION_CRICK_STRAND
        elif incoming == self.head:
            return DIRECTION_WATSON_STRAND
        else:
            raise LookupError(("%s is neither connected to head (%s) nor tail (%s) of " \
                    + "gene %s") %(incoming, self.head, self.tail, self.name))

    def map_reduce(self, mapper, reducer, init=None, incoming=None, stop_at=None):

        next_g = None

        if incoming == self.head:
            next_g = self.tail
        elif incoming == self.tail:
            next_g = self.head
        else:
            raise LookupError(("%s is neither connected to head (%s) nor tail (%s) of " \
                    + "gene %s") %(incoming, self.head, self.tail, self.name))

        if stop_at == None:
            stop_at = self

        if next_g != None and next_g != stop_at:
            return reducer(mapper(self), next_g.map_reduce(mapper, reducer,
                init=init, incoming=self, stop_at=stop_at))
        elif init:
            return reducer(mapper(self), init)
        else:
            return mapper(self)

class Telomere(Gene):

    def __init__(self, next_gene=None):
        Gene.__init__(self, 'o', head=next_gene)

    def orient(self, incoming=None):
        return Gene.orient(self, incoming)


def print_genome_unimog(genome, genomeName, out):

    out.write('>%s\n' %genomeName)
    for initG in genome:
        chro = initG.map_reduce(lambda x: x, lambda x, y: [(x.orient(y[0][1]), x)] + y,
                init=[(DIRECTION_CRICK_STRAND, initG.tail and initG or None)],
                incoming=initG.tail)
        chro.pop()
        if isinstance(chro[0][1], Telomere):
            chro = chro[1:-1]
        out.write(' '.join(map(lambda x: x[0] == DIRECTION_WATSON_STRAND and str(x[1]) or
            '-%s'%str(x[1]), chro)))
        out.write(' %s'% (initG.tail==None and '|' or ')'))
        if initG != genome[-1]:
            out.write('\n')

def print_genome(genome, out, genomeName=None, printLengths=False, short=False):

    if genomeName:
        out.write(genomeName)
        if short:
            out.write('\t')
        else:
            out.write('\n')

    count = 1
    for initG in genome:
        if genomeName and not short:
            out.write('\t')

        if not short:
            out.write('Chr%s' %count)
            if printLengths:
                length = initG.map_reduce(lambda x: 1, add, incoming=initG.tail)
                out.write(' (len %s)' %length)
            out.write('\t')

        out.write(initG.tail==None and '[ ' or '( ')

        chro = initG.map_reduce(lambda x: x, lambda x, y: [(x.orient(y[0][1]), x)] + y,
                init=[(DIRECTION_CRICK_STRAND, initG.tail and initG or None)],
                incoming=initG.tail)
        chro.pop()

#        if len(filter(lambda x: x[0] == DIRECTION_CRICK_STRAND or x[0] == DIRECTION_BOTH_STRAND, chro)) > len(chro)/2:
#            out.write(' '.join(map(lambda x: (x[0] == DIRECTION_CRICK_STRAND or x[0] == DIRECTION_BOTH_STRAND) and str(x[1]) or
#                '-%s'%str(x[1]), chro)))
#        else:
        out.write(' '.join(map(lambda x: x[0] == DIRECTION_WATSON_STRAND and str(x[1]) or
            '-%s'%str(x[1]), chro)))

        out.write(' %s'% (initG.tail==None and ']' or ')'))

        if short and initG != genome[-1]:
            out.write(', ')
        else:
            out.write('\n')

        count+=1

def segment_to_str(startG, endG, prev):

    res = StringIO()
    cur = startG
    while cur != endG:
        res.write('%s ' %cur.name)
        prev, cur = cur, cur.next(prev)

    if endG != None:
        res.write(str(endG.name))
    return res.getvalue()


def copy_segment(startPos, prev=None, stopAt=None, numberOfGenes=0):

    assert numberOfGenes == 0 or stopAt == None

    if isinstance(startPos, Telomere):
        copyChr = Telomere()
    else:
        copyChr = Gene(startPos.name)

    if stopAt == None:
        stopAt = startPos

    if numberOfGenes == 0:
        numberOfGenes = maxint

    if prev == None:
        prev = startPos.tail

    curG = startPos.next(prev)
    curOrient = curG.orient(startPos)

    copyPrev = copyChr
    if curOrient == DIRECTION_CRICK_STRAND or curOrient == DIRECTION_BOTH_STRAND:
        stickyStart = partial(joinTail, end=copyChr)
        stickyEnd = partial(joinHead, end=copyChr)
    else:
        stickyStart = partial(joinHead, end=copyChr)
        stickyEnd = partial(joinTail, end=copyChr)

    i = 1
    while curG != None and curG != stopAt and i < numberOfGenes:
        # copy gene
        if isinstance(curG, Telomere):
            copyG = Telomere()
        else:
            copyG = Gene(name=curG.name)
        # join ends with previous gene
        stickyEnd(copyG)

        if curOrient == DIRECTION_CRICK_STRAND or curOrient == DIRECTION_BOTH_STRAND:
            # join ends with previous gene
            copyG.tail = copyPrev
            # create sticky ends for next iteration
            stickyEnd = partial(joinHead, end=copyG)
            # set up gene and orientation for next iteration
            curOrient = curG.head != None and curG.head.orient(curG) or None
            curG = curG.head
        else:
            # join ends with previous gene
            copyG.head = copyPrev
            # create sticky ends for next iteration
            stickyEnd = partial(joinTail, end=copyG)
            # set up gene and orientation for next iteration
            curOrient = curG.tail != None and curG.tail.orient(curG) or None
            curG = curG.tail

        copyPrev = copyG
        i += 1

    return stickyStart, stickyEnd, copyChr, copyPrev

def copy_genome(genome):

    res = list()

    for chromosome in genome:
        stickyStart, stickyEnd, chrStart, chrEnd = copy_segment(chromosome)
        if chromosome.tail != None:
            stickyEnd(chrStart)
            stickyStart(chrEnd)
        res.append(chrStart)
    return res

def checkIntegrity(genome):

    # integrity check is flawed.. check last connection between chrx.tail and
    # chrx for circular chromosomes?
    for i, chrx in enumerate(genome):
        prev = chrx.tail
        cur = chrx
        cur_orient = DIRECTION_CRICK_STRAND
        prev_orient = prev == None and DIRECTION_CRICK_STRAND
        if prev != None:
            prev_orient = prev.head == cur and DIRECTION_CRICK_STRAND or \
                    DIRECTION_WATSON_STRAND
        c = 0
        while cur != chrx.tail:
            if prev != None:
                if cur_orient == DIRECTION_CRICK_STRAND:
                    if prev_orient == DIRECTION_CRICK_STRAND and \
                        (prev.head != cur or cur.tail != prev):
                        LOG.error(('integrity test of chromosome %s failed ' + \
                                'at index %s: %s.head != %s or %s.tail ' + \
                                '!= %s') %(i+1, c, prev, cur, cur, prev))
                        import pdb; pdb.set_trace()
                        return False
                    if prev_orient == DIRECTION_WATSON_STRAND and \
                        (prev.tail != cur or cur.tail != prev):
                        LOG.error(('integrity test of chromosome %s failed ' + \
                                'at index %s: %s.tail != %s or %s.tail ' + \
                                '!= %s') %(i+1, c, prev, cur, cur, prev))
                        import pdb; pdb.set_trace()
                        return False
                    prev, prev_orient, cur = cur, cur_orient, cur.head
                    if cur != None:
                        cur_orient = cur.tail == prev and \
                                DIRECTION_CRICK_STRAND or \
                                DIRECTION_WATSON_STRAND
                elif cur_orient == DIRECTION_WATSON_STRAND:
                    if prev_orient == DIRECTION_CRICK_STRAND and \
                        (prev.head != cur or cur.head != prev):
                        LOG.error(('integrity test of chromosome %s failed ' + \
                                'at index %s: %s.head != %s or %s.head' + \
                                '!= %s') %(i+1, c, prev, cur, cur, prev))
                        import pdb; pdb.set_trace()
                        return False

                    if prev_orient == DIRECTION_WATSON_STRAND and \
                        (prev.tail != cur or cur.head != prev):
                        LOG.error(('integrity test of chromosome %s failed ' + \
                                'at index %s: %s.tail != %s or %s.head' + \
                                '!= %s') %(i+1, c, prev, cur, cur, prev))
                        import pdb; pdb.set_trace()
                        return False
                    prev, prev_orient, cur = cur, cur_orient, cur.tail,
                    if cur != None:
                        cur_orient = cur.tail == prev and \
                                DIRECTION_CRICK_STRAND or \
                                DIRECTION_WATSON_STRAND
            else:
                prev, prev_orient, cur, cur_orient = cur, cur_orient, \
                        cur.head, cur.head.tail == cur and \
                        DIRECTION_CRICK_STRAND or DIRECTION_WATSON_STRAND
            c += 1
    return True

def generateRootGenome(no_genes, no_chrs, do_circular=False):

    genome = list()

    seq = range(1, no_genes+1)

    cuts = set()
    cuts.add(0)

    while len(cuts) < no_chrs:
        cuts.add(np.random.randint(1, len(seq)-1))

    prev = None
    for i in xrange(len(seq)):
        if i in cuts:
            if prev != None:
                if do_circular:
                    genome[-1].tail = prev
                    prev.head = genome[-1]
                else:
                    t = Telomere(prev)
                    prev.head = t

            g = Gene(seq[i])
            if not do_circular:
                t = Telomere(g)
                g.tail= t
                genome.append(t)
            else:
                genome.append(g)
        else:
            g = Gene(seq[i], tail=prev)
            prev.head = g

        prev = g

    if prev != None:
        if do_circular:
            genome[-1].tail = prev
            prev.head = genome[-1]
        else:
            t = Telomere(prev)
            prev.head = t

    return genome


def goToPos(chromosome, pos, prev=None):
    curG = chromosome
    if prev == None:
        prev = chromosome.tail
    while pos > 0:
        pos -= 1
        tmp = curG
        curG = curG.next(prev)
        prev = tmp
    return curG, prev


def cut(chromosome, pos):

    curG, prev = goToPos(chromosome, pos-1)

    nextG = curG.next(prev)
    if curG.orient(prev) == DIRECTION_CRICK_STRAND or \
            curG.orient(prev) == DIRECTION_BOTH_STRAND:
        stickyEnd1 = partial(joinHead, end=curG)
    else:
        stickyEnd1 = partial(joinTail, end=curG)

    if nextG.orient(curG) == DIRECTION_CRICK_STRAND or \
            nextG.orient(curG) == DIRECTION_BOTH_STRAND:
        stickyEnd2 = partial(joinTail, end=nextG)
    else:
        stickyEnd2 = partial(joinHead, end=nextG)

    return stickyEnd1, stickyEnd2, curG, nextG


def sampleSites(genome, chr_lengths, number_of_sites=1, min_dist=0):
    """ sampling function that works for finding (i) cut points for DCJ
    operations, (ii) donor and receiver sites for duplication events and (iii)
    indel regions """

    assert number_of_sites > 0 and number_of_sites <= 2 and min_dist >= 0

    # for first site, only select sites from chromosomes whose size is larger
    # than min_dist
    feasible_chrs = [x for x in xrange(len(genome)) if chr_lengths[x]
        - (genome[x].tail == None and 2 or 0) >= min_dist]
    assert len(feasible_chrs) > 0

    sites = list()
    while len(sites) < number_of_sites:
        # choose chromosome
        if not sites:
            c = np.random.choice(feasible_chrs)
        else:
            c = np.random.randint(0, len(genome))
        # choose site, do exlude site before first and after last telomere, but
        # if chromosome is circular, allow site betweeen first and last gene
        end_offset = (genome[c].tail == None and 1 or 0)
        #
        # also, impose min_dist constraint only for first site (required for
        # donor region of duplication event or deletion region); that is, if
        # chromosome is linear keep do not allow sites within the end region
        # [-1+min_dist:] of the chromosome
        if not sites:
            end_offset = (genome[c].tail == None and 1+min_dist or  max(0,
                min_dist-1))
        site = (c, np.random.randint(1, chr_lengths[c] + 1 - end_offset))
        # don't allow second site to be equal to first or fall within min_dist
        # region
        if not sites or sites[-1][0] != c or site[1] < sites[0][1] or \
                site[1] >= sites[0][1]+max(1, min_dist):
            sites.append(site)
    return sites


def simulateDCJ(root_genome, tree, gfCount=None, duplication_rate=0,
        duplication_zipf=DEFAULT_ZIPF_S, insertion_rate=0, deletion_rate=0,
        indel_zipf=DEFAULT_ZIPF_S):

    tree = deepcopy(tree)

    stack = list()
    stack.extend(tree.subtree.subtrees)
    setattr(tree.subtree, 'genome', root_genome)

    tree.subtree.label = 'root'

    no_nodes = len(tree.getNodes())
    node_count = 1
    cur_c = 0

    if gfCount == None:
        gfCount = sum(i.map_reduce(lambda x: 1, add, incoming=i.tail) for i in
                root_genome)

    while stack:
        x = stack.pop()

        if type(x) is not Leaf:
            stack.extend(x.subtrees)

        genome = copy_genome(x.getAncestor().genome)

        chr_lengths = [i.map_reduce(lambda x: 1, add, incoming=i.tail) for i in genome]

        LOG.info('mutate genome over evolutionary distance %s PAM: ' %x.length)

        evo_dist = x.length
        while np.random.random() < evo_dist:
            evo_dist -= 1

            # determine whether to perform a duplication
            dr = duplication_rate * min(1, evo_dist)
            while np.random.random() < dr:
                dr -= 1
                segment_size = np.random.zipf(duplication_zipf)
                # randomly sample donor and receiver site
                while segment_size:
                    try:
                        donor, receiver = tuple(sampleSites(genome, chr_lengths,
                            2, segment_size))
                        break
                    except AssertionError, e:
                        segment_size -= 1
                if not segment_size:
                    LOG.fatal('Unable to perform duplication: %s' %e.message)
                    exit(1)

                curG, prev = goToPos(genome[donor[0]], donor[1])
                stickyEnd1, stickyEnd2, dup_end1, dup_end2 = copy_segment(curG,
                        prev=prev, numberOfGenes=segment_size)
                stickyEnd3, stickyEnd4, end3, end4 = cut(genome[receiver[0]],
                        receiver[1])

                LOG.info(('++ duplication of donor region Chr%s:[%s-%s] = ' + \
                        '%s to receiving site Chr%s:%s') %(donor[0]+1, \
                        donor[1], donor[1]+segment_size-1,\
                        segment_to_str(dup_end1, dup_end2, None), \
                        receiver[0]+1, receiver[1]))
                # randomly choose orientation of inserted segment
                if np.random.randint(0, 2):
                    stickyEnd1(end3)
                    stickyEnd3(dup_end1)
                    stickyEnd2(end4)
                    stickyEnd4(dup_end2)
                else:
                    stickyEnd1(end4)
                    stickyEnd4(dup_end1)
                    stickyEnd2(end3)
                    stickyEnd3(dup_end2)
                chr_lengths[receiver[0]] += segment_size

                # Chromosome handle does not need to be updated, because genes
                # are only added, not removed. For linear chromosomes, insertion
                # before telomere is prohibited, so the handle on the telomere
                # is unaffected
                #checkIntegrity(genome)

            # determine whether to perform a deletion
            ir = deletion_rate * min(1, evo_dist)
            while np.random.random() < ir:
                ir -= 1
                segment_size = np.random.zipf(indel_zipf)
                # randomly sample site
                while segment_size:
                    try:
                        (site, ) = tuple(sampleSites(genome, chr_lengths, \
                                min_dist=segment_size))
                        break
                    except AssertionError, e:
                        LOG.warning(('Sampling deletion site with segment ' + \
                                'size %s was unsuccessful, but will try ' + \
                                'again with smaller segment size: %s') %(
                                    segment_size, e.message))
                        segment_size -= 1
                if not segment_size:
                    LOG.fatal('Unable to perform deletion: %s' %e.message)
                    exit(1)
                startG, startPrev = goToPos(genome[site[0]], site[1])
                endG, endPrev = goToPos(startG, segment_size, prev=startPrev)

                LOG.info(('++ deletion of region Chr%s:[%s-%s] = %s') %(
                    site[0]+1, site[1], site[1] + segment_size-1,
                    segment_to_str(startG, endPrev, startPrev)))

                # We aim to orient startPrev using startG. Gene startPrev is the
                # last gene prior to the deletion. In doing so, we come from the
                # WATSON STRAND (-)
                #             <--- incoming --
                # startPrev <cut> startG
                if startPrev.orient(startG) == DIRECTION_WATSON_STRAND or \
                        startPrev.orient(startG) == DIRECTION_BOTH_STRAND:
                    stickyEnd1 = partial(joinHead, end=startPrev)
                else:
                    stickyEnd1 = partial(joinTail, end=startPrev)

                # For the end, everything is as usual:
                # -incoming -->
                #    endPrev <cut> endG
                if endG.orient(endPrev) == DIRECTION_CRICK_STRAND or \
                        endG.orient(endPrev) == DIRECTION_BOTH_STRAND:
                    stickyEnd2 = partial(joinTail, end=endG)
                else:
                    stickyEnd2 = partial(joinHead, end=endG)

                stickyEnd1(endG)
                stickyEnd2(startPrev)

                # update handle for circular genomes
                if genome[site[0]].tail and genome[site[0]] == endPrev:
                    LOG.debug(('handle of chromosome %s was deleted, updating' + \
                        ' to new handle %s') %(site[0]+1, endG))
                    genome[site[0]] = endG

                chr_lengths[site[0]] -= segment_size
                # remove circular chromosome if all has been deleted, but keep
                # telomeres of empty linear chromosomes
                if not chr_lengths[site[0]]:
#                    or (genome[site[0]].tail == None \
#                        and chr_lengths[site[0]] == 2):
                    del chr_lengths[site[0]]
                    del genome[site[0]]

                if not sum(chr_lengths):
                    LOG.fatal('One of the simulated genomes is empty! ' + \
                            'Adjust evolution parameters!  Exiting')
                    exit(1)

                #checkIntegrity(genome)

            # determine whether to perform an insertion
            ir = insertion_rate * min(1, evo_dist)
            while np.random.random() < ir:
                ir -= 1
                segment_size = np.random.zipf(indel_zipf)
                # randomly sample site
                (site, ) = tuple(sampleSites(genome, chr_lengths))

                LOG.info('++ insertion of sequence (%s) at site Chr%s:%s' %( \
                        ' '.join(map(str, xrange(gfCount + 1, gfCount + 1 + \
                        segment_size))), site[0]+1, site[1]))
                gfCount += 1
                startG = Gene(gfCount)
                prev = startG
                i = segment_size -1
                while i > 0:
                    gfCount += 1
                    g = Gene(gfCount, tail=prev)
                    prev.head = g
                    prev = g
                    i -= 1
                stickyEnd1 = partial(joinTail, end=startG)
                stickyEnd2 = partial(joinHead, end=prev)
                end3, end4 = goToPos(genome[site[0]], site[1])

                # same situation as above for the deletion case
                if end3.orient(end4) == DIRECTION_WATSON_STRAND or \
                        end3.orient(end4) == DIRECTION_BOTH_STRAND:
                    stickyEnd3 = partial(joinHead, end=end3)
                else:
                    stickyEnd3 = partial(joinTail, end=end3)
                if end4.orient(end3) == DIRECTION_CRICK_STRAND or \
                        end4.orient(end3) == DIRECTION_BOTH_STRAND:
                    stickyEnd4 = partial(joinTail, end=end4)
                else:
                    stickyEnd4 = partial(joinHead, end=end4)

                if np.random.randint(0, 2):
                    stickyEnd1(end3)
                    stickyEnd3(startG)
                    stickyEnd2(end4)
                    stickyEnd4(prev)
                else:
                    stickyEnd1(end4)
                    stickyEnd4(startG)
                    stickyEnd2(end3)
                    stickyEnd3(prev)

                chr_lengths[site[0]] += segment_size
                # Chromosome handle does not need to be updated, because genes
                # are only added, not removed. For linear chromosomes, insertion
                # before telomere is prohibited, so the handle on the telomere
                # is unaffected
                #checkIntegrity(genome)

            # randomly sample two cuts
            # sorting cut points is important for the case where both cuts
            # reside on the same chromosome
            (chr1, p1), (chr2, p2) = sorted(sampleSites(genome, chr_lengths, 2))

            info_msg = '++ double cut (Chr%s[%s|%s], Chr%s[%s|%s])' %(chr1+1, p1,
                    p1+1, chr2+1, p2, p2+1)
            stickyEnd1, stickyEnd2, end1, end2 = cut(genome[chr1], p1)
            stickyEnd3, stickyEnd4, end3, end4 = cut(genome[chr2], p2)

            # choose join
            join_combination = np.random.randint(0, 2)
            # ---1|2---  ---3|4---
            if join_combination:
                # ---1|4---  ---2|3---
                stickyEnd1(end4)
                stickyEnd4(end1)
                stickyEnd2(end3)
                stickyEnd3(end2)

                info_msg += (' and join (Chr%s:%s|Chr%s:%s, ' \
                        + 'Chr%s:%s|Chr%s:%s)') % (chr1+1, p1, chr2+1, p2+1,
                                chr1+1, p1+1, chr2+1, p2)
            else:
                # ---1|3---  ---2|4---
                stickyEnd1(end3)
                stickyEnd3(end1)
                stickyEnd2(end4)
                stickyEnd4(end2)

                info_msg += (' and join (Chr%s:%s|Chr%s:%s, ' \
                        + 'Chr%s:%s|Chr%s:%s)') % (chr1+1, p1, chr2+1, p2,
                                chr1+1, p1+1, chr2+1, p2+1)

            # update chromosome number and lengths
            if chr1 == chr2 and join_combination == 1:

                # CHROMOSOME FISSION EVENT


                # circular chromosome is created, hence add it to genome
                # list
                genome.append(end2)

                # update genome handle for circular genomes
                if genome[chr1].tail != None:
                    genome[chr1] = end1

                new_len = abs(p2-p1)
                chr_lengths.append(new_len)
                chr_lengths[chr1] -= new_len

                info_msg += (' results in chromosome fission: create new ' \
                        + 'circular chromosome Chr%s of length %s, old ' \
                        + 'chromosome Chr%s has new length %s') %(len(genome),
                                new_len, chr1+1, chr_lengths[chr1])
            elif chr1 != chr2:
                if genome[chr1].tail == None and genome[chr2].tail == None:
                    # TRANSLOCATION EVENT
                    if join_combination == 1:
                        tmp1 = chr_lengths[chr1]
                        chr_lengths[chr1] = p1 + chr_lengths[chr2] - p2
                        chr_lengths[chr2] = p2 + tmp1 - p1
                    else:
                        # assign new telomere end to genome list at position
                        # genome[chr2]
                        curG = end2.head
                        prev = end2
                        while curG != None and curG != end2:
                            tmp = curG
                            curG = curG.next(prev)
                            prev = tmp

                        genome[chr2] = prev

                        tmp1 = chr_lengths[chr1]
                        chr_lengths[chr1] = p1 + p2
                        chr_lengths[chr2] = tmp1 - p1 + chr_lengths[chr2] - p2

                    info_msg += (' results in chromosome translocation: new' \
                            + ' lengths of chromosomes are |Chr%s|=%s and ' \
                            + '|Chr%s|=%s') %(chr1+1, chr_lengths[chr1],
                                    chr2+1, chr_lengths[chr2])

                else:
                    # CHROMSOME FUSION EVENT
                    if genome[chr1].tail == None:
                        chr_lengths[chr1] += chr_lengths[chr2]
                        info_msg += (' results in chromosome fusion: fused' \
                                + ' chromosome Chr%s has length %s') %(chr1+1,
                                        chr_lengths[chr1])
                        del genome[chr2]
                        del chr_lengths[chr2]
                    else:
                        chr_lengths[chr2] += chr_lengths[chr1]
                        info_msg += (' results in chromosome fusion: fused' \
                                + ' chromosome Chr%s has length %s') %(chr2+1,
                                        chr_lengths[chr2])
                        del genome[chr1]
                        del chr_lengths[chr1]

            LOG.info(info_msg)
            #checkIntegrity(genome)

        if not x.label or x.label == '@':
            x.label = 'node_%s' %(node_count)

        tmp = int(round(node_count * 10./no_nodes))
        node_count += 1

        if tmp > cur_c:
            LOG.info('processed %s%% of evolutionary tree' %(tmp * 10))
            cur_c = tmp
        setattr(x, 'genome', genome)

    return tree, gfCount

def write_genomes(tree, out):
    stack = list()
    stack.append(tree.subtree)

    while stack:
        x = stack.pop()
        print_genome([chx for chx in x.genome if chx.tail != None or
            chx.head.tail != None], out, genomeName=x.label, short=True)
        out.write('\n')
        if type(x) is not Leaf:
            stack.extend(x.subtrees)

def write_genomes_unimog(tree, out, leaves_only = False):
    stack = list()
    stack.append(tree.subtree)

    while stack:
        x = stack.pop()
        if not leaves_only or type(x) is Leaf:
            print_genome_unimog([chx for chx in x.genome if chx.tail != None or
                chx.head.tail != None], x.label, out)
            out.write('\n')
        if type(x) is not Leaf:
            stack.extend(x.subtrees)

def toAdjcencyGraph(genome1, genome2):

    G = nx.MultiGraph()

    adjs1 = dict()
    for i, genome in enumerate((genome1, genome2)):
        for chrx in genome:

            if chrx.tail == None and chrx.head == None:
                continue
            out = StringIO()
            print_genome([chrx], out, short=True)
            seq = out.getvalue().strip()

            prev = None
            circFirst = None
            for si in seq.split(' '):

                if si == '[' or si== ']':
                    continue
                if si == 'o' or si == '-o':
                    if prev != None:
                        adj = (i, (prev, ))
                        G.add_node(adj)
                        if i == 0:
                            adjs1[prev] =  adj
                        else:
                            G.add_edge(adj, adjs1[prev])
                            prev = None
                elif si == '(':
                    circFirst = True
                elif si == ')':
                    adj = (i, (prev, circFirst))
                    G.add_node(adj)
                    if i == 0:
                        adjs1[prev] = adj
                        adjs1[circFirst] = adj
                    else:
                        G.add_edge(adj, adjs1[prev])
                        G.add_edge(adjs1[circFirst], adj)
                    prev = None
                elif prev != None:
                    gid = int(si)
                    if gid > 0:
                        adj = (i, (prev, (gid, 't')))
                        G.add_node(adj)
                        if i == 0:
                            adjs1[prev] = adj
                            adjs1[(gid, 't')] = adj
                        else:
                            G.add_edge(adj, adjs1[(gid, 't')])
                            G.add_edge(adj, adjs1[prev])
                        prev = (gid, 'h')

                    else:
                        adj = (i, (prev, (abs(gid), 'h')))
                        G.add_node(adj)
                        if i == 0:
                            adjs1[prev] = adj
                            adjs1[(abs(gid), 'h')] = adj
                        else:
                            G.add_edge(adj, adjs1[(abs(gid), 'h')])
                            G.add_edge(adj, adjs1[prev])
                        prev = (abs(gid), 't')
                elif circFirst:
                    gid = int(si)
                    if gid > 0:
                        circFirst = (gid, 't')
                        prev = (gid, 'h')
                    else:
                        circFirst = (abs(gid), 'h')
                        prev = (abs(gid), 't')
                else:
                    gid = int(si)
                    if gid > 0:
                        adj = (i, ((gid, 't'), ))
                        G.add_node(adj)
                        if i == 0:
                            adjs1[(gid, 't')] = adj
                        else:
                            G.add_edge(adj, adjs1[(gid, 't')])
                        prev = (gid, 'h')
                    else:
                        adj = (i, ((abs(gid), 'h'), ))
                        G.add_node(adj)
                        if i == 0:
                            adjs1[(abs(gid), 'h')] = adj
                        else:
                            G.add_edge(adj, adjs1[(abs(gid), 'h')])

                        prev = (abs(gid), 't')
    return G

def calcDCJ(G):
    c = 0
    i = 0
    for C in nx.connected_components(G):
        degrees = tuple(sorted(set(G.degree(C).values())))
        l = len(G.edges(C))
        if degrees == (2,):
            c += 1
        elif degrees == (1,2) or degrees == (1, ):
            if l %2:
                i += 1
        else:
            import pdb; pdb.set_trace()
            raise Exception, 'There is an error in the program, the adjacency ' + \
                    'graph has components that are neither simple cycles ' + \
                    'nor simple paths'
    return c + i/2

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-g', '--number_of_genes', type=int, default=1000,
            help='number of genes in simulated genomes')
    parser.add_argument('-x', '--number_of_chromosomes', type=int, default=1,
            help='number of chromosomes in root genome')
    parser.add_argument('-c', '--circular_only', action = 'store_true',
            help='only simulate circular chromosomes')
    parser.add_argument('-r', '--reroot', action='store_true', default=False,
            help='reroot input tree to node labeled "@"')
    parser.add_argument('-s', '--rescale', type=int, default=0,
            help='scale tree to the given, absolute number of evolutionary' \
                    + ' events along longest path in tree')
    parser.add_argument('-d', '--duplication_rate', type=float, default=0,
            help='rate of duplication, measured relative to DCJ events')
    parser.add_argument('--duplication_size_zipf', type=float,
            default=DEFAULT_ZIPF_S,
            help='size of duplicated segments, controlled by zipfian ' + \
                    'distribution parameter')
    parser.add_argument('-i', '--insertion_rate', type=float, default=0,
            help='rate of indels, measured relative to DCJ events')
    parser.add_argument('-e', '--deletion_rate', type=float, default=0,
            help='rate of deletions, measured relative to DCJ events')
    parser.add_argument('--indel_size_zipf', type=float, default=DEFAULT_ZIPF_S,
            help='size of indels, controlled by zipfian distribution parameter')
    parser.add_argument('-l', '--leaves_only', action = 'store_true',
            help='only output genomes corresponding to the leaves of the ' + \
                    'phylogeny')
    parser.add_argument('-o', '--out_file', type=FileType('w'),
            help='base filename of output files (genomes and tree)' + \
                    '[default=<basename of input file>]')
    parser.add_argument('newick_tree', type=file, default=0,
            help='rate of indels, measured relative to DCJ events')
    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    #
    # read and manipulate tree
    #

    LOG.info('Parameter of evolution run:\n\t\t%s' %'\n\t\t'.join((
        'path to file containing evolutionary tree: %s' %args.newick_tree.name,
        'number of genes in root genome: %s' %args.number_of_genes,
        'number of chromosomes in root genome: %s' %args.number_of_chromosomes,
        'scale tree to given, absolute number of evolutionary events: ' + \
                '%s' %args.rescale and 'Yes, %s' %args.rescale or 'No',
        'duplication rate: %s' %args.duplication_rate,
        'zipfian distribution shape parameter for duplication size: ' + \
                '%s' %args.duplication_size_zipf,
        'insertion rate: %s' %args.insertion_rate,
        'deletion rate: %s' %args.deletion_rate,
        'zipfian distribution shape parameter for indel size: ' + \
                '%s' %args.indel_size_zipf)))

    tree = parse_tree(args.newick_tree)

    if args.rescale > 0:
        LOG.info(('rescaling tree to %s evolutionary events on longest ' \
                + 'path') %args.rescale)
        tree = rescale_absolute(tree, args.rescale)
        LOG.info('rescaled tree: %s' %tree)
    if args.reroot:
        reroot_br = [x for x in tree.getNodes() if x.label == '@']
        if reroot_br:
            tree = reroot(reroot_br[0])
            LOG.info('rerooting tree, new tree: %s' %tree)
    #
    # generate root genome
    #

    LOG.info(('generate root genome with %s genes partitioned into ' \
            '%s chromosomes') %(args.number_of_genes,
                args.number_of_chromosomes))
    root = generateRootGenome(args.number_of_genes, args.number_of_chromosomes,
            args.circular_only)

    s = StringIO()
    print_genome(root, s, 'root', printLengths=True)
    LOG.info('\n%s\n\n' %s.getvalue())

    #
    # simulate DCJ evolution along the given phylogeny
    #

    LOG.info('start simulating DCJ evolution along phylogeny')
    dcj_tree, _ = simulateDCJ(root, tree, gfCount=args.number_of_genes,
            duplication_rate=args.duplication_rate,
            duplication_zipf=args.duplication_size_zipf,
            insertion_rate=args.insertion_rate,
            deletion_rate=args.deletion_rate, indel_zipf=args.indel_size_zipf)

#    if args.duplication_rate == 0 and args.insertion_rate == 0 and \
#            args.deletion_rate == 0:
#        genomes = [l.genome for l in dcj_tree.getLeaves()]
#        G = toAdjcencyGraph(*genomes)
#        c = calcDCJ(G)
#        if args.rescale < args.number_of_genes-c:
#            LOG.fatal(('DCJ distance between generated genomes is larger than' + \
#                    'the requested one: given: %s, generated: %s') %(args.rescale,
#                    args.number_of_genes-c))
#            exit(1)
#        elif args.rescale > args.number_of_genes-c:
#            LOG.warning(('DCJ distance between generated genomes is smaller than' + \
#                    'the requested one: given: %s, generated: %s') %(args.rescale,
#                    args.number_of_genes-c))
#            exit(1)


    if args.out_file == None:
        write_genomes_unimog(dcj_tree, stdout, leaves_only=args.leaves_only)
        LOG.info('Genomes printed to stdout. Done')
    else:
        if args.out_file == None:
            treeOut = '%s.dcjsim' %args.newick_tree.name
        else:
            treeOut = '%s.tree' %args.out_file.name
        treeOutHandle = open(treeOut, 'w')
        print >> treeOutHandle, dcj_tree
        treeOutHandle.close()
        write_genomes_unimog(dcj_tree, args.out_file,
                leaves_only=args.leaves_only)
        LOG.info('Output written to %s and %s. Done' %(args.out_file.name,
            treeOut))
