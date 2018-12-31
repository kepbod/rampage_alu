#!/usr/bin/env python

'''
Usage: entropy.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
'''


import os.path
import math
import pysam
from collections import defaultdict
from seqlib.path import check_dir
from seqlib.ngs import check_bed
from joblib import Parallel, delayed

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '1.0.0'


def entropy(options):
    '''
    Calculate entropy for each cluster
    '''
    # parse options
    folder = check_dir(options['<rampagedir>'])
    link_f = check_bed(os.path.join(folder, 'rampage_link.bed'),
                       return_handle=False)
    threads = int(options['--thread'])
    with open(os.path.join(folder, 'rampage_peaks.txt'), 'r') as peak:
        result = Parallel(n_jobs=threads)(delayed(cal_entropy)(line, link_f)
                                          for line in peak)
    with open(os.path.join(folder, 'rampage_entropy.txt'), 'w') as out:
        for r in result:
            out.write(r)


def cal_entropy(line, link_f):
    link = pysam.TabixFile(link_f)
    chrom, start, end, _, _, strand = line.split()[:6]
    start = int(start)
    end = int(end)
    link_lst = defaultdict(int)
    total_counts = 0
    # fetch links
    for l in link.fetch(chrom, start, end):
        s, e, _, reads, ss = l.split()[1:6]
        if ss != strand:  # not identical strand
            continue
        if ss == '+':
            read1_pos = int(s)
            read2_pos = e
        else:
            read1_pos = int(e)
            read2_pos = s
        if start <= read1_pos <= end:  # within peak region
            link_lst[read2_pos] += int(reads)
            total_counts += int(reads)
    # calculate entropy
    entropy = 0
    for pos in link_lst:
        p = link_lst[pos] * 1.0 / total_counts
        entropy += p * math.log(p, 2)
    entropy = -entropy
    return line.rstrip() + '\t%f\t' % entropy + '|'.join(link_lst) + '\n'


if __name__ == '__main__':
    from docopt import docopt
    entropy(docopt(__doc__, version=__version__))
