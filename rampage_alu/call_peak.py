#!/usr/bin/env python

'''
Usage: call_peak.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    --wig                          Create Wig files.
    -p PERCENT                     Retained percent of reads in resized peaks.
                                   [default: 0.95]
'''

import sys
import os.path
import tempfile
import glob
import operator
from collections import Counter
from seqlib.path import which, check_dir
from seqlib.helper import run_command
from seqlib.ngs import check_bed

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def fseq(options):
    '''
    Call peaks using F-seq
    '''
    # parse options
    if not which('fseq'):
        sys.exit('Error: No F-seq installed!')
    folder = check_dir(options['<rampagedir>'])
    flength = options['-l']
    wig_flag = options['--wig']
    percent = float(options['-p'])
    with open(os.path.join(folder, 'total_counts.txt'), 'r') as f:
        total = int(f.read().rstrip())
    # run F-seq
    flist = {'+': 'rampage_plus_5end.bed', '-': 'rampage_minus_5end.bed'}
    all_peak_f = os.path.join(folder, 'rampage_peaks.txt')
    with open(all_peak_f, 'w') as out:
        for strand in flist:
            peak_f = run_fseq(folder, flist[strand], strand, flength, wig_flag,
                              percent)
            with open(peak_f, 'r') as f:
                for line in f:
                    if total:  # calculate RPM
                        reads = int(line.rstrip().split()[9])
                        rpm = reads * 1000000.0 / total
                        out.write(line.rstrip() + '\t%f\n' % rpm)
                    else:
                        out.write(line)


def run_fseq(folder, bed, strand, flength, wig_flag, percent):
    prefix = os.path.splitext(bed)[0]
    # create bed files
    temp_dir = tempfile.mkdtemp()
    bed_f = os.path.join(folder, bed)
    # run fseq
    command = 'fseq -f 0 -l %s -of bed -o %s %s' % (flength, temp_dir, bed_f)
    run_command(command, 'Error in F-seq!')
    # cat fseq files
    peak_f = os.path.join(folder, prefix + '_fseq.bed')
    cat_files(temp_dir, peak_f)
    # resize peaks
    resized_peak_f = os.path.join(folder, prefix + '_peak.bed')
    resize_peak(peak_f, bed_f, resized_peak_f, strand, percent)
    # create wig files
    if wig_flag:
        # run fseq
        temp_dir = tempfile.mkdtemp()
        command = 'fseq -f 0 -l %s -o %s %s' % (flength, temp_dir, bed_f)
        run_command(command, 'Error in F-seq!')
        # cat fseq wig files
        wig = prefix + '_fseq.wig'
        wig_f = os.path.join(folder, wig)
        cat_files(temp_dir, wig_f, is_wig=True, name=wig)
    return resized_peak_f


def cat_files(temp_dir, outf, is_wig=False, name=None):
    with open(outf, 'w') as out:
        if is_wig:  # for wig files, add header
            out.write('track type=wiggle_0 name=%s description=%s\n' % (name,
                                                                        name))
        for fname in glob.iglob(os.path.join(temp_dir, 'chr*')):
            with open(fname) as f:
                out.write(f.read())


def resize_peak(peak, bed, resized_peak, strand, percent):
    bed_f = check_bed(bed)
    with open(peak, 'r') as f, open(resized_peak, 'w') as out:
        for line in f:
            chrom, start, end = line.split()[:3]
            start = int(start)
            end = int(end) + 1
            # count total tags
            total = 0
            sites = []
            for read in bed_f.fetch(chrom, start, end):
                sites.append(int(read.split()[1]))
                total += 1
            if total == 0:
                continue
            # fetch peak tag
            sites = Counter(sites)
            loc, height = sites.most_common()[0]
            # count peak region tags
            peak_sites = 0
            for i in range(loc - 2, loc + 2):
                peak_sites += sites.get(i, 0)
            # resize cluster region
            required_sites = int(total * percent)
            sub_sites = 0
            region = []
            for site, num in sorted(sites.items(), key=operator.itemgetter(1),
                                    reverse=True):
                sub_sites += num
                region.append(site)
                if sub_sites >= required_sites:
                    break
            region.sort()
            new_start, new_end = region[0], region[-1] + 1
            # output result
            out_format = '%s\t%d\t%d\tpeak\t0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n'
            out.write(out_format % (chrom, new_start, new_end, strand, loc,
                                    height, peak_sites, total, start, end))


if __name__ == '__main__':
    from docopt import docopt
    fseq(docopt(__doc__, version=__version__))
