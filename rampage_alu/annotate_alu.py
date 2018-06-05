#!/usr/bin/env python

'''
Usage: annotate_alu.py [options] -f ref (-a alu | -r rep) <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -f ref                         Gene annotations.
    -t type                        File type of gene annotations.
                                   [default: ref]
    --promoter region              Promoter region. [default: 250]
    -a alu                         Alu annotations.
    -r rep                         Repeatmasker annotations.
    --extend length                Alu extended length. [default: 50]
    --span span                    Span cutoff. [default: 1000]
    -o out                         Output file. [default: alu_peak.txt]
'''

import sys
import os.path
from seqlib.parse import Annotation
from pybedtools import BedTool
import tempfile
import numpy as np
from collections import defaultdict

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def annotate(options):
    '''
    Annotate expressed Alu elements
    '''
    # check options
    if options['-t'] not in ('ref', 'bed'):
        sys.exit('Error: -t should be "ref" or "bed"!')
    p_r = int(options['--promoter'])
    extend = int(options['--extend'])
    span = int(options['--span'])
    output = options['-o']
    # prepare tmp files
    temp_dir = tempfile.mkdtemp()
    promoter = temp_dir + '/promoter.bed'
    alu = temp_dir + '/alu.bed'
    # parse gene annotations to fetch promoter
    ann = Annotation(options['-f'], ftype=options['-t'])
    with open(promoter, 'w') as out:
        for info in ann:
            if info.strand == '+':
                tss = info.tx_start
            else:
                tss = info.tx_end
            if tss - p_r < 0:
                continue
            out.write('%s\t%d\t%d\tpromoter\t0\t%s\n' % (info.chrom,
                                                         tss - p_r,
                                                         tss + p_r,
                                                         info.strand))
    # parse alu annotations
    with open(alu, 'w') as out:
        for chrom, start, end, name, strand in parse_alu(options):
            if strand == '+':
                start -= extend
                if start < 0:
                    continue
            else:
                end += extend
            out.write('%s\t%d\t%d\t%s\t0\t%s\n' % (chrom, start, end, name,
                                                   strand))
    # fetch Alu peaks
    rampage_bed = BedTool(os.path.join(options['<rampagedir>'],
                                       'rampage_entropy.txt'))
    promoter_bed = BedTool(promoter)
    alu_bed = BedTool(alu)
    alu_peak = rampage_bed.intersect(promoter_bed, v=True,
                                     s=True).intersect(alu_bed, s=True,
                                                       wa=True, wb=True)
    # filter alu peaks
    peak_lst = defaultdict(list)
    for p in alu_peak:
        if float(p[13]) < 2.5:
            continue
        if not check_span(p, span):
            continue
        peak = int(p[6])
        alu_chr = p[15]
        alu_start = int(p[16])
        alu_end = int(p[17])
        alu_name = p[18]
        alu_strand = p[20]
        if alu_strand == '+':
            alu_start += extend
            coverage = alu_end - peak
        else:
            alu_end -= extend
            coverage = peak - alu_start
        alu_len = alu_end - alu_start
        if coverage < alu_len * 0.5:
            continue
        alu_info = '%s\t%d\t%d\t%s\t0\t%s' % (alu_chr, alu_start, alu_end,
                                              alu_name, alu_strand)
        peak_info = '\t'.join(p[:15])
        peak_lst[peak_info].append(alu_info)
    alu_exp = {}
    alu_info = {}
    for peak_info in peak_lst:
        expression = float(peak_info.split()[13])
        if len(peak_lst[peak_info]) == 1:
            alu = peak_lst[peak_info][0]
            if alu in alu_exp:
                if expression > alu_exp[alu]:
                    alu_exp[alu] = expression
                    alu_info[alu] = peak_info
            else:
                alu_exp[alu] = expression
                alu_info[alu] = peak_info
        else:
            strand, peak = peak_info.split()[5:7]
            peak = int(peak)
            for n, alu in enumerate(peak_lst[peak_info]):
                start, end = alu.split()[1:3]
                alu_site = int(start) if strand == '+' else int(end)
                alu_d = abs(alu_site - peak)
                if n == 0:
                    final_alu = alu
                    final_d = alu_d
                else:
                    if alu_d < final_d:
                        final_alu = alu
                        final_d = alu_d
            if final_alu in alu_exp:
                if expression > alu_exp[final_alu]:
                    alu_exp[final_alu] = expression
                    alu_info[final_alu] = peak_info
            else:
                alu_exp[final_alu] = expression
                alu_info[final_alu] = peak_info
    with open(output, 'w') as out:
        for alu in alu_info:
            out.write(alu_info[alu] + '\t' + alu + '\n')


def parse_alu(options):
    if options['-a']:
        alu_f = options['-a']
        index = {'c': 0, 's': 1, 'e': 2, 'n': 3, 't': 5}
        check = False
    else:
        alu_f = options['-r']
        index = {'c': 5, 's': 6, 'e': 7, 'n': 10, 't': 9}
        check = True
    with open(alu_f, 'r') as f:
        for line in f:
            items = line.rstrip().split()
            if check and items[12] != 'Alu':
                    continue
            chrom = items[index['c']]
            start = int(items[index['s']])
            end = int(items[index['e']])
            name = items[index['n']]
            strand = items[index['t']]
            yield (chrom, start, end, name, strand)


def check_span(peak, span_cutoff):
    sites = [int(x) for x in peak[14].split('|')]
    if peak.strand == '+':
        span = np.percentile(sites, 75) - peak.start
    else:
        span = peak.end - np.percentile(sites, 25)
    if span <= span_cutoff:
        return True
    else:
        return False


if __name__ == '__main__':
    from docopt import docopt
    annotate(docopt(__doc__, version=__version__))
