#!/usr/bin/env python

'''
Usage: rm_pcr.py [options] <rampage>...

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
    --min=MIN                      Minimum read counts. [default: 1]
'''

import os
import os.path
from collections import defaultdict
from seqlib.path import create_dir
from seqlib.ngs import check_bam
from seqlib.seq import dna_to_rna

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def rm_pcr(options):
    '''
    Remove PCR duplicates
    '''
    # parse options
    rampage_lst = options['<rampage>']
    thread = int(options['--thread'])
    folder = create_dir(options['--output'])
    min_read = int(options['--min'])
    # remove PCR duplicates
    if thread > 1:
        from multiprocessing import Pool
        p = Pool(thread)
        result = []
        for rampage in rampage_lst:  # for each rampage replicate
            for chrom in check_bam(rampage).references:  # for each chromosome
                if chrom.startswith('chr'):  # only parse normal chromosome
                    result.append(p.apply_async(remove_pcr, args=(rampage,
                                                                  chrom,)))
        p.close()
        p.join()
        collapsed_pairs = []
        for r in result:
            collapsed_pairs.extend(r.get())
    else:
        collapsed_pairs = remove_pcr(rampage_lst)
    total_counts = write_signal(collapsed_pairs, folder, min_read)
    with open(os.path.join(folder, 'total_counts.txt'), 'w') as out:
        out.write('%d\n' % total_counts)


def write_signal(pairs, folder, min_read):
    bed5p = open(os.path.join(folder, 'rampage_plus_5end.bed'), 'w')
    bed5m = open(os.path.join(folder, 'rampage_minus_5end.bed'), 'w')
    bed3p = open(os.path.join(folder, 'rampage_plus_3read.bed'), 'w')
    bed3m = open(os.path.join(folder, 'rampage_minus_3read.bed'), 'w')
    bed = open(os.path.join(folder, 'rampage_link.bed'), 'w')
    uniq_pairs = defaultdict(int)
    for pair in pairs:
        info = pair.rsplit('\t', 1)[0]  # remove barcode info
        uniq_pairs[info] += 1
    total_counts = 0
    for pair in uniq_pairs:
        read_count = uniq_pairs[pair]
        if read_count < min_read:  # not enough reads
            continue
        total_counts += read_count
        pair_info = pair.split()
        r1_chrom, r1_start, r1_end, strand = pair_info[:4]  # read1
        r2_chrom, r2_start, r2_end = pair_info[4:]  # read2
        if strand == '+':
            start = int(r1_start)
            end = int(r2_end)
            bed5p.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, start,
                                                       start + 1, strand) *
                        read_count)
            bed3p.write('%s\t%d\t%d\t3read\t0\t%s\n' % (r2_chrom, end - 1,
                                                        end, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            bed.write('\t'.join([r1_chrom, r1_start, r2_end, 'link',
                                 str(read_count), strand, r1_start, r1_start,
                                 '0,0,0', '2', '1,1', offset]) + '\n')
        else:
            start = int(r2_start)
            end = int(r1_end)
            bed5m.write('%s\t%d\t%d\t5end\t0\t%s\n' % (r1_chrom, end - 1,
                                                       end, strand) *
                        read_count)
            bed3m.write('%s\t%d\t%d\t3read\t0\t%s\n' % (r2_chrom, start,
                                                        start - 1, strand) *
                        read_count)
            offset = '0,' + str(end - start - 1)
            bed.write('\t'.join([r1_chrom, r2_start, r1_end, 'link',
                                 str(read_count), strand, r2_start, r2_start,
                                 '0,0,0', '2', '1,1', offset]) + '\n')
    return total_counts


def remove_pcr(bam_f, chrom=None):
    if type(bam_f) is list:  # multiple bam files
        collapsed_pairs = []
        for f in bam_f:
            bam = check_bam(f)
            read2 = fetch_read2(bam, chrom)  # fetch read2
            # fetch read1
            collapsed_pairs.extend(fetch_read1(bam, read2, chrom))
    else:  # single bam file
        bam = check_bam(bam_f)
        read2 = fetch_read2(bam, chrom)  # fetch read2
        collapsed_pairs = fetch_read1(bam, read2, chrom)  # fetch read1
    return collapsed_pairs


def fetch_read2(bam, chrom):
    read2_lst = {}
    for read in bam.fetch(chrom):
        # not read1 or secondary alignment
        if read.is_read1 or read.is_secondary:
            continue
        if not read.is_proper_pair:  # not proper pair
            continue
        if read.get_tag('NH') != 1:  # not unique read
            continue
        chrom = read.reference_name
        mate_chrom = read.next_reference_name
        if chrom != mate_chrom:  # not same chromosome
            continue
        if not read.is_reverse and read.mate_is_reverse:
            strand = '+'
            mate_strand = '-'
        elif read.is_reverse and not read.mate_is_reverse:
            strand = '-'
            mate_strand = '+'
        else:
            continue
        mate_pos = str(read.next_reference_start)
        name = read.query_name
        start = str(read.reference_start)
        end = str(read.reference_end)
        if strand == '+':
            barcode = dna_to_rna(read.query_sequence[:15])
        else:
            barcode = dna_to_rna(read.query_sequence[-15:],
                                 strand=strand)
        read_id = '\t'.join([name, mate_chrom, mate_pos, mate_strand])
        read2_lst[read_id] = [chrom, start, end, barcode]
    return read2_lst


def fetch_read1(bam, read2, chrom):
    collapsed_pairs = set()
    for read in bam.fetch(chrom):
        # not read2
        if read.is_read2:
            continue
        if not read.is_proper_pair:  # not proper pair
            continue
        name = read.query_name
        chrom = read.reference_name
        start = str(read.reference_start)
        # parse strand info
        strand = '+' if not read.is_reverse else '-'
        read_id = '\t'.join([name, chrom, start, strand])
        if read_id not in read2:
            continue
        end = str(read.reference_end)
        # remove PCR duplicates
        collapsed_pairs.add('\t'.join([chrom, start, end, strand] +
                                      read2[read_id]))
    return collapsed_pairs


if __name__ == '__main__':
    from docopt import docopt
    rm_pcr(docopt(__doc__, version=__version__))
