# !/usr/bin/env python
# -*-coding:utf-8 -*-

import logging
import os.path
import pandas as pd
import pysam
import pyfaidx

fasta = "/work/Database/GATK_db/Other_RefGenome/hg38_rmContig/hg38.fa"

def do_coverage(bed_fname, bam_fname, fasta, by_count=False, min_mapq=0, processes=10, gc=False,norm=False):
    """Calculate coverage in the given regions from BAM read depths."""
    if not ensure_bam_sorted(bam_fname):
        raise RuntimeError("BAM file %s must be sorted by coordinates"
                           % bam_fname)
    cnarr = interval_coverages(bed_fname, bam_fname, fasta, by_count, min_mapq,
                               processes, gc=gc, norm=False)
    return cnarr


def interval_coverages(bed_fname, bam_fname, fasta, by_count, min_mapq, processes, gc, norm=False):
    if by_count:  # only support read count by now
        results = interval_coverages_count(bed_fname, bam_fname, fasta, min_mapq,
                                           gc, processes)
        read_counts, cna_rows = zip(*results)
        read_counts = pd.Series(read_counts)
        cna_rows = list(cna_rows)
    else:
        pass
    tot_mapped_reads = bam_total_reads(bam_fname)

    return cna_rows


def interval_coverages_count(bed_fname, bam_fname, fasta, min_mapq, gc, procs=1):
    """Calculate depth in the BAM file at each interval."""
    if procs==1:
        regions = read_auto(bed_fname)
        bamfile = pysam.Samfile(bam_fname, 'rb')
        for chrom, subregions in regions.items():
            logging.info("Processing chromosome %s of %s",
                            chrom, os.path.basename(bam_fname))
            for count, row in _rdc(bamfile, subregions, fasta, min_mapq,gc):
                yield [count, row]
    else:
        from concurrent import futures
        with futures.ProcessPoolExecutor(procs) as pool:
            args_iter = ((bamfile, subregions, fasta, min_mapq, gc)
                         for chr, subregions in regions.items())
            for chunk in pool.map(_rdc, args_iter):
                for count, row in chunk:
                    yield [count, row]

def _rdc(*args):
    """Wrapper for parallel."""
    return list(_rdc_chunk(*args))


def _rdc_chunk(bamfile, regions, fasta, min_mapq,gc):
    if isinstance(bamfile, str):
        bamfile = pysam.Samfile(bamfile, 'rb')
    for chrom, start, end in regions:
        yield region_depth_count(bamfile, chrom, start, end, fasta, min_mapq, gc)


def region_depth_count(bamfile, chrom, start, end, fasta, min_mapq, gc):
    """Calculate depth of a region via pysam count.
    i.e. counting the number of read starts in a region, then scaling for read
    length and region width to estimate depth.
    Coordinates are 0-based, per pysam.
    """

    def filter_read(read):
        """True if the given read should be counted towards coverage."""
        return not (read.is_duplicate
                    or read.is_secondary
                    or read.is_unmapped
                    or read.is_qcfail
                    or read.mapq < min_mapq)

    count = 0
    bases = 0
    for read in bamfile.fetch(reference=chrom, start=start, end=end):
        if filter_read(read):
            count += 1
            # Only count the bases aligned to the region
            rlen = read.query_alignment_length
            if read.pos < start:
                rlen -= start - read.pos
            if read.pos + read.query_alignment_length > end:
                rlen -= read.pos + read.query_alignment_length - end
            bases += rlen
    depth = bases / (end - start) if end > start else 0
    pos_info = "{}:{}-{}".format(chrom,start,end)
    row = [pos_info, depth]
    if gc:
        interval = "{}:{}-{}".format(chrom, start, end)
        gc_ratio = calculate_gc_lo(fasta_extract_regions(fasta, interval))
        row = [pos_info, gc_ratio, depth]
    return count, row


def read_auto(infile):
    chrom_region = {}
    try:
        with open(infile) as fo:
            for line in fo:
                if line.startswith("Chr"):
                    logging.info("Please check {} format".format(infile))
                line = line.strip().split("\t")
                chrom_region.setdefault(line[0],list()).append([line[0], int(line[1]), int(line[2])])
    except pd.io.common.EmptyDataError:
        logging.info("{} file is empty".format(infile))
    return chrom_region


# bam function to check bam file
def ensure_bam_sorted(bam_fname, span=50):
    from itertools import islice
    
    def out_of_order(read, prev):
        return not (prev is None or
                    read.tid != prev.tid or
                    prev.pos <= read.pos)

    # ENH - repeat at 50%, ~99% through the BAM
    bam = pysam.Samfile(bam_fname, 'rb')
    last_read = None
    for read in islice(bam, span):
        if out_of_order(read, last_read):
            return False
        last_read = read
    bam.close()
    return True


def bam_total_reads(bam_fname):
    # maybe norm by total reads
    table = idxstats(bam_fname, drop_unmapped=True)
    return table.mapped.sum()



def idxstats(bam_fname, drop_unmapped=False):
    """
    Get chromosome names, lengths, and number of mapped/unmapped reads.
    Use the BAM index (.bai) to get the number of reads and size of each
    chromosome.
    """
    from io import StringIO
    handle = StringIO(pysam.idxstats(bam_fname, split_lines=False))
    table = pd.read_csv(handle, sep='\t', header=None,
                        names=['chromosome', 'length', 'mapped', 'unmapped'])
    if drop_unmapped:
        table = table[table.mapped != 0].drop('unmapped', axis=1)
    return table


def calculate_gc_lo(subseq):
    """Calculate the GC and lowercase (RepeatMasked) content of a string."""
    cnt_at_lo = subseq.count('a') + subseq.count('t')
    cnt_at_up = subseq.count('A') + subseq.count('T')
    cnt_gc_lo = subseq.count('g') + subseq.count('c')
    cnt_gc_up = subseq.count('G') + subseq.count('C')
    tot = float(cnt_gc_up + cnt_gc_lo + cnt_at_up + cnt_at_lo)
    if not tot:
        return 0.0
    frac_gc = (cnt_gc_lo + cnt_gc_up) / tot
    return frac_gc


def fasta_extract_regions(fa_fname, intervals):
    with pyfaidx.Fasta(fa_fname, as_raw=True) as fa_file:
        interval = Region(intervals)
        _chrom = interval.chr
        start = interval.start
        end = interval.end
        return fa_file[_chrom][start:end]


class Region:
    def __init__(self, corrdinate):
        chr_list = ["chrX", "chrY", "chrM"] + ["chr{}".format(i) for i in range(1, 23)]
        corrdinate = corrdinate.replace("-", ":").split(":")
        self.chr = corrdinate[0]
        if len(corrdinate) > 1:
            self.start = int(corrdinate[1])
            self.end = int(corrdinate[2])
        else:
            self.start = float("-inf")
            self.end = float("inf")
        if self.chr not in chr_list:
            raise ValueError("please check the input chromosome ID")
        elif self.start > self.end:
            raise ValueError("please check the start position and end position")


def main():
    import sys
    if len(sys.argv) != 6:
        print("Please check your input")
        print("python %s <bam> <bed> <fasta> <sample> <GC N|Y>"% sys.argv[0])
    in_bam = sys.argv[1]
    Sample = sys.argv[4]
    in_bed = sys.argv[2]
    fasta = sys.argv[3]
    GC_status = sys.argv[5]
    if GC_status == "Y":
        header=["Pos","GC",Sample]
        gc_status = True
    else:
        header=["Pos",Sample]
        gc_status = False
    
    df = do_coverage(in_bed, in_bam, fasta, by_count=True, min_mapq=0, processes=1, norm=False,gc=gc_status)
    fo = open("{}.depth.txt".format(Sample),"w")
    fo.write("\t".join(header)+"\n")
    for line in df:
        line = [str(i) for i in line]
        fo.write("\t".join(line)+"\n")
        

if __name__ == '__main__':
    main()
