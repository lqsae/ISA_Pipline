#! /work/Software/Anaconda3-2019.10/bin/python

import sys
from collections import Counter, defaultdict
from functools import partial
from multiprocessing import Pool

import numpy as np
import pysam

doc = f"""
  Version:      v0.0.2 (2021-05-08)
  Authors:      SingChen
  Description:  Statistic Genome Coverage for WGS and WES

  Usage:        mosdepth -t 4 -x <Sample_name> <bam>
                {__file__} <per-base.bed> <RefGenome size> <Sample_name>
                
  Example:      {__file__} MS20041155P1.per-base.bed.gz 3095693983 MS20041155P1.Genome
"""


def calc_contig(contig, per_base, index):
    print("+", contig)
    cell = defaultdict(int)
    bed = pysam.TabixFile(per_base, index=index, parser=pysam.asTuple())
    for record in bed.fetch(contig):
        contig, start, end, depth = str(record).strip().split("\t")
        length = int(end) - int(start)
        cell[depth] += length
    bed.close()
    print("-", contig)
    return cell


def get_contigs(per_base, index):
    bed = pysam.TabixFile(per_base, index=index)
    contigs = bed.contigs
    bed.close()
    return contigs


def get_coverage_stat(counter, genome_size, sample, x1, avg_depth):
    out = sample + ".Coverage.stat"
    header = ["Sample", "TotalBases", "CovBases", "CovRatio(%)", "Ave_Depth"]
    row = [sample, str(genome_size), str(x1), f"{x1/genome_size*100:.2f}", f"{avg_depth:.2f}"]
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")


def get_depth_freq(counter, genome_size, sample, x1):
    out = sample + ".Depth.freq"
    header = ["Depth", "Freq", "Ratio(%)", "Cumulation(%)"]
    cum = 0
    depth_q20 = 0
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        for depth, length in sorted(counter.items(), key=lambda t: int(t[0])):
            if depth != "0":
                freq = length / x1 * 100
                cum += freq
                if not depth_q20 and cum >= 20:
                    depth_q20 = depth
                row = [depth, str(length), f"{freq:.3f}", f"{cum:.3f}"]
                f.write("\t".join(row) + "\n")
    return depth_q20


def get_depth_stat(counter, genome_size, sample):
    out = sample + ".Depth.stat"
    header = ["Sample", "Ave_Depth(X)", ">=1X", ">=10X", ">=20X", ">=30X", ">=50X", ">=100X", ">=150X", ">=200X"]
    cov_length = 0
    x0 = x1 = x10 = x20 = x30 = x50 = x100 = x150 = x200 = 0
    data = 0
    for depth, length in counter.items():
        depth = int(depth)
        if depth >= 200:
            x200 += length
            x150 += length
            x100 += length
            x50 += length
            x30 += length
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 150:
            x150 += length
            x100 += length
            x50 += length
            x30 += length
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 100:
            x100 += length
            x50 += length
            x30 += length
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 50:
            x50 += length
            x30 += length
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 30:
            x30 += length
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 20:
            x20 += length
            x10 += length
            x1 += length
        elif depth >= 10:
            x10 += length
            x1 += length
        elif depth >= 1:
            x1 += length
        else:
            x0 += length
        data += depth * length
    avg_depth = data / x1
    row = [
        sample,
        f"{avg_depth:.2f}",
        f"{x1/genome_size*100:.2f}%",
        f"{x10/genome_size*100:.2f}%",
        f"{x20/genome_size*100:.2f}%",
        f"{x30/genome_size*100:.2f}%",
        f"{x50/genome_size*100:.2f}%",
        f"{x100/genome_size*100:.2f}%",
        f"{x150/genome_size*100:.2f}%",
        f"{x200/genome_size*100:.2f}%",
    ]
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")
    return x1, avg_depth


def get_cv_cov20(per_base, index, avg_depth, sample):
    array = []
    avg_p20 = avg_depth * 0.2
    small = large = 0
    bed = pysam.TabixFile(per_base, index=index, parser=pysam.asTuple())
    for record in bed.fetch():
        contig, start, end, depth = str(record).strip().split("\t")
        length = int(end) - int(start)
        depth = int(depth)
        for i in range(length):
            array.append(depth)
        if depth >= avg_p20:
            large += length
        else:
            small += length
    bed.close()
    cv = np.std(array) / np.mean(array)
    q20 = large / (large + small)
    return cv, q20


def main(per_base, genome_size, sample):
    genome_size = int(genome_size)
    counter = Counter()
    index = per_base + ".csi"
    contigs = get_contigs(per_base, index)
    func = partial(calc_contig, per_base=per_base, index=index)
    with Pool(8) as p:
        cells = p.map(func, contigs)
    for cell in cells:
        counter.update(cell)
    x1, avg_depth = get_depth_stat(counter, genome_size, sample)
    get_coverage_stat(counter, genome_size, sample, x1, avg_depth)
    depth_q20 = get_depth_freq(counter, genome_size, sample, x1)
    fold80 = avg_depth / int(depth_q20)
    cv, q20 = get_cv_cov20(per_base, index, avg_depth, sample)
    stat = sample + ".Depth.stat"
    out = sample + ".Depth-fold80.stat"
    with open(stat) as fi, open(out, "w+") as fo:
        fo.write(next(fi).strip() + "\tFold80\tCV\t>=20%X\n")
        fo.write(next(fi).strip() + f"\t{fold80:.4f}\t{cv:.2f}\t{q20*100:.2f}%\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(doc)
    else:
        main(*sys.argv[1:])
