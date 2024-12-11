#! /work/Software/Anaconda3-2019.10/bin/python

import sys
import time
from pathlib import Path
from subprocess import call

doc = f"""
  Version:         v0.0.1 (2021-04-07)
  Authors:         SingChen
  Description:     Count soft-clip in target regions

  Usage:          {__file__} <Sample.realigned.bam> <Target.bed>
"""


def main(bam, bed):
    soft_bam = Path(f"/tmp/soft.{int(time.time())}.bam")
    soft_bed = Path(f"/tmp/soft.{int(time.time())}.bed")
    call(f"cut -f 1-3 {bed} > {soft_bed}", shell=True)
    call(
        f"samtools view -@ 4 -h {bam} "
        "| awk '$6 ~ /S/ || $1 ~ /@/' "
        "| grep -v M[1-9]S$'\t' "
        "| grep -v $'\t'[1-9]S "
        f"| samtools view -@ 4 -o {soft_bam}",
        shell=True)
    call(f"bedtools coverage -a {soft_bed} -b {soft_bam} | cut -f 1-4",
         shell=True)
    soft_bam.unlink()
    soft_bed.unlink()


if __name__ == "__main__":
    if len(sys.argv) < 3 or not (sys.argv[1].endswith(".bam")
                                 and sys.argv[2].endswith(".bed")):
        sys.exit(doc)
    else:
        main(*sys.argv[1:])
