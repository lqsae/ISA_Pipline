#!/bin/bash

if [ $# != 3 ]; then
    echo ""
    echo "  Version:         0.0.1 (2020-05-08)"
    echo "  Authors:         SingChen"
    echo "  Description:     Get target per-base from mosdepth per-base.bed.gz and target.bed"
    echo ""
    echo "  Requirements:    bedtk, bedtools, bgzip, tabix"
    echo ""
    echo "  Usage:           bash $0 <mosdepth per-base.bed.gz> <target.bed> <output prefix>"
    echo ""
    echo "  Example:         bash $0 demo/MS21011216F1.per-base.bed.gz iWES_D.hg38.bed demo/MS21011216F1"
    echo ""
    exit 1
fi

perbase=$1
target=$2
prefix=$3

if [ ! -f $perbase ]; then
    echo "Error: $perbase not exists"
    exit 1
fi
if [ ! -f $target ]; then
    echo "Error: $target not exists"
    exit 1
fi

export PATH=/work/Software/htslib-1.8/bin/:$PATH
export PATH=/work/Software/bedtools-2.29.0/bin/:$PATH
export PATH=/work/BI/xcwu/tools/bedtk-master/:$PATH

set -e

bedtk flt $target $perbase | bedtools makewindows -b - -w 1 -i src >${prefix}.per-base.split.bed
bedtk flt $target ${prefix}.per-base.split.bed | bgzip -@ 8 >${prefix}.per-base.target.bed.gz
tabix -f -C ${prefix}.per-base.target.bed.gz

set +e
