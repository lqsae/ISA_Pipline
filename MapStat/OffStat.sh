#!/bin/bash
if [ $# != 5 ]
then
        echo ""
        echo "  Version:         v1.0.0 (2020-2-28)"
        echo "  Authors:         Li Min"
        echo "  Description:     stat off target "
        echo ""
        echo "  Usage:           bash $0 <bam file> <BED>  <hg19|hg38>  <work_dir> <Sample>"
        echo ""
        exit 1
fi

Bam_file=$1
Bed_file=$2
Ref_version=$3
Output_dir=$4
work_dir=$Output_dir
Sample=$5

Step_Name="OffTarget"

#source software environment
Pipeline=$(cd `dirname $0`;cd ..; pwd)

source $Pipeline/Source/software.sh
source $Pipeline/Source/log.sh

# creat work directory
mkdir -p $work_dir
cd $work_dir

# creat log file
step_log=$work_dir/$Sample.$Step_Name.steps.log
log_Command $0 $@ > $step_log

detail_log=$work_dir/$Sample.$Step_Name.details.log
exec >${detail_log} 2>&1
set -e -x

# detect Bed file exists
if [ $Bed_file ]
then
	check_file "BED" $Bed_file
fi

# reference databases
if [ $Ref_version != "hg19" ] && [ $Ref_version != "hg38" ]
then
	echo "Error: $Ref_version version not correct, please specified hg19 or hg38"
	exit 1;
fi

source $Pipeline/Source/$Ref_version.env.sh

#-----------------------------------------------------------------------
# 1. Mapping and sorting
#-----------------------------------------------------------------------

log_event auto.1 "Stat begin" >> $step_log

bedtools bamtobed -i $Bam_file |grep -v "chrM"> $work_dir/$Sample.bam.bed
bedtools closest -t first -d -a $work_dir/$Sample.bam.bed -b $Bed_file > $work_dir/$Sample.stat.bed

python $Pipeline/MapStat/stat_off.py $work_dir/$Sample.stat.bed $Sample > $work_dir/$Sample.offTarget.stat

rm -rf $Sample.bam.bed
rm -rf $Sample.stat.bed

set +e

log_event auto.2 "Stat end" >> $step_log
