#!/bin/bash

if [ $# != 7 ]
then
	echo ""
	echo "  Version:         v2.01 (2024-06-06)"
	echo "  Authors:         Qingshan Liu"
	echo "  Description:     Integration Site Analysis from bam file"
	echo ""
	echo "  Usage:           bash $0 <bam> <genome> <ltr_bed>  <gene_exon_bed> <gtf_file>  <Output_dir>  <SampleName>"
        echo ""
	exit 1
fi


bam=$1
genome=$2
ltr_bed=$3
gene_exon_bed=$4
gtf_file=$5
Output_dir=$6
work_dir=$Output_dir
Sample=$7

Step_Name="breakpoint"

#source software environment
Pipeline=$(cd `dirname $0`; pwd)


Genome=$genome


#source $Pipeline/Source/software.sh
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

# cpu number used
nt=6

log_event 1.5 "$Sample $Step_Name start" >> $step_log

#-----------------------------------------------------------------------
# 0.1 获取整合到LTR的reads
#-----------------------------------------------------------------------
#sambamba markdup  $bam  -r  $Sample.rmdup.bam
sambamba view  $bam -L $ltr_bed   -f bam  -o  $Sample.map2LTR.bam

#-----------------------------------------------------------------------
# 0.2 获取 比对到宿主和整合序列的split reads
#-----------------------------------------------------------------------
python3 $Pipeline/src/get_split_reads.py -b  $Sample.map2LTR.bam  -l $ltr_bed   > $Sample.split_event.txt

#-----------------------------------------------------------------------
# 0.3 统计整合位点信息
#-----------------------------------------------------------------------
python3 $Pipeline/src/stat_break_point.py  $Sample.split_event.txt  $Sample  $work_dir $bam


#-----------------------------------------------------------------------
# 0.4 注释 整合位点
#-----------------------------------------------------------------------
cut -f2-14 $Sample.insertion_break_point.txt |grep -v support_reads >  $Sample.insertion_break_point.bed
bedtools intersect -a $Sample.insertion_break_point.bed -b $gene_exon_bed -wao  > $Sample.insertion_break_point.anno.tmp.txt
python3 $Pipeline/src/anno.py $gtf_file  $Sample.insertion_break_point.anno.tmp.txt   > $Sample.insertion_break_point.anno.xls

#-----------------------------------------------------------------------	
# 0.5 获取 整合位点附近的reads
#-----------------------------------------------------------------------	
samtools sort $Sample.insertion_break_point.check.bam -O BAM -o $Sample.insertion_break_point.check.sorted.bam
samtools index $Sample.insertion_break_point.check.sorted.bam

#-----------------------------------------------------------------------
# igv 可视化
create_report $Sample.insertion_break_point.bed \
   --fasta $Genome \
   --tracks  $Sample.insertion_break_point.check.sorted.bam \
   --output $Sample.insertion.igv.html

#-----------------------------------------------------------------------
log_event 1.5 "$Sample $Step_Name finished" >> $step_log
#-----------------------------------------------------------------------