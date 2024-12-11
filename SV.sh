#!/bin/bash
if [ $# != 4 ]
then
        echo ""
        echo "  Version:         v3.0 (2024-02-01)"
        echo "  Authors:         Qingshan Liu"
        echo "  Description:     Call SV from recaled or realigned BAM file"
        echo ""
        echo "  Usage:           bash $0 <bam> <genome> <Output_dir>  <SampleName>"
        echo ""
        echo ""
        exit 1
fi

Bam_file=$1 
genome=$2
Output_dir=$3
work_dir=$Output_dir
Sample=$4
Ref_version=hg38
RefGenome=$genome

Step_Name="SV"

#source software environment
Pipeline=$(cd `dirname $0`; pwd)

gene_bed=$Pipeline/DATABASE/hg38.gene.bed

# source $Pipeline/Source/software.sh
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
nt=10

# check file
check_file "BAM" $Bam_file

#-----------------------------------------------------------------------
#  SV detect
#-----------------------------------------------------------------------

# Delly
log_event auto.2 "Delly begin:" >> $step_log
delly  call -g $RefGenome -o ${Sample}_delly.bcf $Bam_file
bcftools view ${Sample}_delly.bcf >${Sample}_delly.vcf
grep "#" ${Sample}_delly.vcf >${Sample}_delly_filter.vcf
grep -v "#" ${Sample}_delly.vcf| awk '$7=="PASS" {print $0}'  >> ${Sample}_delly_filter.vcf
log_event auto.2 "Delly finished:" >> $step_log

#-----------------------------------------------------------------------
#  SV annotation
#-----------------------------------------------------------------------

log_event auto.3 "Annotation begin:" >> $step_log

# 使用Python脚本进行注释
python3 $Pipeline/src/anno.sv.py \
    -v ${Sample}_delly_filter.vcf \
    -b $gene_bed \
    -o ${Sample}_sv.anno.xlsx \
    -u 2000 \
    -d 2000

log_event auto.3 "Annotation finished:" >> $step_log