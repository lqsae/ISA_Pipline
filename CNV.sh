#!/bin/bash
if [ $# != 6 ]
then
        echo ""
        echo "  Version:         v3.0 (2024-02-01)"
        echo "  Authors:         Qingshan Liu"
        echo "  Description:     Call CNV using cnvkit"
        echo ""
        echo "  Usage:           bash $0 <bam> <genome> <genome_bed> <window_size>   <Output_dir>  <SampleName>"
        echo ""
        echo ""
        exit 1
fi
   
Bam_file=$1
genome=$2
genome_bed=$3
window_size=$4
Output_dir=$5
work_dir=$Output_dir
Sample=$6
Ref_version=hg38
RefGenome=$genome
Step_Name="CNV"

#source software environment
Pipeline=$(cd `dirname $0`; pwd)

Genome=$RefGenome



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
#  CNV detect using cnvkit
#-----------------------------------------------------------------------


log_event auto.0 "CNVkit analysis begin:" >> $step_log
# 滑窗口
bedtools  makewindows -b $genome_bed -w $window_size > $work_dir/genome.window.bed

access_bed=$work_dir/genome.window.bed


log_event auto.1 "CNVkit analysis begin:" >> $step_log

# 1. 运行CNVkit主分析 (flat参考模式)
log_event auto.1.1 "Running CNVkit main analysis:" >> $step_log
cnvkit.py batch $Bam_file \
    --method wgs \
    -n  \
    --fasta $Genome \
    --access $access_bed \
    --processes $nt \
    --scatter --diagram \
    -d ${Sample}_cnvkit_output

# 2. 调用CNV
log_event auto.1.2 "Calling CNVs:" >> $step_log
cnvkit.py call ${Sample}_cnvkit_output/${Sample}.deduped.cns \
    -o ${Sample}_cnvkit_output/${Sample}.call.cns \
    --filter cn \
    --method threshold \
    --ploidy 2 \
    --drop-low-coverage

# 3. 导出分段结果
log_event auto.1.3 "Exporting segments:" >> $step_log
cnvkit.py export seg \
    ${Sample}_cnvkit_output/${Sample}.call.cns \
    -o ${Sample}_cnvkit_output/${Sample}.segments.seg

# 4. 注释CNV区间
log_event auto.1.4 "Annotating CNVs:" >> $step_log
python3 $Pipeline/src/anno.cnv.py \
    -i ${Sample}_cnvkit_output/${Sample}.segments.seg \
    -b $gene_bed \
    -o ${Sample}_cnvkit_output/${Sample}.segments.anno.txt

log_event auto.1.6 "CNVkit analysis finished:" >> $step_log