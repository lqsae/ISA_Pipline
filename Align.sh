#!/bin/bash
2|if [ $# != 5 ]
then
        echo ""
        echo "  Version:         v2.01 (2024-02-01)"
        echo "  Authors:         Qingshan Liu"
        echo "  Description:     Alignment RefGenome (hg19 or hg38) and generate BAM file"
        echo ""
        echo "  Usage:           bash $0  <genome> <trim.R1.fastq.gz>  <trim.R2.fastq.gz>  <Output_dir>  <SampleName> "
        echo ""
        echo "       dedup   :   for WGS, WES, and other capture DNA sequencing"
        echo "       nodedup :   for PCR product sequencing"
        echo ""
        exit 1
fi

R1=$1
R2=$2
genome=$3
Output_dir=$4
work_dir=$Output_dir
Sample=$5

Step_Name="Alignment"

# 获取脚本所在目录
Pipeline=$(cd `dirname $0`; pwd)
data_dir=$Pipeline/DATABASE

# 创建工作目录
mkdir -p $work_dir
cd $work_dir

# 创建日志文件
step_log=$work_dir/$Sample.$Step_Name.steps.log
detail_log=$work_dir/$Sample.$Step_Name.details.log
exec >${detail_log} 2>&1
set -e -x

# 设置线程数
nt=16

# 参考基因组版本和路径
Ref_version="hg38"
ref=$genome

#-----------------------------------------------------------------------
# 1. 准备参考基因组索引
#-----------------------------------------------------------------------
echo "[$(date)] Checking reference genome indexes..."

# 检查并创建参考基因组的序列字典
if [ ! -f ${genome%.fa}.dict ] && [ ! -f ${genome%.fasta}.dict ]; then
    echo "[$(date)] Creating sequence dictionary..."
    picard CreateSequenceDictionary \
        R=$genome \
        O=${genome%.fa}.dict
fi

# 检查并创建BWA索引
if [ ! -f ${genome}.bwt ]; then
    echo "[$(date)] Creating BWA index..."
    bwa index $genome
fi

# 检查并创建samtools索引
if [ ! -f ${genome}.fai ]; then
    echo "[$(date)] Creating FASTA index..."
    samtools faidx $genome
fi

#-----------------------------------------------------------------------
# 2. Mapping and sorting
#-----------------------------------------------------------------------
# 2.1 BWA mem + samtools: alignment and multi-thread sorting
#-----------------------------------------------------------------------
echo "[$(date)] Aligning reads to reference genome..."

bwa mem \
    -t $nt \
    -M \
    -R "@RG\tID:${Sample}\tLB:${Sample}\tSM:${Sample}\tPL:ILLUMINA" \
    $ref $R1 $R2 | \
samtools sort \
    -@ $nt \
    -m 4G \
    -o ${Sample}.sorted.bam \
    -

samtools index -@ $nt ${Sample}.sorted.bam


#-----------------------------------------------------------------------
# 2.2 Picard: mapping metrics
#-----------------------------------------------------------------------
echo "[$(date)] Generating alignment metrics..."

# 计算比对质量指标
# picard CollectAlignmentSummaryMetrics \
#     R=$ref \
#     I=${Sample}.sorted.bam \
#     O=${Sample}.aln_metrics.txt \
#     VALIDATION_STRINGENCY=LENIENT \
#     USE_JDK_DEFLATER=true \
#     USE_JDK_INFLATER=true

# # 计算插入片段大小指标
# picard CollectInsertSizeMetrics \
#     I=${Sample}.sorted.bam \
#     O=${Sample}.insert_metrics.txt \
#     H=${Sample}.insert_histogram.pdf \
#     VALIDATION_STRINGENCY=LENIENT \
#     USE_JDK_DEFLATER=true \
#     USE_JDK_INFLATER=true

# # 计算GC偏好性指标
# picard CollectGcBiasMetrics \
#     I=${Sample}.sorted.bam \
#     O=${Sample}.gc_bias_metrics.txt \
#     CHART=${Sample}.gc_bias_metrics.pdf \
#     S=${Sample}.gc_bias_summary.txt \
#     R=$ref \
#     VALIDATION_STRINGENCY=LENIENT \
#     USE_JDK_DEFLATER=true \
#     USE_JDK_INFLATER=true

#-----------------------------------------------------------------------
# 2.3 Picard: remove PCR duplicates
#-----------------------------------------------------------------------
echo "[$(date)] Marking duplicates..."

picard MarkDuplicates \
    I=${Sample}.sorted.bam \
    O=${Sample}.deduped.bam \
    M=${Sample}.dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true \
    CREATE_INDEX=true

# 生成比对统计信息
perl $Pipeline/MapStat/stat_Map.pl $Sample ${Sample}.dedup_metrics.txt > $Sample.Map.stat

#-----------------------------------------------------------------------
# 2.4 生成深度文件
#-----------------------------------------------------------------------
echo "[$(date)] Generating depth files..."

mosdepth \
    -t $nt \
    -x \
    --d4 \
    $Sample \
    ${Sample}.deduped.bam

# 计算覆盖度统计
bamcov \
    -r $data_dir/hg38.genome.bed \
    -d $Sample.per-base.d4 \
    > $Sample.bamcov.stat

#-----------------------------------------------------------------------
# 2.5 清理临时文件
#-----------------------------------------------------------------------
echo "[$(date)] Cleaning up temporary files..."

rm -f ${Sample}.sorted.bam*

echo "[$(date)] $Sample alignment finished"
echo "[$(date)] Align finished" > $work_dir/$Sample.Align.done