#!/bin/bash
if [ $# != 4 ]
then
    echo ""
    echo "  Version:         v1.0 (2024-02-01)"
    echo "  Authors:         Qingshan Liu"
    echo "  Description:     Call SNV/Indel using GATK4 Docker"
    echo ""
    echo "  Usage:          bash $0 <bam> <genome> <Output_dir> <Sample>"
    echo ""
    echo "  Parameters:"
    echo "      <bam>          Input BAM file"
    echo "      <genome>       Reference genome fasta"
    echo "      <Output_dir>   Output directory"
    echo "      <Sample>       Sample name"
    echo ""
    exit 1
fi

# 参数设置
bam_file=$1
genome=$2
output_dir=$3
Sample=$4

# 获取脚本所在目录
Pipeline=$(cd `dirname $0`; pwd)
data_dir=$Pipeline/DATABASE

# 创建输出目录
mkdir -p $output_dir
cd $output_dir

# 设置日志文件
step_log=$output_dir/$Sample.SNV.steps.log
detail_log=$output_dir/$Sample.SNV.details.log
exec >${detail_log} 2>&1
set -e -x

# 设置线程数和内存
nt=16
memory=32g

data_dir=$Pipeline/DATABASE
# 设置GATK Docker命令
GATK="docker run -v $Pipeline:$Pipeline -v $data_dir:$data_dir -v $output_dir:$output_dir -w $output_dir broadinstitute/gatk:4.4.0.0"

echo "[$(date)] Starting SNV/Indel calling pipeline..."

#-----------------------------------------------------------------------
# 3. SNP和INDEL检测(HaplotypeCaller)
#-----------------------------------------------------------------------
echo "[$(date)] Calling variants..."

# 3.1 使用HaplotypeCaller检测变异
$GATK gatk --java-options "-Xmx${memory}" HaplotypeCaller \
    -R $genome \
    -I $bam_file \
    -O $Sample.raw.vcf \
    --native-pair-hmm-threads $nt

#-----------------------------------------------------------------------
# 4. 变异过滤
#-----------------------------------------------------------------------
echo "[$(date)] Filtering variants..."

# 4.1 分离SNP和INDEL
$GATK gatk --java-options "-Xmx${memory}" SelectVariants \
    -R $genome \
    -V $Sample.raw.vcf \
    --select-type-to-include SNP \
    -O $Sample.raw.snps.vcf

$GATK gatk --java-options "-Xmx${memory}" SelectVariants \
    -R $genome \
    -V $Sample.raw.vcf \
    --select-type-to-include INDEL \
    -O $Sample.raw.indels.vcf

# 4.2 SNP过滤
$GATK gatk --java-options "-Xmx${memory}" VariantFiltration \
    -R $genome \
    -V $Sample.raw.snps.vcf \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "FAIL" \
    -O $Sample.filtered.snps.vcf

# 4.3 INDEL过滤
$GATK gatk --java-options "-Xmx${memory}" VariantFiltration \
    -R $genome \
    -V $Sample.raw.indels.vcf \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "FAIL" \
    -O $Sample.filtered.indels.vcf

# 4.4 合并过滤后的SNP和INDEL
$GATK gatk --java-options "-Xmx${memory}" MergeVcfs \
    -I $Sample.filtered.snps.vcf \
    -I $Sample.filtered.indels.vcf \
    -O $Sample.filtered.vcf

# 对过滤后的vcf文件进行注释
python3 $Pipeline/src/anno.snv.py -i $Sample.filtered.vcf -o $Sample.filtered.anno.xls

#-----------------------------------------------------------------------
# 5. 清理中间文件
#-----------------------------------------------------------------------
echo "[$(date)] Cleaning up..."

rm -f $Sample.raw.snps.vcf* \
      $Sample.raw.indels.vcf*

echo "[$(date)] SNV/Indel calling completed"
echo "[$(date)] SNV finished" > $output_dir/$Sample.SNV.done