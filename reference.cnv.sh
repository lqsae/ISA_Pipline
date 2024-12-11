#!/bin/bash

if [ $# != 5 ]
then
    echo ""
    echo "  Version:         v1.0 (2024-02-01)"
    echo "  Authors:         Qingshan Liu"
    echo "  Description:     Reference CNV analysis"
    echo ""
    echo "  Usage:          bash $0 <d4_file> <reference_bed> <target_bed> <output> <Sample>"
    echo ""
    echo "  Parameters:"
    echo "      <d4_file>        Input d4 file"
    echo "      <reference_bed>  Reference bed file"
    echo "      <target_bed>     Target bed file"
    echo "      <output>         Output file prefix"
    echo "      <Sample>         Sample name"
    echo ""
    exit 1
fi

# 获取脚本所在目录
Pipeline=$(cd `dirname $0`; pwd)


# 参数设置
d4_file=$1
reference_bed=$2
target_bed=$3
output=$4
sample=$5

# 创建输出目录
output_dir=$(dirname $output)
mkdir -p $output_dir

cd $output_dir

# 设置日志文件
log_file=${output}.log
exec 1> >(tee -a "$log_file") 2>&1

echo "[$(date)] Starting reference CNV analysis..."
echo "[$(date)] Parameters:"
echo "  D4 file: $d4_file"
echo "  Reference bed: $reference_bed"
echo "  Target bed: $target_bed"
echo "  Output: $output"
echo "  Sample: $sample"

# 检查输入文件是否存在
for file in "$d4_file" "$reference_bed" "$target_bed"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist!"
        exit 1
    fi
done



# 运行 Python 脚本进行分析
echo "[$(date)] Running CNV analysis..."
python3 $Pipeline/src/reference.cnv.py \
    ${d4_file} \
    ${reference_bed} \
    ${target_bed} \
    > ${output}.cnv.txt

# 检查运行结果
if [ $? -eq 0 ]; then
    echo "[$(date)] Analysis completed successfully"
    echo "[$(date)] Results saved to ${output}.cnv.txt"
else
    echo "[$(date)] Error: Analysis failed!"
    exit 1
fi


# '''
# gene	index	start	end	gc	target_corrected_CN	target_raw_CN	reference_CN
# CLC4B147K562	1	0	50	0.36	3.5078284627940772	3.631240673048385	0.6990616781502836
# CLC4B147K562	2	50	100	0.44	3.9006468915462724	4.105225786324214	0.7669647499133825
# CLC4B147K562	3	100	150	0.62	4.309270910731567	4.0503461855039395	1.0250409193216126
# '''



# 计算每个基因的平均 CNV（转换 log2 值为实际拷贝数）
echo "[$(date)] Calculating gene average CNV..."
echo -e "Gene\tAverage_CN" > ${sample}.gene_average.txt
awk -F'\t' '
    function pow2(x) {
        return 2 ^ x
    }
    NR>1 {
        sum[$1] += pow2($6);  # 将 log2 值转换为实际拷贝数后再累加
        count[$1]++;
    }
    END {
        for (gene in sum) {
            printf "%s\t%.2f\n", gene, sum[gene]/count[gene];
        }
    }
' ${output}.cnv.txt | sort -k1,1 >> ${sample}.gene_average.txt


echo "[$(date)] All done!"