# 生信分析流程说明

## 流程概述
本流程用于处理全基因组测序数据，包含质控、比对、变异检测等多个模块。支持并行分析和断点续跑。


1. 配置数据库路径
编辑config.ini文件，设置必要的参数：
```ini
[paths]
genome = /path/to/reference/genome.fa
gene_bed = /path/to/gene.bed
ltr_bed = /path/to/ltr.bed
gene_exon_bed = /path/to/gene_exon.bed
gtf_file = /path/to/annotation.gtf
```

### 运行方式

#### 1. 运行完整流程
```bash
python3 pipeline.py \
  --fq1 input.R1.fq.gz \
  --fq2 input.R2.fq.gz \
  --sample SAMPLE_ID \
  --output-dir OUTPUT_DIR \
  --platform illumina
```

#### 2. 跳过特定步骤
```bash
python3 pipeline.py \
  --fq1 input.R1.fq.gz \
  --fq2 input.R2.fq.gz \
  --sample SAMPLE_ID \
  --output-dir OUTPUT_DIR \
  --skip QC,Alignment  # 跳过QC和比对步骤
```

#### 3. 强制重新运行
```bash
python3 pipeline.py \
  --fq1 input.R1.fq.gz \
  --fq2 input.R2.fq.gz \
  --sample SAMPLE_ID \
  --output-dir OUTPUT_DIR \
  --force-all  # 忽略已完成状态，重新运行所有步骤
```

### 参数说明
必需参数:
- `--fq1`: Read1 FASTQ文件路径
- `--fq2`: Read2 FASTQ文件路径
- `--sample`: 样本名称，用于命名输出文件
- `--output-dir`: 分析结果输出目录

可选参数:
- `--platform`: 测序平台，可选值：illumina/nextera/bgi，默认：illumina
- `--force-all`: 强制重新运行所有步骤，忽略已完成状态
- `--skip`: 跳过指定步骤，多个步骤用逗号分隔
  - 可选值：QC,Alignment,SNV,IS,SV,CNV,RefCNV

### 运行环境要求
- 操作系统: Linux
- Python: 3.6+
- 内存: ≥32GB
- 存储空间: 样本数据大小的5-10倍

### 输出目录结构
```
OUTPUT_DIR/
├── QC/              # 质控结果
├── Align/           # 比对结果
├── SNV/             # 变异检测结果
├── IS/              # 插入序列分析结果
├── SV/              # 结构变异结果
├── CNV/             # 拷贝数变异结果
├── RefCNV/          # 参考CNV分析结果
├── logs/            # 运行日志
└── pipeline_status/ # 流程状态文件
```

### 运行状态监控
- 查看运行状态：`cat OUTPUT_DIR/logs/[sample].pipeline.log`
- 检查步骤完成情况：`ls OUTPUT_DIR/pipeline_status/`
- 查看错误信息：`cat OUTPUT_DIR/logs/[sample].error.log`

### 常见问题处理

1. 断点续跑
   - 检查pipeline_status目录确认已完成步骤
   - 直接重新运行命令，会自动跳过已完成步骤

2. 强制重新运行某步骤
   - 删除对应的.done文件
   - 使用--force-all参数

### 依赖软件
生信分析软件
- cutadapt
- BWA
- Samtools
- GATK
- Delly
- CNVkit
- Python 3.6+
- picard

python 依赖包
- igv-reports
- pysam
- intervaltree
- pybedtools
- cyvcf2
- pyd4
- scipy
- numpy



## 模块说明

### 1. QC质控分析

**调用命令**: `QC.sh [fq1] [fq2] [genome] [output_dir] [sample]`

**输出目录结构**:
QC/[sample]/
├── [sample].trim.R1.fastq.gz    # 去除接头和低质量序列后的Read1数据
├── [sample].trim.R2.fastq.gz    # 去除接头和低质量序列后的Read2数据
├── [sample].Raw.stat            # 原始数据的碱基质量、GC含量等统计信息
├── [sample].Clean.stat          # 质控后数据的碱基质量、GC含量等统计信息
├── [sample].trim.log            # 质控过程的详细日志，包含过滤统计信息
├── [sample].trim.qual.stat      # 碱基质量分布的详细统计数据
└── [sample].Ratio.stat          # 数据过滤前后的reads数量、碱基数等比例统计


### 2. 序列比对分析
**调用命令**: `Align.sh [fq1] [fq2] [genome] [output_dir] [sample]`

**输出目录结构**:
Align/[sample]/
├── [sample].deduped.bam         # 去除PCR重复后的比对结果文件
├── [sample].deduped.bai         # BAM文件的索引，用于快速访问
├── [sample].per-base.d4         # 每个碱基位置的测序深度统计
├── [sample].Map.stat            # 比对率、覆盖度等统计信息
├── [sample].dedup_metrics.txt   # PCR重复率等去重统计指标
└── [sample].bamcov.stat         # 基因组不同区域的覆盖度分布统计


### 3. SNV分析
**调用命令**: `SNV.sh [bam] [genome] [output_dir] [sample]`

**输出目录结构**:
SNV/[sample]/
├── [sample].raw.vcf             # GATK HaplotypeCaller检测的原始变异结果
├── [sample].filtered.snps.vcf   # 应用质量过滤条件后的SNP结果
├── [sample].filtered.indels.vcf # 应用质量过滤条件后的INDEL结果
└── [sample].filtered.vcf        # 合并后的高质量SNP和INDEL变异结果

**输出文件格式说明**:

`[sample].filtered.vcf`: 最终过滤后的变异结果
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
```
列说明:
- CHROM: 染色体编号
- POS: 变异位置
- ID: 变异ID(若无则为'.')
- REF: 参考基因组序列
- ALT: 变异序列
- QUAL: 变异质量值
- FILTER: 过滤标记(PASS表示通过质控)
- INFO: 变异注释信息
  - VT: 变异类型(SNP/INDEL)
  - QD: 变异质量深度标准化得分
  - FS: Fisher链偏好性检验得分
  - MQ: 平均比对质量值
  - SOR: 链偏好性得分
  - MQRankSum: 变异与参考reads的比对质量比较得分
  - ReadPosRankSum: 变异在reads中的位置偏好性得分
  - IndelType: INDEL类型(仅INDEL有该字段)
  - IndelLength: INDEL长度(仅INDEL有该字段)
- FORMAT: 样本信息格式标签
- SAMPLE: 样本基因型信息
  - GT: 基因型(0/0:参考型;0/1:杂合;1/1:纯合变异)
  - AD: 参考序列,变异序列的支持reads数
  - DP: 该位点总测序深度
  - GQ: 基因型质量值
  - PL: 基因型似然值


### 4. IS(插入序列)分析
**调用命令**: `ISPipline.sh [bam] [genome] [ltr_bed] [gene_exon_bed] [gtf] [output_dir] [sample]`

**输出目录结构**:
IS/[sample]/
├── [sample].map2LTR.bam         # 比对到LTR区域的reads结果
├── [sample].split_event.txt     # 断裂reads事件记录
├── [sample].insertion_break_point.txt    # 整合位点信息
├── [sample].insertion_break_point.anno.xls # 注释后的整合位点
├── [sample].insertion_break_point.check.bam # 整合位点验证数据
└── [sample].insertion.igv.html  # 整合位点的IGV可视化报告

**输出文件格式说明**:

1. `[sample].split_event.txt`: 断裂reads详细信息
```
read_id  sequence  read_number  query_length  cigar  next_ref  next_cigar  ref_name  strand  mapq  ref_start  ref_end  query_start  query_end  chr  pos  ltr_name  ltr_pos  strand
```
列说明:
- read_id: reads ID
- sequence: reads序列
- read_number: read1或read2
- query_length: reads长度
- cigar: 比对的CIGAR字符串
- next_ref: 配对reads比对的参考序列名
- next_cigar: 配对reads的CIGAR字符串
- ref_name: 比对的参考序列名
- strand: 比对链方向(+/-)
- mapq: 比对质量值
- ref_start: 参考基因组起始位置
- ref_end: 参考基因组终止位置
- query_start: reads起始位置
- query_end: reads终止位置
- chr: 染色体编号
- pos: 整合位点位置
- ltr_name: LTR序列名称
- ltr_pos: LTR序列位置
- strand: 整合方向

2. `[sample].insertion_break_point.txt`: 整合位点信息
```
sample  chr  pos  end_pos  ltr_pos  support_reads  all_reads_split_reads  all_map2_pos_reads  ratio  AF  Type  combine_sequence  strand  mapq_info
```
列说明:
- sample: 样本名称
- chr: 染色体编号
- pos: 整合起始位置
- end_pos: 整合终止位置
- ltr_pos: LTR序列位置
- support_reads: 支持整合的reads数量
- all_reads_split_reads: 所有断裂reads数量
- all_map2_pos_reads: 比对到该位置的总reads数量
- ratio: 支持reads占断裂reads的比例
- AF: 支持reads占总reads的比例
- Type: 整合类型(Host@Vector/Vector@Host)
- combine_sequence: 整合位点的序列组合
- strand: 整合方向
- mapq_info: 支持reads的比对质量统计

3. `[sample].insertion_break_point.anno.xls`: 整合位点注释结果
```
chr  start  strand  mapq  support_reads  ratio  AF  Type  sequence  gene  feature  description
```
列说明:
- chr: 染色体编号
- start: 整合位置
- strand: 整合方向
- mapq: 平均比对质量值
- support_reads: 支持reads数量
- ratio: 支持reads比例
- AF: 整合位点频率
- Type: 整合类型
- sequence: 整合序列
- gene: 受影响的基因名称
- feature: 基因组特征(exon/intron/intergenic)
- description: 基因功能描述


### 5. SV(结构变异)分析
**调用命令**: `SV.sh [bam] [genome] [output_dir] [sample]`

**输出目录结构**:
SV/[sample]/
├── [sample]_delly.bcf           # Delly工具输出的原始结构变异结果
├── [sample]_delly.vcf           # 转换为VCF格式的结构变异结果
├── [sample]_delly_filter.vcf    # 经过质量过滤的高可信结构变异
└── [sample]_sv.anno.xlsx        # 结构变异的基因组注释结果

**输出文件格式说明**:

1. `[sample]_delly_filter.vcf`: 过滤后的结构变异结果
```
#CHROM  POS     ID  REF  ALT   QUAL  FILTER  INFO
```
列说明:
- CHROM: 染色体编号
- POS: 变异起始位置
- ID: 变异ID
- REF: 参考基因组序列
- ALT: 变异序列
- QUAL: 变异质量值
- FILTER: 过滤标记(PASS表示通过)
- INFO: 包含以下信息
  - SVTYPE: 变异类型(DEL/DUP/INV/TRA)
  - END: 变异终止位置
  - SVLEN: 变异长度
  - PE: 支持的paired-end reads数量
  - SR: 支持的split-reads数量

2. `[sample]_sv.anno.xlsx`: 结构变异注释结果
```
Chr  Start  End  SV_Type  SV_ID  Gene  Location  Distance
```
列说明:
- Chr: 染色体编号
- Start: 变异起始位置
- End: 变异终止位置
- SV_Type: 变异类型(DEL/DUP/INV/TRA)
- SV_ID: 变异ID
- Gene: 受影响的基因名称
- Location: 相对基因的位置(within_gene/upstream/downstream/overlap_gene)
- Distance: 与基因的距离(bp)

### 6. CNV(拷贝数变异)分析
**调用命令**: `CNV.sh [bam] [genome] [gene_bed] [cnvkit_bed] [output_dir] [sample]`

**输出目录结构**:
CNV/[sample]/[sample]_cnvkit_output/
├── [sample].deduped.cns         # CNVkit分析的原始拷贝数信号
├── [sample].call.cns            # 分析后的拷贝数变异检测结果
├── [sample].segments.seg        # 基因组分段的CNV结果
└── [sample].segments.anno.txt   # CNV区域的基因组注释结果

**输出文件格式说明**:

1. `[sample].call.cns`: CNV检测结果
```
chromosome  start  end  gene  log2  cn  depth  probes  weight
```
列说明:
- chromosome: 染色体编号
- start: 区域起始位置
- end: 区域终止位置
- gene: 基因名称
- log2: log2拷贝数比值
- cn: 预测的拷贝数
- depth: 平均测序深度
- probes: 区域内探针数量
- weight: 区域权重

2. `[sample].segments.anno.txt`: CNV区间注释结果
```
Chr  Start  End  CNV_Type  Log2_Ratio  Num_Probes  Gene  Gene_Start  Gene_End  Location
```
列说明:
- Chr: 染色体编号
- Start: 区域起始位置
- End: 区域终止位置
- CNV_Type: 变异类型(DUP/DEL)
- Log2_Ratio: log2拷贝数比值
- Num_Probes: 区域内探针数量
- Gene: 受影响的基因名称
- Gene_Start: 基因起始位置
- Gene_End: 基因终止位置
- Location: 相对基因的位置

### 7. RefCNV(参考CNV)分析
**调用命令**: `reference.cnv.sh [d4_file] [ref_bed] [target_bed] [output_dir] [sample]`

**输出目录结构**:

RefCNV/[sample]/
├── [sample].cnv.txt             # 基于参考区域的CNV分析结果
├── [sample].gene_average.txt    # 每个基因的平均拷贝数统计
└── [sample].stats.txt           # CNV分析的统计信息，包含区域数量等

**输出文件格式说明**:

1. `[sample].cnv.txt`: 基于参考区域的CNV分析结果
```
gene  index  start  end  gc  target_corrected_CN  target_raw_CN  reference_CN
```
列说明:
- gene: 基因名称
- index: 窗口编号
- start: 窗口起始位置
- end: 窗口终止位置
- gc: GC含量
- target_corrected_CN: GC校正后的拷贝数
- target_raw_CN: 原始拷贝数
- reference_CN: 参考区域的拷贝数

2. `[sample].gene_average.txt`: 基因平均拷贝数统计
```
Gene  Average_CN
```
列说明:
- Gene: 基因名称
- Average_CN: 基因区域的平均拷贝数