#!/usr/bin/env python3
import sys
import argparse
from pybedtools import BedTool

def get_location(cnv_start, cnv_end, gene_start, gene_end, upstream_dist=2000, downstream_dist=2000):
    """
    确定CNV相对于基因的位置关系
    """
    # CNV完全在基因内
    if cnv_start >= gene_start and cnv_end <= gene_end:
        return "within_gene", 0
    
    # CNV在基因上游
    if cnv_end < gene_start:
        distance = gene_start - cnv_end
        if distance <= upstream_dist:
            return "upstream_2k", distance
        return "upstream", distance
    
    # CNV在基因下游
    if cnv_start > gene_end:
        distance = cnv_start - gene_end
        if distance <= downstream_dist:
            return "downstream_2k", distance
        return "downstream", distance
    
    # CNV与基因部分重叠
    return "overlap_gene", 0

def parse_cnv_file(cnv_file):
    """解析CNVkit的segments文件"""
    cnv_list = []
    with open(cnv_file) as f:
        # 跳过头部
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            chr_name = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            num_probes = int(fields[4])
            log2 = float(fields[5])
            
            # 确定CNV类型
            if log2 >= 0.3:
                cnv_type = 'DUP'
            elif log2 <= -0.3:
                cnv_type = 'DEL'
            else:
                continue  # 跳过正常区域
            
            cnv_list.append([
                chr_name,
                start,
                end,
                cnv_type,
                str(log2),
                str(num_probes)
            ])
    
    return BedTool(cnv_list)

def main():
    parser = argparse.ArgumentParser(description='注释CNV区间')
    parser.add_argument('-i', '--input', required=True, help='CNVkit segments文件')
    parser.add_argument('-b', '--bed', required=True, help='基因注释BED文件(4列: chr,start,end,genename)')
    parser.add_argument('-o', '--output', required=True, help='输出文件')
    args = parser.parse_args()

    try:
        # 解析CNV文件并转换为BED格式
        cnv_bed = parse_cnv_file(args.input)
        gene_bed = BedTool(args.bed)
        
        # 使用intersect找到所有重叠的基因
        intersect = cnv_bed.intersect(gene_bed, wa=True, wb=True)
        
        # 输出注释结果
        with open(args.output, 'w') as out:
            # 写入表头
            header = ['Chr', 'Start', 'End', 'CNV_Type', 'Log2_Ratio', 
                     'Num_Probes', 'Gene', 'Gene_Start', 'Gene_End', 
                     'Location']
            out.write('\t'.join(header) + '\n')
            
            # 处理每个CNV-基因对
            current_cnv = None
            for hit in intersect:
                cnv_chr = hit[0]
                cnv_start = int(hit[1])
                cnv_end = int(hit[2])
                cnv_type = hit[3]
                log2_ratio = hit[4]
                num_probes = hit[5]
                
                # 获取基因信息
                gene_name = hit[9]
                gene_start = int(hit[7])
                gene_end = int(hit[8])
                
                # 确定位置关系
                location, _ = get_location(cnv_start, cnv_end, gene_start, gene_end)
                
                # 输出结果
                result = [
                    cnv_chr,
                    str(cnv_start),
                    str(cnv_end),
                    cnv_type,
                    log2_ratio,
                    num_probes,
                    gene_name,
                    str(gene_start),
                    str(gene_end),
                    location
                ]
                out.write('\t'.join(result) + '\n')
            
            # 处理没有重叠基因的CNV
            no_overlap = cnv_bed.intersect(gene_bed, v=True)
            for hit in no_overlap:
                result = [
                    hit[0],  # chr
                    hit[1],  # start
                    hit[2],  # end
                    hit[3],  # cnv_type
                    hit[4],  # log2_ratio
                    hit[5],  # num_probes
                    'NA',    # gene
                    'NA',    # gene_start
                    'NA',    # gene_end
                    'intergenic'  # location
                ]
                out.write('\t'.join(result) + '\n')
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
