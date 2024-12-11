#!/usr/bin/env python3
import sys
import argparse
from cyvcf2 import VCF
from pybedtools import BedTool

def get_location(sv_start, sv_end, gene_start, gene_end, upstream_dist=2000, downstream_dist=2000):
    """
    确定SV相对于基因的位置关系
    """
    # SV完全在基因内
    if sv_start >= gene_start and sv_end <= gene_end:
        return "within_gene", 0
    
    # SV在基因上游
    if sv_end < gene_start:
        distance = gene_start - sv_end
        if distance <= upstream_dist:
            return "upstream_2k", distance
        return "upstream", distance
    
    # SV在基因下游
    if sv_start > gene_end:
        distance = sv_start - gene_end
        if distance <= downstream_dist:
            return "downstream_2k", distance
        return "downstream", distance
    
    # SV与基因部分重叠
    return "overlap_gene", 0

def sv_to_bed(vcf_file):
    """将VCF文件转换为BedTool对象"""
    vcf = VCF(vcf_file)
    sv_list = []
    
    for variant in vcf:
        # 获取SV的基本信息
        chrom = variant.CHROM
        start = variant.POS            
        sv_type = variant.INFO.get('SVTYPE')
        end = variant.INFO.get('END', start)
        if sv_type not in ['DEL', 'DUP', 'INV', 'INS']:
            continue
        
        # 创建BED格式条目
        sv_list.append([
            chrom,
            start - 1,  # BED格式是0-based
            end,
            sv_type,
            variant.ID if variant.ID else '.',
            '.'
        ])
    
    return BedTool(sv_list)

def parse_gene_info(hit):
    """解析基因信息，适配4列gene.bed文件"""
    try:
        # gene.bed格式: chr start end genename
        gene_name = hit[9]  # closest结果中基因名在第10列
        gene_start = int(hit[7])  # 基因起始位置在第8列
        gene_end = int(hit[8])    # 基因终止位置在第9列
        return gene_name, gene_start, gene_end
    except (IndexError, ValueError):
        return '.', 0, 0

def main():
    parser = argparse.ArgumentParser(description='使用cyvcf2和pybedtools注释SV')
    parser.add_argument('-v', '--vcf', required=True, help='输入的VCF文件')
    parser.add_argument('-b', '--bed', required=True, help='基因注释BED文件(4列: chr,start,end,genename)')
    parser.add_argument('-o', '--output', required=True, help='输出文件')
    parser.add_argument('-u', '--upstream', type=int, default=2000, help='上游区域距离阈值')
    parser.add_argument('-d', '--downstream', type=int, default=2000, help='下游区域距离阈值')
    args = parser.parse_args()

    try:
        # 将VCF转换为BED格式
        sv_bed = sv_to_bed(args.vcf)
        
        # 读取基因注释文件
        gene_bed = BedTool(args.bed)
        
        # 使用closest找到最近的基因
        closest = sv_bed.closest(gene_bed, d=True, t="first")
        
        # 输出注释结果
        with open(args.output, 'w') as out:
            # 写入表头
            out.write('#Chr\tStart\tEnd\tSV_Type\tSV_ID\tGene\tLocation\tDistance\n')
            
            # 处理每个SV的注释结果
            for hit in closest:
                chr_name = hit[0]
                sv_start = int(hit[1]) + 1  # 转回1-based
                sv_end = int(hit[2])
                sv_type = hit[3]
                sv_id = hit[4]
                
                # 获取基因信息
                gene_name, gene_start, gene_end = parse_gene_info(hit)
                
                # 如果没有找到最近的基因或距离为-1
                if hit[-1] == '-1' or gene_name == '.':
                    out.write(f'{chr_name}\t{sv_start}\t{sv_end}\t{sv_type}\t{sv_id}\t'
                            f'NA\tintergenic\tNA\n')
                    continue
                
                # 确定位置关系
                location, distance = get_location(
                    sv_start, sv_end,
                    gene_start, gene_end,
                    args.upstream, args.downstream
                )
                
                # 输出结果
                out.write(f'{chr_name}\t{sv_start}\t{sv_end}\t{sv_type}\t{sv_id}\t'
                         f'{gene_name}\t{location}\t{distance}\n')
                
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
