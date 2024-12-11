#!/usr/bin/env python3
import sys
import argparse
from cyvcf2 import VCF
from pybedtools import BedTool

def get_location(var_pos, gene_start, gene_end, upstream_dist=2000, downstream_dist=2000):
    """
    确定变异相对于基因的位置关系
    """
    # 变异在基因内
    if gene_start <= var_pos <= gene_end:
        return "within_gene", 0
    
    # 变异在基因上游
    if var_pos < gene_start:
        distance = gene_start - var_pos
        if distance <= upstream_dist:
            return "upstream_2k", distance
        return "upstream", distance
    
    # 变异在基因下游
    if var_pos > gene_end:
        distance = var_pos - gene_end
        if distance <= downstream_dist:
            return "downstream_2k", distance
        return "downstream", distance

def variants_to_bed(vcf_file):
    """将VCF文件转换为BedTool对象"""
    vcf = VCF(vcf_file)
    variant_list = []
    
    for variant in vcf:
        # 获取变异的基本信息
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = ','.join(variant.ALT)  # 可能���多个替代等位基因
        var_type = 'SNV' if len(ref) == len(alt) == 1 else 'INDEL'
        filter_status = variant.FILTER if variant.FILTER else 'PASS'
        
        # 获取质量相关信息
        qual = variant.QUAL if variant.QUAL else '.'
        dp = variant.format('DP')[0][0] if variant.format('DP') is not None else '.'

        if dp == '.' or  int(dp) < 10:
            continue
        
        # 获取AD和AF信息
        ad = variant.format('AD')
        if ad is not None:
            ref_depth = ad[0][0]  # 参考等位基因深度
            alt_depth = ad[0][1]  # 替代等位基因深度
            af = alt_depth / (ref_depth + alt_depth) if (ref_depth + alt_depth) > 0 else 0
            ad_str = f"{ref_depth},{alt_depth}"
            af_str = f"{af:.3f}"
        else:
            ad_str = "."
            af_str = "."
        
        # 创建BED格式条目 (使用0-based坐标)
        variant_list.append([
            chrom,
            pos - 1,  # 转换为0-based
            pos,
            f"{var_type}|{ref}>{alt}",  # 变异信息
            f"{qual}|{dp}|{filter_status}|{ad_str}|{af_str}",  # 质量和深度信息
            '.'
        ])
    
    return BedTool(variant_list)

def parse_gene_info(hit):
    """解析基因信息"""
    try:
        gene_name = hit[9]  # closest结果中基因名在第10列
        gene_start = int(hit[7])  # 基因起始位置
        gene_end = int(hit[8])    # 基因终止位置
        return gene_name, gene_start, gene_end
    except (IndexError, ValueError):
        return '.', 0, 0

def main():
    parser = argparse.ArgumentParser(description='注释GATK检测的SNV/Indel')
    parser.add_argument('-v', '--vcf', required=True, help='GATK输出的VCF文件')
    parser.add_argument('-b', '--bed', required=True, help='基因注释BED文件(4列: chr,start,end,genename)')
    parser.add_argument('-o', '--output', required=True, help='输出文件')
    parser.add_argument('-u', '--upstream', type=int, default=2000, help='上游区域距离阈值')
    parser.add_argument('-d', '--downstream', type=int, default=2000, help='下游区域距离阈值')
    args = parser.parse_args()

    try:
        # 将VCF转换为BED格式
        var_bed = variants_to_bed(args.vcf)
        
        # 读取基因注释文件
        gene_bed = BedTool(args.bed)
        
        # 使用closest找到最近的基因
        closest = var_bed.closest(gene_bed, d=True, t="first")
        
        # 输出注释结果
        with open(args.output, 'w') as out:
            # 写入表头
            header = ['#Chr', 'Position', 'Variant_Type', 'Ref>Alt', 'Quality', 
                     'Depth', 'Filter', 'AD', 'AF', 'Gene', 'Location', 'Distance']
            out.write('\t'.join(header) + '\n')
            
            # 处理每个变异的注释结果
            for hit in closest:
                chr_name = hit[0]
                pos = int(hit[2])  # 使用1-based位置
                var_info = hit[3].split('|')  # 变异信息
                qual_info = hit[4].split('|')  # 质量信息
                
                var_type = var_info[0]
                ref_alt = var_info[1]
                qual = qual_info[0]
                depth = qual_info[1]
                filter_status = qual_info[2]
                ad = qual_info[3]
                af = qual_info[4]
                
                # 获取基因信息
                gene_name, gene_start, gene_end = parse_gene_info(hit)
                
                # 如果没有找到最近的基因或距离为-1
                if hit[-1] == '-1' or gene_name == '.':
                    out.write(f'{chr_name}\t{pos}\t{var_type}\t{ref_alt}\t{qual}\t'
                            f'{depth}\t{filter_status}\t{ad}\t{af}\t'
                            f'NA\tintergenic\tNA\n')
                    continue
                
                # 确定位置关系
                location, distance = get_location(
                    pos, gene_start, gene_end,
                    args.upstream, args.downstream
                )
                
                # 输出结果
                out.write(f'{chr_name}\t{pos}\t{var_type}\t{ref_alt}\t{qual}\t'
                         f'{depth}\t{filter_status}\t{ad}\t{af}\t'
                         f'{gene_name}\t{location}\t{distance}\n')
                
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
