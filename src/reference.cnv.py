import sys
import collections

from pyd4 import D4File

import numpy as np

from scipy.signal import savgol_filter
from scipy.stats import binned_statistic


def get_reference_mean_depth(reference_bedfile, d4_file):
    d4file = D4File(d4_file)
    all_depth = 0
    all_length = 0
    with open(reference_bedfile) as f:
        for line in f:
            line_split = line.strip().split('\t')
            chr_name, start, end, gene = line_split
            length = int(end) - int(start)
            all_length += length
            data = d4file[(chr_name, int(start), int(end ))]
            all_depth += data.sum()
    mean_data = all_depth / all_length
    
    return mean_data

def get_reference_win_depth(reference_bedfile, d4_file, mean, wins=100):
    d4file = D4File(d4_file)
    gene_win_bed_dict = collections.defaultdict(list)
    with open(reference_bedfile) as f:
        for line in f:
            line_split = line.strip().split('\t')
            chr_name, start, end, gene = line_split
            gene_name = gene.split('_')[0]
            gene_win_bed_dict[gene_name].append((chr_name, start, end))
    # for gene in gene_win_bed_dict:
    for i in range(wins):
        genes_win = [gene_win_bed_dict[gene][i] for gene in gene_win_bed_dict]
        all_depth = 0
        all_length = 0
        for win in genes_win:
            chr_name, start, end = win
            data = d4file[(chr_name, int(start), int(end ))]
            all_length += int(end) - int(start)
            all_depth += data.sum()
        mean_data = all_depth / all_length *2

        yield  mean_data / mean 



def smooth_data(data, window_length=7, polyorder=1):
    """
    使用 Savitzky-Golay 滤波器对数据进行平滑处理
    """

    return savgol_filter(data, window_length, polyorder)



def gc_correction(depths, gc_contents):
    # 将GC含量分成20个bin
    gc_bins = np.linspace(0, 1, 21)
    
    # 计算每个bin的平均深度
    mean_depths, _, _ = binned_statistic(gc_contents, depths, statistic='mean', bins=gc_bins)
    
    # 计算全局平均深度
    global_mean_depth = np.mean(depths)
    
    # 计算每个bin的校正因子
    correction_factors = global_mean_depth / mean_depths
    
    # 对每个数据点进行校正
    corrected_depths = np.zeros_like(depths)
    for i, (depth, gc) in enumerate(zip(depths, gc_contents)):
        bin_index = np.digitize(gc, gc_bins) - 1
        corrected_depths[i] = depth * correction_factors[bin_index]
    
    return corrected_depths


def get_target_depth(d4_file, bed_file, mean):
    d4file = D4File(d4_file)
    gene_win_mean_data = collections.defaultdict(list)
    gene_bed = collections.defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            line_split = line.strip().split('\t')
            chr_name, start, end, gene, gc = line_split
            gene_name = gene.split('_')[0]
            gene_bed[gene_name].append((chr_name, start, end, float(gc)))
    for gene in gene_bed:
        depths = []
        gc_contents = []
        positions = []  # 存储位置信息
        for win in gene_bed[gene]:
            chr_name, start, end, gc = win
            data = d4file[(chr_name, int(start), int(end))]
            depths.append(data.mean())
            gc_contents.append(gc)
            positions.append((start, end))  # 保存位置信息
        corrected_depths = gc_correction(depths, gc_contents)
        for corrected_depth, raw_depth, gc, pos in zip(corrected_depths, depths, gc_contents, positions):
            start, end = pos
            gene_win_mean_data[gene].append((corrected_depth/mean*2, raw_depth/mean*2, gc, start, end))
    return gene_win_mean_data


def main():
    d4_file = sys.argv[1]
    reference_bed_file = sys.argv[2]
    target_bed_file = sys.argv[3]
    mean = get_reference_mean_depth(reference_bed_file, d4_file)
    reference_data = get_reference_win_depth(reference_bed_file, d4_file, mean)
    target_data = get_target_depth(d4_file, target_bed_file, mean)
    reference_data = list(reference_data)
    
    # 修改输出标题，添加 start 和 end 列
    print('gene\tindex\tstart\tend\tgc\ttarget_corrected_CN\ttarget_raw_CN\treference_CN')
    
    # 需要修改 get_target_depth 函数的返回值，使其包含位置信息
    for gene in target_data:
        index = 0
        for target, reference in zip(target_data[gene], reference_data):
            corrected_depth, raw_depth, gc, start, end = target  # 解包增加的位置信息
            index += 1
            print(gene, index, start, end, gc, np.log2(corrected_depth), np.log2(raw_depth), np.log2(reference), sep='\t')


if __name__ == '__main__':
    main()
