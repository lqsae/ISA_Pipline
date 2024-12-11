import sys
import os
import collections
import re
import statistics


import pysam


dict_chr_map = {'chr1':1, "chr2":2, 
               "chr3":3, "chr4":4,
               "chr5":5, "chr6":6,
               "chr7":7, "chr8":8,
               "chr9":9, "chr10":10,
               "chr11":11, "chr12":12,
               "chr13":13, "chr14":14,
               "chr15":15, "chr16":16,
               "chr17":17, "chr18":18,
               "chr19":19, "chr20":20,
               "chr21":21, "chr22":22,
               "chrX":23, "chrY":24,
               "chrM":25}


def stat_map(all_map_list):
    mapq_count = collections.Counter(all_map_list)
    all_mapq_str_list = [f'{mapq}:{count}' for mapq, count in mapq_count.items()]
    return '|'.join(all_mapq_str_list) if len(all_mapq_str_list) <= 15 else ' '


def get_soft_sequce(cigar, sequence):
    all_cigar_value =  re.findall('(\d+)([A-Z]+)', cigar) 
    vector_sequence = ''
    host_sequence = ''
    host = 1
    if 'S' in cigar:
        if len(all_cigar_value) == 2:
            if all_cigar_value[0][1] == 'S':
                soft_clip = int(all_cigar_value[0][0])
                vector_sequence = sequence[:soft_clip]
                host_sequence = sequence[soft_clip:]
                host = 0
            else:
                soft_clip = int(all_cigar_value[1][0])
                vector_sequence = sequence[-soft_clip:]
                host_sequence = sequence[:-soft_clip]
                host = 1
    return vector_sequence, host_sequence, host


def get_break_point(sorted_all_point, samfile, all_reads_split_reads, sample, check_bam):
    all_data_list = []
    out_bam = pysam.AlignmentFile(check_bam, 'wb', template=samfile)
    for point in sorted_all_point:
        support_reads = len(point[1])
        ratio = support_reads / all_reads_split_reads
        chr_name, pos = point[0]
        all_map2_pos_reads = 0
        ltr_pos = [(i[0], i[1]) for i in point[1]]
        strands = [i[4] for i in point[1]]
        strand = max(strands, key=strands.count)
        all_read_id_number = set([(i[2], i[3]) for i in point[1]])
        ltr_pos = get_ltr_pos(ltr_pos)
        ltr_pos_str = f'{ltr_pos[0]}:{ltr_pos[1]}'
        all_host_vector_reads = []
        for i  in samfile.fetch(chr_name, int(pos)-1, int(pos) +1):
            all_map2_pos_reads += 1
        all_mapq_list = []
        for i  in samfile.fetch(chr_name, int(pos)-5, int(pos) +5):
            out_bam.write(i)
            read_id = i.query_name
            cigar = i.cigarstring
            sequence = i.query_sequence
            is_read1 = i.is_read1
            mapq = i.mapping_quality
            read_number = 'read2'
            if is_read1:
                read_number = 'read1'
            read_id_number =  (read_id, read_number)
            if read_id_number in all_read_id_number:
                all_mapq_list.append(mapq)
                vector_sequence, host_sequence, host = get_soft_sequce(cigar, sequence)
                all_host_vector_reads.append((host_sequence, vector_sequence, host, cigar, sequence))
        
        mapq_mean = round(statistics.mean(all_mapq_list),2)
        mapq_str = stat_map(all_mapq_list)
        mapq_str = f'{mapq_str};{mapq_mean}'
        
        host_sequence, vector_sequence, host, cigar, sequence = max(all_host_vector_reads, key=lambda x:len(x[1]))
        if host == 1:
            types = 'Host@Vector'
            combine_sequce =  f'{host_sequence}@{vector_sequence}'
        else:
            combine_sequce =  f'{vector_sequence}@{host_sequence}'
            types = 'Vector@Host'
        if all_map2_pos_reads ==0:
            AF = 0
        else:
            AF = support_reads / all_map2_pos_reads
        data = [sample, chr_name, pos, int(pos) +1, ltr_pos_str, support_reads, 
                all_reads_split_reads, all_map2_pos_reads, ratio, AF , types, combine_sequce, strand, mapq_str]
        all_data_list.append(data)
    return all_data_list


def get_break_point_from_split_reads(infile, sample, bam_file, check_bam):
    samfile = pysam.AlignmentFile(bam_file, "r")
    all_point = collections.defaultdict(list)
    all_reads_split_reads = 0
    with open(infile, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            read_id = line_split[0]
            read_number = line_split[2]
            chrom = line_split[-5]
            point = line_split[-4]
            ltr_pos = line_split[-2]
            ltr_chr = line_split[-3]
            strand = line_split[-1]

            all_point[(chrom, point)].append((ltr_chr, ltr_pos, read_id, read_number, strand))
            all_reads_split_reads += 1
    sorted_all_point = sorted(all_point.items(), key=lambda x: len(x[1]), reverse=True)
    all_data_list = get_break_point(sorted_all_point, samfile, all_reads_split_reads, sample, check_bam)
    for data in all_data_list:
        yield data
    

def get_ltr_pos(all_ltr_pos):
    all_ltr_support_reads = collections.Counter(all_ltr_pos)
    ltr_pos =  max(all_ltr_support_reads, key=lambda x: all_ltr_support_reads[x])
    return ltr_pos


def main():
    infile = sys.argv[1]
    sample = sys.argv[2]
    out_dir = sys.argv[3]
    bam_file = sys.argv[4]
    out_file = os.path.join(out_dir, sample + '.insertion_break_point.txt')
    out_file_raw = os.path.join(out_dir, sample + '.insertion_break_point.raw.txt')
    check_bam = os.path.join(out_dir, sample + '.insertion_break_point.check.bam')
    out_w = open(out_file, 'w')
    out_raw = open(out_file_raw, 'w')
    head = ['sample', 'chr', 'pos', 'end_pos', 'ltr_pos', 'support_reads', 
            'all_reads_split_reads', 'all_map2_pos_reads', 'ratio', 'AF', 'Type', 'combine_sequence', 'strand']
    print('\t'.join(head), file=out_w)
    print('\t'.join(head), file=out_raw)
    all_break_point_data = get_break_point_from_split_reads(infile, sample, bam_file, check_bam)
    sorted_data = sorted(all_break_point_data, key=lambda x: (dict_chr_map[x[1]], int(x[2])))
    for data in sorted_data:
        AF =  float(data[-5])
        support_reads = float(data[5])
        print('\t'.join([str(i) for i in data]), file=out_raw)
        if AF >= 0.05 and support_reads >=5:
            print('\t'.join([str(i) for i in data]), file=out_w)
    out_w.close()
    out_raw.close()


if __name__ == '__main__':
    main()
