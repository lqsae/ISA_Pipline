import re
import collections
import argparse
import pysam

MapRecord = collections.namedtuple('MapRecord', ['reference_name', 'strand', 'map_q', 'reference_start', 
                                                 'reference_end', 'query_start', 'query_end'])

def filter_from_cigar(cigar:str):
    flag = False
    h = cigar.count('H')
    s = cigar.count('S')
    if h + s ==0:
        flag = True
    if h +s >1:
        flag = True
    return flag


def get_breakpoint(all_MapRecords, ltr_list):
    # 判断断点
    sorted_all_MapRecords = sorted(all_MapRecords, key=lambda x: int(x.query_start))
    record1 = sorted_all_MapRecords[0]
    record2 = sorted_all_MapRecords[1]
    strand = '+'
    if record1.strand != record2.strand:
        strand = '-'
    chr_name = record1.reference_name
    # 如果第一段比对到 载体
    if chr_name in ltr_list:
        genome_record = record2
        ltr_record = record1
        ltr_strand = ltr_record.strand
        ltr_name = ltr_record.reference_name
        genome_strand = genome_record.strand
        genome_chr_name = genome_record.reference_name
        if genome_strand == '-':
            break_point = genome_record.reference_end
        else:
            break_point = genome_record.reference_start
        if ltr_strand == '-':
            ltr_break_point = ltr_record.reference_start
        else:
            ltr_break_point = ltr_record.reference_end
    # 如果第一段比对到 宿主
    else:
        genome_record = record1
        ltr_record = record2
        ltr_strand = ltr_record.strand
        ltr_name = ltr_record.reference_name
        genome_strand = genome_record.strand
        genome_chr_name = genome_record.reference_name
        if genome_strand == '-':
            break_point = genome_record.reference_start
        else:
            break_point = genome_record.reference_end
        if ltr_strand == '-':
            ltr_break_point = ltr_record.reference_end
        else:
            ltr_break_point = ltr_record.reference_start
    return [genome_chr_name, break_point, ltr_name, ltr_break_point, strand]


def is_reads_map2_ref_ltr(all_MapRecords, ltr_list):
    #判断reads是否比对到ltr和宿主基因组
    # reads 是否拆分为两段
    if len(all_MapRecords) ==2:
        chr1 = all_MapRecords[0].reference_name
        chr2 = all_MapRecords[1].reference_name
        # 其中一段比对到ltr另一段比对到宿主
        if chr1 in ltr_list and chr2 not  in ltr_list:
            return True
        elif chr1 not in ltr_list and chr2 in ltr_list:
            return True

    return False


def get_optimal_reference_query_pos(read, cigar_string):
    '''get the end and start of the query and start'''
    # 转换为1-based坐标
    reference_start = read.reference_start +1
    reference_end = read.reference_end 
    # 转换为1-based坐标
    query_start = read.query_alignment_start +1
    query_end = read.query_alignment_end
    map_q = int(read.mapping_quality)
    reference_name = read.reference_name
    if read.is_reverse:
        strand = '-'
    else:
        strand = '+'
    query_start, query_end = get_read_pos(cigar_string, strand)
    optimal_info_list = [reference_name, strand, map_q, reference_start, reference_end, query_start, query_end]
    optimal_record = MapRecord._make(optimal_info_list)
    return optimal_record


def get_read_pos(cigar, strand):
    '''get the end and start of query and referencd in SA tag'''
    sum_M = 0
    sum_D = 0
    sum_I = 0
    all_cigar_value =  re.findall('(\d+)([A-Z]+)', cigar)
    if strand == '-':
        all_cigar_value = all_cigar_value[::-1]
    for value, cigar_s in all_cigar_value:
        value = int(value)
        if cigar_s == 'M':
            sum_M += value
        if cigar == 'I':
            sum_I += value
        elif cigar == 'D':
            sum_D += value    
    bias_query = sum_M + sum_I

    if all_cigar_value[0][1] == 'S' or all_cigar_value[0][1] == 'H':
        query_start = int(all_cigar_value[0][0]) +1
    else:
        query_start = 1
    query_end = query_start + bias_query
    return query_start, query_end


def get_query_reference_pos(cigar, reference_start, strand):
    '''get the end and start of query and referencd in SA tag'''
    sum_M = 0
    sum_D = 0
    sum_I = 0
    all_cigar_value =  re.findall('(\d+)([A-Z]+)', cigar)
    if strand == '-':
        all_cigar_value = all_cigar_value[::-1]
    for value, cigar_s in all_cigar_value:
        value = int(value)
        if cigar_s == 'M':
            sum_M += value
        if cigar == 'I':
            sum_I += value
        elif cigar == 'D':
            sum_D += value    
    bias_reference = sum_D + sum_M
    bias_query = sum_M + sum_I
    reference_end = reference_start + bias_reference
    if all_cigar_value[0][1] == 'S' or all_cigar_value[0][1] == 'H':
        query_start = int(all_cigar_value[0][0]) +1
    else:
        query_start = 1
    query_end = query_start + bias_query
    optimal_pos = [reference_start, reference_end, query_start, query_end]
    optimal_pos = [str(i) for i in optimal_pos]
    return optimal_pos


def get_Sub_optimal(sa_tag):
    '''get the SA_tag cigar information'''
    tags = sa_tag.strip(';').split(';')
    sub_optimal_list = []
    for j in tags:
        chr_name, pos, strand, cigar, map_q, *others = j.split(',')
        list_pos = get_query_reference_pos(cigar, int(pos), strand)
        total_list = [chr_name, strand, map_q] + list_pos
        map_record = MapRecord._make(total_list)
        sub_optimal_list.append(map_record)
    return sub_optimal_list


def Lentiviral_Insertion_Site_Analysis(bam, ltr_list):
    samfile = pysam.AlignmentFile(bam, "r", check_sq=False)
    for record in samfile:
        # if record.is_duplicate:
        #     continue
        if record.has_tag('SA'):
            cigar_string = record.cigarstring
            sequence=record.query_sequence
            query_length = record.query_length
            mapq = record.mapping_quality

            all_map_records = []
            SA = record.get_tag('SA')
            if record.has_tag('MC'):
                next_segment_cigar = record.get_tag('MC')
            else:
                next_segment_cigar = '-'
            next_reference_name= record.next_reference_name
            read_name = record.query_name
            is_read1 = record.is_read1
            read = 'read2'
            if is_read1:
                read = 'read1'
            optimal_map_record = get_optimal_reference_query_pos(record, cigar_string)
            read_info = [read_name, sequence, read, query_length, cigar_string,  next_reference_name, next_segment_cigar]
            all_map_records.append(optimal_map_record)
            sub_optimal_list = get_Sub_optimal(SA)
            for Sub_optimal_map_record in sub_optimal_list:
                all_map_records.append(Sub_optimal_map_record)
            if is_reads_map2_ref_ltr(all_map_records, ltr_list):
                break_point = get_breakpoint(all_map_records, ltr_list)
                all_info = []
                for data in all_map_records:
                    all_info.extend(data)
                print('\t'.join(map(str, read_info  + all_info + break_point)))


def mian():
    parser = argparse.ArgumentParser(description='Lentiviral_Insertion_Site_Analysis')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-l', '--ltr', help='ltr bed file', required=True)
    args = parser.parse_args()
    with open(args.ltr, 'r') as f:
        all_lines = f.readlines()
        ltr_list = [i.strip().split('\t')[0] for i in all_lines]
        Lentiviral_Insertion_Site_Analysis(args.bam, ltr_list)


if __name__ == '__main__':
    mian()
