import sys
import os
from collections import defaultdict
from collections import Counter
_usage = ''' For stat off target reads ratio'''


def flag(list):
    list=set(list)
    if "0" in list:
        type = "0"
    elif "1" in list and "0" not in list:
        type = "1"
    else:
        type = "2"
    return type


def identify_type(input):
    input = int(input)
    if int(input) == 0:
        type = "0"
    elif int(input) <= 150 and int(input) != 0 and int(input) != (-1):
        type = "1"
    else:
        type = "2"

    return type


def find_dupReads(infile, id_col=3):
    read_records = set()
    dup_record=set()
    with open(infile) as fo:
        for line in fo:
            line = line.strip().split("\t")
            readID=line[id_col]
            if readID in read_records:
                dup_record.add(readID)
            else:
                read_records.add(readID)
    del read_records
    return dup_record


def process_file(infile, dup_records, id_col=3):
    type_reads = {"0": 0, "1":0, "2":0}
    dup_reads = defaultdict(list)
    with open(infile) as fo:
        for line in fo:
            line = line.strip().split("\t")
            if line[id_col] not in dup_records:
                type_reads[identify_type(line[-1])] = type_reads[identify_type(line[-1])] +1
            else:
                dup_reads[line[id_col]].append(identify_type((line[-1])))

    return type_reads, dup_reads


def main():
    if len(sys.argv) != 3:
        print("python {} <stat file> <sample>".format(sys.argv[0]))
        print("<stat file>: bamtobed and bedtools closest -d -t first -a *.bam.bed -b probe.bed -wao > stat.file")
        exit(1)
    infile = sys.argv[1]
    sample = sys.argv[2]
    dup_records = find_dupReads(infile)
    type_reads, dup_reads = process_file(infile, dup_records)
    final_dict = [flag(i) for i in dup_reads.values()]
    count = Counter(final_dict)
    del final_dict
    header=["Sample", "TotalReads", "OverlapReads", "OverlapRatio(%)","Dist<=150 Reads","Dist<=150 Ratio(%)", "Dist>150 Reads", "Dist>150 Ratio(%)"]
    # total reads
    TotalReads = count["0"]+type_reads["0"]+count["1"]+type_reads["1"]+count["2"]+type_reads["2"]
    # overlap , in 150 ratio and off ratio
    outratio = float((count["2"]+type_reads["2"])/TotalReads * 100)
    overlap_ratio = float((count["0"]+type_reads["0"])/TotalReads * 100)
    in150_ratio = float((count["1"]+type_reads["1"])/TotalReads * 100)
    content = [sample, TotalReads, count["0"]+type_reads["0"], "%.2f"%overlap_ratio, \
               count["1"]+type_reads["1"], "%.2f"%in150_ratio, \
               count["2"]+type_reads["2"], "%.2f"%outratio]
    print("\t".join(header))
    print("\t".join([str(i) for i in content]))


if __name__ == '__main__':
    main()
