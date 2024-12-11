import sys

def gtf2bed(gtf_file):
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            gene_name = fields[8].split(';')[0].split('"')[1]
            tag = fields[2]
            if tag == 'gene':
                print(chrom, start, end, gene_name, sep='\t')


def main():
    gtf_file = sys.argv[1]
    gtf2bed(gtf_file)


if __name__ == '__main__':
    main()