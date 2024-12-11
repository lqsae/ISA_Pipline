import sys


def get_gene_list(infile):
    with open(infile) as f:
        for line in f:
            gene = line.strip()
            yield gene


def anno_gene_des(gtf_file):
    gene_des = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene':
                if 'description' in fields[8]:
                    gene_name = fields[8].split(';')[0].split('"')[1]
                    product = fields[8].split(';')[5].split('"')[1]
                    gene_des[gene_name] = product
                else:
                    gene_name = fields[8].split(';')[0].split('"')[1]
                    product = '-'
                    gene_des[gene_name] = product
    return gene_des


def main():
    gtf_file = sys.argv[1]
    gene_list_file = sys.argv[2]
    gene_des = anno_gene_des(gtf_file)
    for gene in get_gene_list(gene_list_file):
        print(gene, gene_des.get(gene, '-'), sep='\t')


if __name__ == '__main__':
    main()