import sys
import collections


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


def anno(infile, gene_des_dict):
    dict_anno = collections.defaultdict(list)
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                line_split = line.strip().split('\t')
                chr_= line_split[0]
                start = line_split[1]
                end = line_split[2]
                anno = line_split[-2]
                dict_anno[(chr_, start, end)].append(line_split)
    for key, value_list in dict_anno.items():
        use_list = [i for i in value_list if i[-2] != '.']
        tag = '-'
        if len(use_list) > 0:
            tag_list = [i[-2].split(',')[1] for i in use_list]
            gene = [i[-2].split(',')[0] for i in use_list][0]
            gene_des = gene_des_dict.get(gene, '-').replace('Gene', '-')
            if 'exon' in tag_list:
                tag = 'exon'
            else:
                tag = 'intron'
            data = use_list[0]
            strand = data[11]
            map_q = data[12]
        else:
            gene = '-'
            tag = 'intergenic'
            data = value_list[0]
            map_q = data[12]
            strand = data[11]
            gene_des = gene_des_dict.get(gene, '-')
        print(data[0],data[1], strand, map_q,   *data[3:11],   gene, tag, gene_des,  sep='\t')


if __name__ == '__main__':
    gtf_file = sys.argv[1]
    anno_bed  = sys.argv[2]
    gene_des = anno_gene_des(gtf_file)
    anno(anno_bed, gene_des)
