import argparse
from pysam import VariantFile
import pandas as pd
from bits import Bits


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
parser.add_argument('-workdir', help='working directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')
parser.add_argument('-out' , help='output tsv file file name', default='cut_sites_GV_matrix.tsv')


args = parser.parse_args()

cut_sites = pd.read_csv('hg38.cuts.tabs.bed', sep='\t', header=None)

# read samples 

vcf_in = VariantFile(f'{args.workdir}/all_formal_intersected.vcf.gz')

header = str(vcf_in.header.samples.header)
header_samples_line_offset =header.index('#CHROM')
header_samples_line = header[header_samples_line_offset:]
samples = header_samples_line.rstrip().split('\t')[9:]
nsamples = len(samples)


cut_sites_vectors = {}
sites2index = {}
for i in cut_sites.index:
    cut_sites_vectors[i] = Bits(nsamples)
    site = str(cut_sites.at[i, 0]) + '_' + str(cut_sites.at[i, 1])
    sites2index[site] = i


# read from combined bed file information for a cpra2sites dict 

cpra2sites = {}
with open(f'{args.workdir}/all_intersect_affected_sorted.bed') as bfh:
    for line in bfh:
        data = line.rstrip().split('\t')
        cut_site = data[0] + '_' + data[1]
        cpra = data[7]
        if cpra not in cpra2sites:
            cpra2sites[cpra] = [cut_site]
        else:
            cpra2sites[cpra].append(cut_site)
        print(cut_site, cpra)


# Go over combined vcf and read sample-wise information to fill the bit vectors


vcf_in = VariantFile(f'{args.workdir}/all_formal_intersected.vcf.gz')
for rec in vcf_in:
    cpra = '_'.join([rec.chrom, str(rec.pos), rec.ref, rec.alts[0]])
    print(cpra)
    gts = [s['GT'] for s in rec.samples.values()]
    for i in range(len(gts)):
        if gts[i] == (0, 1) or gts[i] == (1, 1): 
            # turn on item i in the vector of tis cutsite:
            if cpra in cpra2sites:
                for s in cpra2sites[cpra]:
                    si = sites2index[s]
                    cut_sites_vectors[si].on(i)


# go over the cutsites and print matrix

with open(args.out, 'w') as ofh:
    # print header 
    print('\t' + '\t'.join(samples), file=ofh)
    for i in cut_sites.index:
        cut_site = str(cut_sites.at[i, 0]) + '_' + str(cut_sites.at[i, 1] + 1) + '_' + str(cut_sites.at[i, 2]) + '_' + str(cut_sites.at[i, 3]) 
        print(cut_site + '\t' + '\t'.join([str(cut_sites_vectors[i].get(s)) for s in range(nsamples)]), file=ofh)
