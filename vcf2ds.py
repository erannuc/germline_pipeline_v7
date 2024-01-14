import os
from Bio import SeqIO
from pysam import VariantFile
import numpy as np
import pandas as pd
import argparse
# from bits import Bits
import gzip

tiny= 1e-100

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="convert from vcf to DS table format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='input standard formal vcf file which include format field for each variant. could be gzipped or not', required=True)
    parser.add_argument('-o', help='output tsv file for DS', required=True)
    parser.add_argument('-control', help='file with control samples', default='control_samples_v7.txt')
    parser.add_argument('-case', help='file with case samples', default='case_samples_v7.txt')

    parser.add_argument('-local', help='local server or aws (default)', action="store_true")


    # this script convert from out pseudo vcf format to standard VCF format compatible with bcftools csq
    # v2 add more info tags added by the filter_variants_freqs_reduce_more.py script and also allows filtering by repetative nature
    # of the region

    global args
    args = parser.parse_args()

    with VariantFile(args.i) as vcf_in, open(args.o, 'w') as ofh:  # auto-detect input format

        header = str(vcf_in.header.samples.header)
        header_samples_line_offset = header.index('#CHROM')
        header_samples_line = header[header_samples_line_offset:]
        samples = header_samples_line.rstrip().split('\t')[9:]
        nsamples = len(samples)

        print('\t'.join(['hg38_chr', 'hg38_pos', 'nucleix_pos', 'ref', 'alt', 'gv', 'annotation',
                         'p_control_gnomad', 'p_case_gnomad', 'p_case_control_binomial', 'p_case_control_prop',
                         'gv_alleles_control', 'total_alleles_control', 'af_control', 'gv_alleles_case',
                         'total_alleles_case', 'af_case', 'odds_ratio', 'gnomad_v3.1_af'] + samples), file=ofh)
        for rec in vcf_in:
            chr = rec.chrom
            pos = rec.pos
            ref = rec.ref
            alt = rec.alts[0]
            op = rec.info.get('OP')[0]
            if rec.info.get('BCSQ') is not None:
                annotation = rec.info.get('BCSQ')[0]
            else:
                annotation = ''
            p_control_gnomad = rec.info.get('PC')[0]
            p_case_gnomad = rec.info.get('PT')[0]
            p_case_control_binomial = rec.info.get('PB')[0]
            p_case_control_proportions = rec.info.get('PR')[0]
            gv = alt
            if rec.info.get('AF_CONT')[0] > 0.5:
                gv = ref
            if gv == alt:
                gv_alleles_control = rec.info.get('AC_CONT')[0]
                af_control = rec.info.get('AF_CONT')[0]
            else:
                gv_alleles_control = rec.info.get('AN_CONT')[0] - rec.info.get('AC_CONT')[0]
                af_control = 1 - rec.info.get('AF_CONT')[0]
            total_alleles_control = rec.info.get('AN_CONT')[0]
            if gv == alt:
                gv_alleles_case = rec.info.get('AC_CASE')[0]
                af_case = rec.info.get('AF_CASE')[0]
            else:
                gv_alleles_case = rec.info.get('AN_CASE')[0] - rec.info.get('AC_CASE')[0]
                af_case = 1 - rec.info.get('AF_CASE')[0]
            total_alleles_case = rec.info.get('AN_CASE')[0]

            odds_ratio = af_case / (af_control + tiny)
            if gv == alt:
                gnomad_af = rec.info.get('GA')[0]
            elif gv == ref:
                gnomad_af = 1 - rec.info.get('GA')[0]

            print(chr, ref, alt, op, annotation)
            gts = [s['GT'] for s in rec.samples.values()]

            genotypes = []
            for i in range(len(gts)):
                if gts[i] == (0, 1):
                    genotypes.append(1)
                elif gts[i] == (1, 1):
                    if gv == alt:
                        genotypes.append(2)
                    elif gv == ref:
                        genotypes.append(0)
                elif gts[i] == (0, 0):
                    if gv == alt:
                        genotypes.append(0)
                    elif gv == ref:
                        genotypes.append(2)
                elif gts[i] == (None, None):
                    genotypes.append(-1)
            genotypes = [str(i) for i in genotypes]
            print('\t'.join([chr, str(pos), str(op), ref, alt, gv, annotation, str(p_control_gnomad), str(p_case_gnomad),
                             str(p_case_control_binomial), str(p_case_control_proportions), str(gv_alleles_control), str(total_alleles_control),
                             str(af_control), str(gv_alleles_case), str(total_alleles_case), str(af_case), str(odds_ratio),
                             str(gnomad_af)] + genotypes), file=ofh)



