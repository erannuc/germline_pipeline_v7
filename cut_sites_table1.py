import argparse
from pysam import VariantFile
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
    parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
    parser.add_argument('-workdir', help='output directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')
    parser.add_argument('-chr', help='first chunk to sbatch', default=0)
    parser.add_argument('-start', help='last chunk to sbatch', type=int)
    parser.add_argument('-end', help='last chunk to sbatch', type=int)
    parser.add_argument('-nsamples', help='overall number of samples', type=int)

    args = parser.parse_args()
    chunk_str = chunk_string = '_'.join([args.chr, str(args.start), str(args.end)])

    vcf_in = VariantFile(f'{args.workdir}/{chunk_str}_formal.vcf.bgz')  # auto-detect input format

    affected_regions = {}

    for rec in vcf_in:
        chr = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]
        # print(chr, pos, ref, alt)
        # now consider to cases according to mutation type
        # the regions stored are in bed like inclusive-exclusive [) format
        # and are defined as ALL positions that if methylation motif START there it
        # will be affected by the variant
        if len(ref) == len(alt): # substitution
            affected_regions[(pos-3, pos +1)] = (chr, pos, ref, alt)
        if len(ref) > len(alt):  # deletion
            affected_regions[(pos - 2, pos + 1 + len(ref) - len(alt))] = (chr, pos, ref, alt)
        if len(ref) < len(alt):  # insertion
            affected_regions[(pos - 2, pos + 1)] = (chr, pos, ref, alt)

    with open(f'{args.workdir}/{chunk_str}_affected_regions.bed', 'w') as ofh:
        for i in affected_regions:
            name = affected_regions[i][0] + '_' + str(affected_regions[i][1]) + '_' + affected_regions[i][2] + '_' + affected_regions[i][3]
            # name = '_'.join(affected_regions[i])
            print('\t'.join([chr, str(i[0]), str(i[1]), name, '0', '+']), file=ofh)

    # create bed file of only cut sites
    with open(f'{args.workdir}/{chunk_str}.bed', 'w') as bfh:
        new_start = args.start
        if new_start != 0:
            new_start -= 100
        print('\t'.join([args.chr, str(new_start), str(args.end + 100)]), file=bfh)
    os.system(f'bedtools intersect -a hg38.cuts.tabs.firstbase.bed -b {args.workdir}/{chunk_str}.bed > {args.workdir}/{chunk_str}_cuts_firstbase.bed')

    os.system(f'bedtools intersect -wo -a {args.workdir}/{chunk_str}_cuts_firstbase.bed -b {args.workdir}/{chunk_str}_affected_regions.bed > {args.workdir}/{chunk_str}_intersect_affected.bed')

