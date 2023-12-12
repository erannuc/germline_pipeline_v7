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

    args = parser.parse_args()
    chunk_str = chunk_string = '_'.join([args.chr, str(args.start), str(args.end)])

    bed_in = f'{args.workdir}/all_intersect_affected_sorted.bed'
    cpra_dict = {}
    print("*")
    with open(bed_in) as bfh:
        # read only variants on this chromosome
        for line in bfh:
            data = line.strip().split('\t')
            chr = data[0]
            if chr != args.chr:
                continue
            cpra_dict[data[7]] = 1
    print("**")

    # Now read the vcf file and create a new one only for the variants which cover cut sites

    vcf_in = VariantFile(f'{args.workdir}/{chunk_str}_formal.vcf.bgz')
    vcf_out = VariantFile(f'{args.workdir}/{chunk_str}_formal_intersected.vcf.gz', 'w', header=vcf_in.header)
    for rec in vcf_in:
        chr = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]
        cpra = '_'.join([chr, str(pos), ref, alt])
        print(cpra)
        if cpra in cpra_dict:
            vcf_out.write(rec)
    print("***")
