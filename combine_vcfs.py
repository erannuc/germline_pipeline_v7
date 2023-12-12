import argparse
from pysam import VariantFile


def read_chr_sizes(genome):
    chr_sizes = {}
    with open(genome + '.fai', 'r') as gfh:
        for line in gfh:
            chr, size, *rest = line.split('\t')
            chr_sizes[chr] = int(size)
    return chr_sizes


def create_chunks(chr_sizes, chunk_size):
    # create chunks of size chunk_size, but restricted to each chromosome, last chunk will hold the rest and will
    # typically be of different size
    chunks = []   # array of tuples of chunk interval. The interval coordinates are 0 based [..) python style

    for c in chr_sizes:
        pos = 0
        while pos < chr_sizes[c]:
            end = pos + chunk_size
            if end > chr_sizes[c]:
                end = chr_sizes[c]

            chunks.append((c, pos, end))
            pos += chunk_size
    return chunks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
    parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
    parser.add_argument('-workdir', help='output directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')

    args = parser.parse_args()

    chr_sizes = read_chr_sizes(args.g)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}
    all_chunks = create_chunks(chr_sizes, args.chunk)

    vcf_out = VariantFile(f'{args.workdir}/all_formal_intersected.vcf.gz', 'w')
    header_printed = False
    for chunk in all_chunks:
        chunk_str = '_'.join([str(i) for i in chunk])
        # print once header
        vcf_in = VariantFile(f'{args.workdir}/{chunk_str}_formal_intersected.vcf.gz')
        print(dir(vcf_in))
        if not header_printed:
            vcf_out.header = vcf_in.header
            header_printed = True
        for rec in vcf_in:
            vcf_out.write(rec)
