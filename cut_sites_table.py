import argparse

import pysam





if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
    parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
    parser.add_argument('-workdir', help='output directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')
    parser.add_argument('-chr', help='first chunk to sbatch', type=int, default=0)
    parser.add_argument('-start', help='last chunk to sbatch', type=int)
    parser.add_argument('-end', help='last chunk to sbatch', type=int)
    parser.add_argument('-nsamples', help='overall number of samples', type=int)

