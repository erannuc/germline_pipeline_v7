import argparse
from pysam import VariantFile
import pandas as pd



parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
parser.add_argument('-workdir', help='output directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')

args = parser.parse_args()

cut_sites = pd.read_csv('hg38.cuts.tabs.bed', sep='\t')

