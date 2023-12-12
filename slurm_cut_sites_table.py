import pandas as pd
import argparse
import os
import stat

def create_bash(c, outdir, datadir):
    chr, start, end = [str(i) for i in c]
    chunk_string = '_'.join([chr, start, end])
    sh_file = f'{outdir}/{chunk_string}.sh'
    with open(sh_file, 'w') as ofh:
        print('#!/bin/sh\n', file=ofh)
        print(f'/data/users/erane/germline/venv_germline/bin/python {os.getcwd()}/convert2vcf_v7_before_qc.py -i {datadir}/{chunk_string}_stat_v7.vcf.bgz ' \
              f'-control {os.getcwd()}/control_samples_v7.txt -case {os.getcwd()}/case_samples_v7.txt -o {outdir}/{chunk_string}_formal.vcf -mode full', file=ofh)
    os.chmod(sh_file, stat.S_IRWXU)
    return sh_file

def slurm_jobs(c, outdir, datadir):
    chunk_string = '_'.join([str(i) for i in c])
    if not os.path.exists(f'{outdir}/{chunk_string}_formal.vcf.bgz') or not os.path.exists(f'{outdir}/{chunk_string}_formal.vcf.bgz.csi'):
        bash_file = create_bash(c, outdir, datadir)
        print(f'sbatch {bash_file}')
        os.system(f'sbatch -c 2 --nice=5000 -J vcf{chunk_string} {bash_file}')


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

def create_sample_chunks(nsamples, chunksize):
    sample_chunks = []
    total = 0
    while total < nsamples:
        if total + chunksize < nsamples:
            end = total + chunksize
        else:
            end = nsamples
        sample_chunks.append((total, end))
        total += chunksize
    return sample_chunks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
    parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
    parser.add_argument('-outdir', help='output directory. full path', default='/data/users/erane/germline/variants_full_formal_vcf_v7')
    parser.add_argument('-datadir', help='directory witg variant file. full path', default='/data/users/erane/germline/variants_filtered_v7')
    parser.add_argument('-cstart', help='first chunk to sbatch', type=int, default=0)
    parser.add_argument('-cend', help='last chunk to sbatch', type=int)
    parser.add_argument('-nsamples', help='overall number of samples', type=int)


    global args
    args = parser.parse_args()

    chr_sizes = read_chr_sizes(args.g)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}
    all_chunks = create_chunks(chr_sizes, args.chunk)

    if not args.cend:
        args.cend = len(all_chunks)

    for c in all_chunks[args.cstart:args.cend]:
        slurm_jobs(c, args.outdir, args.datadir)
