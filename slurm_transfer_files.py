import pandas as pd
import argparse
import os
import stat

# take files in pseudo vcf format and convert to formal vcf format

def check_local_exists(fd, sample, genome):
    # check that all files exist locally and 'ok' file exists

    # get chromsome sizes:
    chr_sizes = read_chr_sizes(genome)
    print(f'{fd}/{sample}')
    if not os.path.exists(f'{fd}/{sample}'):
        print(f'{sample} is missing')
        return False
    for chr in range(1, 23):
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs'):
            print(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs')
            return False
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs.ok'):
            print(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs.ok')
            return False
        if os.path.getsize(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs') != chr_sizes[str(chr)] * 16:
            print(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Nucs size') 
            return False
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Ins'):
            print(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Ins')
            return False
        if os.path.getsize(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Ins') == 0:
            print(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Ins size')
            return False
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Ins.ok'):
            return False
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Del'):
            return False
        if not os.path.exists(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Del.ok'):
            return False
        if os.path.getsize(f'{fd}/{sample}/{sample}.Chr{str(chr)}.Consensus.Del') == 0:
            return False
    return True


def get_sub_dir(sample, db_file):
    db_data = pd.read_csv(db_file, index_col=0)
    if sample not in db_data.index:
        raise Exception(f'Sorry, sample {sample} can not be found in the db file {db_file}')
    return db_data.at[sample, 'Path'].split('\\')[4].lower()

def create_bash(fd, subdir, sample, c):
    sh_file = f'{fd}/{sample}_{c}_transfer.sh'
    bucket = "nucleix.features"
    experiment, subject, *rest = sample.split('.')
    with open(sh_file, 'w') as ofh:
        print('#!/bin/sh\n', file=ofh)

        com = ' '.join(['aws', 's3', 'cp',
                  f's3://{bucket}/{subdir}/{experiment}/{sample}/Genomic/Nucs/{sample}.Chr{c}.Consensus.Nucs.zip', f'{fd}',
                  '&&', '7z', 'x', f'{fd}/{sample}.Chr{c}.Consensus.Nucs.zip', '-o' + f'{fd}/{sample}',
                  '&&', 'rm', '-f', f'{fd}/{sample}.Chr{c}.Consensus.Nucs.zip', '&&', 'touch',
                  f'{fd}/{sample}/{sample}.Chr{c}.Consensus.Nucs.ok'])
        print(com,  file=ofh)

        com = ' '.join(['aws', 's3', 'cp',
                  f's3://{bucket}/{subdir}/{experiment}/{sample}/Genomic/Ins/{sample}.Chr{c}.Consensus.Ins.zip', f'{fd}',
                  '&&', '7z', 'x', f'{fd}/{sample}.Chr{c}.Consensus.Ins.zip', '-o' + f'{fd}/{sample}',
                  '&&', 'rm', '-f', f'{fd}/{sample}.Chr{c}.Consensus.Ins.zip', '&&', 'touch',
                  f'{fd}/{sample}/{sample}.Chr{c}.Consensus.Ins.ok'])
        print(com,  file=ofh)

        com = ' '.join(['aws', 's3', 'cp',
                  f's3://{bucket}/{subdir}/{experiment}/{sample}/Genomic/Del/{sample}.Chr{c}.Consensus.Del.zip', f'{fd}',
                  '&&', '7z', 'x', f'{fd}/{sample}.Chr{c}.Consensus.Del.zip', '-o' + f'{fd}/{sample}',
                  '&&', 'rm', '-f', f'{fd}/{sample}.Chr{c}.Consensus.Del.zip', '&&', 'touch',
                  f'{fd}/{sample}/{sample}.Chr{c}.Consensus.Del.ok'])
        print(com,  file=ofh)

    os.chmod(sh_file, stat.S_IRWXU)
    return sh_file

# python convert2vcf.py -i 22_40000000_50000000_stat.vcf.bgz  -o 22_40000000_50000000_stat_formal.vcf   -control control_samples_orna.txt -case case_samples_orna.txt

def slurm_jobs(fd, subdir, sample, c):
    bash_file = create_bash(fd, subdir, sample, c)
    print(f'sbatch {bash_file}')
    os.system(f'sbatch -c 2 --nice=10000 -J file_transfer{sample}.{c} {bash_file}')


def read_chr_sizes(genome):
    chr_sizes = {}
    with open(genome + '.fai', 'r') as gfh:
        for line in gfh:
            chr, size, *rest = line.split('\t')
            chr_sizes[chr] = int(size)
    return chr_sizes



if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', help='list of samples to transfer. one per line, no header', default='/data/users/erane/germline/list_of_samples_to_transfer.txt')
    parser.add_argument('-g', help='genome fasta file. index should be in the same place',
                        default='/data/db/hg38.fa')
    parser.add_argument('-b', help='s3 bucket', default='nucleix.features')
    parser.add_argument('-env', help='environment', default='aws')

    parser.add_argument('-fd', help='output directory. full path', default='/data/users/erane/germline/features/')


    global args
    args = parser.parse_args()

    genome_path = '/data/db/hg38.fa'
    chr_sizes = read_chr_sizes(args.g)

    if args.env == 'eran':
        db_file = '/mnt/V/Nucleix Public/Nucleix in-house software/Wgrs/db.csv'
    else:
        db_file = '/data/users/erane/germline/db.csv'

    samples = []
    with open(args.s, 'r') as ifh:
        for i in ifh:
            if i.startswith("#"):
                continue
            sample = i.rstrip().split("\t")[0]
            if check_local_exists(args.fd, sample, genome_path):
                continue
            samples.append(sample)
    print(samples)

    for s in samples:
        subdir = get_sub_dir(s, db_file)
        for c in range(1, 23):
            slurm_jobs(args.fd, subdir, s, c)
