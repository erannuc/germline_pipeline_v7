from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import argparse
import subprocess
from Bio import SeqIO
import io
import multiprocessing
from multiprocessing import set_start_method
import json
import os
import struct

from bits import Bits

global tiny
tiny= 0.0001
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}
global dir_map

def get_dir_map(db_file):
    db_data = pd.read_csv(db_file, index_col=0)
    samples2path = {}
    exp2dir = {}
    for i in db_data.index:
        path = str(db_data.at[i, 'Path'])
        print(path.split('\\'))
        samples2path[i] = path
        path_split = path.split('\\')
        experiment, *rest = i.split('.')
        exp2dir[experiment] = path_split[3].lower()
    return samples2path, exp2dir

def read_genome(genome):
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    return genome_dict

def read_chr_sizes(genome):
    chr_sizes = {}
    with open(genome + '.fai', 'r') as gfh:
        for line in gfh:
            chr, size, *rest = line.split('\t')
            chr_sizes[chr] = int(size)
    return chr_sizes

def get_genome_seq(genome, chr, pos, strand):
    params = ['samtools', 'faidx']
    if strand == '-':
        params += ['--reverse-complement']
    params += [genome, chr + ':' + str(pos - 29) + '-' + str(pos + 30)]
    proc = subprocess.Popen(params, stdout=subprocess.PIPE)
    proc.wait()
    with io.TextIOWrapper(proc.stdout, encoding="utf-8") as stdout:
        next(stdout)  # skip header
        seq = next(stdout).rstrip()
    return seq


def create_chunks(chr_sizes, chunk_size, outdir):
    # create chunks of size chunk_size, but restricted to each chromosome, last chunk will hold the rest and will
    # typically be of different size
    chunks = []   # array of tuples of chunk interval. The interval coordinates are 0 based [..) python style

    # afh = open(f'{outdir}/pool_variants_sbatches.sh', 'w')

    for c in chr_sizes:
        pos = 0
        while pos < chr_sizes[c]:
            end = pos + chunk_size
            if end > chr_sizes[c]:
                end = chr_sizes[c]

            '''
            # create bash script for sbatch
            with open(f'{outdir}/pool_variants_{c}_{str(pos)}_{str(end)}.sh', 'w') as ofh:
                print(f'#!/bin/sh\n\n/data/users/erane/germline/venv_germline/bin/python /data/users/erane/germline/germline_with_indel_coverage_eco_aws.py -chr ' \
                       f'{c} -start {str(pos)} -end {str(end)}', file=ofh)
            print(f'sbatch -c 2 --nice=5000 {outdir}/pool_variants_{c}_{str(pos)}_{str(end)}.sh', file=afh)
            '''

            chunks.append((c, pos, end))
            pos += chunk_size

    # afh.close()

    return chunks

def get_done_chunks(outdir):
    ok_files = [f for f in os.listdir(outdir) if isfile(join(outdir, f)) and f.endswith(".ok")]
    done_chunks = []
    for f in ok_files:
        # make sure vcf file also exists
        tsv_file = f.replace('.ok', '.tsv')
        if os.path.exists(join(outdir, tsv_file)):
            file_tokens = f.replace('.', '_').split('_')
            done_chunks.append((file_tokens[0], int(file_tokens[1]), int(file_tokens[2])))
    return done_chunks


def chunk_analysis(fargs):

    chunk, chunk_seq, control, case, outdir = fargs
    chr = chunk[0]
    start = int(chunk[1])
    end = int(chunk[2])
    chunk_size = end - start
    chunk_string = '_'.join([str(chr), str(start), str(end)])
    print(chunk_string)
    # foreach chunk go over all nucleotides in the chunk and do the following:
    snv_data = {}
    indel_data = {}
    genotypes = {}
    positions_to_report = {}  # positions with changed genotype
    for c in list(control.items()) + list(case.items()):
        snv_data, indel_data[c[0]] = sample_analysis(c[0], chunk)

        # indel data is stored in the memory
        for p in range(start, end):
            if p in positions_to_report:
                continue
            if p in indel_data[c[0]]:
                indel_data_p = indel_data[c[0]][p]
            else:
                indel_data_p = {}
            estimated_alleles = estimate_secondary_alleles(c[1], snv_data[p - start], chunk_seq[p - start], indel_data_p)
            if estimated_alleles > 0:
                positions_to_report[p] = 1


    # start here to filter variants. return only filtered variants

    outfile = f'{outdir}/{chunk_string}.tsv'
    with open(outfile, 'w') as ofh:
        for p in positions_to_report:
            print(chr + "\t" + str(p), file=ofh)

    okfile = f'{outdir}/{chunk_string}.ok'
    with open(okfile, 'w') as ofh:
        pass
    return 0


def position_coverage(counts, indel_data):
    return sum(counts) + sum(indel_data.values())


def sum_non_ref(position, start, refrence_base, samples, snv_data, indel_data):
    snr = 0
    for s in samples:
        for i in range(0, 4):
            if base_map[i] != refrence_base:
                snr += snv_data[s][position-start][i]
        if position in indel_data[s]:
            for i in indel_data[s][position]:
                snr += indel_data[s][position][i]
    return snr


def format_print_header(samples):
    line_array = ['#chromosome', 'position', 'refbase']
    for s in samples:
        for nuc in ['T', 'C', 'G', 'A']:
            line_array.append('s' + ':' + nuc)
        line_array.append(s + ':indels')
        line_array.append(s + ':genotype')
        line_array.append(s + ':coverage')
    return '\t'.join(line_array)

def format_print_line(chromosome, position, refbase, start_chunk, samples, snv_data, indel_data, pgenotypes):
    line_array = [str(chromosome), str(position + 1), refbase]
    for s in samples:
        coverage = 0
        for nuc in range(4):
            line_array.append(str(snv_data[s][position-start_chunk][nuc]))
            coverage += snv_data[s][position - start_chunk][nuc]
        if position in indel_data[s]:
            for i in indel_data[s][position]:
                coverage += indel_data[s][position][i]
            indel_str = ';'.join([str(i) + ':' + str(indel_data[s][position][i]) for i in indel_data[s][position]])
        else:
            indel_str = ''
        line_array.append(indel_str)
        line_array.append(str(pgenotypes[s][0]) + '|' + str(pgenotypes[s][1]))
        line_array.append(str(coverage))

    return '\t'.join(line_array)

def estimate_secondary_alleles(n, counts, reference_base, indel_data):
    # estimate the number of alleles in the pool by dividing the rare allele counts by coverage/(2*nsamples)

    indel_coverage = sum(indel_data.values())
    sumreads = sum(counts) + indel_coverage + tiny
    if sumreads < 20:
        return 0
    single_allele_coverage = sumreads / (2 * n)
    allele_types = {}
    for i in range(4):
        if counts[i] < single_allele_coverage / 2 :
            allele_types[base_map[i]] = 0
        elif counts[i] < single_allele_coverage :
            allele_types[base_map[i]] = 1
        else:
            allele_types[base_map[i]] = int(counts[i] / single_allele_coverage)

    for i in indel_data:
        if indel_data[i] < single_allele_coverage / 2 :
            allele_types[i] = 0
        elif indel_data[i] < single_allele_coverage :
            allele_types[i] = 1
        else:
            allele_types[i] = int(indel_data[i] / single_allele_coverage)

    # sum how many non-dominant alleles exits here
    secondary_alleles = sorted(allele_types.values()) # Sort the values

    # In this version also check positions with homozygous non reference allele and return 2
    if reference_base != 'N' and allele_types[reference_base] == 0:
        return 2
    return sum(secondary_alleles[:-1])


def sample_analysis(sample, chunk):
    # get sample and chunk interval and collect stat for all nucleotides within the chunk
    print("in sample analysis " + sample + ' chunk: ' + str(chunk[0]) + ' ' + str(chunk[1]) + ' ' + str(chunk[2]))
    chromosome = chunk[0]
    start = int(chunk[1])
    end = int(chunk[2])

    # Read from binary file the nucleotide counts
    experiment, subject, source, treatment = sample.split('.')

    nuc_file = f'/data/users/erane/germline/features/{sample}/{sample}.Chr{chromosome}.Consensus.Nucs'
    ins_file = f'/data/users/erane/germline/features/{sample}/{sample}.Chr{chromosome}.Consensus.Ins'
    dels_file = f'/data/users/erane/germline/features/{sample}/{sample}.Chr{chromosome}.Consensus.Del'

    with open(ins_file, 'r') as ifh, open(dels_file, 'r') as dfh:
        # first read the insertions and deletions
        # for each indel event, seek Nucs coverage in the right place, return the frequency of the event
        # bfh.seek(chunk[1]*16)
        indels = {}
        coverage = {}
        for line in ifh:
            line = line.rstrip()
            if line.startswith('@'):
                position = int(line[1:])
                if position < start:
                    continue
                if position >= end:
                    break
                indels[position] = {}
            else:
                if position < start:
                    continue
                seq, count = line.split(',')
                indels[position][seq] = int(count)
                # print(position, seq, count)
        for line in dfh:
            line = line.rstrip()
            if line.startswith('@'):
                position = int(line[1:])
                if position < start:
                    continue
                if position >= end + 150:
                    break
                if position not in indels:
                    indels[position] = {}
            else:
                if position < start:
                    continue
                length, count = line.split(',')
                for i in range(0, int(length)):
                    if position + i in indels:
                        indels[position + i]['_'.join([str(-abs(i)), length])] = int(count)
                    else:
                        indels[position + i] = {'_'.join([str(-abs(i)), length]): int(count)}
        counter = 0

    snv_data = np.fromfile(nuc_file, dtype='i', count=4*(end-start), offset=start*16).reshape(end-start, 4)

    return snv_data, indels


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', help='file with control pool samples names', default='/data/users/erane/germline/control_pool_samples_v7.txt', required=False)
    parser.add_argument('-t', help='file with tumor pool samples names', default='/data/users/erane/germline/case_pool_samples_v7.txt', required=False)
    parser.add_argument('-chr', help='restrict analysis to specific chromosome', required=False)
    parser.add_argument('-start', help='start position for analysis - zero based', type=int)
    parser.add_argument('-end', help='end position for analysis - zero based non inclusive. chr, start and end overhide chunk', type=int)
    parser.add_argument('-g', help='genome fasta file. index should be in the same place', default='/data/db/hg38.fa')
    parser.add_argument('-chunk', help='chunk size for analysis of one process', type=int, default=10000000)
    parser.add_argument('-n', help='number of processes', type=int, default=2)
    parser.add_argument('-outdir', help='output directory', default='/data/users/erane/germline/variants_pools_v7')
    parser.add_argument('-samples_db', help='mapping file', default='/data/users/erane/germline/db.csv')
    parser.add_argument('-auto', help='analysis on autosomes only', action='store_true')


    global args
    args = parser.parse_args()

    # samples_map, dir_map = get_dir_map(args.samples_db)

    control_samples = {}
    if args.c is not None:
        with open(args.c, 'r') as ifh:
            for line in ifh:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                sample, n = line.split('\t')
                control_samples[sample] = int(n)
    tumor_samples = {}
    if args.t is not None:
        with open(args.t, 'r') as ifh:
            for line in ifh:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                sample, n = line.split('\t')
                tumor_samples[sample] = int(n)

    chr_sizes = read_chr_sizes(args.g)
    if args.auto:
        chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}

    if args.chr is not None:
        chr_sizes = {c: chr_sizes[c] for c in chr_sizes if c in args.chr}
    print(chr_sizes)

    all_chunks = create_chunks(chr_sizes, args.chunk, args.outdir)
    print(len(all_chunks))

    done_chunks = get_done_chunks(args.outdir)   # done chunks are those chunks with available .vcf file and .ok file in the output dir
    all_chunks = [i for i in all_chunks if i not in done_chunks] 
    genome_dict = read_genome(args.g)

    if args.start is not None and args.end is not None and args.chr is not None:
        # do analysis only for this chunk
        if (args.chr, args.start, args.end) not in done_chunks:
            chunk_seq = str(genome_dict[args.chr].seq)[args.start:args.end].upper()
            status = chunk_analysis([(args.chr, args.start, args.end), chunk_seq, control_samples, tumor_samples, args.outdir])
        exit(0)
    exit(0)

    pool_obj = multiprocessing.Pool(processes=args.n)
    
    args_lists = []

    print(all_chunks)

    for ch in all_chunks[:1]:
        chr, start, end = ch
        chunks_seq = str(genome_dict[chr].seq)[start:end].upper()
        args_lists.append([ch, chunks_seq, control_samples, tumor_samples, args.outdir])
    out = pool_obj.map(chunk_analysis, args_lists)

    for process in range(len(out)):
        if out[process] != 0:
            print("None zero exit for ", args_lists[process])
