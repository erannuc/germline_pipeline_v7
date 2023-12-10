import os
from Bio import SeqIO
import numpy as np
import argparse
from memory_profiler import profile
from bits import Bits

global tiny
tiny= 0.0001
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}

# for a given chunk of 1M positions - open all files and extract genotype.

def read_positions_in_chunk(file):
    # read positions with numpy
    data = np.fromfile(file, dtype=np.int64, sep='\t').reshape(-1, 2)
    return data

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

def read_samples(control_file, case_file, pools_control_file, pools_case_file):
    control_samples = {}
    with open(control_file, 'r') as ifh:
        for line in ifh:
            if line.startswith("#"):
                continue
            sample = line.rstrip()
            control_samples[sample] = 1
    case_samples = {}
    with open(case_file, 'r') as ifh:
        for line in ifh:
            if line.startswith("#"):
                continue
            sample = line.rstrip()
            case_samples[sample] = 1
    pools_control_samples = {}
    with open(pools_control_file, 'r') as ifh:
        for line in ifh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            sample, n = line.split('\t')
            pools_control_samples[sample] = int(n)
    pools_case_samples = {}
    with open(pools_case_file, 'r') as ifh:
        for line in ifh:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            sample, n = line.split('\t')
            pools_case_samples[sample] = int(n)
    return control_samples, case_samples, pools_control_samples, pools_case_samples


def chunk_analysis(positions, chr, start, end, samples, chunk_seq, sample_start, sample_end):
    nucs = {}
    indels = {}
    chunk_string = '_'.join([str(chr), str(start), str(end)])
    print(samples)

    # for each sample, read only the chunk from the Consensus nucs.
    # subtrat chunk start from the positions and extract only information for these sites

    positions_in_chunk = positions[:,1] - start

    genotypes = {}
    for s in samples:
        print(s)
        genotypes[s] = Bits(2 * len(positions_in_chunk))  # each position need two bits for genotype
        nucs[s], indels[s] = sample_analysis(s, (chr, start, end), positions_in_chunk)
        for p in range(len(positions_in_chunk)):
            pos = positions_in_chunk[p]
            if pos + start in indels[s]:
                pindels = indels[s][pos + start]
            else:
                pindels = {}
            genotype = assign_genotype(pos + start, nucs[s][p], pindels, chunk_seq[pos])
            if genotype[0] == 1:
                genotypes[s].on((p) * 2)
            if genotype[1] == 1:
                genotypes[s].on((p) * 2 + 1)

    # print to file
    outfile = f'{args.outdir}/{chunk_string}_{str(sample_start)}_{str(sample_end)}.vcf'
    with open(outfile, 'w') as ofh:
        headerline = format_print_header(samples)
        print(headerline, file=ofh)
        for p in range(len(positions_in_chunk)):
            pos = positions_in_chunk[p]
            outline = format_print_line(chr, p, pos, chunk_seq[pos], start, samples, nucs, indels, genotypes)
            print(outline, file=ofh)

    os.system(' '.join(['bgzip', '-f', outfile]))
    os.rename(outfile + '.gz', outfile + '.bgz')
    os.system(' '.join(['tabix -pvcf', '-fC', outfile + '.bgz']))


def sample_analysis(sample, chunk, requested_positions):
    # get sample and chunk interval and collect stat for all nucleotides within the chunk
    print(sample)
    print(requested_positions)
    print("in sample analysis " + sample + ' chunk: ' + str(chunk[0]) + ' ' + str(chunk[1]) + ' ' + str(chunk[2]))
    chromosome = chunk[0]
    start = int(chunk[1])
    end = int(chunk[2])

    # Read from binary file the nucleotide counts
    experiment, subject, source, treatment = sample.split('.')[:4]

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

    # snv_data = np.fromfile(nuc_file, dtype='i', count=4*(end-start), offset=start*16).reshape(end-start, 4)
    snv_data = np.take(np.fromfile(nuc_file, dtype='i', count=4*(end-start), offset=start*16).reshape(end-start, 4), requested_positions, axis=0)

    # leave only indel requested positions
    requested_positions_set  = set(requested_positions + start)
    for i in list(indels.keys()):
        if i not in requested_positions_set:
            del indels[i]

    return snv_data, indels


def assign_genotype(pos, counts, indel_data, reference_base):
    indel_coverage = sum(indel_data.values())
    sumreads = sum(counts) + indel_coverage
    if sumreads==0:
        sumreads = tiny
    nalleles = [0, 0, 0, 0]
    genotype = [0, 0]
    for i in range(4):
        if counts[i] / sumreads > 0.8:
            nalleles[i] = 2
            if base_map[i] != reference_base:
                genotype = [1, 1]
            else:
                genotype = [0, 0]
        elif counts[i] / sumreads >= 0.2:
            nalleles[i] = 1
            if base_map[i] != reference_base:
                genotype = [0, 1]
        else:
            nalleles[i] = 0
    for i in indel_data:
        # check if one of the indels can be called:
        if indel_data[i] / sumreads > 0.8:
            genotype = [1, 1]
        elif indel_data[i] / sumreads >= 0.2:
            genotype = [0, 1]
    return genotype

def format_print_header(samples):
    line_array = ['#chromosome', 'position', 'refbase']
    for s in samples:
        for nuc in ['T', 'C', 'G', 'A']:
            line_array.append(s + ':' + nuc)
        line_array.append(s + ':indels')
        line_array.append(s + ':genotype')
        line_array.append(s + ':coverage')
    return '\t'.join(line_array)


def format_print_line(chromosome, p, position, refbase, start_chunk, samples, snv_data, indel_data, genotypes):
    line_array = [str(chromosome), str(position + start_chunk + 1), refbase]
    for s in samples:
        coverage = 0
        for nuc in range(4):
            line_array.append(str(snv_data[s][p][nuc]))
            coverage += snv_data[s][p][nuc]
        if position + start_chunk in indel_data[s]:
            for i in indel_data[s][position + start_chunk]:
                coverage += indel_data[s][position + start_chunk][i]
            indel_str = ';'.join([str(i) + ':' + str(indel_data[s][position + start_chunk][i]) for i in indel_data[s][position + start_chunk]])
        else:
            indel_str = ''
        line_array.append(indel_str)
        line_array.append(str(genotypes[s].get(p * 2)) + '|' + str(genotypes[s].get(p * 2 + 1)))
        line_array.append(str(coverage))

    return '\t'.join(line_array)


def analysis(fargs):
    job = fargs
    nucs_file = r'W:\\Exp40W\Exp40W.UBC_CON015.WB.digested\Genomic\Nucs\Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Nucs'
    ins_file = r'W:\\Exp40W\Exp40W.UBC_CON015.WB.digested\Genomic\Ins\Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Ins'
    del_file = r'W:\\Exp40W\Exp40W.UBC_CON015.WB.digested\Genomic\Del\Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Del'
    readtext(ins_file, del_file, job)
    readfile(nucs_file, job)
    readnumpy(nucs_file, job)
    return 0


def readtext(insfile, delfile, job):
    print(shared_x, shared_y)
    with open(insfile, 'r') as ifh, open(delfile, 'r') as dfh:
        counter = 0
        for line in ifh:
            counter += 1
            if counter%10000 == 0:
                print(job , " i:", counter)
        counter = 0
        for line in dfh:
            counter += 1
            if counter%10000 == 0:
                print(job , " d:", counter)


def readfile(file, job):
    with open(file, 'br') as bfh:
        print("*")
        for p in range(0, 10000000):
            TB = bfh.read(4)
            T = int.from_bytes(TB, "little")
            CB = bfh.read(4)
            C = int.from_bytes(CB, "little")
            GB = bfh.read(4)
            G = int.from_bytes(GB, "little")
            AB = bfh.read(4)
            A = int.from_bytes(AB, "little")

            if p%1000000 == 0:
                print(job , ":", p)
        print("**")

def readnumpy(file, job):
    print("*")
    data = np.fromfile(file, dtype='i', count=800000000).reshape((200000000, 4))


    # print(data)
    print("**")
    print(data[2][3])


def check_read_file_numpy(file, chrsize, n):
    print(file, chrsize)
    print("*")
    positions = [random.randrange(chrsize) for i in range(n)]
    positions.sort()
    print("**")
    data = np.fromfile(file, dtype='i', count=4*chrsize).reshape((chrsize, 4))
    print("***")
    printnparray = np.array(positions).reshape((n, 1))
    print("****")
    # np.savetxt("testnp.tsv",  printnparray, fmt='%d', delimiter='\t')
    print("*****")
    np.savetxt("testnp.tsv",  np.concatenate((printnparray, np.take(data, positions, axis=0)), axis=1), fmt='%d', delimiter='\t')
    print("******")

    with open("test.tsv", 'w') as ofh:
        # print(np.take(data, positions, axis=0), file=ofh)
        for p in positions:
            # print('\t'.join([str(p), str(data[p][0]), str(data[p][1]), str(data[p][2]), str(data[p][3])]), file=ofh)
            print(str(p) + np.array2string(data[p], separator='\t'), file=ofh)
    print("*******")

def search_with_seek(file, chrsize, n):
    positions = [random.randrange(chrsize) for i in range(n)]
    positions.sort()
    with open(file, 'br') as bfh:
        print("*")
        for p in range(n):
            # bfh.seek(positions[p] * 16)
            TB = bfh.read(4)
            T = int.from_bytes(TB, "little")
            CB = bfh.read(4)
            C = int.from_bytes(CB, "little")
            GB = bfh.read(4)
            G = int.from_bytes(GB, "little")
            AB = bfh.read(4)
            A = int.from_bytes(AB, "little")

            if p%10000 == 0:
                print(p)
        print("**")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-chr', help='chromosome', required=False)
    parser.add_argument('-start', help='start position for analysis - zero based', type=int)
    parser.add_argument('-end', help='end position for analysis - zero based non inclusive', type=int)
    parser.add_argument('-case', help='file with list of case samples', default='/data/users/erane/germline/case_samples_v7.txt')
    parser.add_argument('-control', help='file with list of control samples', default='/data/users/erane/germline/control_samples_v7.txt')
    parser.add_argument('-pools_case', help='file with list of pool of case samples', default='/data/users/erane/germline/case_pool_samples_v7.txt')
    parser.add_argument('-pools_control', help='file with list of pool of control samples', default='/data/users/erane/germline/control_pool_samples_v7.txt')
    parser.add_argument('-posdir', help='input directory with positions to scan', default='/data/users/erane/germline/variants_pools_v7/')
    parser.add_argument('-outdir', help='output directory', default='/data/users/erane/germline/variants_genotypes_v7/')
    parser.add_argument('-sample_start', help='index of first sample to analyze', type=int)
    parser.add_argument('-sample_end', help='index of last sample to analyze', type=int)

    global args
    args = parser.parse_args()

    genome = '/data/db/hg38.fa'
    genome_dict = read_genome(genome)
    chr_sizes = read_chr_sizes(genome)

    control, case, control_pools, case_pools = read_samples(args.case, args.control, args.pools_case, args.pools_control)


    # nucs_file = '/data/users/erane/germline/features/Exp40W.UBC_CON015.WB.digested/Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Nucs'
    # ins_file = '/data/users/erane/germline/features/Exp40W.UBC_CON015.WB.digested/Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Ins'
    # del_file = '/data/users/erane/germline/features/Exp40W.UBC_CON015.WB.digested/Exp40W.UBC_CON015.WB.digested.Chr1.Consensus.Del'
    pos_file =  f'{args.posdir}/{args.chr}_{args.start}_{args.end}_sorted.tsv'

    positions = read_positions_in_chunk(pos_file)

    # read the Consensus.Nucs of the required positions only from each sample

    samples = list(control.keys()) + list(case.keys()) + list(control_pools.keys()) + list(case_pools.keys())
    chunks_seq = str(genome_dict[args.chr].seq)[args.start:args.end].upper()

    nucs = chunk_analysis(positions, args.chr, args.start, args.end, samples[args.sample_start:args.sample_end], chunks_seq, args.sample_start, args.sample_end)
