
import os
from Bio import SeqIO
import numpy as np
import argparse
from memory_profiler import profile
from bits import Bits
import gzip
from scipy.stats import binomtest

global tiny
tiny= 1e-321
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}



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


def chunk_analysis(chr, start, end, control_samples, case_samples, chunk_seq, workdir):


    nucs = {}
    indels = {}
    chunk_string = '_'.join([str(chr), str(start), str(end)])

    num_alleles_dict = {}

    # Get gnomeAD positions in this chunk

    variants_file = f'{workdir}{chunk_string}_stat_v7.vcf.bgz'
    reduced_variants_file = f'{workdir}{chunk_string}_stat_reduced_v7.vcf'

    # Go over file and filter variants based on the results of the stat tests

    # read each line and check the results of teh stat tests.
    with gzip.open(variants_file, 'rt') as ifh, open(reduced_variants_file, 'w') as ofh:
        header = ifh.readline().strip()
        print(header, file=ofh)
        header_fields = header.split('\t')
        pvalue_case_twoside_index = header_fields.index('pvalue_case_twoside')
        pvalue_control_twoside_index = header_fields.index('pvalue_control_twoside')
        pvalue_case_control_binomial_index = header_fields.index('pvalue_case_control_twoside_binomial')
        pvalue_case_control_prop_index = header_fields.index('pvalue_case_control_twoside_prop')
        af_index = header_fields.index('gnomad_af')
        alt_alleles_control_index = header_fields.index('alt_alleles_control')
        tot_alleles_control_index = header_fields.index('tot_alleles_control')
        alt_alleles_case_index = header_fields.index('alt_alleles_case')
        tot_alleles_case_index = header_fields.index('tot_alleles_case')

        control_coverage_indices = get_coverage_indices(control_samples, header_fields)
        case_coverage_indices = get_coverage_indices(case_samples, header_fields)

        total_variants = 0
        pass_variants = 0
        case_population_test_fail = 0
        case_control_test_fail = 0
        coverage_control_test_fail = 0
        coverage_case_test_fail = 0
        af_order_test_fail = 0

        for line in ifh:
            line_data = line.strip().split('\t')
            # print(line_data[pvalue_case_twoside_index], line_data[pvalue_case_control_twoside_index])
            control_mean_coverage = mean([int(line_data[i]) for i in control_coverage_indices])
            case_mean_coverage = mean([int(line_data[i]) for i in case_coverage_indices])
            population_af = float(line_data[af_index])
            control_af = int(line_data[alt_alleles_control_index]) / (int(line_data[tot_alleles_control_index]) + tiny)
            case_af = int(line_data[alt_alleles_case_index]) / (int(line_data[tot_alleles_case_index]) + tiny)
            # print(control_mean_coverage, case_mean_coverage, population_af, control_af, case_af)
            pass_filters = True
            total_variants += 1
            if float(line_data[pvalue_case_twoside_index]) > 1e-10:
                case_population_test_fail += 1
                pass_filters = False
            if float(line_data[pvalue_case_control_prop_index]) > 1e-3:
                case_control_test_fail += 1
                pass_filters = False
            if control_mean_coverage < 10:
                coverage_control_test_fail += 1
                pass_filters = False
            if case_mean_coverage < 10:
                coverage_case_test_fail += 1
                pass_filters = False
            if not pass_freq_order(population_af, control_af, case_af):
                af_order_test_fail += 1
                pass_filters = False

            # print line to the output file if passed all filters
            if pass_filters:
                pass_variants += 1
                print(line, end='', file=ofh)

    print("\t".join(["chunk", "total variants", "case_population_test_fail", "case_control_test_fail", "coverage_control_test_fail",
                     "coverage_case_test_fail", "af_order_test_fail", "pass variants"]))
    print("\t".join([chunk_string, str(total_variants), str(case_population_test_fail), str(case_control_test_fail), str(coverage_control_test_fail),
                     str(coverage_case_test_fail), str(af_order_test_fail), str(pass_variants)]))
    '''
    print("case_population_test_fail\t" + str(case_population_test_fail))
    print("case_control_test_fail\t" + str(case_control_test_fail))
    print("coverage_control_test_fail\t" + str(coverage_control_test_fail))
    print("coverage_case_test_fail\t" + str(coverage_case_test_fail))
    print("af_order_test_fail\t" + str(af_order_test_fail))
    print("pass variants\t" + str(pass_variants))
    '''
    # bgzip and tabix

    # os.system(' '.join(['bgzip', '-f', positions_variants_stat_file]))
    # os.rename(positions_variants_stat_file + '.gz', positions_variants_stat_file + '.bgz')
    # os.system(' '.join(['tabix -pvcf', '-fC', positions_variants_stat_file + '.bgz']))

def get_coverage_indices(samples, header_fields):
    return [i for i in range(len(header_fields)) if header_fields[i].endswith('coverage') and header_fields[i][:-9] in samples]

def mean(list):
    return sum(list)/len(list)

def pass_freq_order(population_af, control_af, case_af):
    # note the change in test in v7, the last condition was added
    if case_af > population_af and population_af > control_af and abs(case_af - population_af) > 2 * abs(control_af - population_af):
    # if case_af > population_af and case_af > control_af:
        return True
    if case_af < population_af and population_af < control_af and abs(case_af - population_af) > 2 * abs(control_af - population_af):
        return True
    if case_af > control_af and control_af >= population_af:
        return True
    if case_af < control_af and control_af <= population_af:
       return True
    return False

def pass_freq_order_old(population_af, control_af, case_af):
    if case_af > population_af and case_af > control_af:
        return True
    if case_af < population_af and case_af < control_af:
        return True
    return False

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-chr', help='chromosome', required=False)
    parser.add_argument('-start', help='start position for analysis - zero based', type=int)
    parser.add_argument('-end', help='end position for analysis - zero based non inclusive', type=int)
    parser.add_argument('-case', help='file with list of case samples', default='/data/users/erane/germline/case_samples_v7.txt')
    parser.add_argument('-control', help='file with list of control samples', default='/data/users/erane/germline/control_samples_v7.txt')
    parser.add_argument('-pools_case', help='file with list of pool of case samples', default='/data/users/erane/germline/case_pool_samples_v7.txt')
    parser.add_argument('-pools_control', help='file with list of pool of control samples', default='/data/users/erane/germline/control_pool_samples_v7.txt')
    parser.add_argument('-workdir', help='output directory', default='/data/users/erane/germline/variants_filtered_v7/')
    parser.add_argument('-sample_start', help='index of first sample to analyze', type=int)
    parser.add_argument('-sample_end', help='index of last sample to analyze', type=int)

    global args
    args = parser.parse_args()

    genome = '/data/db/hg38.fa'
    genome_dict = read_genome(genome)
    chr_sizes = read_chr_sizes(genome)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}

    control, case, control_pools, case_pools = read_samples(args.control, args.case, args.pools_control, args.pools_case)

    samples = list(control.keys()) + list(case.keys()) + list(control_pools.keys()) + list(case_pools.keys())
    control_samples = list(control.keys())
    case_samples = list(case.keys())

    chunks_seq = str(genome_dict[args.chr].seq)[args.start:args.end].upper()

    nucs = chunk_analysis(args.chr, args.start, args.end, control_samples, case_samples, chunks_seq, args.workdir)
