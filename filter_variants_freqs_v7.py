import os
from Bio import SeqIO
import numpy as np
import argparse
from memory_profiler import profile
from bits import Bits
import gzip
from scipy.stats import binomtest
from statsmodels.stats.proportion import proportions_ztest
import math

global tiny
tiny= 1e-321
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}


def read_positions(file):
    # read positions with numpy
    data = np.fromfile(file, dtype=np.int64, sep='\t').reshape(-1, 2)
    positions_dict = {}
    for i in range(len(data)):
        if str(data[i][0]) not in positions_dict:
            positions_dict[str(data[i][0])] = {}
        positions_dict[str(data[i][0])][data[i][1]] = set()


    return positions_dict

def extract_af(af_vcf, positions_dict, outtsv):
    with open(af_vcf, 'r') as vfh:
        linecounter = 0
        for line in vfh:
            linecounter += 1
            if line.startswith('#'):
                continue
            chr, pos, id, ref, alt, qual, filter, info, *rest = line.strip().split('\t')
            alts = alt.split(',')
            info_data = info.split(';')
            for d in info_data:
                k, v = d.split('=')
                if k == 'AF':
                    afs = v.split(',')
                    break

            # now for each ref-alt combination - insert all potential affecting variants to positions_dict
            # for alt-ref of the same length - insert only the original position. if the length is different insert also
            # to the positions of deleted bases

            for i in range(len(alts)):
                if len(alts[i]) >= len(ref):
                    # pos is in the 1-based vcf coordinates while the dict is in our zero based coors
                    if int(pos) - 1 in positions_dict[chr]:
                        var_data_string = ":".join([chr, pos, id, ref, alts[i], afs[i]])
                        if var_data_string not in positions_dict[chr][int(pos) - 1]:
                            # print("*", var_data_string)
                            positions_dict[chr][int(pos) - 1].add(var_data_string)
                else:
                    # in case of deletion the position of missing bases themselves are considered (like in the features)
                    # not that there is no '-1' corrections to the coordinates because in the HGVS
                    # the deletion ref position is of the base prior to the actuall missing bases
                    # We consider the real missing positions. THis effect cancels the 1 base difference due to the
                    # difference between one-based and zero-based coordinates
                    for p in range(len(ref) - len(alts[i])):
                        if int(pos) + p  in positions_dict[chr]:
                            var_data_string = ":".join([chr, pos, id, ref, alts[i], afs[i]])
                            # print("**", var_data_string)
                            positions_dict[chr][int(pos) + p].add(var_data_string)

        # create the output file:

        with open(outtsv, 'w') as ofh:
            for i in positions_dict:
                for j in positions_dict[i]:
                    covering_variants = ','.join(positions_dict[i][j])
                    print('\t'.join([i, str(j), covering_variants]), file=ofh)


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

def read_samples(control_file, case_file):
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

    return control_samples, case_samples


def chunk_analysis(chr, start, end, control_samples, case_samples, chunk_seq, workdir, min_allele_count, min_coverage):
    nucs = {}
    indels = {}
    chunk_string = '_'.join([str(chr), str(start), str(end)])

    num_alleles_dict = {}

    # Get gnomeAD positions in this chunk

    positions_file = f'/data/users/erane/germline/variants_pools_v7/{chunk_string}_sorted.tsv'
    positions_dict = read_positions(positions_file)
    # positions_in_chunk = [int(i) - start for i in list(positions_dict[chr].keys())]

    # extract variants in the range from the af gnomeAD vcf file

    os.system(f'bcftools view /data/db/hg38.gnomad.v3.vcf.gz {chr}:{start-100}-{end} > {workdir}/{chunk_string}.af.vcf')

    # now extract from the region af vcf file the frequencies
    extract_af(f'{workdir}/{chunk_string}.af.vcf', positions_dict, f'{workdir}/{chunk_string}.af.tsv')

    positions_af = {}
    with open(f'{workdir}/{chunk_string}.af.tsv', 'r') as ifh:
        for line in ifh:
            afchr, afpos, afdata = line.strip('\n').split('\t')
            if afdata:
                alleles = afdata.split(',')
                nonref_af = 0
                for a in alleles:
                    nonref_af += float(a.split(':')[5])
            else:
                nonref_af = 0.000008
            positions_af[int(afpos)] = nonref_af

    # Now count genotypes
    # For each variant with more than 1 allele get allele frequency and do binomial test for the two populations (control/case) vs gnomad
    # and between themselves
    # Create output file with information on the two tests, the tests ratio, number of individuals and genotypes in case and control

    # read from the combined vcf file

    positions_variants_info_file = f'/data/users/erane/germline/variants_genotypes_v7/{chunk_string}.vcf.bgz'
    positions_variants_stat_file = f'{workdir}/{chunk_string}_stat_v7.vcf'

    with gzip.open(positions_variants_info_file, 'rt') as ifh, open(positions_variants_stat_file, 'w') as ofh:
        header_line = ifh.readline().strip('\n')
        header_fields = header_line.split('\t')
        print('\t'.join(header_fields[:3]) + '\tgnomad_af\talt_alleles_control\ttot_alleles_control\tcar_control\ttot_control\t' +
              'alt_alleles_case\ttot_alleles_case\tcar_case\ttot_case\tpvalue_case_oneside\t' +
              'pvalue_control_twoside\tpvalue_case_twoside\tpvalue_case_control_twoside_binomial\tpvalue_case_control_twoside_prop\t' + '\t'.join(header_fields[3:]), file=ofh)
        control_genotype_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in control_samples]
        case_genotype_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in case_samples]

        linecounter = 0
        for line in ifh:
            linedata = line.strip('\n').split('\t')
            control_genotypes = [linedata[i] for i in control_genotype_indices]
            case_genotypes = [linedata[i] for i in case_genotype_indices]
            sum_alt_control_alleles = 0
            sum_alt_control_carriers = 0
            low_coverage_samples_control = 0
            for i in control_genotype_indices:
                if int(linedata[i+1]) < min_coverage: # coverage test. the coverage is the next index to the genotype
                    linedata[i] = '.|.'
                    low_coverage_samples_control += 1 
                elif linedata[i] == '0|1':
                    sum_alt_control_alleles += 1
                    sum_alt_control_carriers +=1
                elif linedata[i] == '1|1':
                    sum_alt_control_alleles += 2
                    sum_alt_control_carriers += 1

            sum_alt_case_alleles = 0
            sum_alt_case_carriers = 0
            low_coverage_samples_case = 0
            for i in case_genotype_indices:
                if int(linedata[i+1]) < min_coverage: # coverage test. the coverage is the next index to the genotype
                    linedata[i] = '.|.' 
                    low_coverage_samples_case += 1 
                elif linedata[i] == '0|1':
                    sum_alt_case_alleles += 1
                    sum_alt_case_carriers += 1
                elif linedata[i] == '1|1':
                    sum_alt_case_alleles += 2
                    sum_alt_case_carriers += 1

            # perform test and report only positions with more than one varied allele:

            if sum_alt_case_alleles + sum_alt_control_alleles < min_allele_count:
                continue
            # in version 7 NOT filtering out variants with uniformally 100% non-reference allele
            # if sum_alt_case_alleles + sum_alt_control_alleles > 2*(len(control_samples) + len(case_samples)) - min_allele_count:
            #    continue

            # stat binomial test

            pos = int(linedata[1])
            if pos-1 in positions_af:
                af = positions_af[pos-1]  # the coordinates in the tsv af file are zero-based while in positions_variants_info_file are 1-based
            else:
                af = 0.000008
            if af >= 1:
                print(chr, pos, af)
                af = 0.999992  # edge cases of accumulated alleles with overall greater than 1.0 frequency. also need maximal freq not being exactly 1 for the stat test
            if af < 0.5:
                alternative = 'greater'
            else:
                alternative = 'less'
            total_control_alleles = (len(control_samples) - low_coverage_samples_control) * 2
            if total_control_alleles > 0:
                control_binom_test = binomtest(sum_alt_control_alleles, n=total_control_alleles, p=af, alternative='two-sided')
                pval_control_twoside = control_binom_test.pvalue
                if pval_control_twoside == 0:
                    pval_control_twoside = 1e-321
            else:
                pval_control_oneside = 1.0
                pval_control_twoside = 1.0

            total_case_alleles = (len(case_samples) - low_coverage_samples_case) * 2
            if total_case_alleles > 0:
                total_case_alleles = (len(case_samples) - low_coverage_samples_case) * 2
                case_binom_test = binomtest(sum_alt_case_alleles, n=total_case_alleles, p=af, alternative=alternative)
                pval_case_oneside = case_binom_test.pvalue
                if pval_case_oneside==0:
                    pval_case_oneside = 1e-321

                case_binom_test = binomtest(sum_alt_case_alleles, n=total_case_alleles, p=af, alternative='two-sided')
                pval_case_twoside = case_binom_test.pvalue
                if pval_case_twoside==0:
                    pval_case_twoside = 1e-321
            else:
                pval_case_oneside = 1.0
                pval_case_twoside = 1.0


            control_af = sum_alt_control_alleles / (total_control_alleles + tiny)
            detection_freq = 1/(1.5 * total_control_alleles + tiny)

            if control_af < 0.0000000000001:
                control_af = detection_freq
            if control_af > 0.9999999999999:
                control_af = 1 - detection_freq
            if total_case_alleles > 0 and total_control_alleles > 0:
                case_control_binom_test = binomtest(sum_alt_case_alleles, n=total_case_alleles, p=control_af, alternative='two-sided')
                pval_case_control_twoside = case_control_binom_test.pvalue
                if pval_case_control_twoside==0:
                    pval_case_control_twoside = 1e-321
                # add proportaions test for case vs control
                stat_case_control_prop, pval_case_control_prop = proportions_ztest([sum_alt_control_alleles, sum_alt_case_alleles], [total_control_alleles, total_case_alleles])
                if pval_case_control_prop==0:
                    pval_case_control_prop = 1e-321
            else:
                pval_case_control_twoside = 1.0
                pval_case_control_prop = 1.0

            linecounter += 1

            # print details in output table

            print('\t'.join(linedata[:3]) + '\t' + '\t'.join([str(af), str(sum_alt_control_alleles), str(total_control_alleles), str(sum_alt_control_carriers), str(len(control_samples)-low_coverage_samples_control), \
                  str(sum_alt_case_alleles), str(total_case_alleles), str(sum_alt_case_carriers), str(len(case_samples)-low_coverage_samples_case), str(pval_case_oneside), \
                  str(pval_control_twoside), str(pval_case_twoside)]) + '\t' + str(pval_case_control_twoside) + '\t' + str(pval_case_control_prop) + '\t' + '\t'.join(linedata[3:]), file=ofh)

    # bgzip and tabix

    os.system(' '.join(['bgzip', '-f', positions_variants_stat_file]))
    os.rename(positions_variants_stat_file + '.gz', positions_variants_stat_file + '.bgz')
    os.system(' '.join(['tabix -pvcf', '-fC', positions_variants_stat_file + '.bgz']))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-chr', help='chromosome', required=False)
    parser.add_argument('-start', help='start position for analysis - zero based', type=int)
    parser.add_argument('-end', help='end position for analysis - zero based non inclusive', type=int)
    parser.add_argument('-case', help='file with list of case samples', default='/data/users/erane/germline/case_samples_v7.txt')
    parser.add_argument('-control', help='file with list of control samples', default='/data/users/erane/germline/control_samples_v7.txt')
    parser.add_argument('-workdir', help='output directory', default='/data/users/erane/germline/variants_filtered_v7/')
    parser.add_argument('-min_allele_count', help='min number of minor alleles to be considered', type=int, default=1)
    parser.add_argument('-min_coverage', help='min coverage of sample to call gentype', type=int, default=10)

    global args
    args = parser.parse_args()

    genome = '/data/db/hg38.fa'
    genome_dict = read_genome(genome)
    chr_sizes = read_chr_sizes(genome)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}

    control, case  = read_samples(args.control, args.case)

    control_samples = list(control.keys())
    case_samples = list(case.keys())

    chunks_seq = str(genome_dict[args.chr].seq)[args.start:args.end].upper()

    nucs = chunk_analysis(args.chr, args.start, args.end, control_samples, case_samples, chunks_seq, args.workdir, args.min_allele_count, args.min_coverage)

