import argparse
from scipy.stats import chi2_contingency
from Bio import SeqIO
import numpy as np
import pandas as pd
import gzip

from nucleix_features import FeaturesAnalysis


tiny= 1e-300
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}

def get_overlapping_chits(samples):
    ochits_dict = {}
    for s in samples:
        features_obj = FeaturesAnalysis(s, 'eran', genome, '?')
        ochits_dict[s] = features_obj.get_wgrs_median('Overlapping.CHits')
    return ochits_dict

def find_atomic_rep(indel):
    rep_dict = {len(indel): 1}
    max_repeat_size = 5
    max_repeat_size = min(len(indel), max_repeat_size)
    print(max_repeat_size)
    for i in range(1, max_repeat_size):
        print(i)
        s = indel[:i]
        tmpindel = indel
        counter = 0
        while len(tmpindel) >= i:
            print(tmpindel)
            if tmpindel.startswith(s):
                counter += 1
                tmpindel = tmpindel[i:]
            else:
                break
        if len(tmpindel) == 0:
            rep_dict[i] = counter
    # take the repeat with the maximum number of appearances

    return sorted(rep_dict.items(), key=lambda x: (x[1]), reverse=True)[0]



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
            if len(line.rstrip().split('\t'))==1:
                sample = line.rstrip()
                control_samples[sample] = 'None'
            else:
                sample, busample = line.rstrip().split('\t')[:2]
                control_samples[sample] = busample
    case_samples = {}
    with open(case_file, 'r') as ifh:
        for line in ifh:
            if line.startswith("#"):
                continue
            if len(line.rstrip().split('\t'))==1:
                sample = line.rstrip()
                case_samples[sample] = 'None'
            else:
                sample, busample = line.rstrip().split('\t')[:2]
                case_samples[sample] = busample
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


def get_alt_reads(nucs, indels_str, ref_base):
    if indels_str:
        indels = indels_str.split(';')
    else:
        indels = []
    indels_sample = {}
    for i in indels:
        k, v = i.split(':')
        indels_sample[k] = int(v)
    indel_coverage = sum(indels_sample.values())
    alt_nucs_coverage = sum([int(nucs[i]) for i in range(4) if base_map[i] != ref_base])
    return alt_nucs_coverage + indel_coverage


def find_genotype_indices(header_line, control_samples, case_samples):
    header_fields = header_line.split('\t')

    genotype_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype')]
    control_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in control_samples]
    case_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in case_samples]
    index2sample = {i: header_fields[i][:-9] for i in range(len(header_fields)) if header_fields[i].endswith('genotype')}
    sample2index = {header_fields[i][:-9]: i for i in range(len(header_fields)) if header_fields[i].endswith('genotype')}

    return genotype_indices, control_indices, case_indices, index2sample, sample2index


def calculate_mean_non_called_reads(line_data, control_indices, case_indices):
    fracs_non_called_control = []
    for i in control_indices:  # the indices are of the genotype column
        sum_non_called_reads = non_called_reads([int(n) for n in line_data[i-5:i-1]], line_data[i-1])
        # print(i, [int(n) for n in line_data[i-5:i-1]], sum_non_called_reads)
        frac_non_called_control = sum_non_called_reads / (int(line_data[i+1]) + tiny)
        fracs_non_called_control.append(frac_non_called_control)
    fracs_non_called_case = []
    for i in case_indices:  # the indices are of the genotype column
        sum_non_called_reads = non_called_reads([int(n) for n in line_data[i-5:i-1]], line_data[i-1])
        frac_non_called_case = sum_non_called_reads / (int(line_data[i+1]) + tiny)
        fracs_non_called_case.append(frac_non_called_case)

    return sum(fracs_non_called_control) / len(fracs_non_called_control), sum(fracs_non_called_case) / len(fracs_non_called_case)

def non_called_reads(nucs, indels_str):
    if indels_str:
        indels = indels_str.split(';')
    else:
        indels = []
    indels_sample = {}
    for i in indels:
        k, v = i.split(':')
        indels_sample[k] = int(v)
    indel_coverage = sum(indels_sample.values())
    sumreads = sum(nucs) + indel_coverage + tiny
    sum_non_called_reads = 0
    for i in range(4):
        if nucs[i] / sumreads < 0.2:
            sum_non_called_reads += nucs[i]
    for i in indels_sample:
        if indels_sample[i] / sumreads < 0.2:
            sum_non_called_reads += indels_sample[i]

    return sum_non_called_reads

def calculate_alt_coverage_ratio(line_data, control_indices, case_indices, ochits_index_dict):
    ngen1 = 0
    sum1 = 0
    ngen0 = 0
    sum0 = 0
    for i in control_indices + case_indices:  # the indices are of the genotype column
        if ochits_index_dict[i] == 0:
            continue
        if line_data[i] == '0|1' or line_data[i] == '1|1':
            ngen1 += 1
            sum1 += int(line_data[i+1]) / ochits_index_dict[i]
        elif line_data[i] == '0|0':
            ngen0 += 1
            sum0 += int(line_data[i+1]) / ochits_index_dict[i]
    if ngen0 == 0 or ngen1 == 0:
        alt_coverage_ratio = 1.00
    else:
        alt_coverage_ratio = (sum1 / (ngen1 + tiny)) / (sum0 / (ngen0 + tiny) + tiny)

    return alt_coverage_ratio

def genotype_at_position(sample, chr, position, reference_base, has_indels):

    position = int(position) - 1   # move to zero based coordinates

    experiment, subject, source, treatment, *rest = sample.split('.')

    consensus_nucs_path = '/'.join(['/data', 'features', sample, f'{sample}.Chr{chr}.Consensus.Nucs'])
    snv_data = np.fromfile(consensus_nucs_path, dtype='i', count=4, offset=position*16)
    consensus_del_path = '/'.join(['/data', 'features', sample, f'{sample}.Chr{chr}.Consensus.Del'])
    consensus_ins_path = '/'.join(['/data', 'features', sample, f'{sample}.Chr{chr}.Consensus.Ins'])

    indels = {}
    indels[position] = {}

    if not has_indels:   # if the original sample does not have indels - no check in the backup
        return sample_genotype(snv_data, indels[position], reference_base)

    with open(consensus_ins_path, 'r') as ifh, open(consensus_del_path, 'r') as dfh:
        # read the insertions and deletions
        # for each indel event, seek Nucs coverage in the right place, return the frequency of the event
        # bfh.seek(chunk[1]*16)
        coverage = {}
        for line in ifh:
            line = line.rstrip()
            if line.startswith('@'):
                pos = int(line[1:])
                if pos < position:
                    continue
                if pos > position:
                    break
                indels[pos] = {}
            else:
                if pos < position:
                    continue
                seq, count = line.split(',')
                indels[pos][seq] = int(count)
        for line in dfh:
            line = line.rstrip()
            if line.startswith('@'):
                pos = int(line[1:])
                if pos < position - 150:
                    continue
                if pos > position + 5:
                    break
                if pos not in indels:
                    indels[pos] = {}
            else:
                if pos < position - 150:
                    continue
                if pos > position + 5:
                    continue
                length, count = line.split(',')
                for i in range(0, int(length)):
                    if pos + i in indels:
                        indels[pos + i]['_'.join([str(-abs(i)), length])] = int(count)
                    else:
                        indels[pos + i] = {'_'.join([str(-abs(i)), length]): int(count)}

    return sample_genotype(snv_data, indels[position], reference_base)

def sample_genotype(nucs, indels, reference_base):
    indel_coverage = sum(indels.values())
    nucs = [int(n) for n in nucs]
    nucs_coverage = sum(nucs)
    sumreads = nucs_coverage + indel_coverage + tiny
    genotypes = {}
    for i in range(4):
        if nucs[i] / sumreads > 0.80:
            if base_map[i] != reference_base:
                if sumreads < 10:
                    genotypes[base_map[i]] = '.|.'
                else:
                    genotypes[base_map[i]] = '1|1'
        elif nucs[i] / sumreads >= 0.2:
           if base_map[i] != reference_base:
                if sumreads < 10:
                    genotypes[base_map[i]] = '.|.'
                else:
                    genotypes[base_map[i]] = '0|1'
    for i in indels:
        # check if one of the indels can be called:
        if indels[i] / sumreads > 0.80:
            if sumreads < 10:
                genotypes[i] = '.|.'
            else:
                genotypes[i] = '1|1'
        elif indels[i] / sumreads >= 0.2:
            if sumreads < 10:
                genotypes[i] = '.|.'
            else:
                genotypes[i] = '0|1'

    return genotypes


def list_analysis(input_quazi_vcf, output_quazi_vcf, control_samples, case_samples, ochits_dict):
    with (gzip.open if input_quazi_vcf.endswith("gz") else open)(input_quazi_vcf, 'tr') as ifh, open(args.o, 'w') as ofh:
    # with open(input_quazi_vcf, 'r') as ifh, open(output_quazi_vcf, 'w') as ofh:
        header = ifh.readline().strip()
        header_fields = header.split('\t')

        new_columns = ['af_controls', 'rf_controls', 'afrf_controls', 'af_case', 'rf_case', 'afrf_case', 'mean_non_called_reads_frac_controls',
                       'mean_non_called_reads_frac_case', 'alt_coverage_ratio', 'fraction_alt_genotype_validated', 'number_alt_genotype_checked']

        output_header_fields = header_fields[:17] + new_columns + header_fields[17:]
        print('\t'.join(output_header_fields), file=ofh)

        total_variants = 0
        pass_variants = 0
        alt_alleles_control_index = header_fields.index('alt_alleles_control')
        tot_alleles_control_index = header_fields.index('tot_alleles_control')
        alt_alleles_case_index = header_fields.index('alt_alleles_case')
        tot_alleles_case_index = header_fields.index('tot_alleles_case')
        genotype_indices, control_indices, case_indices, index2sample, sample2index = find_genotype_indices(header, list(control_samples.keys()), list(case_samples.keys()))

        ochits_index_dict = {sample2index[s]: ochits_dict[s] for s in ochits_dict}

        Exp900Vcontrols_first_index = header_fields.index('Exp900V.controls.plasma.digested:T')
        Exp900Vcancer_first_index = header_fields.index('Exp900V.cancer.plasma.digested:T')

        for line in ifh:
            line_data = line.strip().split('\t')
            # print(line_data[:5])
            chr = line_data[0]
            position = line_data[1]
            ref = line_data[2]

            pass_filters = True
            total_variants += 1

            # check the 900V samples and make chi-square statistical test for the sample

            Exp900Vcontrols_alt_reads = get_alt_reads(line_data[Exp900Vcontrols_first_index:Exp900Vcontrols_first_index+4],
                                                      line_data[Exp900Vcontrols_first_index+4], ref)
            Exp900Vcontrols_tot_reads = int(line_data[Exp900Vcontrols_first_index+6])
            Exp900Vcontrols_ref_reads = Exp900Vcontrols_tot_reads - Exp900Vcontrols_alt_reads
            Exp900Vcase_alt_reads = get_alt_reads(line_data[Exp900Vcancer_first_index:Exp900Vcancer_first_index+4],
                                                      line_data[Exp900Vcancer_first_index+4], ref)
            Exp900Vcase_tot_reads = int(line_data[Exp900Vcancer_first_index+6])
            Exp900Vcase_ref_reads = Exp900Vcase_tot_reads - Exp900Vcase_alt_reads

            # stat, p, dof, expected = chi2_contingency([[Exp900Vcontrols_alt_reads, Exp900Vcontrols_ref_reads],
            #                                            [Exp900Vcases_alt_reads, Exp900Vcases_ref_reads]])


            af_controls = int(line_data[alt_alleles_control_index]) / (int(line_data[tot_alleles_control_index])+ tiny)
            rf_controls = Exp900Vcontrols_alt_reads / (Exp900Vcontrols_tot_reads + tiny)
            afrf_controls = af_controls / (rf_controls + tiny)
            if afrf_controls > 2:
                afrf_controls = 1.000

            af_case = int(line_data[alt_alleles_case_index]) / (int(line_data[tot_alleles_case_index]) + tiny)
            rf_case = Exp900Vcase_alt_reads / (Exp900Vcase_tot_reads + tiny)
            afrf_case = af_case / (rf_case + tiny)
            if afrf_case > 2:
                afrf_case = 1.000

            # print(af_controls, rf_controls, afrf_controls, af_case, rf_case, afrf_case)

            mean_non_called_reads_controls, mean_non_called_reads_case = calculate_mean_non_called_reads(line_data, control_indices, case_indices)

            alt_coverage_ratio = calculate_alt_coverage_ratio(line_data, control_indices, case_indices, ochits_index_dict)

            # print(mean_non_called_reads_controls, mean_non_called_reads_case)

            # now search for genotype agreement percentage in the backup samples

            # for each positive sample with back-up: calculate the genotype at the backup. Calculate overall positive allele agreement

            alt_genotypes_supported = 0
            alt_genotypes_checked = 0

            for i in case_samples:

                # check only hetero/homzygous genotype and only sample with back-up

                if case_samples[i] == 'None':
                    continue
                if line_data[sample2index[i]] == '0|0' or line_data[sample2index[i]] == '.|.' or line_data[sample2index[i]] == '?|?':
                    continue
                if line_data[sample2index[i] - 1]: # if there are indels, othewise flag it and don't check indels also in the backup sample
                    has_indels = True
                else:
                    has_indels = False
                # case_samples[i] is the back-up sample for sample i
                genotype_bu = genotype_at_position(case_samples[i], chr, position, ref, has_indels)
                # print(case_samples[i], genotype_bu)
                case_genotype = line_data[sample2index[i]]
                # print(case_genotype, genotype_bu)
                bu_any_genotype_called = False
                bu_match_genotype_called = False
                for i in genotype_bu:
                    if genotype_bu[i].count('1') or genotype_bu[i].count('0'):
                        bu_any_genotype_called = True
                    if genotype_bu[i] == case_genotype:
                        bu_match_genotype_called = True

                if bu_any_genotype_called:
                    alt_genotypes_checked += 1
                if bu_match_genotype_called:
                    alt_genotypes_supported += 1


            fraction_alt_genotype_validated = alt_genotypes_supported / (alt_genotypes_checked + tiny)
            print("#######", total_variants, fraction_alt_genotype_validated, alt_genotypes_checked, alt_coverage_ratio)

            new_data_columns = [str(af_controls), str(rf_controls), str(afrf_controls), str(af_case), str(rf_case), str(afrf_case), str(mean_non_called_reads_controls),
                        str(mean_non_called_reads_case), str(alt_coverage_ratio), str(fraction_alt_genotype_validated), str(alt_genotypes_checked)]

            output_data_fields = line_data[:17] + new_data_columns + line_data[17:]
            print('\t'.join(output_data_fields), file=ofh)
            if total_variants % 1000 == 0:
                print(total_variants)

# additional filtering step to account for:
#   genotype calling errors, by going back to V900 pools and making chi-square between the ref/non-ref counts in al and control
#   For each input variant with called genotype check existence of variant in another case sample, preferably WB
#
#   Repetative regions. Mark variants with repetition flag according to Dani's somatic pipeline definition
#


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='file in our internal vcf format with input variants. could be gzipped or not')
    parser.add_argument('-case', help='file with list of case samples', default='case_samples_with_backup_v7.txt')
    parser.add_argument('-control', help='file with list of control samples and bachup samples', default='control_samples_v7.txt')
    parser.add_argument('-pools_case', help='file with list of pool of case samples and backup samples', default='case_pool_samples_v7.txt')
    parser.add_argument('-pools_control', help='file with list of pool of control samples', default='control_pool_samples_v7.txt')
    parser.add_argument('-o', help='file in our internal vcf format with output variants', default='variants_reduced_v7.vcf')

    global args
    args = parser.parse_args()

    genome = '/home/eraneyal/Genomes/hg38.fa'
    genome_dict = read_genome(genome)
    chr_sizes = read_chr_sizes(genome)

    '''
    features_obj = FeaturesAnalysis('TESExp19N.SL_CON358.plasma.digested.nodup', 'eran', genome, '?')
    overlapping_chits_median_coverage = features_obj.median_coverage('Overlapping.CHits', '/home/eraneyal/Projects/targeted/hg38.AloofSane.bed')
    print(overlapping_chits_median_coverage)
    exit(0)
    '''
    # ochits_median = features_obj.get_wgrs_median('Overlapping.CHits')
    # def median_coverage(self, feature, bed_file, mode='nuc'):

    db = pd.read_csv('/mnt/Data/Wgrs/db.csv', sep=',', index_col='Name')

    control, case, control_pools, case_pools = read_samples(args.control, args.case, args.pools_control, args.pools_case)

    samples = list(control.keys()) + list(case.keys())
    ochits_dict = get_overlapping_chits(samples)
    # for i in ochits_dict:
    #    print(i, ochits_dict[i])
    # ochits_dict['TESExp19N.SL_CON2744.plasma.digested.nodup'] = 1

    list_analysis(args.i, args.o, control, case, ochits_dict)
