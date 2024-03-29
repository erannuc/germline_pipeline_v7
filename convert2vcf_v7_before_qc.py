import os
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import argparse
from memory_profiler import profile
from bits import Bits
import gzip
from scipy.stats import binomtest
from datetime import datetime

global tiny
tiny= 1e-321
global base_map
base_map = { 0: 'T', 1: 'C', 2: 'G', 3: 'A'}


def find_genotype_indices(header_line, control_samples, case_samples):
    header_fields = header_line.split('\t')

    genotype_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype')]
    control_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in control_samples]
    case_indices = [i for i in range(len(header_fields)) if header_fields[i].endswith('genotype') and header_fields[i][:-9] in case_samples]
    index2sample = {i: header_fields[i][:-9] for i in range(len(header_fields)) if header_fields[i].endswith('genotype')}
    sample2index = {header_fields[i][:-9]: i for i in range(len(header_fields)) if header_fields[i].endswith('genotype')}

    return genotype_indices, control_indices, case_indices, index2sample, sample2index


def infer_genotypes(line, genotype_indices, reference_base, mode):

    line_data = line.strip().split('\t')
    # get all samples with alt genotype
    alt_genotypes_indices = [i for i in genotype_indices if '1' in line_data[i]]

    lineg = {}   # dict with genotyp
    linecov = {}
    samples_indices_to_check = alt_genotypes_indices

    if mode == 'full':
        samples_indices_to_check = genotype_indices

    for a in samples_indices_to_check:
        # make dict from the indels line
        indels_sample = {}
        # parse the indels field and populate
        indels_str = line_data[a-1]
        if indels_str:
            indels = indels_str.split(';')
            # print(indels)
            for i in indels:
                k, v = i.split(':')
                indels_sample[k] = int(v)
        sampleg, samplecov = sample_genotype(line_data[a-5:a-1], indels_sample, reference_base)
        for g in samplecov:
            if g not in linecov:
                linecov[g] = {a: samplecov[g]}
            else:
                linecov[g].update({a: samplecov[g]})
        for g in sampleg:
            if g not in lineg:
                lineg[g] = {a: sampleg[g]}
            else:
                lineg[g].update({a: sampleg[g]})
        # print(line_data[0], line_data[1], reference_base, line_data[a-5:a-1], indels_sample, sampleg)

    return lineg, linecov

def sample_genotype(nucs, indels, reference_base):
    indel_coverage = sum(indels.values())
    nucs = [int(n) for n in nucs]
    nucs_coverage = sum(nucs)
    sumreads = nucs_coverage + indel_coverage
    if sumreads==0:
        sumreads = tiny
    genotypes = {}
    coverages = {}
    for i in range(4):
        coverages[base_map[i]] = nucs[i]
        if nucs[i] / sumreads > 0.8:
            if base_map[i] != reference_base:
                if sumreads < 10:
                    genotypes[base_map[i]] = ('.', '.')
                else:
                    genotypes[base_map[i]] = (1, 1)
        elif nucs[i] / sumreads >= 0.2:
           if base_map[i] != reference_base:
                if sumreads < 10:
                    genotypes[base_map[i]] = ('.', '.')
                else:
                    genotypes[base_map[i]] = (0, 1)
    for i in indels:
        coverages[i] = indels[i]
        # check if one of the indels can be called:
        if indels[i] / sumreads > 0.8:
            if sumreads < 10:
                genotypes[i] = ('.', '.')
            else:
                genotypes[i] = (1, 1)
        elif indels[i] / sumreads >= 0.2:
            if sumreads < 10:
                genotypes[i] = ('.', '.')
            else:
                genotypes[i] = (0, 1)

    return genotypes, coverages


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


def left_normalization(chrseq, chr, pos, ref, alt):

    # check that the given variant is valid (consistent with the chromosome sequence)

    if not ref_valid(chrseq, pos, ref):
        # raise Exception(f"left_normalization: ref allele is not consistent with chromosome {chr} {pos} sequence")
        print(f"left_normalization: ref allele is not consistent with chromosome {chr} {pos} sequence")
        return chr, pos, ref, alt
    # first check that $ref and $alt are not empty

    if ref == '' or alt == '':
        pos -= 1
        previous_base = get_subseq(chrseq, pos, pos)
        ref = previous_base + ref
        alt = previous_base + alt

    while ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
        if ref == '' or alt == '':
            pos -= 1
            previous_base = get_subseq(chrseq, pos, pos)
            ref = previous_base + ref
            alt = previous_base + alt

    while ref[0] == alt[0] and len(ref) >= 2 and len(alt) >= 2:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    # ref = short_seq_annoation(ref)
    # alt = short_seq_annoation(alt)

    return chr, pos, ref, alt


def right_normalization(chrseq, chr, pos, ref, alt):

    # check that the given variant is valid (consistent with the chromosome sequence)

    if not ref_valid(chrseq, pos, ref):
        raise Exception(f"right_normalization: ref allele is not consistent with chromosome {chr} sequence")

    # first check that $ref and $alt are not empty

    if ref == '' or alt == '':
        pos -= 1
        previous_base = get_subseq(chrseq, pos, pos)
        ref = previous_base + ref
        alt = previous_base + alt

    if len(ref) != len(alt):
        ref = ref[1:] + get_subseq(chrseq, pos + len(ref), pos + len(ref))
        alt = alt[1:] + get_subseq(chrseq, pos + len(ref), pos + len(ref))
        pos += len(ref)
    else:
        pos += len(ref) - 1    # empirically justified

    while ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        if ref == '' or alt == '':
            pos += 1
            next_base = get_subseq(chrseq, pos, pos)
            ref = ref + next_base
            alt = alt + next_base

    while ref[-1] == alt[-1] and len(ref) >= 2 and len(alt) >= 2:
        ref = ref[:-1]
        alt = alt[:-1]
        pos -= 1

    # ref = short_seq_annoation(ref)
    # alt = short_seq_annoation(alt)

    return chr, pos, ref, alt

def short_seq_annoation(seq):
    if len(seq) > 50:
        seq = seq[0] + str(len(seq) - 2) + 'bases' + seq[-1]
    return seq


def ref_valid(chrseq, pos, ref):
    return str(chrseq[pos - 1:pos + len(ref) - 1]).upper() == ref


def get_subseq(chrseq, start, end):
    # will return genomic fragments - 1-based inclusive e.g get_subseq(3, 4, 7) will return string of 4 bases -
    # located in the forth, fifth, sixth and seventh position in chr3
    return str(chrseq[start - 1:end])


def print_vcf_header(ofh, samples, mode):
    print("##fileformat=VCFv4.3", file=ofh)
    print("##date=" + datetime.now().strftime("%d/%m/%Y:%H:%M:%S"), file=ofh)
    print("##genome=hg38", file=ofh)
    for i in range(1, 26):
        print(f"##contig=<ID={str(i)}>", file=ofh)
    print("##INFO=<ID=GA,Number=G,Type=Float,Description=\"gnomAD v2.1 frequencies for all variants in that position combined\">", file=ofh)
    print("##INFO=<ID=AC_CONT,Number=G,Type=Integer,Description=\"Total number of alternative alleles in control samples called genotypes\">", file=ofh)
    print("##INFO=<ID=AN_CONT,Number=G,Type=Integer,Description=\"Total number of alleles in control samples called genotypes\">", file=ofh)
    print("##INFO=<ID=AF_CONT,Number=G,Type=Float,Description=\"Alternative allele frequency in control samples\">", file=ofh)
    print("##INFO=<ID=IC_CONT,Number=G,Type=Integer,Description=\"Total number of alternative allele carriers in control samples\">", file=ofh)
    print("##INFO=<ID=IN_CONT,Number=G,Type=Float,Description=\"Total number of control samples\">", file=ofh)
    print("##INFO=<ID=AC_CASE,Number=G,Type=Integer,Description=\"Total number of alternative alleles in case samples called genotypes\">", file=ofh)
    print("##INFO=<ID=AN_CASE,Number=G,Type=Integer,Description=\"Total number of alleles in case samples called genotypes\">", file=ofh)
    print("##INFO=<ID=AF_CASE,Number=G,Type=Float,Description=\"Alternative allele frequency in case samples\">", file=ofh)
    print("##INFO=<ID=IC_CASE,Number=G,Type=Integer,Description=\"Total number of alternative allele carriers in case samples\">", file=ofh)
    print("##INFO=<ID=IN_CASE,Number=G,Type=Float,Description=\"Total number of case samples\">", file=ofh)
    print("##INFO=<ID=PC,Number=G,Type=Float,Description=\"Binomial test (two-sided) between control alleles and gnomAD frequency\">", file=ofh)
    print("##INFO=<ID=PT,Number=G,Type=Float,Description=\"Binomial test (two-sided) between case (tumor) alleles and gnomAD frequency\">", file=ofh)
    print("##INFO=<ID=PR,Number=G,Type=Float,Description=\"Proportions test (two-sided) between cases alleles and control alleles\">", file=ofh)
    if mode=='full':
        print("##FORMAT=<ID=GT,Number=G,Type=String,Description=\"Unphased genotype\"", file=ofh)
        print("##FORMAT=<ID=AP,Number=G,Type=Integer,Description=\"Number of covering reads with pair-end consensus supporting the allele call at this position\"", file=ofh)
        print("##FORMAT=<ID=DP,Number=G,Type=Integer,Description=\"Number of covering reads with pair-end consensus call at this position\"", file=ofh)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(samples), file=ofh)
    else:
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=ofh)


def missing_genotypes(line, indices):
    line_data = line.strip().split('\t')
    # search for unresolved genotype
    missing_genotypes_indices = [i for i in indices if '?' in line_data[i] or '.' in line_data[i]]
    return missing_genotypes_indices

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='input our vcf', required=True)
    parser.add_argument('-o', help='output standard vcf file, with no .bgz suffix which will be added (the output file is bgzipped)', required=True)
    parser.add_argument('-control', help='file with control samples', required=True)
    parser.add_argument('-case', help='file with case samples', required=True)

    parser.add_argument('-mode', help='vcf mode compact (only info) or full (with samples info)', choices=['compact', 'full'], default='compact')

    # this script convert from out pseudo vcf format to standard VCF format compatible with bcftools csq, before adding additional QC fields 

    global args
    args = parser.parse_args()

    control_samples = list(pd.read_csv(args.control, sep='\t', header=None)[0])
    control_samples = [i for i in control_samples if not i.startswith('#')]
    case_samples = list(pd.read_csv(args.case, sep='\t', header=None)[0])
    case_samples = [i for i in case_samples if not i.startswith('#')]

    genome = '/data/db/hg38.fa'
    genome_dict = read_genome(genome)

    chr_sizes = read_chr_sizes(genome)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}

    written_variants = set()
    with (gzip.open if args.i.endswith("gz") else open)(args.i, 'tr') as ifh, open(args.o + '.tmp.vcf', 'w') as ofh:
        print_vcf_header(ofh, control_samples + case_samples, args.mode)
        header_line = ifh.readline()
        genotype_indices, control_indices, case_indices, index2sample, sample2index = find_genotype_indices(header_line.strip(), control_samples, case_samples)

        linecount = 0
        for line in ifh:
            line_data = line.strip().split('\t')
            chr = line_data[0]
            pos = int(line_data[1])
            ref = line_data[2]
            gnomad = line_data[3]
            p_control = line_data[13]
            p_case = line_data[14]
            p_case_control = line_data[16]
            line_genotypes, line_coverages = infer_genotypes(line, genotype_indices, genome_dict[chr][pos-1], args.mode)
            control_samples_missing_genotype = len(missing_genotypes(line, control_indices))
            case_samples_missing_genotype = len(missing_genotypes(line, case_indices))
            # if control_samples_missing_genotype > 0:
            #    print(control_samples_missing_genotype, case_samples_missing_genotype)
            # print(line_genotypes)
            # print(line_coverages)
            refs = []
            alts = []
            poss = []
            ac_cont = []
            an_cont = []
            af_cont = []
            ic_cont = []
            in_cont = []
            ac_case = []
            an_case = []
            af_case = []
            ic_case = []
            in_case = []
            formats = []
            gnomads = []            
            p_controls = []
            p_cases = []
            p_case_controls = []

            for i in line_genotypes:

                # calculate allele frequencies to case and control

                ac_cont.append(sum([line_genotypes[i][s].count(1) for s in line_genotypes[i] if s in control_indices]))
                an_cont.append((len(control_indices) - control_samples_missing_genotype) * 2)   # fix later for sex chromosomes
                af_cont.append(ac_cont[-1] / (an_cont[-1] + tiny))
                ic_cont.append(len([1 for s in line_genotypes[i] if s in control_indices]))
                in_cont.append(len(control_indices) - control_samples_missing_genotype)
                ac_case.append(sum([line_genotypes[i][s].count(1) for s in line_genotypes[i] if s in case_indices]))
                an_case.append((len(case_indices) - case_samples_missing_genotype) * 2)   # fix later for sex chromosomes
                af_case.append(ac_case[-1] / (an_case[-1] + tiny))
                ic_case.append(sum([1 for s in line_genotypes[i] if s in case_indices]))
                in_case.append(len(case_indices) - case_samples_missing_genotype)
                gnomads.append(gnomad)
                p_controls.append(p_control)
                p_cases.append(p_case)
                p_case_controls.append(p_case_control)

                # print(i, line_genotypes[i])
                # print(ac_cont[-1], an_cont[-1], af_cont[-1], ic_cont[-1], in_cont[-1], ac_case[-1], an_case[-1], af_case[-1], ic_case[-1], in_case[-1])

                if len(i) == 1:  # SNV
                    poss.append(pos)
                    refs.append(genome_dict[chr][pos - 1])
                    alts.append(i)
                elif '_' not in i: # insertion
                    inref = i[0]
                    inalt = i
                    # left normalization here
                    chr, invcfpos, inref, inalt = left_normalization(genome_dict[chr].seq, chr, pos, inref, inalt)
                    poss.append(pos)
                    refs.append(inref)
                    alts.append(inalt)
                else:
                    # deletion. Nucleix features report deleted positions. to convert to HGVS (where the position is before the actual deletion we need the previous base. This also
                    # changes the position to the base before the reported deletion.
                    delstart, dellength = i.split('_')
                    delvcfpos = pos - 1
                    delvcfpos += int(delstart)  # delstart can be a negative number or zero
                    delref = get_subseq(genome_dict[chr].seq, delvcfpos, delvcfpos + int(dellength))
                    delalt = get_subseq(genome_dict[chr].seq, delvcfpos, delvcfpos)
                    # left normalization here
                    chr, delvcfpos, delref, delalt = left_normalization(genome_dict[chr].seq, chr, delvcfpos, delref, delalt)
                    poss.append(delvcfpos)
                    refs.append(delref)
                    alts.append(delalt)
                if args.mode == 'full':
                    # create string line with info data for each alternative allele
                    format_str = "\tGT:AP:DP"
                    for s in control_samples + case_samples:
                        dp = line_data[sample2index[s] + 1]
                        format_str += '\t'
                        if sample2index[s] not in line_genotypes[i]:
                            if sample2index[s] not in line_coverages[i]:
                                line_coverages[i][sample2index[s]] = 0
                            format_str += ':'.join(['0/0', str(line_coverages[i][sample2index[s]]), dp])
                        else:
                            geno_str = '/'.join([str(a) for a in line_genotypes[i][sample2index[s]]])
                            format_str += ':'.join([geno_str, str(line_coverages[i][sample2index[s]]), dp])

                    formats.append(format_str)
            # now write all genotypes to the output vcf file. each genotype to separate line and remember the variants to eliminate redundancy.

            for i in range(len(poss)):
                variant_key = '_'.join([chr, str(poss[i]), str(refs[i]), str(alts[i])])
                if variant_key not in written_variants:
                    info_str = f'GA={gnomads[i]};AC_CONT={ac_cont[i]};AN_CONT={an_cont[i]};AF_CONT={round(af_cont[i], 3)};IC_CONT={ic_cont[i]};IN_CONT={in_cont[i]};' \
                               f'AC_CASE={ac_case[i]};AN_CASE={an_case[i]};AF_CASE={round(af_case[i], 3)};IC_CASE={ic_case[i]};IN_CASE={in_case[i]};' \
                               f'PC={format(float(p_controls[i]), ".2E")};PT={format(float(p_cases[i]), ".2E")};PR={format(float(p_case_controls[i]), ".2E")}'
                    print(f'{chr}\t{str(poss[i])}\t.\t{refs[i]}\t{alts[i]}\t.\t.\t{info_str}', end='', file=ofh)
                    if args.mode == 'full':
                        print(formats[i], file=ofh)
                    else:
                        print(file=ofh)
                    written_variants.add(variant_key)
            linecount += 1
            if linecount % 1000 ==0:
                print(linecount)

    # sort the file

    os.system(f"bcftools sort {args.o}.tmp.vcf -o {args.o}.tmp.sorted.vcf")

    # bgzip and tabix

    # pysam.tabix_compress(args.o + '.tmp.sorted.vcf', args.o, force=False)
    pysam.tabix_index(args.o + '.tmp.sorted.vcf', force=True, preset='vcf', keep_original=False, csi=True)

    # rename and remove tmp files

    os.rename(args.o + '.tmp.sorted.vcf.gz', args.o + '.bgz')
    os.rename(args.o + '.tmp.sorted.vcf.gz.csi', args.o + '.bgz.csi')
    os.remove(args.o + '.tmp.vcf')

    # Steps:
    #   1. get information. For all called genotypes infer the variant allele
    #        if it is an SNV: write the SNV as alt
    #        if insertion or deletion => left align the variant and report in the vcf file.
    #        also create info file with various tags
