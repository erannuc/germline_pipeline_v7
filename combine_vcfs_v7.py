import pandas as pd
import argparse
import os
import stat
import gzip

import multiprocessing


def chunk_combine(args):
    c, outdir = args
    print(c)
    chr, start, end = c

    vcf1 = f'{outdir}/{chr}_{start}_{end}_{0}_{50}.vcf.bgz'
    vcf2 = f'{outdir}/{chr}_{start}_{end}_{50}_{100}.vcf.bgz'
    vcf3 = f'{outdir}/{chr}_{start}_{end}_{100}_{150}.vcf.bgz'
    vcf4 = f'{outdir}/{chr}_{start}_{end}_{150}_{200}.vcf.bgz'
    vcf5 = f'{outdir}/{chr}_{start}_{end}_{200}_{250}.vcf.bgz'
    vcf6 = f'{outdir}/{chr}_{start}_{end}_{250}_{300}.vcf.bgz'
    vcf7 = f'{outdir}/{chr}_{start}_{end}_{300}_{350}.vcf.bgz'
    vcf8 = f'{outdir}/{chr}_{start}_{end}_{350}_{400}.vcf.bgz'
    vcf9 = f'{outdir}/{chr}_{start}_{end}_{400}_{450}.vcf.bgz'
    vcf10 = f'{outdir}/{chr}_{start}_{end}_{450}_{500}.vcf.bgz'
    vcf11 = f'{outdir}/{chr}_{start}_{end}_{500}_{550}.vcf.bgz'
    vcf12 = f'{outdir}/{chr}_{start}_{end}_{550}_{600}.vcf.bgz'
    vcf13 = f'{outdir}/{chr}_{start}_{end}_{600}_{650}.vcf.bgz'
    vcf14 = f'{outdir}/{chr}_{start}_{end}_{650}_{700}.vcf.bgz'
    vcf15 = f'{outdir}/{chr}_{start}_{end}_{700}_{750}.vcf.bgz'
    vcf16 = f'{outdir}/{chr}_{start}_{end}_{750}_{764}.vcf.bgz'


    out_vcf = f'{outdir}/{chr}_{start}_{end}.vcf'

    with gzip.open(vcf1, 'rt') as ifh1, gzip.open(vcf2, 'rt') as ifh2, gzip.open(vcf3, 'rt') as ifh3, \
        gzip.open(vcf4, 'rt') as ifh4, gzip.open(vcf5, 'rt') as ifh5, gzip.open(vcf6, 'rt') as ifh6, \
        gzip.open(vcf7, 'rt') as ifh7, gzip.open(vcf8, 'rt') as ifh8, gzip.open(vcf9, 'rt') as ifh9, \
        gzip.open(vcf10, 'rt') as ifh10, gzip.open(vcf11, 'rt') as ifh11, gzip.open(vcf12, 'rt') as ifh12, \
        gzip.open(vcf13, 'rt') as ifh13, gzip.open(vcf14, 'rt') as ifh14, gzip.open(vcf15, 'rt') as ifh15, \
        gzip.open(vcf16, 'rt') as ifh16, open(out_vcf, 'w') as ofh:
        linenum = 0
        for line1 in ifh1:
            line1 = line1.strip()
            line2 = ifh2.readline().strip()
            line3 = ifh3.readline().strip()
            line4 = ifh4.readline().strip()
            line5 = ifh5.readline().strip()
            line6 = ifh6.readline().strip()
            line7 = ifh7.readline().strip()
            line8 = ifh8.readline().strip()
            line9 = ifh9.readline().strip()
            line10 = ifh10.readline().strip()
            line11 = ifh11.readline().strip()
            line12 = ifh12.readline().strip()
            line13 = ifh13.readline().strip()
            line14 = ifh14.readline().strip()
            line15 = ifh15.readline().strip()
            line16 = ifh16.readline().strip()


            linenum += 1
            data1 = line1.split('\t')
            data2 = line2.split('\t')
            data3 = line3.split('\t')
            data4 = line4.split('\t')
            data5 = line5.split('\t')
            data6 = line6.split('\t')
            data7 = line7.split('\t')
            data8 = line8.split('\t')
            data9 = line9.split('\t')
            data10 = line10.split('\t')
            data11 = line11.split('\t')
            data12 = line12.split('\t')
            data13 = line13.split('\t')
            data14 = line14.split('\t')
            data15 = line15.split('\t')
            data16 = line16.split('\t')

            # ensure first field is identical other wise exit with error

            if not identical_first_field([data1[0], data2[0], data3[0], data4[0], data5[0], data6[0], data7[0], data8[0], data9[0], data10[0], 
                                         data11[0], data12[0], data13[0], data14[0], data15[0], data16[0]]):
                raise Exception("Order error in line " + str(linenum))
            else:
                print('\t'.join(data1[0:3]), end='', file=ofh)
                print('\t' + '\t'.join(data1[3:] + data2[3:] + data3[3:] + data4[3:] + data5[3:] + data6[3:] + data7[3:] + data8[3:] + 
                                       data9[3:] + data10[3:] + data11[3:] + data12[3:] + data13[3:] + data14[3:] + data15[3:] +  data16[3:]), file=ofh)

    os.system(' '.join(['bgzip', '-f', out_vcf]))
    os.rename(out_vcf + '.gz', out_vcf + '.bgz')
    os.system(' '.join(['tabix -pvcf', '-fC', out_vcf + '.bgz']))

    return 0


def identical_first_field(first_fields):
    if len(set(first_fields)) == 1:
        return True
    else:
        return False


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
    parser.add_argument('-outdir', help='output directory. full path', required=True)
    parser.add_argument('-cstart', help='first chunk to sbatch', type=int, default=0)
    parser.add_argument('-cend', help='last chunk to sbatch', type=int)
    parser.add_argument('-nsamples', help='overall number of samples', type=int, default=764)
    parser.add_argument('-samples_chunk', help='samples chunk', type=int, default=50)
    parser.add_argument('-n', help='number of processes', type=int, default=5)


    global args
    args = parser.parse_args()

    chr_sizes = read_chr_sizes(args.g)
    chr_sizes = {c: chr_sizes[c] for c in chr_sizes if int(c) <=22}
    all_chunks = create_chunks(chr_sizes, args.chunk)

    # currently not in use. hard coded samples_chunk of 50
    # sample_chunks = create_sample_chunks(args.nsamples, args.samples_chunk)

    pool_obj = multiprocessing.Pool(processes=args.n)
    args_lists = []

    for ch in all_chunks:
        # check that dont exists:
        chunk_string = '_'.join([str(i) for i in ch])
        if (not os.path.exists(f'{args.outdir}/{chunk_string}.vcf.bgz')) or (not os.path.exists(f'{args.outdir}/{chunk_string}.vcf.bgz.csi')):
            args_lists.append([ch, args.outdir])

    out = pool_obj.map(chunk_combine, args_lists)
