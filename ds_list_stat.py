import pandas as pd
import argparse


def read_samples(control_file, case_file):
    control_samples = {}
    with open(control_file, 'r') as ifh:
        for line in ifh:
            sample = line.rstrip().split('\t')[0]
            control_samples[sample] = 1
    case_samples = {}
    with open(case_file, 'r') as ifh:
        for line in ifh:
            sample = line.rstrip().split('\t')[0]
            case_samples[sample] = 1
    return control_samples, case_samples



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="cluster variants on same haplotype", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='input file in the v7 ds format')
    parser.add_argument('-o', help='output file')
    parser.add_argument('-case', help='case samples', default='/home/eraneyal/Projects/germline/germline_pipeline_v7/case_samples_v7p.txt')
    parser.add_argument('-control', help='control samples', default='/home/eraneyal/Projects/germline/germline_pipeline_v7/control_samples_v7p.txt' )
    args = parser.parse_args()

    data = pd.read_csv(args.i, sep='\t')

    control_samples, case_samples = read_samples(args.control, args.case)
    # print(list(data.columns))
    snv_counter = 0
    ins_counter = 0
    del_counter = 0
    for i in data.index:
        if len(data.at[i, 'ref']) == 1 and len(data.at[i, 'alt']) == 1:
            snv_counter += 1
        elif len(data.at[i, 'ref']) == 1 and len(data.at[i, 'alt']) > 1:
            ins_counter += 1
        else:
            del_counter += 1
        if data.at[i, 'odds_ratio'] < 20:
            print(str(data.at[i, 'gnomad_v3.1_af']) + '\t' + str(data.at[i, 'odds_ratio']))
        # print(data.at[i, 'gnomad_v3.1_af'])


    # print(snv_counter, del_counter, ins_counter)