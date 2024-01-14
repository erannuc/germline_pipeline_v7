import argparse

import pandas as pd

tiny = 1E-10

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='input ds format', required=True)
    parser.add_argument('-o', help='results text file', required=True)

    global args
    args = parser.parse_args()

    data = pd.read_csv(args.i, sep='\t')
    with open(args.o, 'w') as ofh:
        for i in data.index:
            # print(i)
            case_counter_total = 0
            case_counter_alt = 0
            control_counter_total = 0
            control_counter_alt = 0
            af_control = data.at[i, 'af_control']
            af_case = data.at[i, 'af_case']
            for c in data.columns:
                if c.startswith('WGS') and 'LNG' in c and data.at[i, c] != -1:
                    case_counter_total += 2
                    case_counter_alt += data.at[i, c]
                if c.startswith('WGS') and 'CON' in c and not '.EG_' in c and data.at[i, c] != -1:
                    control_counter_total += 2
                    control_counter_alt += data.at[i, c]
            af_control_test = control_counter_alt / control_counter_total
            af_case_test = case_counter_alt / case_counter_total
            print(str(af_case/(af_control + tiny)) + '\t' + str(af_case_test/(af_control_test + tiny)), file=ofh)
