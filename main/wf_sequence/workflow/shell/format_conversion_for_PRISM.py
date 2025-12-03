import pandas as pd
import numpy as np
import re
import argparse

def clear_seq(s):
    s_2 = re.sub(r"[^A-Z]", "", s)
    return s_2

def main():
    parser = argparse.ArgumentParser(description='Convert Casanovo output to PRISM input format')
    parser.add_argument('--input_mztab', required=True, help='Path to Casanovo output .mztab')
    parser.add_argument('--output_csv', required=True, help='Path to output .csv')
    args = parser.parse_args()


    with open(args.input_mztab,"r") as fi:
        k=0
        for ln in fi:
                k+=1
                if ln.startswith("PSH"):
                        starts=k-1
                        break
    df = pd.read_csv(args.input_mztab, sep='\t', header=starts)
    df = df[df['search_engine_score[1]']>0]
    df['Peptide'] = df['sequence'].apply(clear_seq)
    df = df.rename(columns={'spectra_ref': 'Fraction', 'search_engine_score[1]': 'ALC (%)'})
    df['Scan'] = df['Fraction'].apply(lambda x: x.split(':')[1])
    df['Fraction'] = df['Fraction'].apply(lambda x: x.split(':')[0])
    df['ALC (%)'] = df['ALC (%)'].apply(lambda x: round(x*100))
    df['-10lgP'] = df['ALC (%)'].apply(lambda x: round(-10*np.log10((100-x)/100)))
    df = df.dropna(axis=1, how='all')
    df = df.drop(columns=['search_engine', 'opt_ms_run[1]_aa_scores'])
    df.to_csv(args.output_csv, index=False)




if __name__ == "__main__":
    main()