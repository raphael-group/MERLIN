#!~/miniconda3/bin/python

'''
Prep the output of simulation for input to pairtree
'''

import argparse
import numpy as np
import pandas as pd
import json

    
def main():
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')

    parser.add_argument('-t', '--total', type=str, required=True, help='total count matrix')
    parser.add_argument("-v",'--variant', type=str, required=True, help='variant count matrix')
    parser.add_argument("-o",'--out', type=str, required=True, help='Output file prefix')
    args = parser.parse_args()
    df_variant = pd.read_csv(args.variant,header=0,index_col=0).T.astype(int)
    df_total = pd.read_csv(args.total,header=0,index_col=0).T.astype(int)
    # print(df_variant.index)
    with open(f'{args.out}/pairtree_input.ssm', 'w') as out:
        out.write('\t'.join(['id', 'name', 'var_reads', 'total_reads', 'var_read_prob']) + '\n')
        for mut_idx, mutation in enumerate(df_variant.index):
            out.write(f's{mut_idx}\tS_{mut_idx}\t')
            out.write(','.join(map(str, df_variant.loc[mutation].values)) + '\t')
            out.write(','.join(map(str, df_total.loc[mutation].values)) + '\t')
            out.write(','.join(map(str, [1]*len(df_variant.columns))))
            out.write('\n')
        
    pairtree_dict = {'samples': list(map(str, df_variant.columns)),
                    'clusters': [[f's{mut_idx}'] for mut_idx in range(len(df_variant))],
                    'garbage': []}
    with open(f'{args.out}/params.json', 'w') as out:
        out.write(json.dumps(pairtree_dict))
    
if __name__ == "__main__":
    main()