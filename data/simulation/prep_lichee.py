#!~/miniconda3/bin/python

'''prep input for lichee'''

import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='convert output of simulation into proper input format for SCITE')
    parser.add_argument('-v','--variant',type=str,required=True,help='variant read count')
    parser.add_argument('-t', '--total',type=str,required=True, help='total read count')
    parser.add_argument('-o',"--output",type=str, required=True,help='output file name')
    # mutation by cell matrix
    # just make a variant count 0/1 matrix??
    args = parser.parse_args()
    df_variant = pd.read_csv(args.variant,header=0, index_col=0)
    df_total = pd.read_csv(args.total,header=0, index_col=0)
    
    
    df_vaf = df_variant/ df_total
    df = df_vaf.T
    cell_list = list(df.columns)
    nmut = len(df.index)
    df['#chr'] = [(i+1) for i in range(nmut)]
    df['position'] =  [(i+1) for i in range(nmut)]
    df['description'] =  df.index
    df['Normal'] = 0.0
    new_col = ["#chr",'position','description','Normal'] + cell_list
    df = df[new_col]
    # print(df)
    
    df.to_csv(args.output, index=None,sep='\t')
    
if __name__ == "__main__":
    main()
