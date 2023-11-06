#!~/miniconda3/bin/python

'''prep input for phiscs and phiscs-bnb'''

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
    # print(df_var)
    # df = (df_var > 0.05).astype(int)
    # print(df.columns)
    # df.columns.name = 'cell/mut'
    # df.to_csv(args.output,sep="\t")
    df_vaf = df_variant/ df_total
    df_vaf = df_vaf.T
    df_thresholded = (df_vaf >= 0.05).astype(int)
    df_thresholded[df_total == 0] = 3
    thresholded_str_array = np.char.replace(df_thresholded.values.astype(str), '3', '?').T
    snv_mat_phiscs = np.vstack((np.hstack((np.array([['cell_idx/mut_idx']]), np.array(df_vaf.T.columns)[np.newaxis,:])),
                                np.hstack((np.array(df_vaf.T.index)[:,np.newaxis], thresholded_str_array))))
    np.savetxt(args.output, snv_mat_phiscs, delimiter="\t", fmt='%s')

    
if __name__ == "__main__":
    main()
