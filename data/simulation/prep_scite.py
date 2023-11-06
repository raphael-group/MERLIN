#!~/miniconda3/bin/python

'''prep output for scite'''

import argparse
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='convert output of simulation into proper input format for SCITE')
    parser.add_argument('-v','--variant',type=str,required=True,help='variant read count')
    parser.add_argument('-o',"--output",type=str, required=True,help='output file name')
    # mutation by cell matrix
    # just make a variant count 0/1 matrix??
    args = parser.parse_args()
    df_var = pd.read_csv(args.variant,header=0, index_col=0).T
    df = (df_var > 0.05).astype(int)
    data = df.values
    np.savetxt(args.output, data, delimiter=" ",fmt="%d")
    
if __name__ == "__main__":
    main()
