
'''
Prep the output of simulation for input to calder
'''

import argparse
import numpy as np
import pandas as pd

def make_matrix(total, variant):
  ''' total, variant are str of file for total and variant matrix
    returns a dataframe of the combined matrix
  '''
  df_total = pd.read_csv(total, header=0, index_col=0)
  df_var = pd.read_csv(variant, header=0,index_col=0)
  # print(df_total)
  # print(df_var)
  
  assert df_total.shape == df_var.shape
  num_cols = df_total.shape[1]
  num_rows = df_total.shape[0]
  
  df = pd.DataFrame()
  for col in df_total.columns:
      df[str(col)] = df_total[col] - df_var[col]
      df[str(col)+"_"] = df_var[col]
  cols = []
  for idx in df_total.columns:
    cols.append(idx)
    cols.append(idx)
  df.columns = cols
  df.columns.name = 'gene_id'
  # print("df",df)
  return df

def main():
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')

    parser.add_argument('--total', type=str, required=True, help='total count matrix')
    parser.add_argument('--variant', type=str, required=True, help='variant count matrix')
    parser.add_argument('--out', type=str, required=True, help='Output file prefix')
    args = parser.parse_args()
    
    df = make_matrix(args.total, args.variant)
    df.to_csv(f'{args.out}', sep="\t", index=True,header=True)
    # stacked_df = 
    # np.savetxt(args.out, df, delimiter="\t", fmt='%s')
    
if __name__ == "__main__":
    main() 