import pandas as pd
import argparse

def is_perfect_phylogeny(matrix):
    for i in range(matrix.shape[1]):
        for j in range(i + 1, matrix.shape[1]):
            group1 = set(matrix.index[matrix.iloc[:, i] == 1])
            group2 = set(matrix.index[matrix.iloc[:, j] == 1])
            if group1 & group2 and group1 - group2 and group2 - group1:
                return False
    return True

def build_newick(matrix, mutations):
    if matrix.shape[0] == 1:
        return matrix.index[0]
    if matrix.empty or not mutations:
        short =  ",".join(matrix.index)
        return short
    for mutation in mutations:
        positive_branch = matrix[matrix[mutation] == 1]
        negative_branch = matrix[matrix[mutation] == 0]
        if len(positive_branch) > 0 and len(negative_branch) > 0:
            break
    if len(positive_branch) == 0:
        return build_newick(negative_branch, mutations[1:])
    if len(negative_branch) == 0:
        return build_newick(positive_branch, mutations[1:])
    positive_tree = build_newick(positive_branch, mutations[1:])
    negative_tree = build_newick(negative_branch, mutations[1:])
    return f"({positive_tree},{negative_tree})"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--csv", required=True, help="binarized U matrix matrix")
    parser.add_argument("-o", "--output", required=True, help="Path to output Newick tree file")
    args = parser.parse_args()
    
    matrix = pd.read_csv(args.csv, index_col=0)
    if not is_perfect_phylogeny(matrix):
        raise ValueError("The input matrix is not a perfect phylogeny")
    
    mutations = list(matrix.columns)
    newick_tree = build_newick(matrix, mutations) + ";"
    
    with open(args.output, "w") as f:
        f.write(newick_tree)