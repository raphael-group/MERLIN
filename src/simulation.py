import argparse,json
import numpy as np
import pandas as pd
import networkx as nx
from datetime import datetime
import random

from scipy.stats import poisson, binom, betabinom

"""
Simulate a clonal tree with n nodes and m mutations,
by assigning each of the m mutations to a random node.

Output:
    - tree: a networkx DiGraph object representing the clonal tree
      where the root is labeled by 0, and the other nodes are labeled
      by 1, 2, ..., n - 1. The attached mutations are stored in the 
      'mutation' attribute of the nodes.
"""

# n = genomes
# m = mutations
def simulate_clonal_tree(m, n):
    assert m >= n # for now, to keep it simple
    tree = nx.DiGraph()
    tree.add_node(0) 

    for i in range(1,n):
        parent = np.random.choice(np.arange(i))
        tree.add_node(i)
        tree.add_edge(parent, i)
    # Ensures that every node has at least one mutation

    mutation_to_clone_mapping = {i:i for i in range(m)}
    return tree, mutation_to_clone_mapping

"""
Constructs the clonal matrix from a clonal tree.
"""
# i don't quite follow here
# B[j,i] = 1 if i is j's parent
def construct_clonal_matrix(tree):
    n = len(tree.nodes)

    B = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if nx.has_path(tree, i, j):
                B[j, i] = 1

    return B

def make_cell_lineage_tree(clone_tree):
    clones = list(nx.topological_sort(clone_tree))
    # print(clone_tree.edges)
    # print(clones)
    cell_tree = nx.DiGraph()
    cell_tree.add_node(0) # by our construction 0 is always an ancestor
    for clone in clones:
        if clone == 0: continue
        parent_in_clone_tree = list(clone_tree.predecessors(clone))[0]
        # print(parent_in_clone_tree)
        assert parent_in_clone_tree in cell_tree.nodes
        subtree = nx.descendants(cell_tree, parent_in_clone_tree) # set of nodes
        subtree.add(parent_in_clone_tree)
        parent_in_cell_tree = random.choice(tuple(subtree))
        cell_tree.add_edge(parent_in_cell_tree, clone)
    # print(cell_tree.edges)
    return cell_tree


"""
Simulate a usage matrix with s samples and n clones,
by randomly selecting k of n clones and assigning each of them
a random usage probability.

Output:
    - matrix: a s x n matrix where each row represents a sample 
    and the sum of the entries is 1.
"""
def simulate_usage_matrix(cell_tree, s, clone_tree, threshold=0.05):
    matrix = np.zeros((s, len(clone_tree)))
    cell_to_clone = dict()
    for i in range(s):
        clone = np.random.choice(cell_tree.nodes)
        cell_to_clone[i] = clone
        # print(clone)
        all_clones = list(nx.ancestors(cell_tree, clone))
        all_clones.append(clone)
        # print(all_clones)
        if len(all_clones) > 20:
            all_clones = np.random.choice(all_clones,size=20,replace = False)
        usages = np.random.dirichlet(np.ones(len(all_clones)), size=1)[0]
        while (usages < threshold).any():
            usages = np.random.dirichlet(np.ones(len(all_clones)), size=1)[0]
        for clone, usage in zip(all_clones, usages):
            matrix[i, clone] = usage
    return matrix, pd.DataFrame.from_dict(cell_to_clone,orient='index').reset_index()

"""
Simulate a mutation read counts from a clonal matrix, usage matrix and 
a mapping from mutations to clone.
"""
def simulate_read_counts(usage_matrix, clonal_matrix, mutation_to_clone_mapping, num_mutations, coverage):
    F = usage_matrix @ clonal_matrix

    # Clamp to [0, 1] for floating point errors
    F[F < 0] = 0
    F[F > 1] = 1
    # print("U:", usage_matrix)
    # print("clonal matrix", clonal_matrix)
    # print("map", mutation_to_clone_mapping)
    variant_count_matrix = np.zeros((usage_matrix.shape[0], num_mutations),dtype=int)
    total_count_matrix   = poisson.rvs(coverage, size=variant_count_matrix.shape)
    # print(F)
    # for mutation in range(num_mutations):
    #     for s in range(usage_matrix.shape[0]):
    #         p = F[s, mutation]
    #         epsilon = 0.01
    #         if p < epsilon: p = epsilon
    #         if (1-p) < epsilon: p = 1-epsilon
    #         precision = 15
    #         ado_alpha = p * precision
    #         ado_beta = precision * (1 - p)
    #         f = betabinom.rvs(
    #                 total_count_matrix[s, mutation],
    #                 ado_alpha,ado_beta
    #         )
    #         # there could be a bug here...
    #         variant_count_matrix[s, mutation] = f
    for mutation in range(num_mutations):
        for s in range(usage_matrix.shape[0]):
            f = binom.rvs(
                    total_count_matrix[s, mutation],
                    F[s, mutation_to_clone_mapping[mutation]]
            )
            # there could be a bug here...
            variant_count_matrix[s, mutation] = f

    assert np.all(variant_count_matrix <= total_count_matrix)
    return variant_count_matrix, total_count_matrix

def observe_frequency_matrix(variant_count_matrix, total_count_matrix, mutation_to_clone_mapping, num_clones):
    clone_to_mutation_mapping = {}
    for clone in range(num_clones+1):
        clone_to_mutation_mapping[clone] = []

    for mutation in mutation_to_clone_mapping:
        clone = mutation_to_clone_mapping[mutation]
        clone_to_mutation_mapping[clone].append(mutation)

    clone_mut = lambda c: clone_to_mutation_mapping[c]

    obs_frequency_matrix = np.zeros((variant_count_matrix.shape[0], len(clone_to_mutation_mapping.keys())))
    collapsed_variant_matrix = np.zeros((variant_count_matrix.shape[0], len(clone_to_mutation_mapping.keys())))
    collapsed_total_matrix   = np.zeros((variant_count_matrix.shape[0], len(clone_to_mutation_mapping.keys())))
    for s in range(obs_frequency_matrix.shape[0]):
        for clone in range(obs_frequency_matrix.shape[1]):
            variant_reads = sum([variant_count_matrix[s, m] for m in clone_mut(clone)])
            total_reads   = sum([total_count_matrix[s, m] for m in clone_mut(clone)])
            if total_reads > 0:
                obs_frequency_matrix[s, clone] = variant_reads / total_reads
            collapsed_variant_matrix[s, clone] = variant_reads
            collapsed_total_matrix[s, clone]   = total_reads
    
    return obs_frequency_matrix, collapsed_variant_matrix, collapsed_total_matrix,clone_to_mutation_mapping

def tree_to_newick(T, root=None):
    ''' 
    takes in a nx.DiGraph object T (need to be a tree)\\
    returns the Newick string format
    '''
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"    


def main():
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')

    parser.add_argument('-m', '--mutations', type=int, required=True, help='Number of mutations.')
    parser.add_argument('-g', '--genomes', type=int, required=True, help='Number of genomes.')
    parser.add_argument('-n', '--cells', type=int, required=True, help='Number of cells.')
    parser.add_argument('-c', '--coverage', type=float, required=True, help='Expected sequencing coverage.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output prefix.')
    parser.add_argument('-t', '--threshold', type=float, required=False, help="minimum allele frequency allowed",default=0.05)
    # parser.add_argument('-h', required=False, action='store_true')
    args = parser.parse_args()
   
    # help_message = '''usage:
    # '''
    ending = (args.output.split("rep")[-1])
    if ending.isnumeric():
        seed = 40 + int(ending)
    else: seed = 42
    np.random.seed(seed)
    
    assert args.mutations >= args.genomes, 'Number of mutations must be greater than or equal to number of genomes.'
    tree, mutation_to_clone_mapping = simulate_clonal_tree(args.mutations, args.genomes)
    
    clonal_matrix = construct_clonal_matrix(tree)
    
    cell_tree = make_cell_lineage_tree(tree)
    
    usage_matrix, cell_to_clone = simulate_usage_matrix(cell_tree, args.cells, tree, args.threshold)
    
    variant_matrix, total_matrix = simulate_read_counts(
            usage_matrix, clonal_matrix, mutation_to_clone_mapping, 
            args.mutations, args.coverage
    )


    f_hat, collapsed_variant_matrix, collapsed_total_matrix, clone2mut = observe_frequency_matrix(variant_matrix, total_matrix, mutation_to_clone_mapping, args.genomes)
    for k in clone2mut:
        if len(clone2mut[k]) == 0: clone2mut[k] = 0
        elif clone2mut[k] == [0]: clone2mut[k] = len(clone2mut)-1
        else: clone2mut[k] = clone2mut[k][0] 
    mut_tree = nx.relabel_nodes(tree, clone2mut)
    
    mut_names = [f"mut{args.mutations}"] + [f'mut{i}' for i in range(1, args.mutations) ]
    cell_names = [f'cell{i+1}' for i in range(args.cells) ]
    variant_matrix,total_matrix = pd.DataFrame(variant_matrix),pd.DataFrame(total_matrix)
    variant_matrix.index = cell_names
    total_matrix.index = cell_names
    variant_matrix.columns = mut_names
    total_matrix.columns = mut_names
    np.savetxt(f'{args.output}/clonal_matrix.txt', clonal_matrix, fmt='%d')
    np.savetxt(f'{args.output}/usage_matrix.txt', usage_matrix, fmt='%.4f')
    
    variant_matrix.to_csv(f'{args.output}/variant_matrix.csv')
    total_matrix.to_csv(f'{args.output}/total_matrix.csv')
    np.savetxt(f'{args.output}/obs_frequency_matrix.txt', f_hat, fmt='%.4f')
    np.savetxt(f'{args.output}/collapsed_variant_matrix.txt', collapsed_variant_matrix, fmt='%d')
    np.savetxt(f'{args.output}/collapsed_total_matrix.txt', collapsed_total_matrix, fmt='%d')
    nx.write_adjlist(mut_tree, f'{args.output}/tree.txt')
    nx.write_adjlist(cell_tree, f'{args.output}/cell_tree.txt')
    mutation_to_clone_mapping[args.mutations] = mutation_to_clone_mapping[0]
    del mutation_to_clone_mapping[0]
    df = pd.DataFrame.from_dict(mutation_to_clone_mapping, orient='index').reset_index()
    df.columns = ['mutation', 'clone']
    cell_to_clone.to_csv(f'{args.output}/cell_to_clone_mapping.txt', sep=',', index=False)
    df.to_csv(f'{args.output}/mutation_to_clone_mapping.txt', sep=',', index=False)

if __name__ == "__main__":
    main()



