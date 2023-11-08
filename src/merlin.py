import numpy as np
import pandas as pd
import itertools,argparse
import gurobipy as gp
import networkx as nx
from scipy.stats import linregress, binom
# potentially unused
from networkx.algorithms import tournament
from networkx.drawing.nx_pydot import graphviz_layout
import seaborn as sns
import matplotlib.pyplot as plt
import pickle as pkl

np.set_printoptions(formatter={'float': lambda x: format(x, '.5f')})

def set_cell_names(df_total,df_variant,namelist):
    if namelist != None:
        cell_list = []
        with open(namelist, 'r') as inp:
            for line in inp:
                cell_list.append(line.rstrip('\n'))
    else:
        cell_list = [f'cell_{i}' for i in range(len(df_total.columns))]
    df_variant.columns = cell_list
    df_total.columns = cell_list
    return df_total,df_variant,cell_list
  # build ancestry graph
def get_ancestry_clustered_data(df_variant, df_total, read_threshold = 25, vaf_threshold = 0.05):
    '''
    @params:
      - df_variant: mutation by cell containing variant read count
      - df_total: mutation by cell df containing total read count
      
    @return:
      - a networkx DiGraph containing all the ancestral relationships
    '''
    G = nx.DiGraph()
    for mut in df_total.index:
        G.add_node(mut)
    
    mutation_list = list(df_total.index)
    forbidden_disjoint_pairs = []
    ncells = len(df_total.columns)
    
    for mut1_idx, mut2_idx in itertools.combinations(range(len(df_total)), 2):
        mut1 = mutation_list[mut1_idx]
        mut2 = mutation_list[mut2_idx]
        
        series1 = df_total.loc[mut1].to_numpy()
        series2 = df_total.loc[mut2].to_numpy()
        
        mut1_vars = df_variant.loc[mut1].to_numpy()
        mut2_vars = df_variant.loc[mut2].to_numpy()
        mut1_tot = df_total.loc[mut1].to_numpy()
        mut2_tot = df_total.loc[mut2].to_numpy()
        
        vaf1 = mut1_vars / mut1_tot
        vaf2 = mut2_vars / mut2_tot
        # print(vaf1,vaf2)
       
        #print(vaf1, vaf2)

        vaf1[np.isnan(vaf1)] = 0
        vaf2[np.isnan(vaf2)] = 0
        
        good1 = (mut1_tot >= read_threshold) & (vaf1 >= vaf_threshold)
        good2 = (mut2_tot >= read_threshold) & (vaf2 >= vaf_threshold)
        good = np.logical_and(good1, good2)
        
        vaf1 = vaf1[good]
        vaf2 = vaf2[good]
        
        vaf1[np.isnan(vaf1)] = 0
        vaf2[np.isnan(vaf2)] = 0
        
        slope = np.dot(vaf1, vaf2) / np.dot(vaf1, vaf1)

        #print(mut1, mut2, slope)

        p = 0
        
        # print(mut1, mut2, slope, sep='\t')
        # slope, _, _, p, _ = linregress(vaf1, vaf2)
        
        if (np.tan(45 * np.pi / 180) > slope) and (slope > 0) and p <= 0.01:
            G.add_edge(mut1, mut2)
            # print(mut1, mut2, slope, p, sep='\t')
        elif np.tan(45 * np.pi / 180) < slope and p <= 0.01:
            G.add_edge(mut2, mut1)
            # print(mut1, mut2, slope, p, sep='\t')
        
    
    G.add_node('root')
    for node in G.nodes:
        if node != 'root':
            if len(G.in_edges(node)) == 0:
                G.add_edge('root', node)
    
    return G

# clustering cells
  # something something about cloen assignment???

# ILP formulation
def set_ILP_formuation(G_ancestry, cell2group, df_cluster_variant, df_cluster_total, prefix, minU=0.05):
    
    nmutations = len(G_ancestry.nodes) - 1
    nclones = nmutations
    ngroups = len(set(cell2group))
    ncells = len(df_cluster_total.columns)

    mutation_list = list(df_cluster_total.index)
    mutation2index = {mutation:idx for idx, mutation in enumerate(mutation_list)}

    freq = (df_cluster_variant / df_cluster_total).values.T

    ancestry_edge_list = list(G_ancestry.edges)
    nedges = len(ancestry_edge_list)
    edge2index = {edge:idx for idx, edge in enumerate(ancestry_edge_list)}

   
    model = gp.Model('phyloMito')
    model.setParam('TimeLimit', 2*60*60)
    x = model.addVars(nedges, vtype=gp.GRB.BINARY, name='x')
    # support of mixture matrix
    a = model.addVars(ngroups, nmutations, vtype=gp.GRB.BINARY, name='a')
    # cell-level frequencies and proportions
    u = model.addVars(ncells, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='u')
    f = model.addVars(ncells, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='f')
    # edge contribution of frequency
    fe = model.addVars(ncells, nedges, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='fe')
    # perfect phylogeny variables
    ya = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='ya')
    za = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='za')
    # L1 objective
    error = model.addVars(ncells, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='error')


    # edge sum
    #def add_edge_sum_constraints():
    for node in G_ancestry.nodes:
        curr_in_edges = G_ancestry.in_edges(node)
        edge_sum = gp.LinExpr()
        for edge in curr_in_edges:
            edge_sum += x[edge2index[edge]]
        
        if len(curr_in_edges) > 0:
            model.addConstr(edge_sum == 1)


    # U matrix
    #def add_other_constraints():
    for cell_idx in range(ncells):
        for mut_idx in range(nmutations):
            group_idx = cell2group[cell_idx]
            model.addConstr(u[cell_idx, mut_idx] >= minU * a[group_idx, mut_idx])
            model.addConstr(u[cell_idx, mut_idx] <= a[group_idx, mut_idx])

    for cell_idx in range(ncells):
        usum = gp.LinExpr()
        for mut_idx in range(nmutations):
            usum += u[cell_idx, mut_idx]
        model.addConstr(usum <= 1)
    

    # concordance between clonal tree and mixture matrix
    #def concordance_constraint():
    for mut_idx in range(nmutations):
        curr_out_edges = G_ancestry.out_edges(mutation_list[mut_idx])
        for edge in curr_out_edges:
            target_mut = edge[1]
            target_mut_idx = mutation2index[target_mut]
            edge_idx = edge2index[edge]
            for group_idx in range(ngroups):
                model.addConstr(a[group_idx, mut_idx] >= x[edge_idx] + a[group_idx, target_mut_idx] - 1)
    

    #def frequency_matrix():
    for cell_idx in range(ncells):
        for mut_idx in range(nmutations):
            curr_mut = mutation_list[mut_idx]
            
            fusum = gp.LinExpr()
            fusum += u[cell_idx, mut_idx]
            
            curr_out_edges = G_ancestry.out_edges(curr_mut)            
            for edge in curr_out_edges:
                edge_idx = edge2index[edge]
                target_mut_idx = mutation2index[edge[1]]
                
                model.addConstr(fe[cell_idx, edge_idx] <= f[cell_idx, target_mut_idx])
                model.addConstr(fe[cell_idx, edge_idx] <= x[edge_idx])
                model.addConstr(fe[cell_idx, edge_idx] >= x[edge_idx] + f[cell_idx, target_mut_idx] - 1)
            
                fusum += fe[cell_idx, edge_idx]
            
            model.addConstr(f[cell_idx, mut_idx] == fusum)
          

    #def perfect_phylogeny_MM():
    for mut1_idx, mut2_idx in itertools.combinations(range(nmutations), 2):
        
        for group_idx in range(ngroups):
            model.addConstr(a[group_idx, mut1_idx] + a[group_idx, mut2_idx] <= 2 - za[mut1_idx, mut2_idx])
        
        for group_idx in range(ngroups):
            model.addConstr(a[group_idx, mut1_idx] <= a[group_idx, mut2_idx] + (1 - ya[mut1_idx, mut2_idx]))    
        
        for group_idx in range(ngroups):
            model.addConstr(a[group_idx, mut2_idx] <= a[group_idx, mut1_idx] + (1 - ya[mut2_idx, mut1_idx]))
            
        for group_idx in range(ngroups):
            model.addConstr(za[mut1_idx, mut2_idx] + ya[mut1_idx, mut2_idx] + ya[mut2_idx, mut1_idx] >= 1)
    

    #def set_objective_fn():
    obj_sum = gp.LinExpr()
    for cell_idx in range(ncells):    
        for mut_idx in range(nmutations):
            if not np.isnan(freq[cell_idx][mut_idx]): # freq = some vaf stuff
                model.addConstr(error[cell_idx, mut_idx] >= freq[cell_idx][mut_idx] - f[cell_idx, mut_idx])
                model.addConstr(error[cell_idx, mut_idx] >= f[cell_idx, mut_idx] - freq[cell_idx][mut_idx])
                obj_sum += error[cell_idx, mut_idx]
    
    model.setObjective(obj_sum,gp.GRB.MINIMIZE)
    

    model.setParam(gp.GRB.Param.Threads, 1)
    model.setParam(gp.GRB.Param.Method, 4)
    model.optimize()

    # write output

    ancestry_edge_list = list(G_ancestry.edges)
    with open(f'{prefix}_clone_tree_edge_list.txt', 'w') as out:
        for edge_idx, val in model.getAttr('x', x).items():
            if val > 0.5:
                out.write(f'{ancestry_edge_list[edge_idx][0]}, {ancestry_edge_list[edge_idx][1]}\n')

    df_a = pd.DataFrame(np.reshape(model.getAttr('x', a).values(), (ngroups, nmutations)).astype(int),
                        columns=df_cluster_total.index, index=df_cluster_total.columns)
    df_a.to_csv(f'{prefix}_Amatrix.csv')
    df_u = pd.DataFrame(np.reshape(model.getAttr('x', u).values(), (ncells, nmutations)),
                        columns=df_cluster_total.index,index=df_cluster_total.columns)
    df_u.to_csv(f'{prefix}_Umatrix.csv')
    with open(f'{prefix}_ancestry_edge_list.txt', 'w') as out:
        for edge in G_ancestry.edges:
            out.write(f"{edge[0]}, {edge[1]}\n")
    
  
def main():
    
    parser = argparse.ArgumentParser(description='Simulate a clonal matrix, usage matrix and read counts.')
    parser.add_argument('-t','--total', type=str, required=True, help='mutation x cell csv file for total read count')
    parser.add_argument('-v','--variant', type=str, required=True, help='mutation x cell csv file for variant read count')
    parser.add_argument('-o','--out', type=str, required=True, help='prefix for output files')
    parser.add_argument('--cell_list', type=str, required=False, help='list of cell name (optional)')
    parser.add_argument('--cell_groups',type=str,required=False, help='clustering of cells, ')
    args = parser.parse_args()

    df_variant = pd.read_csv(args.variant, header=0, index_col = 0).T.astype(int)
    df_total = pd.read_csv(args.total, header=0, index_col = 0).T.astype(int)

    G_ancestry = get_ancestry_clustered_data(df_variant, df_total)
    
    cell2group = {i:i for i in range(len(df_total.columns))}
    set_ILP_formuation(G_ancestry, cell2group, df_variant, df_total, args.out)
if __name__ == "__main__":
    main()
