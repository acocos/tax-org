###############################
# maxtransforest.py
###############################
#
# Max-Trans-Forest method for taxonomy induction from the paper:
# @article{berant2015efficient,
#   title={Efficient global learning of entailment graphs},
#   author={Berant, Jonathan and Alon, Noga and Dagan, Ido and Goldberger, Jacob},
#   journal={Computational Linguistics},
#   year={2015},
#   publisher={MIT Press}
# }

import util

import time
import sys
import os
import itertools
import networkx as nx
import numpy as np
import argparse

from gurobipy import *


def flatten(l):
    return [item for sublist in l for item in sublist]

def cc(dgraph):
    '''
    Find connected components of directed graph, when transformed to undirected
    '''
    ugraph = nx.DiGraph()
    ugraph.add_nodes_from(dgraph)
    ugraph.add_edges_from([e for e in dgraph.edges(data=True) if e[-1]['weight'] > 0])
    ugraph = ugraph.to_undirected()
    return [cc for cc in nx.connected_components(ugraph)]


def mtf(entities, relfile, reltype='hypernym', minscore=1.5, lmda=1.0,
        noise=0., wn=False):
    '''
    :param lmda : Lambda hyperparameter (set to negative value to use raw scores)
    '''
    ## Retrieve nodes, and reltype weight between each pair
    d = util.get_relspecific_graph(entities, relfile, minscore=minscore, 
                                   reltypes=[reltype], equivrel=reltype)
    
    ## Convert weights to logodds and subtract lambda
    for lnk in d['links']:
        lnk['weight'] = lnk['weight'] - lmda

    if wn:
        d = util.add_wn_relations(d)
    
    if noise > 0:
        d = util.add_noise(d, noise)
    
    g = util.dct2nx(d)
    
    ## Initialize solution
    g_ = nx.DiGraph()
    g_.add_nodes_from(g)
    stats = {'node_cnt': 0,
             'num_vars': 0,
             'num_constrs': 0,
             'runtime': 0,
             'keptedge_cnt': 0,
             'possedge_cnt': 0,
             'gt0edge_cnt': 0,
             'timeout': 0
            }
    
    ## Decompose
    ccomps = cc(g)
    for v_i in ccomps:
        ## Get subgraph
        g_i = g.subgraph(v_i)
        
        ## Solve ILP problem restricted to g_i
        keptedges, g_i_stats = solve_exact(g_i)
        
        g_.add_edges_from([(i,j,{'relation': reltype, 'weight': g[i].get(j,{}).get('weight',0)}) 
                          for i,j in keptedges])
        stats = util.update_stats(stats, g_i_stats)
    
    pruned = util.nx2dct(g_)
    
    # checks
    assert len(g_.edges()) == stats['keptedge_cnt']
    assert len(g_.nodes()) == stats['node_cnt']
    
    return pruned, stats
    

def solve_exact(g_i):
    # Generate list of edge variables and weights
    possedges = [(i,j) for (i,j) in itertools.product(g_i.nodes(), g_i.nodes()) if i!=j]
    weights = {(i,j): g_i[i].get(j,{}).get('weight',0.) for i,j in possedges}
    
    # Initialize model
    m = Model('taxonomy')
    timelimit = 4*60
    overalllimit = 5*60
    m.setParam('TimeLimit', timelimit)
    x = m.addVars(possedges, vtype=GRB.BINARY, name='x')
    m.setObjective(x.prod(weights), GRB.MAXIMIZE)
    
    all_triples = itertools.permutations(g_i.nodes(), 3)
    trans_constr = m.addConstrs((x[i,j]+x[j,k]-x[i,k] <= 1 
                                 for i,j,k in all_triples),
                                 name='trans_c')
    
    # Solve ILP using active constraints
    iter = 0
    timedout = False
    start_overall = time.time()
    while True:
        iter += 1
        if iter >= 10000:
            sys.stderr.write('10000 iterations reached\n')
            break
        start_iter = time.time()
        m.optimize()
        end_iter = time.time()
        elapsed = end_iter - start_iter
        if elapsed >= timelimit:
            timedout = True
        totaltime = end_iter - start_overall
        if totaltime >= overalllimit:
            timedout = True
            break
        # Get violated constraints
        violated_trans_constraints = get_viol_trans_constraints(x)
        violated_tree_constraints = get_viol_tree_constraints(x)
        if len(violated_trans_constraints | violated_tree_constraints) == 0:
            break
        trans_constr = m.addConstrs((x[i,j]+x[j,k]-x[i,k] <= 1 
                                     for i,j,k in violated_trans_constraints), 
                                    name='trans_c')
        tree_constr = m.addConstrs((x[i,j]+x[i,k]-x[j,k]-x[k,j] <= 1 
                                    for i,j,k in violated_tree_constraints), 
                                   name='tree_c')
    end_overall = time.time()
    elapsed_overall = end_overall - start_overall
    
    # Return result
    keptedges = [(i,j) for i,j in possedges if x[i,j].x > 0]
    sys.stderr.write('Kept %d of %d possible edges in %d iter\n' 
                     % (len(keptedges), len(possedges), iter))
    # Log
    comp_stats = {'node_cnt': len(g_i.nodes()),
                  'num_vars': m.numVars,
                  'num_constrs': m.numConstrs,
                  'runtime': elapsed_overall,
                  'keptedge_cnt': len(keptedges),
                  'possedge_cnt': len(possedges),
                  'gt0edge_cnt': len([k for k,v in weights.items() if v > 0]),
                  'timeout': int(timedout)
                 }
    return keptedges, comp_stats

def get_viol_trans_constraints(x):
    '''
    '''
    viol = set([])
    edges = {k: v for k,v in x.items() if v.x > 0}
    
    for i,j in edges: 
        for k in [kk for (jj,kk) in edges if jj==j]:
            if k==i:
                continue
            if (i,k) not in edges:
                viol.add((i,j,k))
    return viol

def get_viol_tree_constraints(x):
    '''
    '''
    viol = set([])
    edges = {k: v for k,v in x.items() if v.x > 0}
    nodes = set(flatten(edges.keys()))
    
    for i,j in edges:
        for k in [kk for (ii,kk) in edges if ii==i]:
            if k==j:
                continue
            if (k,j) not in edges:
                if (j,k) not in edges:
                    viol.add((i,j,k))
    return viol
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Max-Trans-Forest with a single relationtype.')
    parser.add_argument('-f', '--f', type=str, dest='filename',
                        help='File containing terms to organize')
    parser.add_argument('-r', '--r', type=str, dest='relationfile',
                        help='File containing predicted relations')
    parser.add_argument('-t', '--t', type=str, default='hypernym', dest='relationtype',
                        help='Relation type')
    parser.add_argument('-o', '--o', type=str, default='./', dest='outdir',
                        help='Output directory')
    parser.add_argument('-m', '--minscore', type=float, default=0., dest='minscore',
                        help='Minimum PPDB score to establish link between terms')
    parser.add_argument('-l', '--lambda', type=float, default=1.0, dest='p_lambda',
                        help='Lambda hyperparameter')
    parser.add_argument('-L', '--logdir', type=str, default='./logs', dest='logdir',
                        help='Directory for writing logs')
    parser.add_argument('-n', '--noise', type=float, default=0., dest='noise',
                        help='Likelihood of adding noise to each weight')
    parser.add_argument('-w', '--wn', action="store_true", help="Set to add WordNet edges=1")
    args = parser.parse_args()
    
    termsfile = args.filename
    relationfile = args.relationfile
    relationtype = args.relationtype
    outdir = args.outdir
    minscore = args.minscore
    p_lambda = args.p_lambda
    logdir = args.logdir
    noise = args.noise
    wn = args.wn
    
    entities = util.read_terms(termsfile)
    
    pruned, stats = mtf(entities, relationfile, reltype=relationtype, wn=wn,
                        minscore=minscore, lmda=p_lambda, noise=noise)
    
    bn = os.path.basename(termsfile)
    util.write_taxo(pruned, fname=os.path.join(outdir, bn.replace('.terms','.taxo')))
    util.write_graph(pruned, fname=os.path.join(outdir, bn.replace('.terms','.graph')))
    util.write_log(pruned, stats, fname=os.path.join(logdir, 'mtf.terms-%s.rel-%s.minsc%0.1f.lambda-%0.1f' % (bn, relationtype, minscore, p_lambda)))
  

    