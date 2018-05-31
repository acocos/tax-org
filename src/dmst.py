#!/usr/bin/env python

###############################
# dmst.py
###############################
#
# A directed Maximum-spanning-tree taxonomy induction, operates over
# nodes which are clusters of synonyms

import util

import sys
import time
import os
import itertools
import networkx as nx
import numpy as np
import argparse

from dependency_decoding import chu_liu_edmonds

def flatten(l):
    return [item for sublist in l for item in sublist]

def cc(dgraph):
    '''
    Find connected components of directed graph, when transformed to undirected
    '''
    ugraph = dgraph.to_undirected()
    return [cc for cc in nx.connected_components(ugraph)]

def tax_mst(entities, relfile, clusfile=None, minscore=0., noise=0., wn=False):
    '''
    :param lmda : Lambda hyperparameter (set to negative value to use raw scores)
    '''
    ## Retrieve nodes, and IS-A relation weight between each pair
    d = util.get_relspecific_graph(entities, relfile, minscore=minscore)
    if wn:
        d = util.add_wn_relations(d)
    
    if noise > 0:
        d = util.add_noise(d, noise)
    
    g = util.dct2nx(d)
    
    ## Consolidate same-cluster nodes
    if clusfile is not None:
        g = util.consolidate_clusters(g, clusfile, resolveweights=np.mean)
    
    ## Initialize solution
    g_ = nx.DiGraph()
    g_.add_nodes_from(g)
    stats = {'node_cnt': 0,
             'runtime': 0,
             'keptedge_cnt': 0
            }
    
    ## Decompose
    ccomps = cc(g)
    for v_i in ccomps:
        ## Get subgraph
        g_i = g.subgraph(v_i)
        
        ## Solve MST problem restricted to g_i
        keptedges, g_i_stats = solve_mst(g_i)
        
        g_.add_edges_from(keptedges)
#         g_.add_edges_from([(i,j,{'relation': reltype, 'weight': g[i].get(j,{}).get('weight',-1000)}) 
#                           for i,j in keptedges])
        stats = util.update_stats(stats, g_i_stats)
    
    pruned = util.nx2dct(g_)
    
    # checks
    assert len(g_.edges()) == stats['keptedge_cnt']
    assert len(g_.nodes()) == stats['node_cnt']
    
    return pruned, stats
    
def solve_mst(g):
    start_overall = time.time()
    # Create dummy root node
    nodes = g.nodes()
    gc = g.copy()
    gc.add_node('ROOT')
    for w in g.nodes():
        gc.add_edge(w, 'ROOT', {'weight': 0.1})
    
    # Create matrix
    N = len(gc)
    gcmat = np.zeros((N,N))
    newnodes = ['ROOT'] + nodes
    for i, w1 in enumerate(newnodes):
        for j, w2 in enumerate(newnodes):
            if w1==w2: continue
            gcmat[i][j] = gc[w1].get(w2, {}).get('weight', 0.)
    
    # Solve MST
    heads, score = chu_liu_edmonds(gcmat)
    
    # Keep only edges in MST
    keptedges = []  # [(i,j,{'relation':'hypernym','weight':X.X})]
    for i, headidx in enumerate(heads[1:]):  # skip ROOT
        if headidx == 0:
            continue
        w1 = nodes[i]
        w2 = nodes[headidx-1]
        keptedges.append((w1, w2, g[w1][w2]))
    
    # Calculate stats
    end_overall = time.time()
    stats = {'node_cnt': len(nodes),
             'runtime': end_overall - start_overall, 
             'keptedge_cnt': len(keptedges)
            }
    return keptedges, stats

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='MST Taxonomy Induction.')
    parser.add_argument('-f', '--f', type=str, dest='filename',
                        help='File containing terms to organize')
    parser.add_argument('-r', '--r', type=str, dest='relationfile',
                        help='File containing predicted relations')
    parser.add_argument('-o', '--o', type=str, default='./', dest='outdir',
                        help='Output directory')
    parser.add_argument('-m', '--minscore', type=float, default=0., dest='minscore',
                        help='Minimum PPDB score to establish link between terms')
    parser.add_argument('-L', '--logdir', type=str, default='./logs', dest='logdir',
                        help='Directory for writing logs')
    parser.add_argument('-c', '--clusfile', type=str, default=None, dest='clusfile',
                        help='Cluster file')
    parser.add_argument('-n', '--noise', type=float, default=0.,
                        help='Likelihood of adding noise')
    parser.add_argument('-w', '--wn', action="store_true", help="Set to add WordNet edges=1")
    args = parser.parse_args()
    
    termsfile = args.filename
    relationfile = args.relationfile
    outdir = args.outdir
    minscore = args.minscore
    logdir = args.logdir
    clusfile = args.clusfile
    noise = args.noise
    wn = args.wn
    
    entities = util.read_terms(termsfile)
    
    pruned, stats = tax_mst(entities, relationfile, clusfile=clusfile, minscore=minscore,
                            noise=noise, wn=wn)
    
    bn = os.path.basename(termsfile)
    util.write_taxo(pruned, fname=os.path.join(outdir, bn.replace('.terms','.taxo')))
    util.write_graph(pruned, fname=os.path.join(outdir, bn.replace('.terms','.graph')))
    util.write_log(pruned, stats, fname=os.path.join(logdir, 'mst.terms-%s.minsc%0.1f' % (bn, minscore)))