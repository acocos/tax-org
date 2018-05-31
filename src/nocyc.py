#!/usr/bin/env python

###############################
# tax_nocyc.py
###############################
#
# A cycle-breaking method for taxonomy organization, operates over
# nodes which are clusters of synonyms

import util

import sys
import time
import os
import itertools
import networkx as nx
import numpy as np
import argparse

def flatten(l):
    return [item for sublist in l for item in sublist]

def prune_cycle(g, cyc):
    minweight = 1000
    minedge = None
    for n1, n2 in itertools.product(cyc, cyc):
        if (n1, n2) in g.edges():
            wgt = g.edge[n1][n2]['weight']
            if wgt < minweight:
                minweight = wgt
                minedge = (n1, n2)
    g.remove_edge(*minedge)
    return g

def tax_nocyc(entities, relfile, clusfile=None, minscore=0., noise=0., lmda=0.5, wn=False):
    '''
    :param lmda : Lambda hyperparameter (set to negative value to use raw scores)
    '''
    ## Retrieve nodes, and IS-A relation weight between each pair
    d = util.get_relspecific_graph(entities, relfile, minscore=minscore)
    
    if wn:
        d = util.add_wn_relations(d)
    
    if noise > 0:
        d = util.add_noise(d, noise)
    
    for lnk in d['links']:
        lnk['weight'] = lnk['weight'] - lmda
    
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
    
    keptedges, g_stats = solve_nocyc(g)
    g_.add_edges_from(keptedges)
    
    stats = util.update_stats(stats, g_stats)
    
    pruned = util.nx2dct(g_)
    
    # checks
    assert len(g_.edges()) == stats['keptedge_cnt']
    assert len(g_.nodes()) == stats['node_cnt']
    
    return pruned, stats
    
def solve_nocyc(g):
    start_overall = time.time()
    
    # Prune negative edges
    nodes = g.nodes()
    gc = nx.DiGraph()
    gc.add_nodes_from(nodes)
    gc.add_edges_from([e for e in g.edges(data=True) if e[-1]['weight'] > 0])
    
    # Iteratively break cycles, longest first
    cnt = 0
    cycles = [c for c in nx.strongly_connected_components(gc) if len(c) > 1]
    while len(cycles) > 0:
        orderedcyc = sorted(cycles, key=len)
        cyc = orderedcyc[-1]
        gc = prune_cycle(gc, cyc)
        cycles = [c for c in nx.strongly_connected_components(gc) if len(c) > 1] # Tarjan
        cnt += 1
        if cnt >= 5000:
            sys.stderr.write('More than 5000 pruning iterations\n')
            break
    
    keptedges = gc.edges(data=True)
    
    # Calculate stats
    end_overall = time.time()
    stats = {'node_cnt': len(nodes),
             'runtime': end_overall - start_overall, 
             'keptedge_cnt': len(keptedges)
            }
    return keptedges, stats

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='No-Cycle Taxonomy Induction.')
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
    parser.add_argument('-l', '--lambda', type=float, default=0.5, dest='p_lambda',
                        help='Lambda hyperparamter')
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
    p_lambda = args.p_lambda
    wn = args.wn
    
    entities = util.read_terms(termsfile)
    
    pruned, stats = tax_nocyc(entities, relationfile, clusfile=clusfile, minscore=minscore,
                              noise=noise, lmda=p_lambda, wn=wn)
    
    bn = os.path.basename(termsfile)
    util.write_taxo(pruned, fname=os.path.join(outdir, bn.replace('.terms','.taxo')))
    util.write_graph(pruned, fname=os.path.join(outdir, bn.replace('.terms','.graph')))
    util.write_log(pruned, stats, fname=os.path.join(logdir, 'nocyc.terms-%s.minsc%0.1f' % (bn, minscore)))