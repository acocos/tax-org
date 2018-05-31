'''
convert_ti_output_to_formatted.py

Convert raw taxonomies as output by each algorithm into one suitable
for evaluation. This includes:
-- For all methods, exploding clusters (i.e. word1|word2|word3)
   into individual nodes with equiv edges
-- For ILP methods: turning cycles into equiv edges, and non-cycle
   edges into hypernym edges
-- For NoCyc methods: identifying synonyms as nodes sharing
   the same direct hypernyms and hyponyms
-- For all methods, taking the transitive closure of the result
'''
import sys, os
import matplotlib
matplotlib.use("Agg")
from matplotlib import pylab as plt
import networkx as nx

import util

def flatten(l):
    return [item for sublist in l for item in sublist]

def samehh2equiv(g):
    '''
    Consolidate nodes that have the same direct hypernyms and hyponyms
    '''
    gtr = util.transitive_reduction(g)
    newedges = []
    for w1 in gtr.nodes():
        in1 = set([n1 for (n1, n2, d) in gtr.in_edges(w1, data=True) if d['relation']=='hypernym'])
        out1 = set([n2 for (n1, n2, d) in gtr.out_edges(w1, data=True) if d['relation']=='hypernym'])
        for w2 in gtr.nodes():
            if w1>=w2: continue
            in2 = set([n1 for (n1, n2, d) in gtr.in_edges(w2, data=True) if d['relation']=='hypernym'])
            out2 = set([n2 for (n1, n2, d) in gtr.out_edges(w2, data=True) if d['relation']=='hypernym'])
            if in1==in2 and out1==out2:  # equiv!
                if g[w1].get(w2,{}).get('relation', '') == 'equiv': continue  # already there
                newedges.append((w1, w2, {'relation': 'equiv', 'weight': 1.0}))
                newedges.append((w2, w1, {'relation': 'equiv', 'weight': 1.0}))
    g.add_edges_from(newedges)
    return g

def entmax2hyp(g):
    '''
    Simply change the name of 'entmax' edges to 'hypernym'
    '''
    for n1, n2, d in g.edges(data=True):
        if d['relation']=='entmax':
            d['relation'] = 'hypernym'
    return g
    
    
def scc2equiv(g):
    '''
    Change the relation type of any edges in a strongly connected component to 'equiv'
    :param g: nx.MultiDiGraph
    '''
    for scc in nx.strongly_connected_components(g):
        if len(scc) < 2:
            continue
        for n1, n2 in g.edges():
            if n1 in scc and n2 in scc:
                try:
                    g[n1][n2][0]['relation'] = 'equiv' # nx.MultiDiGraph
                except KeyError:
                    g[n1][n2]['relation'] = 'equiv' # nx.DiGraph
    return g

def filter_nodes(G, filt):
    '''
    If the list filt has node names in it, filter G to be only
    the subgraph including those nodes. If the filt list is empty,
    do nothing
    '''
    if len(filt)==0:
        return G
    filt = set(filt)
    for n in G.nodes():
        if n not in filt:
            G.remove_node(n)
    return G

def hypernym_transitive_closure(G):
    hyperedges = [e for e in G.edges(data=True) if e[-1]['relation']=='hypernym']
    hyperedgenodes = set(flatten([(n1,n2) for n1,n2,__ in hyperedges]))
    hypergraph = nx.DiGraph(G.subgraph(hyperedgenodes))
    tc = nx.transitive_closure(hypergraph)
    for n1, n2 in set(tc.edges()) - set(hypergraph.edges()):
        G.add_edges_from([(n1, n2, {'relation': 'hypernym'})])
    return G

def convert_graph(graphfile, clusmethod, filt=[]):
    assert clusmethod in ['filteronly','mtg','mtf','nocyc','dmst']
    
    G = util.dct2nx(util.read_graph(graphfile))
    
    # Step 0: If just filtering, do that
    if clusmethod=='filteronly':
        G = filter_nodes(G, filt)
        return G

    # Step 1: Find equivalent terms based on method
    if clusmethod in ['mtg','mtf']:
        G = entmax2hyp(G)
        G = scc2equiv(G)
    
    if clusmethod=='nocyc':
        # Find nodes with same direct hypernyms and hyponyms
        G = samehh2equiv(G)
    
    # Step 2: Collapse equivalent nodes
    G = util.consolidate_equiv(G)
    
    # Step 3: Find transitive closure
    G = hypernym_transitive_closure(G)
    
    # Step 4: Explode equiv nodes again
    G = util.expand_equiv_nodes(G)
    
    # Step 5: Remove nodes not in filter
    G = filter_nodes(G, filt)
    
    return G


if __name__=="__main__":
    graphfile = sys.argv[1]
    clusmethod = sys.argv[2]
    outfile = sys.argv[3]
    if len(sys.argv) > 4:
        filtterms = util.read_terms(sys.argv[4])
    else:
        filtterms = []
    
    newgraph = convert_graph(graphfile, clusmethod, filt=filtterms)
    newtaxo = util.nx2dct(newgraph)
    
    util.write_taxo(newtaxo, outfile)
    
