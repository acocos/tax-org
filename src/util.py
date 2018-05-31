##############################
# util.py
##############################
#
# Utilities for loading pairwise relation data from .tsv file (as output by 
# binary relation classification model)

import sys
import json
import pandas as pd
import numpy as np
import networkx as nx
from nltk.corpus import wordnet as wn
from nltk.tag import pos_tag
from nltk.stem import WordNetLemmatizer
import operator
from copy import deepcopy
import itertools
import sqlite3
import pydot
import matplotlib
matplotlib.use("Agg")
from matplotlib import pylab as plt


ltzr = WordNetLemmatizer()

def flatten(l):
    return [item for sublist in l for item in sublist]


def read_terms(f):
    terms = set([line.strip().split('\t')[-1] 
                 for line in open(f, 'rU').readlines()])
    return terms

def update_stats(sold, supdate):
    snew = {}
    for k,v in sold.items():
        snew[k] = v + supdate.get(k, 0)
    return snew

def write_taxo(d, fname):
    if type(fname) == str:
        fout = open(fname,'w')
    elif type(fname) == file:
        fout = fname
    for i, lnk in enumerate(d['links']):
        print >> fout, '\t'.join((str(i), lnk['source'], lnk['target'], lnk['relation']))

def write_subtree(g, node, depth):
    children = g.predecessors(node)
    print '    '*depth+node
    for child in children:
        write_subtree(g, child, depth+1)

def print_hier(g):
    roots = [n for n,d in g.out_degree().items() if d==0]
    __, rts = zip(*sorted([(d, n) for n,d in g.in_degree().items() if n in roots], reverse=True))
    for node in rts:
        write_subtree(g, node, 0)

def draw_hier(g, fn):
#     if not nx.is_directed_acyclic_graph(g):
#         sys.stderr.write('Cannot draw graphs if not a DAG\n')
#         return None
    def draw_subtree(g, node, gph, seen=set()):
        children = g.predecessors(node)
        if len(children) == 0:
            return
        if node in seen:
            return 
        seen.add(node)
        for child in children:
            edge = pydot.Edge(node, child)
            gph.add_edge(edge)
            draw_subtree(g, child, gph, seen)
    gph = pydot.Dot(graph_type='graph')
    roots = [n for n,d in g.out_degree().items() if d==0]
    __, rts = zip(*sorted([(d, n) for n,d in g.in_degree().items() if n in roots], reverse=True))
    for node in rts:
        # check if singleton
        if len(g.predecessors(node)) == 0:
            singleton = pydot.Node(node)
            gph.add_node(singleton)
            continue
        draw_subtree(g, node, gph)
    gph.write_png(fn+'.png')
    gph.write_pdf(fn+'.pdf')

def easy_print(gfile):
    print_hier(consolidate_equiv(transitive_reduction(dct2nx(read_graph(gfile)))))


def draw_graph(g, fname):
    plt.clf()
    nx.draw(g, with_labels=True)
    plt.savefig(fname)

def write_graph(d, fname):
    with open(fname, 'w') as fout:
        print >> fout, json.dumps(d, indent=2)

def write_log(d, stats, fname):  ## TODO: Update as necessary 
    with open(fname, 'w') as fout:
        print >> fout, json.dumps(stats, indent=2)

def read_graph(fname):
    with open(fname, 'rU') as fin:
        d = json.loads(open(fname,'rU').read())
    return d

def read_taxo(fname):
    '''
    Reads .taxo file into JSON graph
    '''
    nodes = set()
    links = set()
    with open(fname,'rU') as fin:
        for line in fin:
            __, w1, w2, rel = line.strip().split('\t')
            nodes.add(w1)
            nodes.add(w2)
            links.add((w1, w2, rel))
    d = {'nodes': [{'name': n} for n in nodes], 
         'links': [{'source': w1, 'target': w2, 'relation': rel, 'weight': 1.0} for (w1,w2,rel) in links]}
    return d

def nx2dct(g):
  '''
  Convert networkx graph to node/link dict
  :param g: networkx Graph or DiGraph
  :returns: dict
  '''
  d = {'nodes': [], 'links': []}
  d['nodes'] = [{'name': n} for n in g.nodes()]
  d['links'] = [{'source': s, 'target': t, 'relation': dat['relation'], 'weight': dat['weight']} for s,t,dat in g.edges(data=True)]
  return d


def dct2nx(d, multi=False):
  '''
  Convert node/link dict to networkx graph
  :param d: dict
  :returns: networkx DiGraph
  '''
  if multi:
      g = nx.MultiDiGraph()
  else:
      g = nx.DiGraph()
  nodenames = [n['name'] for n in d['nodes']]
  g.add_nodes_from(nodenames)
  edges = [(l['source'], l['target'], {'relation': l['relation'], 'weight': l['weight']}) for l in d['links']]
  g.add_edges_from(edges)
  return g

def transitive_reduction(g):
    if type(g) == dict:
        g = dct2nx(g)
    consol = consolidate_equiv(g)
    if not nx.is_directed_acyclic_graph(consol):
        sys.stderr.write('Transitive reduction only uniquely defined on directed acyclic graphs\n')
#         return None
    TR = nx.DiGraph()
    TR.add_nodes_from(consol.nodes())
    for u in consol:
        u_edges = set(consol[u])
        for v in consol[u]:
            u_edges -= {y for x,y in nx.dfs_edges(consol,v)}
        TR.add_edges_from((u,v) for v in u_edges)
    rels = {(w1, w2): d['relation'] for w1, w2, d in consol.edges(data=True) if (w1,w2) in TR.edges()}
    wgts = {(w1, w2): d['weight'] for w1, w2, d in consol.edges(data=True) if (w1,w2) in TR.edges()}
    nx.set_edge_attributes(TR, 'relation', rels)
    nx.set_edge_attributes(TR, 'weight', wgts)
    return TR
        
    
def norm(term):
    def haspos(t) : return t[0] in ['N', 'V', 'J', 'R']
    def mappos(t) : 
        if not haspos(t) : return wn.NOUN
        else : 
            if t.startswith('N') : return wn.NOUN   
            if t.startswith('V') : return wn.VERB   
            if t.startswith('R') : return wn.ADV    
            if t.startswith('J') : return wn.ADJ
    toks = term.split()
    if len(toks)==1:  ## TODO: Update when we do other-than-nouns
        pos = 'NN'
        postags=['NN']
    else:
        postags = [t[-1][:2] for t in pos_tag(toks)] ## only take first-two letters (i.e. NN instead of NNS)
        pos = ' '.join(postags)
    lem = [ltzr.lemmatize(tok.decode('utf8'),mappos(tag)).encode('utf8') 
           for tok,tag in zip(toks, postags)]
    lems = ' '.join(lem)
    return lems

def direct_hypernyms(w, pos):
    syns = wn.synsets(w,pos)
    h = set([])
    for syn in syns:
        h |= set(flatten([[l.name() for l in hyp.lemmas()] for hyp in syn.hypernyms()]))
    h -= set([w])
    return h
    
def all_hypernyms(w, pos):
    h = set([])
    hyper = lambda s: s.hypernyms()
    for s in wn.synsets(w,pos):
        h |= set(flatten([[l.name() for l in hyp.lemmas()] for hyp in s.closure(hyper)]))
    h -= set(w)
    return h

def all_equiv(w, pos):
    equiv = set(flatten([[l.name() for l in syn.lemmas()] for syn in wn.synsets(w,pos)]))
    equiv -= set(w)
    return equiv

def build_wn_graph(intsct, pos):
    W = nx.DiGraph()
    W.add_nodes_from(intsct)
    nodes = set(W.nodes())
    for n1 in intsct:
        hyp = all_hypernyms(n1, pos)
        for n2 in hyp & nodes: 
            W.add_edge(n1, n2, {'relation': 'hypernym'})
        equiv = all_equiv(n1, pos)
        for n2 in equiv & nodes:
            W.add_edge(n1, n2, {'relation': 'equiv'})
    return W

def get_relspecific_graph(terms, relfile, minscore=0., reltypes=['hypernym'], 
                          ppdbcol='PPDB20Score', normalize=True, equivrel='equiv'):
    '''
    Load nodes and edges from .tsv or .db file corresponding to specified target and pos.
    All edges of the specified types between any pair of nodes are returned.
    '''
    d = {'nodes': [], 'links': []} 
    d['nodes'] = [{'name': p} for p in terms]
    if normalize:
        normterms = {}
        for t in terms:
            if norm(t) not in normterms:
                normterms[norm(t)] = []
            normterms[norm(t)].append(t)
    else:
        normterms = {t: [t] for t in terms}
    for nt, tlist in normterms.items():
        if len(tlist) > 1:
            for t1 in tlist:
                for t2 in tlist:
                    if t1==t2: continue
                    d['links'].append({'source': t1, 'target': t2, 
                                       'relation': equivrel, 'weight': 1.})
    format = relfile.split('.')[-1]
    if format =='csv':
        df = pd.read_csv(relfile)
        def match(row):
            return row['w1'] in normterms and row['w2'] in normterms
        df['match'] = df.apply(lambda row: match(row), axis=1)
        df_filt = df.loc[df['match']==True]
        if minscore > 0:
            df_filt = df_filt.loc[df_filt['%s' % ppdbcol] >= minscore]
        for row in df_filt.iterrows():
            dta = row[-1]
            for rel in reltypes:
                wgt = dta[rel]
                for t1 in normterms[dta['w1']]:
                    for t2 in normterms[dta['w2']]:
                        d['links'].append({'source': t, 'target': t2,
                                           'weight': wgt, 'relation': rel})
    elif format=='db':
        conn = sqlite3.connect(relfile)
        c = conn.cursor()
        colnames = [col[1] for col in c.execute("PRAGMA table_info(relations)")] 
        sel = 'SELECT * FROM relations ' \
              'WHERE (w1 IN ("%s")) ' \
              'AND (w2 IN ("%s")) ' \
              'AND ([%s] >= %0.2f)' % ('","'.join(normterms.keys()), '","'.join(normterms.keys()), ppdbcol, minscore)
        for row in c.execute(sel):
            rowdct = dict(zip(colnames, row))
            for rel in reltypes:
                wgt = rowdct[rel]
                for t1 in normterms[rowdct['w1']]:
                    for t2 in normterms[rowdct['w2']]:
                        d['links'].append({'source': t1, 'target': t2,
                                           'weight': wgt, 'relation': rel})
    return d

def add_wn_relations(d, pos='n'):
    '''
    Augment existing graph of predicted relations by adding edges
    for any nodes that appear as direct or transitive relations in WordNet
    '''
    newlinks = []
    nodes = [n['name'] for n in d['nodes']]
    alledges = set([(l['source'],l['target']) for l in d['links']])
    hyp = {w: all_hypernyms(w, pos) for w in nodes}
    eqv = {w: all_equiv(w,pos) for w in nodes}
    for linkdct in d['links']:
        w1 = linkdct['source']
        w2 = linkdct['target']
        rel = linkdct['relation']
        if rel=='hypernym': 
            if w2 in hyp[w1]:
                newlinks.append({'source': w1, 'target': w2, 'relation': rel, 'weight': 1.0})
            else:
                newlinks.append(linkdct)
        elif rel=='equiv':
            if w2 in eqv[w1]:
                newlinks.append({'source': w1, 'target': w2, 'relation': rel, 'weight': 1.0})
            else:
                newlinks.append(linkdct)
        elif rel in ['entailment', 'entmax', 'ppdb2.0scoreAdj']:
            if w2 in hyp[w1] or w2 in eqv[w1]:
                newlinks.append({'source': w1, 'target': w2, 'relation': rel, 'weight': 1.0})
            else:
                newlinks.append(linkdct)
        elif rel=='ppdb2.0score':
            if w2 in hyp[w1] or w2 in eqv[w1]:
                newlinks.append({'source': w1, 'target': w2, 'relation': rel, 'weight': 5.0})
            else:
                newlinks.append(linkdct)
    assert len(newlinks)==len(d['links'])
    d['links'] = newlinks
    return d

def add_noise(d, noise):
    '''Add noise to links in d'''
    if noise <= 0: 
        return d
    newlinks = []
    jitter = 0.
    for linkdct in d['links']:
        wgt = linkdct['weight']
        if np.random.uniform(0,1) <= noise:
            wgt += np.random.normal(0, 0.2)
            wgt = max(0., wgt)
            wgt = min(1., wgt)
            jitter += abs(linkdct['weight']-wgt)
        linkdct['weight'] = wgt
        newlinks.append(linkdct)
    d['links'] = newlinks
    sys.stderr.write('Added noise to links with total jitter %0.3f\n' % jitter)
    return d

def add_root_node(d, reltype='hypernym', wgt=0.1):
    dc = deepcopy(d)
    newlinks = d['links']
    for n in d['nodes']:
        newlinks.append({'source': n['name'], 'relation': reltype, 'target': 'ROOT', 'weight': wgt})
    newnodes = d['nodes']
    newnodes.append({'name': 'ROOT'})
    dc['links'] = newlinks
    dc['nodes'] = newnodes
    return dc

def consolidate_clusters(g, clusfile, resolveweights=np.mean):
    '''
    Consolidate all same-cluster entities into a single node
    :param g: nx DiGraph
    :param clusfile: str
    :param resolveweights: function (e.g. max, min, or mean)
    '''
    clusedges = set()
    with open(clusfile,'rU') as fin:
        for line in fin:    #  w1 \t w2 \t w3 ...
            cluster = line.strip().split('\t')
            cluster = [c for c in cluster if c in g.nodes()]
            if len(cluster) < 2: continue
            orignodes = g.nodes()
            n = len(orignodes)
            m = len(cluster)
            incoming = np.zeros((m,n))  # j->i for i,j in (clus x nodes)
            outgoing = np.zeros((m,n))  # i->j for i,j in (clus x nodes)
            for i,w1 in enumerate(cluster):
                for j,w2 in enumerate(orignodes):
                    if w1==w2: 
                        incoming[i][j] = 1.
                        outgoing[i][j] = 1.
                        continue
                    incoming_wt = g[w2].get(w1, {}).get('weight',-1)
                    if incoming_wt > 0: incoming[i][j] = incoming_wt
                    outgoing_wt = g[w1].get(w2, {}).get('weight',-1)
                    if outgoing_wt > 0: outgoing[i][j] = outgoing_wt
            incoming_agg = resolveweights(incoming, axis=0)
            outgoing_agg = resolveweights(outgoing, axis=0)
            newnode = '|'.join(cluster)
            g.add_node(newnode)
            for w1 in cluster:
                g.remove_node(w1)
            for i, w2 in enumerate(orignodes):
                if w2 in cluster:
                    continue
                if outgoing_agg[i] > 0:
                    g.add_edge(newnode, w2, {'relation': 'hypernym', 'weight': outgoing_agg[i]})
                if incoming_agg[i] > 0:
                    g.add_edge(w2, newnode, {'relation': 'hypernym', 'weight': incoming_agg[i]})
    return g

def add_cluster_edgeweights(d, clusfile, clusweight=0.75):
    '''
    Augment existing graph of predicted relations by adding weight <clusweight> to any
    edge between words belonging to the same cluster
    '''
    clusedges = set()
    with open(clusfile,'rU') as fin:
        for line in fin:    #  w1 \t w2 \t w3 ...
            cluster = line.strip().split('\t')
            for w1 in cluster:
                for w2 in cluster:
                    if w1==w2: continue
                    clusedges.add((w1,w2))
    for lnk in d['links']:
        if (lnk['source'], lnk['target']) in clusedges:
            ## HARD 
            if lnk['relation']=='equiv' or lnk['relation']=='entmax':
                lnk['weight'] = 1.
            ## SOFT
#             lnk['weight'] += clusweight
    return d

def logodds(p):
    if p > 0.9999999:
        return 17.
    if p < 0.0000001:
        return -17
    return np.log(p/(1-p))
    
def consolidate_equiv(g):
    ## Consolidate equivalent nodes
    reltypes = set([d['relation'] for s,t,d in g.edges(data=True)])
    while 'equiv' in reltypes:
        equiv = [(s,t) for s,t,d in g.edges(data=True) if d['relation']=='equiv']
        if len(equiv) == 0:
            break
        mergeme = equiv[0]
        g = merge_nodes(g, mergeme)
        reltypes = set([d['relation'] for s,t,d in g.edges(data=True)])
    return g

def merge_nodes(g, e):
    assert len(e) == 2
    newnode = '|'.join(e)
    g.add_node(newnode)
    newedges = {}
    for n1, n2, d in g.edges(data=True):
        # outgoing edges
        if n1 in e: 
            existingedges = newedges.get((newnode, n2), {d['relation']: 0.})
            if d['relation'] in existingedges:
                existingedges[d['relation']] += d['weight']
            else:
                existingedges[d['relation']] = d['weight']
            newedges[(newnode, n2)] = existingedges
        # incoming edges
        if n2 in e:
            existingedges = newedges.get((n1, newnode), {d['relation']: 0.})
            if d['relation'] in existingedges:
                existingedges[d['relation']] += d['weight']
            else:
                existingedges[d['relation']] = d['weight']
            newedges[(n1, newnode)] = existingedges
    # Take max-weighted relation for each new edge
    allnewedges = []
    for edge, reldict in newedges.items():
        relmax, maxwgt = max(reldict.iteritems(), key=operator.itemgetter(1))
        allnewedges.append((edge[0], edge[1], {'relation': relmax, 'weight': maxwgt}))
    g.add_edges_from(allnewedges)
    for n in e:
        g.remove_node(n)
    return g
            

def expand_equiv_node(g, n, equivweight=0.1):
  '''
  '''
  gnew = g.copy()
  newnodes = set(n.split('|'))
  gnew.add_nodes_from(newnodes)
  for n1, n2, dat in gnew.edges(data=True):
    if n1==n:
      for nnew in newnodes:
        gnew.add_edges_from([(nnew, n2, dat)])
    if n2==n:
      for nnew in newnodes:
        gnew.add_edges_from([(n1, nnew, dat)])
  for n1, n2 in itertools.product(newnodes, newnodes):
    if n1 == n2:
      continue
    gnew.add_edges_from([(n1, n2, {'relation': 'equiv', 'weight': equivweight})])
  gnew.remove_node(n)
  return gnew
  
def expand_equiv_nodes(g):
  '''
  Expand graph by breaking up equivalence nodes (those of type 'word1|word2|word3'), 
  adding 'equiv' edges between all of them, and projecting edges from the original
  condensed node to the expanded ones
  '''
  for n in g.nodes():
    if len(n.split('|')) > 1:
      g = expand_equiv_node(g, n)
  return g

  
