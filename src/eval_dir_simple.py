'''
eval_dir_simple.py

Run simple P/R/F evaluation on converted, filtered graph files
in directory, compared to WordNet gs_taxo files
'''
import os, sys
import numpy as np
import json
from sklearn.metrics import precision_recall_fscore_support

def read_edges(f):
    return [tuple(l.strip().split('\t')[1:]) for l in open(f,'rU').readlines()]

def eval_file(predfile, goldfile):
    pred_edges = set(read_edges(predfile))
    gs_edges = set(read_edges(goldfile))
    
    N = len(gs_edges)
    
    # Combined eval
    if len(pred_edges) > 0:
        p = float(len(pred_edges & gs_edges)) / len(pred_edges)
    else:
        p = 1.
    if len(gs_edges) > 0:
        r = float(len(pred_edges & gs_edges)) / len(gs_edges)
    else:
        r = 1.
    if p + r > 0:
        f = 2 * p * r / (p + r)
    else:
        f = 0.
    
    # Equiv eval
    eq_pred_edges = set([e for e in pred_edges if e[-1]=='equiv'])
    eq_gs_edges = set([e for e in gs_edges if e[-1]=='equiv'])
    if len(eq_pred_edges) > 0:
        p_eq = float(len(eq_pred_edges & eq_gs_edges)) / len(eq_pred_edges)
    else:
        p_eq = 1.
    if len(eq_gs_edges) > 0:
        r_eq = float(len(eq_pred_edges & eq_gs_edges)) / len(eq_gs_edges)
    else:
        r_eq = 1.
    if p_eq + r_eq > 0:
        f_eq = 2 * p_eq * r_eq / (p_eq + r_eq)
    else:
        f_eq = 0.
    
    # Hypernym eval
    hyp_pred_edges = set([e for e in pred_edges if e[-1]=='hypernym'])
    hyp_gs_edges = set([e for e in gs_edges if e[-1]=='hypernym'])
    if len(hyp_pred_edges) > 0:
        p_hyp = float(len(hyp_pred_edges & hyp_gs_edges)) / len(hyp_pred_edges)
    else:
        p_hyp = 1.
    if len(hyp_gs_edges) > 0:
        r_hyp = float(len(hyp_pred_edges & hyp_gs_edges)) / len(hyp_gs_edges)
    else:
        r_hyp = 1.
    if p_hyp + r_hyp > 0:
        f_hyp = 2 * p_hyp * r_hyp / (p_hyp + r_hyp)
    else:
        f_hyp = 0.
    
    res = {'num_gs_edges': N,
           'numTAXedges': len(pred_edges),
           'numIntsct': len(pred_edges & gs_edges),
           'p_hyp': p_hyp,
           'r_hyp': r_hyp,
           'f_hyp': f_hyp,
           'p_eqv': p_eq,
           'r_eqv': r_eq,
           'f_eqv': f_eq,
           'p': p, 
           'r': r,
           'f': f
          }
    print 'RESULTS:', predfile
    print json.dumps(res, indent=2)
    return res


if __name__=="__main__":
    dirname = sys.argv[1]
    gsdirname = sys.argv[2]
    
    predterms = set(['.'.join(f.split('.')[0:2]) for f in os.listdir(dirname)])
    gsterms = set(['.'.join(f.split('.')[1:3]) for f in os.listdir(gsdirname)])
    
    results = {}
    for fname in os.listdir(dirname):
        if not fname.endswith('.converted.filt.taxo'):
            continue
        term = '.'.join(fname.split('.')[0:2])
        if term not in gsterms & predterms:
            continue
        results[term] = eval_file(os.path.join(dirname, fname), 
                                  os.path.join(gsdirname, 
                                               'wn.'+fname.replace('.converted.filt','')))
    
    wghtavg_prf = sum([d['num_gs_edges']*np.array([d['p'],d['r'],d['f']]) for d in results.values()]) / sum([d['num_gs_edges'] for d in results.values()])
    print 'OVERALL: Avg P/R/F:', wghtavg_prf
    eqv_wghtavg_prf = sum([d['num_gs_edges']*np.array([d['p_eqv'],d['r_eqv'],d['f_eqv']]) for d in results.values()]) / sum([d['num_gs_edges'] for d in results.values()])
    print 'EQUIV: Avg P/R/F:', eqv_wghtavg_prf
    hyp_wghtavg_prf = sum([d['num_gs_edges']*np.array([d['p_hyp'],d['r_hyp'],d['f_hyp']]) for d in results.values()]) / sum([d['num_gs_edges'] for d in results.values()])
    print 'HYPER: Avg P/R/F:', hyp_wghtavg_prf
    wghtavg_f = sum([d['num_gs_edges']*d['f'] for d in results.values()]) / sum([d['num_gs_edges'] for d in results.values()])
    print 'Weighted average F-score: %0.6f' % wghtavg_f
