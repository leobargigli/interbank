from scipy.stats import binom , kendalltau, pearsonr, norm
import networkx as nx
from scipy.sparse import csc_matrix, extract, linalg, isspmatrix_csc,dia_matrix
from scipy.cluster.hierarchy import linkage,fcluster,distance
from itertools import cycle
from cPickle import load
#from sparsesvd import sparsesvd
import numpy as np
import pylab as plt
import os

    
     
def spectral_partition(W,q,method = 'complete', metric = 'cosine'):

    n,m = W.shape
    K = Kmatrix(W).todense()
    r = min(n,m)

    if n == m:
        e,v = linalg.eigen(K, r - 2)

    else:
        u,e,v = linalg.svd(K, r - 1)
        v = np.concatenate((u, v.T), 0)
                
    max_index = e.argmax()
    v = np.delete(v,max_index,1)
    Obs = np.real(v)
    D = distance.pdist(Obs,metric = metric)
    D = np.multiply(D >= 0, D)
    Z = linkage(D, method = method, metric = metric)
    cluster = fcluster(Z, q, criterion = 'maxclust')
            
    cluster += - 1
    cluster = {'spectral' : cluster}

    return cluster

def n_of_modules_norm(W , mean , var, normed = False):

    n, m = W.shape
    r = min(n, m)
    M = { True: Kmatrix(W).todense(), False: W }
    svs =  np.linalg.svd(M[normed],0, 0) 
    eigs = np.power(svs,2)
    eigs = eigs[1:] # remove the unitary sv
    p = np.zeros((r - 1,))

    for j in xrange(r - 1):
        sigma = eigs.sum() 
        p[j]  = 1 - norm.cdf(sigma,mean,np.sqrt(var)) #min(exp_value / sigma,1.)
        eigs = eigs[1:] # remove the largest sv

    pdiff = p[1:] - p[:-1]

    try: 
        delta  = pdiff.argmax() 
       
    except ValueError:
        delta = -1
    
    if normed is False:
        q = 2 + delta # the first addend is explained as follows: there is at least one module & python index starts at 0

    return q, pdiff,p, svs

def n_of_modules_markov(W):

    K = Kmatrix(W).todense()
    v = W.sum()
    n, m = K.shape
    r = min(n, m)
    exp_value = (n * m) / float(v)
    svs =   np.linalg.svd(K,0, 0) 
    eigs = np.power(svs,2)
    eigs = eigs[1:] # remove the unitary sv
    p = np.zeros((r - 1,))

    for j in xrange(r - 1):
        sigma = eigs.sum() 
        p[j]  = min(exp_value / sigma,1.)
        eigs = eigs[1:] # remove the largest sv

    pdiff = p[1:] - p[:-1]

    try: 
        delta  = pdiff.argmax() 
       
    except ValueError:
        delta = -1
    
    q = 2 + delta # the first addend is explained as follows: there is at least one module & python index starts at 0

    return q, pdiff,p, svs


def plot_svs(svs,pdiff,q,filename, loc = (0.6,0.7) ,adj = 0.1, fmt = 'png'):

    eigs = svs ** 2
    ymax = eigs.max()
    rank = len(svs)

    fig = plt.figure(figsize = (8.,6.))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    line_labels =('$\sigma^2$','$\Delta$ p-values')

    ax1.set_ylabel(line_labels[1])
    ax2.set_ylabel(line_labels[0], rotation = 'horizontal')
    ax2.set_ylim((10**-3,ymax))

                    
    line1 = ax2.plot(eigs, 'red' , drawstyle = 'steps-mid')
    line2 = ax1.plot(range(2,rank), pdiff, 'blue' , drawstyle = 'steps-mid')
    line3 = plt.Line2D([q,q],[0.,ymax], color='black' , linewidth=1.)
    ax1.add_line(line3)
    ax2.set_yscale('log')
    ax1.set_ylim((0,pdiff.max() + 0.01))
    fig.legend((line1,line2),line_labels,loc = loc)
    filename = filename + 'n_of_mod'
    plt.savefig(filename + '.' + fmt , format = fmt )



def Kmatrix(W):

    n, m = W.shape
    out_strength = W.sum(1)
    empty = np.where(out_strength == 0)
    out_strength[empty] = 1.

    in_strength = W.sum(0)
    empty = np.where(in_strength == 0)
    in_strength[empty] = 1.

    out_strength = 1.0 / np.sqrt(out_strength)
    W_out = dia_matrix((out_strength.T,0),shape = (n,n))

    in_strength = 1.0 / np.sqrt(in_strength)
    W_in = dia_matrix((in_strength,0),shape = (m,m))

    K = W_out * W * W_in
    K = csc_matrix(K)

    return K

def community_matrix(cluster):

    n = len(cluster)
    cluster.resize((n, 1))
    clu_unique = np.unique(cluster)
    c = len(clu_unique)
    clu_unique.resize((c,1))
    M = np.repeat(cluster,c,1) == np.repeat(clu_unique.T,n,0)
    M = 1. * np.asmatrix(M)
    
    return M

def makecover(svnet, M):
    old = 0 
    new = M.sum()
    n, m = svnet.shape
    pvalue = 0.01
    while old <> new:
        old = M.sum()
        NTC,  CTC = community_matrices(svnet, M)
        P = prob_node_vs_community(CTC, NTC, M)
        M,pvalue = extend_community(M,P,pvalue = pvalue)
        new = M.sum()
        print new
    return M
    
  
def prob_node_vs_community(CTC,NTC, M):

    I = np.diag(np.ones((len(CTC), )))
    CTC = np.multiply(CTC, 1 - I)
    NTC = np.multiply(NTC, 1 - M)
    ext_in_degree = CTC.sum(0)
    v = NTC.sum()
    ext_out_degree = NTC.sum(1)
    P = ext_out_degree * ext_in_degree / v**2
    P = 1 - binom.cdf(NTC-1,v,P)      

    return P

def community_matrices(A, M):

    NTC = A * M 
    CTC = M.T * A * M

    return NTC,  CTC

def extend_community (M, P, pvalue = 0.01):

    while 1:
        r = P.min()
        i,j = np.where(P == r)
        n = len(i)

        if r * n < pvalue:
            M[i,j] = 1.
            pvalue = pvalue - (n * r)
            P[i,j] = 1
        else:
            break

    return M,pvalue
        

def plot_community_distr(covers, labels , logx = False, logy = False):

    if isinstance(covers,list) is False :
        covers = [ covers ]

    if isinstance(labels,list) is False :
        labels = [ labels ]
    

    distr = { 'Community size distr.' : dict(zip(labels,[j.sum(0).T for j in covers])), 'Node part. distr.' : dict(zip(labels,[j.sum(1) for j in covers]))} 

    for i in distr:
        
        data = distr[i]
        fig = plt.figure()
        ax = fig.add_subplot(111,ylabel='ccdf') 
        markers = cycle(['o','^','v','d','s','<','>','+','x','1','2','3','4','h','H','p','|','_'])
        for j in data:
            x = data[j]
            n = len(x)
            bins = np.unique(np.asarray(x))
            minbin = max ( bins.min() , 1 )
            bins = np.concatenate(([minbin - 1], bins))
            h, bins = np.histogram(x, bins)
            cdf = np.cumsum(h)
            ccdf = [1-k/float(n) for k in cdf]
            ax.plot(bins[1:], ccdf,label = str(j) ,marker = markers.next(), ls = 'None', mfc='w')
            

        ax.legend(loc = 1, numpoints = 1)
        
        if logx:
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')
        xmin, xmax = plt.xlim()
        ax.set_xlim((-0.2, xmax))
        ymin, ymax = plt.ylim()
        ax.set_ylim((ymin, 1.1))
        ax.set_title(str(i))
    
    return distr
    
def community_stats(filename,method):
    
    filename = filename.split('.')[0]
    
    M = filename + 'M_' + method + '.pkl'
    M = open(M,'r')
    M = load(M).todense()
    svnet = filename + 'SVGraph.pkl'
    svnet = open(svnet,'r')
    svnet = load(svnet)
    net = filename + 'Graph.pkl'
    net = open(net,'r')
    net = load(net)
    
    W = nx.to_numpy_matrix(net)
    
    stats = dict()
    
    comps = nx.weakly_connected_components(svnet)
    stats ['svnet: largest connected component:'] = len(comps[0])
    inlinks = ( np.multiply( (M * M.T) > 0, W > 0 ) ).sum() 
    totlinks = (W > 0).sum()
    stats ['share of intracommunity links:'] = inlinks / float(totlinks)
    multinodes = (M.sum(1)>1).sum()
    n = float(len(M))
    stats ['share of multicommunity nodes:'] = multinodes / n
    v = W.sum()
    D = M.T * W * M
    inweight  = np.trace(D)
    stats ['share of intracommunity strength'] = inweight / v
    EW = Ematrix(W)
    W = W - EW
    D = M.T * W * M
    inweight  = np.trace(D)
    stats ['average intracommunity discrepancy:'] = inweight/ inlinks
    outweight = D.sum() - inweight
    stats ['average extracommunity discrepancy:'] = outweight/ (totlinks - inlinks)
    stats ['# of nodes in svnet'] = len(svnet)
    stats ['# of valid self-links'] = len(svnet.selfloop_edges())


    stats = stats.items()
    stats = np.array(stats,dtype = [('stat','S50'),('value',np.float32)])
    np.savetxt(filename +'_'+ method + '.comstats',stats,fmt = ['%10s','%10.10f'])
    return stats
    
    
def bip_community_stats(filename):
    
    filename = filename.split('.')[0]
    
    M = filename + 'M.pkl'
    M = open(M,'r')
    M = load(M).todense()
    svnet = filename + 'SVGraph.pkl'
    svnet = open(svnet,'r')
    svnet = load(svnet)
    W = filename + 'BiAdj.pkl'
    W = open(W,'r')
    W = load(W).todense()
    
    n,m = W.shape
    R = M[:n]
    C = M[n:]
    
    stats = dict()
    
    comps = nx.weakly_connected_components(svnet)
    stats ['svnet: largest connected component:'] = len(comps[0])
    inlinks = ( np.multiply( (R * C.T) > 0, W > 0 ) ).sum() 
    totlinks = (W > 0).sum()
    stats ['share of intracommunity links:'] = inlinks / float(totlinks)
    multinodes = (M.sum(1)>1).sum()
    n = float(len(M))
    stats ['share of multicommunity nodes:'] = multinodes / n
    v = W.sum()
    D = R.T * W * C
    inweight  = np.trace(D)
    stats ['share of intracommunity strength'] = inweight / v
    EW = Ematrix(W)
    W = W - EW
    D = R.T * W * C
    inweight  = np.trace(D)
    stats ['average intracommunity discrepancy:'] = inweight/ inlinks
    outweight = D.sum() - inweight
    stats ['average extracommunity discrepancy:'] = outweight/ (totlinks - inlinks)

    stats = stats.items()
    stats = np.array(stats,dtype = [('stat','S50'),('value',np.float32)])
    np.savetxt(filename + '.comstats',stats,fmt = ['%10s','%10.10f'])
    
    return stats

def Ematrix(M):
    out_degree = M.sum(1)
    in_degree = M.sum(0)
    v = M.sum()
    E = out_degree * in_degree / v
    return E


