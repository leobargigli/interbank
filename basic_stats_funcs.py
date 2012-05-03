from scipy.stats import binom , kendalltau, pearsonr, norm
from scipy.stats.mstats import mquantiles
import networkx as nx
from scipy.sparse import csc_matrix, extract, linalg, isspmatrix_csc,dia_matrix
#from scipy.cluster.hierarchy import linkage,fcluster,distance
#from sparsesvd import sparsesvd
import numpy as np
import pylab as plt
import os
from math import ceil



def density(G, nbunch = None):
    A = nx.to_numpy_matrix(G , weight = None)
    indices = np.diag_indices_from(A)
    A[indices] = 0.
    links = A.sum()
    if nbunch is not None:
        n = len(nbunch)
        m = G.number_of_nodes() - n
        max_links = n * (n-1) + n * m
    else:
        n = G.number_of_nodes()
        max_links = n * (n-1)
    return links / float(max_links)
                            
def reciprocity(G, nbunch = None,  weight = None):
    
    if nbunch is not None:
        nodes = set(G.nodes())
        fnodes = np.array(list(nodes - set(nbunch)))
        nbunch = np.array(list(nbunch))
        np.sort(nbunch)
        np.sort(fnodes)
        nodelist = np.append(nbunch, fnodes)
    else:
        nodelist = None
    
    W = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(W)
    W[indices] = 0.
    A = np.matrix(1. * (W > 0))
    M = {None: A,True: W}
    
    M = M[weight]
    l = M.sum()

    F = np.zeros(A.shape)

    if nbunch is None:
        n = M.shape[0]
        a = l / float( ( n - 1 ) *  n )
    else:
        n = len(nbunch)
        m = len(fnodes)
        a = l / float( n * ( n - 1 ) + n * m)
        #F[n:, n:] = a

    #F[indices] = a * np.ones((len(M), ))
    Msq = np.power(M, 2)
    omega = {None: 1., True: Msq.sum() / M.sum()}

    rho = np.multiply(M,M.T).sum() / M.sum()
    rho = ( rho - a ) / (omega[weight] - a)
    
    #~ E = a * np.ones(A.shape)
    #~ M = M - E + F 
    #~ Msq = np.power(M, 2)
    #~ MT = M.T - E + F
    #~ M = np.multiply(MT, M)
    #~ rho = M.sum() / Msq.sum()
#~ 
    return rho

def dists(G, nbunch = None):
    
    G = G.copy()
    
    if nbunch is None:
        nbunch = G.nodes()
    
    try:
        out_degree = G.out_degree(nbunch = nbunch)
        in_degree = G.in_degree(nbunch = nbunch)
        gross_out_weight = G.out_degree(weighted = True, nbunch = nbunch)
        gross_in_weight = G.in_degree(weighted = True, nbunch = nbunch)

    except TypeError:
        out_degree = G.out_degree(nbunch = nbunch)
        in_degree = G.in_degree(nbunch = nbunch)
        gross_out_weight = G.out_degree(weight = 'weight', nbunch = nbunch)
        gross_in_weight = G.in_degree(weight = 'weight', nbunch = nbunch)

        
    A = nx.to_scipy_sparse_matrix(G, nodelist = nbunch)
    i, j, grosscells = extract.find(A)

    selfloops = G.selfloop_edges(data = True)
    G.remove_edges_from(selfloops)
    
    
    try:
        net_out_weight = G.out_degree(weighted = True, nbunch = nbunch)
        net_in_weight = G.in_degree(weighted = True, nbunch = nbunch)

    except TypeError:
        net_out_weight = G.out_degree(weight = 'weight', nbunch = nbunch)
        net_in_weight = G.in_degree(weight = 'weight', nbunch = nbunch)


    A = nx.to_scipy_sparse_matrix(G, nodelist = nbunch)
    i, j, netcells = extract.find(A)

    dists = {
    'out-degree': np.array([out_degree[i] for i in nbunch],dtype = np.float32), 
    'in-degree': np.array([in_degree[i] for i in nbunch],dtype = np.float32), 
    'gross out-weight': np.array([gross_out_weight[i] for i in nbunch],dtype = np.float32), 
    'gross in-weight': np.array([gross_in_weight[i] for i in nbunch],dtype = np.float32),  
    'net out-weight': np.array([net_out_weight[i] for i in nbunch],dtype = np.float32), 
    'net in-weight': np.array([net_in_weight[i] for i in nbunch],dtype = np.float32),  
    'gross cells': grosscells,
    'net cells': netcells
    }
    
    return dists
    
def scatter(x, y, labels, filename, fmt ='png', diag = False):
    xlabel, ylabel, title = labels
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y)
    
    ax.set_title(title, fontsize = 16)
    ax.set_xlabel(xlabel, fontsize = 16)
    ax.set_ylabel(ylabel, fontsize = 16)
    
    if diag == True:
        line = plt.Line2D([0,1],[0,1])
        ax.add_line(line)
        ax.set_xlim((0.,1. ))
        ax.set_ylim((0.,1. ))
    plt.savefig(title + '_' + filename + '.' + fmt, fmt = fmt)
    plt.close()

def k_vs_nnk(G, label, nbunch = None):
    
    selfloops = G.selfloop_edges()
    netG = G.copy()
    netG.remove_edges_from(selfloops)
    
    functions = {
    'out-degree': G.out_degree,
     'gross out-weight': G.out_degree,
     'in-degree': G.in_degree, 
     'gross in-weight': G.in_degree,
    'net out-weight': netG.out_degree,
     'net in-weight': netG.in_degree,
     }
    weight = {
    'out-degree': None,
     'gross out-weight': 'weight',
     'in-degree': None, 
     'gross in-weight': 'weight',
     'net in-weight': 'weight',
     'net out-weight': 'weight',
     }
    weighted = {
    'out-degree': 0, 
    'net out-weight': 1,
    'gross out-weight': 1,
    'in-degree': 0, 
    'net in-weight': 1,
    'gross in-weight': 1
    }
    neighbors = {
    'out-degree':  G.successors, 
    'gross out-weight':  G.successors, 
    'net out-weight':  netG.successors, 
    'gross in-weight':  G.predecessors, 
    'net in-weight':  netG.predecessors, 
    'in-degree':  G.predecessors
    }
    try:
        dist = functions[label](weight = weight[label])
    except TypeError:
        dist = functions[label](weighted = weighted[label])
    dist = np.array(dist.items(),dtype = ([('nodes','|S10'),('degree',np.float32)]))
    values = dist['degree']
    nodes = dist['nodes']
    uniquevalues = np.unique(values)
    indices = np.where(uniquevalues == 0)
    uniquevalues = np.delete(uniquevalues,indices)
    knnk = dict()
    for k in uniquevalues:
        indices = np.where(values == k)
        knodes = np.array(nodes)[indices]
        nbrdeg = list()
        
        for h in knodes:
            nn = set(neighbors[label](h))
            nn.discard(h)
            if nbunch is not None:
                nn = list(nn.intersection(set(nbunch)))
            try:
                nnd = functions[label](nbunch = nn, weighted = weighted[label]).values()
            except TypeError:
                nnd = functions[label](nbunch = nn, weight = weight[label]).values()

            for i in nnd:
                nbrdeg.append(i)
                
        knnk[k] = np.average(nbrdeg)
    
    return knnk

def  participation_ratio(G, direction, nbunch = None, quant = 50, degree = True):
    
    
    functions = {'out': G.out_degree, 'in': G.in_degree}
    
    if nbunch is not None:
        nodes = set(G.nodes())
        fnodes = np.array(list(nodes - set(nbunch)))
        nbunch = np.array(list(nbunch))
        np.sort(nbunch)
        np.sort(fnodes)
        nodelist = np.append(nbunch, fnodes)
        n, m = (len(nbunch), len(fnodes))
    else:
        nodelist = G.nodes()
        nbunch = G.nodes()
        n = len(G)
    
    A = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(A)
    A[indices] = 0.
    
    if direction == 'in':
        A = A.T
    
    # this is the observed participatio ratio
    
    w = np.asarray(A.sum(1))
    indices = np.where(w == 0)
    w[indices] = 1.
    w.resize((len(w),))
    w = 1. / w
    w = np.asmatrix(np.diag(w))
    S = w * A

    S = np.power(S, 2)
    part_r = np.asarray(S.sum(1))
    part_r.resize((len(part_r),))
    part_r = part_r[:n]
    

    if degree is False:

    #this is the expectation under an expected weight model
        
        ew = np.asarray(A.sum(0))
        ew = np.repeat(ew,n,0)
        ew =  ew * np.array(A > 0)
        ewsum = np.array(ew.sum(1))
        indices = np.where(ewsum == 0)
        ewsum[indices] = 1.
        ewsum.resize((len(ewsum),))
        ewsum = 1. / ewsum
        ewsum = np.asmatrix(np.diag(ewsum))
        eS = ewsum * np.asmatrix(ew)
        
        eS = np.power(eS, 2)
        epart_r = np.asarray(eS.sum(1))
        epart_r.resize((len(epart_r),))
        epart_r = epart_r[:n]

        bins = np.unique(epart_r)
        t = len(bins)
        step = max( t / quant , 1)
        indices = range(t - 1,0,-step)[::-1]
        bins = bins[indices]
        t = len(bins)
    
        part_dict = dict()

        for k in range(t - 1): 
            lower = np.ma.masked_less_equal(epart_r,bins[k + 1])
            greater = np.ma.masked_greater(epart_r,bins[k])
            indices = lower.mask * greater.mask
            indices = np.where(indices > 0)
            part_k = part_r[indices]
            part_dict[ bins[k + 1] ] = np.average(part_k)


    #this is the uniform expectation: 1 / degree(i)

    else:
        
        dist = functions[direction](nbunch = nbunch)
        dtype = ([('node','|S10'),('degree',np.float32)])
        epart_r = np.array([(i,dist[i]) for i in nbunch],dtype = dtype)
        epart_r = epart_r['degree']
        indices = np.where(epart_r == 0)
        epart_r = np.delete(epart_r, indices)
        part_r = np.delete(part_r,indices)
        epart_r.resize((len(epart_r),))
        epart_r = 1. / epart_r[:n]
        bins = np.unique(epart_r)

        part_dict = dict()

        for k in bins: 
        
            indices = np.where(epart_r == k)
            part_k = part_r[indices]
            part_dict[k] = np.average(part_k)


    return part_dict
        

def clust_vs_degree(G, filename, weight = None, nbunch = None,  format = None, directed = False):
    if format is None:    format = 'png'
    if directed is False:
        C, degree = clustering_coefficient(G, weight = weight, nbunch = nbunch)
    else: 
        C, degree = dir_clustering_coefficient(G, weight = weight, nbunch = nbunch)
    
    dirdict  = {True: 'directed', False: ''}
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(degree, C)
    avg_C = np.average(C)
    xmin, xmax = ax.get_xlim()
    line = plt.Line2D([xmin,xmax],[avg_C,avg_C])
    ax.add_line(line)
    ax.set_ylabel('Clustering coefficient')
        
    wlabel = {None: 'degree', 'weight':'weight'}
    if nbunch is not None:
        nlabel = '_adj'
    else:
        nlabel = '_unadj'
    if filename.find('dom') <> -1:
        nlabel = ''
    ax.set_xlabel(wlabel[weight])
    ax.set_title(dirdict[directed] + ' cluster vs ' + wlabel[weight])
    plt.savefig(dirdict[directed] + 'cluster_vs_' + wlabel[weight] + '_' + filename+nlabel+'.' + format, format = format)
    plt.close()


def dir_clustering_coefficient(G, weight = None, nbunch = None):
    
    if nbunch is not None:
        nodes = set(G.nodes())
        fnodes = np.array(list(nodes - set(nbunch)))
        nbunch = np.array(list(nbunch))
        np.sort(nbunch)
        np.sort(fnodes)
        nodelist = np.append(nbunch, fnodes)
        n, m = (len(nbunch), len(fnodes))
    else:
        nodelist = None
        nbunch = G.nodes()
        n = len(G)
        
    W = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(W)
    W[indices] = 0.
    A  =  np.matrix(1. * (W > 0))

    out_degree = G.out_degree(nbunch = nbunch)
    in_degree = G.in_degree(nbunch = nbunch)
    degree = np.array([out_degree[i]+in_degree[i] for i in nbunch])
    bil_degree = np.diag(A**2)
    if nbunch is not None:
        bil_degree = bil_degree[0:n]
    Tmax = 2 * (degree * (degree - 1) - 2 * bil_degree)
    indices = np.where(Tmax == 0)
    Tmax[indices] = 1.
    
    if weight is not None:
        A = W / W.max()

    A = A + A.T
    triangles = np.diag(A**3)
    if nbunch is not None:
        triangles = triangles [:n]
    C = triangles / Tmax
    C = C.tolist()
    return C, degree
     
def clustering_coefficient(G, weight = None,  nbunch = None):
    
    G = G.to_undirected()
    
    if nbunch is not None:
        nodes = set(G.nodes())
        fnodes = np.array(list(nodes - set(nbunch)))
        nbunch = np.array(list(nbunch))
        np.sort(nbunch)
        np.sort(fnodes)
        nodelist = np.append(nbunch, fnodes)
        n, m = (len(nbunch), len(fnodes))
    else:
        nodelist = None
        nbunch = G.nodes()
        n = len(G)
        
        
    W = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(W)
    W[indices] = 0.
    A  =  np.matrix(1. * (W > 0))
    degree = G.degree(nbunch = nbunch)
    degree = np.array([degree[i] for i in nbunch])
    Td = degree * (degree - 1) 
    indices = np.where(Td == 0)
    Td[indices] = 1.

    if weight is not None:
        A = W / W.max()
    triangles = np.diag(A**3)
    if nbunch is not None:
        triangles = triangles [:n]

    C = triangles / Td
    C = C.tolist()
    
    return C, degree
    
def node_degree_xy(G, x='out', y='out', nbunch = None,  weighted = False):

    direction = {'out':G.out_degree_iter, 'in':G.in_degree_iter}
    weighted_dict = {True: 'weight', False: None}
    try:
        xdeg = direction[x](weighted = weighted, nbunch = nbunch)
    except TypeError:
       xdeg = direction[x](weight = weighted_dict[weighted], nbunch = nbunch)
    for u,degu in xdeg:
        neighbors = set([nbr for _,nbr in G.edges_iter(u)])
        if nbunch is not None:
            xset = set(nbunch)
            neighbors = neighbors.intersection(xset)
        try:
            neighbors = list(neighbors.discard(u))
        except TypeError:
            pass

        try:
           ydeg = direction[y](neighbors, weighted = weighted)
        except TypeError:
           ydeg = direction[y](neighbors, weight = weighted_dict[weighted])
   
        for v,degv in ydeg:
           yield degu,degv
 
def assortativity(G, x='out', y='out', nbunch = None, weighted = False):
    xy = node_degree_xy(G, x = x, y = y, nbunch = nbunch,  weighted = weighted)
    x , y = ([], [])
    for i, j in xy:
        x.append(i)
        y.append(j)
    x = np.array(x)
    y = np.array(y)
    rho, p = pearsonr(x,y)
    return rho, p
    
     
