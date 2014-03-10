from scipy.stats import binom , kendalltau, pearsonr, norm, betai
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
        max_links = n * (n - 1 + 2 * m)
    else:
        n = G.number_of_nodes()
        max_links = n * (n-1)
    return links / float(max_links)
                            
def reciprocity(G, nbunch = None,  weight = None):
    
    if nbunch is not None:
        nodes = np.sort(G.nodes())
        nbunch = np.sort(nbunch)
        fnodes = np.setdiff1d(nodes,nbunch)
        nodelist = np.append(nbunch, fnodes)
    else:
        nbunch = G.nodes()
        nodelist = G.nodes()
        fnodes = list()
    
    W = np.array(nx.to_numpy_matrix(G, nodelist = nodelist))
    indices = np.diag_indices_from(W)
    W[indices] = 0.
    
    if weight is None:
        
        W = 1. * (W > 0)

    l = float(W.sum())
    n = len(nbunch)
    m = len(fnodes)
    df = n  * (n - 1  + 2 * m)
    a = l / df # this is to take into account that the maximal number of observable links is lower than (n+m)*(n+m-1)
    W[n:,n:] = a # this is to set the terms corresponding to unobserved exposures to zero

    l2 = (W**2).sum()
    omega =  l2 / l
        
    rho = (W * W.T).sum() / l
    rho = (rho - a) / (omega - a)
    
    rho = max(min(rho, 1.0), - 1.0)

    if abs(rho) == 1.0:
        prob = 0.0
    else:
        t_squared = rho * rho * (df / ((1.0 - rho) * (1.0 + rho)))
        prob = betai(0.5*df, 0.5, df / (df + t_squared))    
    
    return rho,prob

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
    'out-degree': 
    np.array([out_degree[i] for i in nbunch],dtype = np.float32), 
    'in-degree': 
    np.array([in_degree[i] for i in nbunch],dtype = np.float32), 
    'gross_out-weight': 
    np.array([gross_out_weight[i] for i in nbunch],dtype = np.float32), 
    'gross_in-weight': 
    np.array([gross_in_weight[i] for i in nbunch],dtype = np.float32),  
    'net_out-weight': 
    np.array([net_out_weight[i] for i in nbunch],dtype = np.float32), 
    'net_in-weight': 
    np.array([net_in_weight[i] for i in nbunch],dtype = np.float32),  
    'gross_cells': grosscells,
    'net_cells': netcells
    }
    
    return dists
    
def scatter_save(x, y, filename, fmt ='pdf', diag = False, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(x, y, **kwargs)
    
    if diag == True:
        line = plt.Line2D([0,1],[0,1])
        ax.add_line(line)
        ax.set_xlim((0.,1. ))
        ax.set_ylim((0.,1. ))
        
    plt.savefig(filename + '.' + fmt, fmt = fmt)
    plt.close()

def k_vs_nnk(G, label, nbunch = None):
    
    selfloops = G.selfloop_edges()
    netG = G.copy()
    G.remove_edges_from(selfloops)
    
    functions = {
    'out-degree': G.out_degree,
     'gross_out-weight': G.out_degree,
     'in-degree': G.in_degree, 
     'gross_in-weight': G.in_degree,
    'net_out-weight': netG.out_degree,
     'net_in-weight': netG.in_degree,
     }
    weight = {
    'out-degree': None,
     'gross_out-weight': 'weight',
     'in-degree': None, 
     'gross_in-weight': 'weight',
     'net_in-weight': 'weight',
     'net_out-weight': 'weight',
     }
    weighted = {
    'out-degree': 0, 
    'net_out-weight': 1,
    'gross_out-weight': 1,
    'in-degree': 0, 
    'net_in-weight': 1,
    'gross_in-weight': 1
    }
    neighbors = {
    'out-degree':  G.successors, 
    'gross_out-weight':  G.successors, 
    'net_out-weight':  netG.successors, 
    'gross_in-weight':  G.predecessors, 
    'net_in-weight':  netG.predecessors, 
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


def matrix_from_graph(G,nbunch = None):
    
    if nbunch is not None:
        nodes = G.nodes()
        fnodes = np.setdiff1d(nodes,nbunch)
        np.sort(nbunch)
        np.sort(fnodes)
        nodelist = np.append(nbunch, fnodes)
    else:
        nodelist = G.nodes()
        nbunch = G.nodes()
        
    A = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(A)
    A[indices] = 0.
    A = np.array(A)
    
    return A

def participation_ratio(A,weight, nbunch = None):
    
    if nbunch is not None:
        n = len(nbunch)
    else:
        n = len(A)
    
    # this is the observed participatio ratio
    
    indices = np.where(weight == 0.)
    weight[indices] = 1.
    q = len(weight)
    weight.resize((q,1))
    W = np.repeat(weight,q,1)
    S = A / W
    S = S ** 2
    part_r = S.sum(1)
    part_r = part_r[:n]

    return part_r
    
    
def exp_part_ratio(G,direction, nbunch = None, quant = 50, degree = True):

    functions = {
    'out': G.out_degree,
    'in': G.in_degree    
    }

    A = matrix_from_graph(G, nbunch = nbunch)
    if direction is 'in':
        A = A.T
    weight = A.sum(1)
    part_r = participation_ratio(A, weight, nbunch = nbunch)

    if nbunch is not None:
        nodes = G.nodes()
        fnodes = np.setdiff1d(nodes,nbunch)
        np.sort(nbunch)
        np.sort(fnodes)
        n,m = len(nbunch),len(fnodes)
    else:
        nodes = G.nodes()     
        nbunch = G.nodes()
        n = len(nodes)

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
        
  
def dir_clustering_coefficient(G, weight = None, nbunch = None):
    
    G = G.copy()    
    selfloops = G.selfloop_edges()
    G.remove_edges_from(selfloops) # exclude selfloops
    
    if nbunch is not None:
        nodes = np.sort(G.nodes())
        nbunch = np.sort(nbunch)
        fnodes = np.setdiff1d(nodes,nbunch)
        nodelist = np.append(nbunch, fnodes)

    else:
        nbunch = G.nodes()
        nodelist = G.nodes()

    n = len(nbunch)
    W = nx.to_numpy_matrix(G, nodelist = nodelist)
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
    
    A = A + A.T
    triangles = np.diag(A**3)[:n]
    cc = triangles / Tmax
    cc = cc.tolist()
    
    return cc,degree
     
def clustering_coefficient(G,  nbunch = None):
    
    G = G.to_undirected()
    selfloops = G.selfloop_edges()
    G.remove_edges_from(selfloops) # exclude selfloops
    
    if nbunch is not None:
        nodes = np.sort(G.nodes())
        nbunch = np.sort(nbunch)
        fnodes = np.setdiff1d(nodes,nbunch)
        nodelist = np.append(nbunch, fnodes)

    else:
        nbunch = G.nodes()
        nodelist = G.nodes()
        
    n = len(nbunch)
    W = nx.to_numpy_matrix(G, nodelist = nodelist) # we take into account all links
    A  =  np.matrix(1. * (W > 0))
    degree = G.degree(nbunch = nbunch)
    degree = np.array([degree[i] for i in nbunch])
    Td = degree * (degree - 1) 
    indices = np.where(Td == 0)
    Td[indices] = 1.
    triangles = np.diag(A**3)[:n]
    cc = triangles / Td
    cc = cc.tolist()
    
    return cc,degree
    
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
    
     
def core_vector(A,tiering = False, 
                selfloops = False, directed = False):
    
    outk = A.sum(1)
    ink = A.sum(0)
    deg = outk + ink
    n = len(deg)
    indices = np.argsort(deg)[::-1]
    Z = np.array(range(1,n + 1)) - 0.5 * (1 + deg[indices])
    core = indices[Z < 0]
    core_v = np.zeros((n,))
    core_v[core] = 1
    c = core_v.sum()
    if tiering:
        indices = np.argsort(core_v)[::-1]
        M = A[indices][:,indices] # this is the sorted adjacency matrix
        CC = M[:c,:c]
        PP = M[c:,c:]
        if selfloops:
            Z = c**2 - CC.sum() + PP.sum()
        else:
            Z = c*(c-1) - CC.sum() + PP.sum()
        out_indices = np.where(outk <= c)
        out_indices = np.intersect1d(out_indices[0],core)
        tier_v = core_v.copy()
        #print out_indices

        for i in out_indices:
            h = np.where(indices == i)[0]
         #   print h
            if M[h,c:].sum() == 0:

                M = np.insert(M,[c],M[h],0) # put the h^th node in the c^th row
                M = np.delete(M,h,0)
                M = np.insert(M,[c],M[:,h],1)# put the h^th node in the c^th col
                M = np.delete(M,h,1)

                if selfloops:
                    Znew = (c-1)**2 - M[:c-1,:c-1].sum() + M[c-1:,c-1:].sum()
                else:
                    Znew = (c-1)*(c-2) - M[:c-1,:c-1].sum() + M[c-1:,c-1:].sum()
               
          #      print Znew
                dZ = c - n + Znew - Z
                if Znew < Z:
                    return 'something wrong'
           #     print dZ
                if dZ < 0:
                    tier_v[i] = 0
                    c -= 1 
                    Z = Znew
            
        in_indices = np.where(ink <= c)
        in_indices = np.intersect1d(in_indices[0],core)
            
        for i in in_indices:
            h = np.where(indices == i)[0]
            if M[c:,h].sum() == 0:
                
                M = np.insert(M,[c],M[h],0)
                M = np.delete(M,h,0)
                M = np.insert(M,[c],M[:,h],1)
                M = np.delete(M,h,1)
                #return M,c
                
                if selfloops:
                    Znew = (c-1)**2 - M[:c-1,:c-1].sum() + M[c-1:,c-1:].sum()
                else:
                    Znew = (c-1)*(c-2) - M[:c-1,:c-1].sum() + M[c-1:,c-1:].sum()
                
                dZ = c - n + Znew - Z
                if Znew < Z:
                    return 'something wrong'
                if dZ < 0:
                    tier_v[i] = 0
                    c -= 1 
                    Z = Znew
        
        return core_v,tier_v,pearsonr(core_v,tier_v)
    else:
        return core_v