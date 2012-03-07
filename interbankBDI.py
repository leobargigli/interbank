from scipy.stats import binom , kendalltau, pearsonr
import networkx as nx
from scipy.sparse import csc_matrix, extract, linalg, isspmatrix_csc
import numpy as np
import pylab as plt
import os

class Year:

    def __init__(self, filename, delimiter = ','): #delimiter = ','
        
        try:
            edgelist = np.loadtxt(filename) 
        except ValueError:
            edgelist = np.loadtxt(filename, delimiter = delimiter) 
        #empty_rows = np.where(edgelist[:, 2] == 0)[0]
        #edgelist = np.delete(edgelist,empty_rows,0)
        #source = np.int32(edgelist[:, 0])
        #destination = np.int32(edgelist[:, 1])
        weight = edgelist[:, 2].copy()
        edgelist = [(str(i[0]),str(i[1]), i[2]) for i in edgelist]
        self.edgetype = np.dtype([('source','|S10'),('dest','S10'),('weight',np.float64)])
        edgelist = np.array(edgelist, dtype = self.edgetype)
        
        G = nx.DiGraph()
        G.add_weighted_edges_from(edgelist)
        self.Graph = G
        self.Adj = nx.to_scipy_sparse_matrix(G) 
        filename = os.path.splitext(filename)[0]
        self.filename = filename
        self.edgelist = edgelist
        self.weights = edgelist['weight']
    
    def stats(self, nbunch = None):
        G = self.Graph
        nodes = G.number_of_nodes()
        size = G.size(weight = 'weight')
        if nbunch is not None:
            nodes = len(nbunch)
        edges = G.number_of_edges()
        selfloops = G.selfloop_edges()
        edges = edges - len(selfloops)
        density = self.density(nbunch = nbunch)
        volume = self.Adj.sum()
        comps = nx.weakly_connected_component_subgraphs(G)
        comp_size = np.array([i.number_of_nodes() for i in comps], dtype = np.int32)
    
        
        assortativity_out = assortativity(G, x='out', y='out', nbunch = nbunch)
        assortativity_in = assortativity(G, x='in', y='in', nbunch = nbunch)

        try:
            w_assortativity_out = assortativity(G, x='out', y='out', weighted = True, nbunch = nbunch)
            w_assortativity_in = assortativity(G, x='in', y='in', weighted = True, nbunch = nbunch)
        except TypeError:
            w_assortativity_out = 'na'
            w_assortativity_in  = 'na'

        recip = reciprocity(G, nbunch = nbunch)
        w_recip = reciprocity(G, weight = True, nbunch = nbunch)
        distG = dists(G, nbunch = nbunch)
        out_degree = distG['out-degree']
        out_weight = distG['out-weight']
        in_degree = distG['in-degree']
        in_weight = distG['in-weight']

        weight = out_weight + in_weight
        degree = out_degree + in_degree

        indices = np.where(out_degree == 0.)
        out_degree = np.delete(out_degree, indices)
        out_weight = np.delete(out_weight, indices)
        
        indices = np.where(in_degree == 0.)
        in_degree = np.delete(in_degree, indices)
        in_weight = np.delete(in_weight, indices)
        
        indices = np.where(degree == 0.)
        degree = np.delete(degree, indices)
        weight = np.delete(weight, indices)

        mean_ow = out_weight / out_degree
        mean_iw = in_weight / in_degree
        mean_w = weight / degree


        try:
            out_tau, out_p = kendalltau(mean_ow, out_degree)
            in_tau, in_p = kendalltau(mean_iw, in_degree)
            tot_tau, tot_p = kendalltau(mean_w, degree) 
        except TypeError:
            out_tau, out_p = (0., 0.)
            in_tau, in_p = (0., 0.)
            tot_tau, tot_p = (0., 0.)
        dC = np.average(clustering_coefficient(G, nbunch = nbunch)[0])
        dWC = np.average(clustering_coefficient(G, weight = 'weight', nbunch = nbunch)[0])

        uC = np.average(und_clustering_coefficient(G, nbunch = nbunch)[0])
        uWC = np.average(und_clustering_coefficient(G, weight = 'weight', nbunch = nbunch)[0])

        
        G.remove_edges_from(selfloops)

        if len(comps) == 1:
            avg_path_length = nx.average_shortest_path_length(G)
            try:
                avg_weight_path_length = nx.average_shortest_path_length(G,weighted = True)
            except TypeError:
                avg_weight_path_length = nx.average_shortest_path_length(G,weight = 'weight')
                
        else:
            avg_path_length = nx.average_shortest_path_length(comps[0])
            try:
                avg_weight_path_length = nx.average_shortest_path_length(comps[0],weighted = True)
            except TypeError:
                avg_weight_path_length = nx.average_shortest_path_length(comps[0],weight = 'weight')
        

        
        if nbunch is not None:
            label = '_adj'
        else:
            label = '_unadj'
        if  self.filename.find('dom') <>-1:
            label = ''
        
        filename = self.filename + label + '_stats.dat'
        output = open(filename, 'wb')
        output.write('# of nodes: %i\n'%(nodes)) 
        output.write('# of edges: %i\n'%(edges)) 
        output.write('Density: %f\n'%(density)) 
        output.write('Volume: %f\n'%(volume)) 
        output.write('Sum of weights: %f\n'%(self.weights.sum())) 
        output.write('Volume(2): %f\n'%(size) )
        output.write('Weak components size distribution:\n')
        np.savetxt(output, comp_size, fmt='%1i')
        output.write('Average path length: %f\n'%( avg_path_length )) 
        output.write('Average weighted path length: %f\n'%( avg_weight_path_length )) 
        output.write('Out-degree assortativity (p-value): %f (%f)\n'%assortativity_out)
        output.write('In-degree assortativity (p-value): %f (%f)\n'%assortativity_in) 
        if isinstance(w_assortativity_out, str):
            output.write('Out-weight assortativity: %s\n'%( w_assortativity_out)) 
            output.write('In-weight assortativity: %s\n'%( w_assortativity_in)) 
        else:
            output.write('Out-weight assortativity (p-value): %f (%f)\n'%w_assortativity_out) 
            output.write('In-weight assortativity (p-value): %f (%f)\n'%w_assortativity_in) 
        output.write('Degree reciprocity: %f\n'%(recip)) 
        output.write('Weight reciprocity: %f\n'%(w_recip))
        output.write('Average directed clustering: %f\n'%(dC)) 
        output.write('Average dir. & weighted clustering: %f\n'%(dWC))
        output.write('Average undirected clustering: %f\n'%(uC)) 
        output.write('Average undir. & weighted clustering: %f\n'%(uWC))
        output.write('Kendall tau w_out / d_out vs d_out (p-value): %f (%f)\n' % (out_tau, out_p))
        output.write('Kendall tau w_in / d_in vs d_in (p-value): %f (%f)\n'% (in_tau, in_p))
        output.write('Kendall tau w / d vs d (p-value): %f (%f)\n'% (tot_tau, tot_p))
        
        output.flush()

    
    def density(self, nbunch = None):
        W = self.Adj.todense()
        A = np.matrix(1. * (W>0))
        indices = np.diag_indices_from(A)
        A[indices] = 0.
        links = A.sum()
        if nbunch is not None:
            n = len(nbunch)
            m = self.Graph.number_of_nodes() - n
            max_links = n * (n-1) + n * m
        else:
            n = self.Graph.number_of_nodes()
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
    try:
        out_degree = np.asarray(G.out_degree(nbunch = nbunch).values(), dtype = np.float32)
        in_degree = np.asarray(G.in_degree(nbunch = nbunch).values(), dtype = np.float32)
        out_weight = np.asarray(G.out_degree(weighted = True, nbunch = nbunch).values(), dtype = np.float32)
        in_weight = np.asarray(G.in_degree(weighted = True, nbunch = nbunch).values(), dtype = np.float32)
    except TypeError:
        out_degree = np.asarray(G.out_degree(nbunch = nbunch).values(), dtype = np.float32)
        in_degree = np.asarray(G.in_degree(nbunch = nbunch).values(), dtype = np.float32)
        out_weight = np.asarray(G.out_degree(weight = 'weight', nbunch = nbunch).values(), dtype = np.float32)
        in_weight = np.asarray(G.in_degree(weight = 'weight', nbunch = nbunch).values(), dtype = np.float32)
        
    A = nx.to_scipy_sparse_matrix(G)
    i, j, cells = extract.find(A)

    dists = {'out-degree': out_degree, 'in-degree': in_degree, 'out-weight': out_weight, 'in-weight': in_weight,  'cells': cells}
    return dists
    
def scatter(x, y, labels, filename, format='png', diag = False):
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
    plt.savefig(title+'_'+filename+'.'+format, format = format)
    plt.close()

def k_vs_nnk(G, label, nbunch = None):
    functions = {'out-degree': G.out_degree, 'out-weight': G.out_degree,'in-degree': G.in_degree, 'in-weight': G.in_degree}
    weight = {'out-degree': None, 'out-weight': 'weight','in-degree': None, 'in-weight': 'weight'}
    weighted = {'out-degree': 0, 'out-weight': 1,'in-degree': 0, 'in-weight': 1}
    neighbors = {'out-degree':  G.successors, 'out-weight':  G.successors, 'in-weight':  G.predecessors, 'in-degree':  G.predecessors}
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

def  participation_ratio(G, direction, nbunch = None):
    
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
        nodelist = None
        nbunch = G.nodes()
        n = len(G)
    
    A = nx.to_numpy_matrix(G, nodelist = nodelist)
    indices = np.diag_indices_from(A)
    A[indices] = 0.
    
    if direction == 'in':
        A = A.T
    w = np.asarray(A.sum(1))
    indices = np.where(w==0)
    w[indices] = 1.
    w.resize((len(w), ))
    w = 1/w
    W = np.asmatrix(np.diag(w))
    S = W * A

    S = np.power(S, 2)
    part_r = np.asarray(S.sum(1))
    part_r.resize((len(part_r, )))
    if nbunch is not None:
        part_r = part_r[:n]
    
    dist = functions[direction](nbunch = nbunch)
    dist = np.array([(i,dist[i]) for i in nbunch],dtype = ([('node','|S10'),('degree',np.float32)]))
    values = dist['degree']
    uniquevalues = np.float32(np.unique(values))
    indices = np.where(uniquevalues == 0)
    uniquevalues = np.delete(uniquevalues, indices)
    part_degree = dict()
    for k in uniquevalues:
        indices = np.where(values == k)
        part_k = part_r[indices]
        part_degree[1/k] = np.average(part_k)

    return part_degree
        

def clust_vs_degree(G, filename, weight = None, nbunch = None,  format = None):
    if format is None:    format = 'png'
    C, degree = clustering_coefficient(G, weight = weight, nbunch = nbunch)
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
    ax.set_title('cluster_vs_'+wlabel[weight])
    plt.savefig('cluster_vs_'+wlabel[weight]+'_'+filename+nlabel+'.' + format, format = format)
    plt.close()


def clustering_coefficient(G, weight = None, nbunch = None):
    
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
     
def und_clustering_coefficient(G, weight = None,  nbunch = None):
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
    
     

    
