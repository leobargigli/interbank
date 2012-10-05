from scipy.stats import binom
import networkx as nx
from scipy.sparse import csc_matrix, extract, linalg, isspmatrix_csc, dia_matrix
from scipy.cluster.hierarchy import linkage,fcluster,distance
import numpy as np
import pylab as plt
#from sparsesvd import sparsesvd
import os
from cPickle import load, dump
from comm_detect_funcs import *
#import interbank_classes as IB

class Year():#IB.Year):

    def __init__(self, filename, delimiter = '\t', dtype = int, divide_factor = 1E0):

        edgelist = np.loadtxt(filename,delimiter,dtype)
        empty_rows = np.where(edgelist[:, 2] == 0)[0]
        edgelist = np.delete(edgelist,empty_rows,0)
        edgelist[:,2] = map(lambda x : np.ceil(x / divide_factor), edgelist[:,2])
            
        edgelist = [('F'+str(i[0]),'B'+str(i[1]), i[2]) for i in edgelist]
        self.edgetype = np.dtype([('source','|S10'),('dest','S10'),('weight',np.float64)])
        edgelist = np.array(edgelist, dtype = self.edgetype)
        self.edgelist = edgelist.copy()
        self.firms = np.unique(edgelist['source'])
        self.banks = np.unique(edgelist['dest'])
        self.filename = filename.split('.')[0]

        G = nx.DiGraph()
        G.add_nodes_from(self.firms)
        G.add_nodes_from(self.banks)
        G.add_weighted_edges_from(edgelist)
        self.Net = G
        
        row_indices = dict(zip(self.firms, range(len(self.firms))))
        col_indices = dict(zip(self.banks, range(len(self.banks))))
        n_edges = len(edgelist)

        ij = np.zeros((2, n_edges))
        for i in xrange(n_edges):
            source = edgelist['source'][i]
            dest = edgelist['dest'][i]
            ij[0, i] = row_indices[source]
            ij[1, i] = col_indices[dest]

        data = edgelist['weight']
        A = csc_matrix((data,ij))
        self.Adj = A
        self.descr = 'original network'
        
    def saveDigraph(self):
        G = self.Adj
        outfile = open(self.filename.split('.')[0] + 'Graph.pkl','wb')
        dump(G,outfile)

class SVnet(Year):#,IB.SVnet):
    
    def __init__(self, Year, pvalue = 0.01):

        A =Year.Adj
        v = float(A.sum())

        n, m = A.shape
        self.sets = (n, m)
        alpha = pvalue / float(n * m)
        in_degree = A.sum(0)
        out_degree = A.sum(1)
        i, j, aij = extract.find(A)
        
        nonzero = len(i)
        pij = np.zeros((nonzero, ))
        for h in xrange(nonzero):                      
            pij[h] = out_degree[i[h]] * in_degree[0,j[h]] / v**2
        P = 1-binom.cdf(aij - 1,v,pij)
        data = 1. * (P<= alpha)
        zero_entries = np.where(data == 0) 
        data = np.delete(data, zero_entries)
        i = np.delete(i, zero_entries)
        j = np.delete(j, zero_entries)
        aij = np.delete(aij,zero_entries)
        ij = np.asarray(zip(i,j)).T
        self.svnet = csc_matrix((data, ij))
        self.Adj = csc_matrix((aij,ij))
        self.filename = Year.filename
        self.edgetype = Year.edgetype
        self.banks = Year.banks
        self.firms = Year.firms
        self.descr = 'valid network'

        
    def edgelist(self, weighted = True):

        i, j, data = extract.find(self.Adj)
        n_edges = len(data)
     
        firms = dict(zip(range(len(self.firms)), self.firms))
        banks = dict(zip(range(len(self.banks)), self.banks))
        edgelist = np.zeros((n_edges, ), dtype = self.edgetype)

        for k in xrange(n_edges):
            row = i[k]
            col = j[k]
            edgelist['source'][k] = firms[row]
            edgelist['dest'][k] = banks[col]
            edgelist['weight'][k] = data[k]

        return  edgelist
        
    def make_Digraph(self):

        G = nx.DiGraph()
        G.add_nodes_from(self.firms)
        G.add_nodes_from(self.banks)
        edgelist = self.edgelist()
        G.add_weighted_edges_from(edgelist)
        filename = self.filename.split('.')[0]
        outfile = open( filename + 'SVGraph.pkl','wb')
        dump(G,outfile)
        
        return G
        
    def makecover(self, partition, q = '',method = 'spectral'):
        svnet = self.svnet
        svnet = adj_from_biadj(svnet)
        M = community_matrix(partition)
        M = makecover(svnet, M)
        M = csc_matrix(M)
        #i, j, mij = extract.find(M)
        #nodelist = np.append(self.firms,self.banks)
        #cover = np.asarray(zip(nodelist, j),dtype = [('nodes','S10'),('community',np.int)])
        #filename = self.filename.split('.')[0]
        #np.savetxt(filename +'_' + str(q) + '.cover', cover, fmt = ['%10s','%10i'])
        #outfile = open(filename +'M_'+ method + '_' + str(q) + '.pkl','wb')
        #dump(M,outfile)
        return M
        
        
#def n_of_modules(A):
#
#    K = Kmatrix(A)
#    v = A.sum()
#    n, m = K.shape
#    r = min(n, m)
#    pvalue_diff = np.zeros((r- 2,1))
#    exp_value = (n * m) / float(v)
#    U, svs, V = sparsesvd(K, r-1)
#    eigs = np.power(svs,2)
#    eigs = eigs[1:] # remove the unitary sv
#    sigma = eigs.sum() 
#    pvalue = np.zeros((r- 1,1))
#    pvalue[0] = exp_value/sigma
#
#    for j in xrange(1, r - 1):
#        eigs = eigs[1:] # remove the largest sv
#        sigma = eigs.sum() 
#        pvalue[j]  = min(exp_value / sigma,1.)
#        pvalue_diff[j - 1] = pvalue[j] - pvalue[j-1]
#        #pvalue_old = pvalue_new
#
#    try: 
#        delta  = pvalue_diff.argmax() 
#       
#    except ValueError:
#        delta = -1
    
    
 #   n_of_modules = 2 + delta # the first addend is explained as follows: a) there is at least one module (the component itself) b) python index starts at 0

  #  return n_of_modules,  pvalue_diff, svs

#def svnet(A, pvalue = 0.01):
#    v = A.sum()
#    n, m = A.shape
#    alpha = pvalue / (float(n) * float(m))
#    in_degree = A.sum(0)
#    out_degree = A.sum(1)
#    i, j, aij = extract.find(A)
#    nonzero = len(i)
#    pij = np.zeros((nonzero, ))
#    for h in xrange(nonzero):                      
#        pij[h] = (out_degree[i[h]] * in_degree[0,j[h]] )/ v**2
#    P = 1-binom.cdf(aij - 1,v,pij)
#    data = 1.*(P<= alpha)
#    zero_entries = np.where(data == 0) [0]
#    data = np.delete(data, zero_entries)
#    i = np.delete(i, zero_entries)
#    j = np.delete(j, zero_entries)
#    ij = np.asarray(zip(i,j)).T
#    svnet = csc_matrix((data, ij))
#    aij = np.delete(aij, zero_entries)
#    A = csc_matrix((aij, ij))    
#    return A

def adj_from_biadj(A):
    n, m = A.shape
    try:
        A = A.todense()
    except AttributeError:
        pass
    UL = np.asmatrix(np.zeros((n, n)))
    LR = np.asmatrix(np.zeros((m, m)))
    U = np.concatenate((UL, A), 1)
    L = np.concatenate((A.T, LR), 1)
    A = np.concatenate((U, L), 0)
    return A

def bipartite_adj(G):
    G = G.to_undirected()
    if not nx.is_bipartite(G):
        return 'Error: the network is not bipartite'
    set1, set2 = nx.bipartite.sets(G)
    if len(set1) > len(set2):
        firms = np.array(list(set1))
        banks = np.array(list(set2))
    else:
        firms = np.array(list(set2))
        banks = np.array(list(set1))
    
    n = len(firms)
    nodelist = np.append(firms, banks)
    
    A = nx.to_scipy_sparse_matrix(G, nodelist = nodelist)
    biAdj = A[0:n, n:] 
    return biAdj

