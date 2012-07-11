from scipy.stats import binom , kendalltau
import networkx as nx
from scipy.sparse import csc_matrix, extract
#from scipy.cluster.hierarchy import linkage,fcluster,distance
#from sparsesvd import sparsesvd
import numpy as np
import os
from comm_detect_funcs import *
from basic_stats_funcs import *
from cPickle import dump

class Year:

    def __init__(self, filename, delimiter = ',',directed = True): #delimiter = ','
        
        try:
            edgelist = np.loadtxt(filename) 
        except ValueError:
            edgelist = np.loadtxt(filename, delimiter = delimiter) 
        weights = edgelist[:, 2].copy()
        edges = np.array([str(int(i[0])) +'*'+ str(int(i[1])) for i in edgelist])
        uni_edges = np.unique(edges)
        q = len(uni_edges)
        edgelist = list()

        for i in range(q):
            indices = np.where(edges == uni_edges[i])
            uni_weight = weights[indices].sum() 
            source,dest = tuple(uni_edges[i].split('*'))
            edgelist.append(tuple((source,dest,uni_weight)))
            
            
        self.edgetype = np.dtype([('source','|S10'),('dest','S10'),('weight',np.float64)])
        edgelist = np.array(edgelist, dtype = self.edgetype)
        if directed is True:
            G = nx.DiGraph()
        else:
            G = nx.Graph()
        G.add_weighted_edges_from(edgelist)
        self.Net = G
        self.nodes = list( np.sort ( G.nodes() ))
        #self.Adj = nx.to_scipy_sparse_matrix(G, nodelist = self.nodes , format = 'csc') 
        self.Adj = csc_matrix(nx.to_numpy_matrix(G, nodelist = self.nodes))
        filename = os.path.splitext(filename)[0]
        self.filename = filename
        self.edgelist = edgelist
        self.weights = edgelist['weight']

    def saveDigraph(self):
        G = self.Net
        outfile = open(self.filename.split('.')[0] + 'Graph.pkl','wb')
        dump(G,outfile)


    
    def stats(self, distG,nbunch = None):
        
        G = self.Net
        nodes = G.number_of_nodes()
        if nbunch is not None:
            nodes = len(nbunch)

        edges = G.number_of_edges()
        selfloops = G.selfloop_edges()
        edges = edges - len(selfloops)
        d = density(G , nbunch = nbunch)
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
        #~ distG = dists(G, nbunch = nbunch)
        out_degree = distG['out-degree']
        out_weight = distG['gross out-weight']
        in_degree = distG['in-degree']
        in_weight = distG['gross in-weight']

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
        except RuntimeError:
            out_tau, out_p = (0., 0.)
            in_tau, in_p = (0., 0.)
            tot_tau, tot_p = (0., 0.)
    
        dC = np.average(dir_clustering_coefficient(G, nbunch = nbunch)[0])
        dWC = np.average(dir_clustering_coefficient(G, weight = 'weight', nbunch = nbunch)[0])

        uC = np.average(clustering_coefficient(G, nbunch = nbunch)[0])
        uWC = np.average(clustering_coefficient(G, weight = 'weight', nbunch = nbunch)[0])

        
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
        output.write('Density: %f\n'%(d)) 
        output.write('Volume: %f\n'%(volume)) 
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


class SVnet(Year):

    def __init__(self, Year, pvalue = 0.01):

        A = Year.Adj
        
        n, m = A.shape
        alpha = pvalue / n / m
        in_degree = A.sum(0)
        out_degree = A.sum(1)
        i, j, wij = extract.find(A)
        indices = np.where(wij > 0)
        eps = max(wij[indices].min(),0.1)
        v = A.sum() / eps
        
        nonzero = len(i)
        pij = np.zeros((nonzero, ))
        for h in xrange(nonzero):                      
            pij[h] = out_degree[i[h]] * in_degree[0,j[h]] / v**2
        P = 1 - binom.cdf(wij - 1,v,pij)
        data = P <= alpha
        zero_entries = np.where(data == 0)
        data = np.delete(data, zero_entries)
        i = np.delete(i, zero_entries)
        j = np.delete(j, zero_entries)
        wij = np.delete(wij,zero_entries)
        ij = np.asarray(zip(i,j)).T
        self.svnet = csc_matrix((data, ij), shape = (n,m) )
        self.Adj = csc_matrix((wij, ij), shape = (n,m) )
        self.nodes = Year.nodes
        self.filename = Year.filename
        self.edgetype = Year.edgetype


    def edgelist(self , weighted = True):

        nodes = dict(zip(range(len(self.nodes)), self.nodes))
        i, j, data = extract.find(self.Adj)

        n_edges = len(data)
        edgelist = np.zeros((n_edges, ), dtype = self.edgetype)

        for k in xrange(n_edges):
            row = i[k]
            col = j[k]
            edgelist['source'][k] = nodes[row]
            edgelist['dest'][k] = nodes[col]
            edgelist['weight'][k] = data[k]
        
        filename = self.filename.split('.')[0]
        np.savetxt(filename + '.svnet', edgelist, fmt = ['%10s','%10s','%10.10f'])
        
        return edgelist

    def make_Digraph(self):

        G = nx.DiGraph()
        G.add_nodes_from(self.nodes)
        edgelist = self.edgelist() 
        G.add_weighted_edges_from(edgelist)
        filename = self.filename.split('.')[0]
        outfile = open( filename + 'SVGraph.pkl','wb')
        dump(G,outfile)
        
        return G
        

    def from_pajek(self):

        os.chdir('Infomap')
        filename = self.filename.split('.')[0]
        filename = '%s.clu' % (filename)
        try:
            cluster = open(filename,'r')
        except IOError:
            os.chdir('..')
            return 'Run the Infomap algorithm for %s first' %(filename.split('.')[0])
                
        cluster = cluster.readlines()
        cluster = cluster[1:] # deletes the header
        cluster = np.array(cluster,dtype='float') - 1 
        os.chdir('..')
        cluster = {'infomap': cluster}

        return cluster

    def to_pajek(self):
        
        os.chdir('Infomap')
    
        W = csc_matrix(self.Adj)
        row,col,data = extract.find(W)
        n = len(self.nodes)
        nodes = np.array(range(n))
        nodes = np.resize(nodes,(n,1))
        row = np.resize(row,(n,1))
        col = np.resize(col,(n,1))
        data = np.resize(data,(n,1))
        row = 1 + row # NB this is for pajek
        col = 1 + col # idem
        edges = np.concatenate((row,col,data),axis = 1)
        nodes = np.concatenate((nodes + 1,nodes),axis=1)

        filename = self.filename.split('.')[0]
        filename = '%s.net' % (filename)
        output = open(filename,'wb')
        output.write("*vertices %i\n"%(n))
        np.savetxt(output,nodes,fmt = ('%1.1i','%1.10s'))
        output.write("*arcs\n")
        np.savetxt(output,edges,fmt = ('%1.1i','%1.1i','%1.18e'))
    
    
        os.chdir('..')

    def Infomap(self, niter = 10,seed = 345234 , selfloops = False): 

        selfdict = {True: 'selflink' , False: ''}
        os.chdir('Infomap')
        filename = self.filename.split('.')[0]
        filename = '%s.net' %(filename)
        command = './infomap %i %s %i %s' % (seed,filename,niter,selfdict[selfloops])
        os.system(command)
        os.chdir('..')

    def makecover(self, partition,selfloops = False, q = ''):
        svnet = self.svnet.todense()

        if selfloops is False:
            indices = np.diag_indices_from(svnet)
            svnet[indices] = 0.
        
        svnet = csc_matrix(svnet)
            
        label,partition = partition.items()[0]
        M = community_matrix(partition)
        M = makecover(svnet, M)
        M = csc_matrix(M)
        i, j, mij = extract.find(M)
        cover = np.asarray(zip(i, j),dtype = [('nodes','S10'),('community',np.int)])
        filename = self.filename.split('.')[0]
        np.savetxt(filename + '_' + label + '_' + str(q) + '.cover', cover, fmt = ['%10s','%10i'])
        outfile = open(filename +'M' + '_' + label + '_' + str(q) + '.pkl','wb')
        dump(M,outfile)
        return cover,M
        


