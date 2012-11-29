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

    def __init__(self, filename, delimiter = ','):

        edgetype = np.dtype([('source','|S10'),
                                 ('dest','|S10' ), 
                                 ('weight',np.float32),
                                 ('location','|S10')])
            
#        edgetype = np.dtype([('source','|S10'),
#                                 ('dest','|S10' ), 
#                                 ('weight',np.float32)])
        
        try:
            edgelist = np.loadtxt(filename,
                                  dtype = edgetype) 
        except IndexError:
            edgelist = np.loadtxt(filename, 
                                  delimiter = ',',
                                  dtype = edgetype) 

        reporters = np.unique(edgelist['source'])        
        counterparts = np.unique(edgelist['dest'])
        non_reporters = np.setdiff1d(counterparts,reporters)

         # this is to clean fake domestic relations
        
        try:
            
            dom_indices = np.where(edgelist['location'] == 'domest')[0]
            
            for i in non_reporters:
            
                nr_indices = np.where(edgelist['dest'] == i)[0]
                dom_nr = np.intersect1d(dom_indices,nr_indices)
                edgelist = np.delete(edgelist,dom_nr,0)
        
        except ValueError:
            pass
            
            
            # this is to treat foreign subsidiaries as a separate node
            
        foreign_links = np.where(edgelist['location'] == 'estero')[0]
        
        for i in foreign_links:
            if edgelist[i]['source'] == edgelist[i]['dest']:
                edgelist[i]['dest'] += 'F'
            
        # this is to sum weights across different link types
        
        weights = edgelist['weight']
        edges = np.array([i['source'] +'*'+ i['dest'] for i in edgelist])
        uni_edges = np.unique(edges)
        q = len(uni_edges)

        self.edgetype = np.dtype([('source','|S10'),
                             ('dest','|S10' ), 
                             ('weight',np.float64)])

        
        edgelist = np.zeros((q,), dtype = self.edgetype)

        for i in range(q):
            indices = np.where(edges == uni_edges[i])
            uni_weight = weights[indices].sum() 
            source,dest = tuple(uni_edges[i].split('*'))
            edgelist[i] = (source,dest,uni_weight)

        G = nx.DiGraph()
        G.add_weighted_edges_from(edgelist)

        self.Net = G
        self.nodes = np.array((np.sort(G.nodes())),dtype = str)
        self.Adj = csc_matrix(nx.to_numpy_matrix(G, nodelist = self.nodes))
        self.filename = os.path.splitext(filename)[0]
        self.edgelist = edgelist
        


    def saveDigraph(self):
        G = self.Net
        outfile = open(self.filename.split('.')[0] + 'Graph.pkl','wb')
        dump(G,outfile)

   
    def stats(self, distG,nbunch = None):
        
        if nbunch is not None:
            label = '_adj'
        else:
            label = '_unadj'
            
        if  self.filename.find('dom') <>-1:
            label = ''
            
        G = self.Net
        
        if nbunch is None:
            nbunch = G.nodes()
            
        nodes = len(nbunch)
        edges = G.subgraph(nbunch).number_of_edges()
        selfloops = G.subgraph(nbunch).selfloop_edges()
        edges = edges - len(selfloops)
        d = density(G , nbunch = nbunch)
        
        try:
            volume = G.subgraph(nbunch).size(weight = 'weight') #self.Adj.sum()
        except TypeError:
            volume = G.subgraph(nbunch).size(weighted = True)
        
        comps = nx.weakly_connected_component_subgraphs(G.subgraph(nbunch))
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
            out_tau, out_p = (None,None)
            in_tau, in_p = (None, None)
            tot_tau, tot_p = (None,None)

        except TypeError:
            out_tau, out_p = (None, None)
            in_tau, in_p = (None, None)
            tot_tau, tot_p = (None,None)
    
        dC = np.average(dir_clustering_coefficient(G, nbunch = nbunch)[0])
        uC = np.average(clustering_coefficient(G, nbunch = nbunch)[0])
        
        G.remove_edges_from(selfloops)

        if len(comps) == 1:
            avg_path_length = nx.average_shortest_path_length(G)
            try:
                avg_weight_path_length = nx.average_shortest_path_length(G,weighted = True)
            except TypeError:
                avg_weight_path_length = nx.average_shortest_path_length(G,weight = 'weight')
                
        elif len(comps[0]) > 1:
            avg_path_length = nx.average_shortest_path_length(comps[0])
            try:
                avg_weight_path_length = nx.average_shortest_path_length(comps[0],weighted = True)
            except TypeError:
                avg_weight_path_length = nx.average_shortest_path_length(comps[0],weight = 'weight')
        
        else:
            avg_path_length = 0
            avg_weight_path_length = 0
            
        

        filename = self.filename + label + '_stats.dat'
        output = open(filename, 'wb')
        output.write('# of nodes: %i\n'%(nodes)) 
        output.write('# of edges (net): %i\n'%(edges)) 
        output.write('# of selfloops: %i\n'%(len(selfloops))) 
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
        output.write('Degree reciprocity (p-value): %f (%f)\n'%(recip)) 
        output.write('Weight reciprocity (p-value): %f (%f)\n'%(w_recip))
        output.write('Average directed clustering: %f\n'%(dC)) 
        output.write('Average undirected clustering: %f\n'%(uC)) 
        if out_tau is not None and out_p is not None:
            output.write('Kendall tau w_out / d_out vs d_out (p-value): %f (%f)\n' % (out_tau, out_p))
        if in_tau is not None and in_p is not None:
            output.write('Kendall tau w_in / d_in vs d_in (p-value): %f (%f)\n'% (in_tau, in_p))
        if tot_tau is not None and tot_p is not None:
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
        #i, j, mij = extract.find(M)
        #cover = np.asarray(zip(i, j),dtype = [('nodes','S10'),('community',np.int)])
        #filename = self.filename.split('.')[0]
        #np.savetxt(filename + '_' + label + '_' + str(q) + '.cover', cover, fmt = ['%10s','%10i'])
        #outfile = open(filename +'M' + '_' + label + '_' + str(q) + '.pkl','wb')
        #dump(M,outfile)
        return M#cover,M
        


