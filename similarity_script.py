# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 11:42:34 2012

@author: leonardo
"""

import interbank_classes as bdi
from optparse import OptionParser
from numpy import intersect1d,savetxt,zeros,loadtxt,dot,array,diag_indices_from,tril_indices_from
from scipy.linalg import svd
from numpy.random import shuffle
from numpy.linalg import norm
from networkx import to_numpy_matrix,weakly_connected_component_subgraphs
import os
from comm_detect_funcs import Kmatrix

USAGE = "%prog (date) (location) (rapporto) (maturity)"
USE_DESCRIPTION = "date is YYYYMMDD. \n location is 'dom', 'tot_adj' or 'tot_unadj'.\n rapporto is 'SECURED', 'UNSECURED' or 'TOT'. \n maturity is 'overnight', 'medium', longterm', 'nonsignif' or 'TOT'."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity = None)


def main():
    
    args = parser.parse_args()[1]
    
    rapporto = [
    'SECURED',
    'UNSECURED',
    'TOT']
        
    maturity = [
    'overnight',
    'longterm',
    'nonsignif',
    'medium',
    'TOT']
    
    location = [
    'dom',
    'tot_adj',
    'tot_unadj']
    
    dates = open('dates.csv', 'r')
    dates = dates.readlines()
    #dates.pop(0)
    dates = [j.replace('"','') for j in dates]
    dates = [j.replace('\n','') for j in dates]
    #dates = [j.replace('\r','') for j in dates]
    #print dates

    foldername = ''

    if len(args) > 0:
        for i in args:
            foldername += i
    else:
        foldername += 'total'
    
    
    for i in args:

        select_dates = list()
        try:
            j = dates.index(i)
            select_dates.append(dates[j])
        except ValueError:
            pass
        if len(select_dates) > 0:
            dates = select_dates

        select_rapporto = list()
        try:
            j = rapporto.index(i)
            select_rapporto.append(rapporto[j])
        except ValueError:
            pass
        if len(select_rapporto) > 0:
            rapporto = select_rapporto
        if len(rapporto) > 1:
            rapporto = [j for j in rapporto if j <> 'TOT']            

        select_maturity = list()
        try:
            j = maturity.index(i)
            select_maturity.append(maturity[j])
        except ValueError:
            pass
        if len(select_maturity) > 0:
            maturity = select_maturity
        if len(maturity) > 1:
            maturity = [j for j in maturity if j <> 'TOT']            

        select_location = list()
        try:
            j = location.index(i)
            select_location.append(location[j])
        except ValueError:
            pass
        if len(select_location) > 0:
            location = select_location



    G = {}
    labels = list()

    for k in dates:
        for h in location:
            for i in rapporto:
                for j in maturity:
                    filename = k + '_' + h.split('_')[0] + '.npy'
                    try:
                        if h.find('tot_adj') is not -1:
                            reporters = k + '_reporters.csv'
                            nodelist = loadtxt(reporters,dtype = str, delimiter = ',')
                            G[k + h + i + j] =  bdi.Year(filename, rapporto = i, maturity = j).Net.subgraph(nodelist)
                        else:
                            G[k + h + i + j] =  bdi.Year(filename,rapporto = i, maturity = j).Net
                            #comps = weakly_connected_component_subgraphs(G[k + h + i + j])
                            #if len(comps) > 0:
                            #    W = to_numpy_matrix(comps[0])
                            #    K = Kmatrix(W)
#                                sv = svd(K,0,0)
#                                try:
#                                    os.chdir('similarity_results' + foldername)
#                                except OSError:
#                                    os.mkdir('similarity_results' + foldername)
#                                    os.chdir('similarity_results' + foldername)
#
#                                savetxt(k + h + i + j + '_svs',sv)
#                                os.chdir('..')
                            
                        labels.append(k + h + i + j)
                    except AttributeError:
                        pass
    
    to_remove = [
    'SECUREDovernight',
    'SECUREDlongterm',
    ]
    
    for i in to_remove:
        labels = [j for j in labels if j.find(i) == -1 or j.find('UNSECURED') <>-1]
        
    n = len(labels) 
   
    print 'number of combinations: ' + str(n)
    
    print 'selected labels:'
    print labels
    
    I = zeros((n,n))
    J = zeros((n,n))
    PJ = zeros((n,n))
    JT = zeros((n,n))
    PJT = zeros((n,n))
    C = zeros((n,n))
    PC = zeros((n,n))
    CT = zeros((n,n))
    PCT = zeros((n,n))    
    
    links = zeros((n,))
    nodes = zeros((n,))
    
    r = n * (n - 1) / 2. + n 
    for k in range(n):
        try:
            x = labels[k]
            row_nodes = G[x].nodes()
            nodes[k] = len(row_nodes)
            links[k] = G[x].number_of_edges()
            for h in range(k,n):
                try:
                    y = labels[h]
                    col_nodes = G[y].nodes()
                    nodelist = intersect1d(row_nodes,col_nodes)
#                    A = to_numpy_matrix(G[x], weight = 'weight')
#                    B = to_numpy_matrix(G[y], weight = 'weight')
#                    indices = diag_indices_from(A)
#                    A[indices] = 0
#                    indices = diag_indices_from(B)
#                    B[indices] = 0

#                    FN = norm(B) * norm(A)
#                    try:
#                        BN = (1. * (A > 0)).sum() + (1. * (B > 0)).sum()
#                    except TypeError:
#                        BN = 0

                    A = to_numpy_matrix(G[x], nodelist = nodelist, weight = 'weight')
                    B = to_numpy_matrix(G[y], nodelist = nodelist, weight = 'weight')
                   
                    indices = diag_indices_from(A)
                    A[indices] = 0
                    B[indices] = 0


                    FN = norm(B) * norm(A)
                    try:
                        BN = (1. * (A > 0)).sum() + (1. * (B > 0)).sum()
                    except TypeError:
                        BN = 0

                    
                    A = array(A).flatten()
                    BT = array(B.T).flatten()
                    B = array(B).flatten()

                    if FN > 0:
                        C[k,h] = dot(A,B) / FN
                        CT[k,h] = dot(A,BT) / FN
                    
                    if BN > 0:
                        BTN = BN - ((A>0)*(BT>0)).sum()
                        try:
                            JT[k,h] = ((A>0)*(BT>0)).sum() / BTN
                        except TypeError:
                            pass

                        BN -= ((A>0)*(B>0)).sum()
                        
                        try:
                            J[k,h] = ((A>0) * (B>0)).sum() / BN
                        except TypeError:
                            pass

                        z = zeros((10**3,))
                        t = zeros((10**3,))

                        for q in range(10**3):
                            if divmod(q,10**2)[1] == 0:
                                print q
                            shuffle(B)
                            t[q] = ((A>0) * (B>0)).sum() / BN
                            
                            if FN > 0:
                                z[q] = dot(A,B) / FN
                            
                        PJ[k,h] = 1. - sum(t <= J[k,h]) / 10.**3
                        PJT[k,h] = 1. - sum(t <= JT[k,h]) / 10.**3
                        
                        if FN > 0:
                            
                            PC[k,h] = 1. - sum(z <= C[k,h]) / 10.**3
                            PCT[k,h] = 1. - sum(z <= CT[k,h]) / 10.**3

                            
                    I[k,h] = len(nodelist)

                        
                except AttributeError:
                    pass
                r -= 1
                print 'steps left: ' + str(r)


        except AttributeError:
            pass
        

    
    indices = tril_indices_from(C)
    C[indices] = C.T[indices]
    J[indices] = J.T[indices]
    JT[indices] = JT.T[indices]
    I[indices] = I.T[indices]
    PJ[indices] = PJ.T[indices]
    PJT[indices] = PJT.T[indices]
    CT[indices] = CT.T[indices]
    PC[indices] = PC.T[indices]
    PCT[indices] = PCT.T[indices]
    
    try:
        os.chdir('similarity_results' + foldername)
    except OSError:
        os.mkdir('similarity_results' + foldername)
        os.chdir('similarity_results' + foldername)

        
    savetxt(foldername + '.wmatrix',C,fmt = '%.4f')
    savetxt(foldername + '.wTmatrix',CT,fmt = '%.4f')
    savetxt(foldername + '.jTmatrix',JT,fmt = '%.4f')
    savetxt(foldername + '.jmatrix',J,fmt = '%.4f')
    savetxt(foldername + '.intersection',I,fmt = '%.4f')
    savetxt(foldername + '.jpvalues',PJ,fmt = '%.4f')
    savetxt(foldername + '.wpvalues',PC,fmt = '%.4f')
    savetxt(foldername + '.jTpvalues',PJT,fmt = '%.4f')
    savetxt(foldername + '.wTpvalues',PCT,fmt = '%.4f')
    savetxt(foldername + '.nodes',nodes,fmt = '%.4f')
    savetxt(foldername + '.links',links,fmt = '%.4f')
    
    output = open(foldername + '.labels','wb')
    output.writelines([i +'\n'  for i in labels])
    
                                    
    os.chdir('..')
                    
    
    
if __name__ == "__main__":
    main()