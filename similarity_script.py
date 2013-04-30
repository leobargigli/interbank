# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 11:42:34 2012

@author: leonardo
"""

import interbank_classes as bdi
from optparse import OptionParser
from numpy import intersect1d,minimum,maximum,delete,savetxt,zeros,loadtxt,dot,array,diag_indices_from
from numpy.random import shuffle
from numpy.linalg import norm
from networkx import to_numpy_matrix
import os

USAGE = "%prog (date) (location) (rapporto) (maturity)"
USE_DESCRIPTION = "date is YYYYMMDD. \n location is 'dom', 'tot_adj' or 'tot_unadj'.\n rapporto is 'SECURED', 'UNSECURED' or 'TOT'. \n maturity is 'overnight', 'longterm', 'nonsignif' or 'TOT'."
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
    print dates
    
    
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

        select_maturity = list()
        try:
            j = maturity.index(i)
            select_maturity.append(maturity[j])
        except ValueError:
            pass
        if len(select_maturity) > 0:
            maturity = select_maturity

        select_location = list()
        try:
            j = location.index(i)
            select_location.append(location[j])
        except ValueError:
            pass
        if len(select_location) > 0:
            location = select_location



    print 'selected labels:'
    print dates, rapporto,maturity, location
    
    
    G = {}
    n = 0    
    for k in dates:
        G[k] = {}
        for h in location:
            G[k][h] = {}
            for i in rapporto:
                G[k][h][i] = {}
                for j in maturity:
                    filename = k + '_' + h.split('_')[0] + '.npy'
                    try:
                        if h.find('tot_adj') is not -1:
                            reporters = k + '_reporters.csv'
                            nodelist = loadtxt(reporters,dtype = str, delimiter = ',')
                            G[k][h][i][j] =  bdi.Year(filename, rapporto = i, maturity = j).Net.subgraph(nodelist)
                        else:
                            G[k][h][i][j] =  bdi.Year(filename,rapporto = i, maturity = j).Net
                    except AttributeError:
                        G[k][h][i][j] =  None
                    n += 1
                    
    
    print 'number of combinations: ' + str(n)

    size = len(dates) * len(location) * len(rapporto) * len(maturity) 
    J = zeros((size,size))
    JT = zeros((size,size))
    I = zeros((size,size))
    PJ = zeros((size,size))
    C = zeros((size,size))
    CT = zeros((size,size))
    PC = zeros((size,size))
    
    links = zeros((size,))
    nodes = zeros((size,))
        
        #print 'values:'
    
    n = 0
    labels = list()
    row_to_delete = list()
    col_to_delete = list()
    for k in dates:
        for h in location:
            for i in rapporto:
                for j in maturity:
                    try:
                        row_nodes = G[k][h][i][j].nodes()
                        labels.append(k + h + i + j + '\n')
                        m = 0
                        nodes[n] = len(row_nodes)
                        links[n] = G[k][h][i][j].number_of_edges()
                        for kb in dates:
                            for hb in location:
                                for ib in rapporto:
                                    for jb in maturity:
                                        try:
                                            col_nodes = G[kb][hb][ib][jb].nodes()
                                            nodelist = intersect1d(row_nodes,col_nodes)
                                            #print len(row_nodes),len(col_nodes),len(nodelist)
                                            A = to_numpy_matrix(G[k][h][i][j], weight = 'weight')
                                            B = to_numpy_matrix(G[kb][hb][ib][jb], weight = 'weight')
                                            indices = diag_indices_from(A)
                                            A[indices] = 0
                                            indices = diag_indices_from(B)
                                            B[indices] = 0

                                            FN = norm(B) * norm(A)
                                            try:
                                                BN = (1. * (A > 0)).sum() + (1. * (B > 0)).sum()
                                            except TypeError:
                                                BN = 0

                                            A = to_numpy_matrix(G[k][h][i][j], nodelist = nodelist, weight = 'weight')
                                            B = to_numpy_matrix(G[kb][hb][ib][jb], nodelist = nodelist, weight = 'weight')
                                            #print 'B is symmetric: ' + str((B == B.T).min()) 
                                            
                                            indices = diag_indices_from(A)
                                            A[indices] = 0
                                            B[indices] = 0
                                            
                                            A = array(A).flatten()
                                            BT = array(B.T).flatten()
                                            B = array(B).flatten()
                                            
                                            

                                            #print dot(A,B),N
                                            if FN > 0:
                                                C[n,m] = dot(A,B) / FN
                                                CT[n,m] = dot(A,BT) / FN
                                            else:
                                                C[n,m] = 0.
                                                CT[n,m] = 0.
                                                
                                            #print C[n,m]
                                            #x = zeros((10**3,))
                                            #for q in range(10**3):
                                            #    shuffle(B)
                                            #    x[q] = dot(A,B) / norm(A) / norm(B)
                                            #PC[n,m] = 1. - sum(x <= C[n,m]) / 10.**3
                                            
                                            A = 1. * (A > 0)
                                            BT = 1. * (BT > 0)
                                            B = 1. * (B > 0)
                                            
                                            if BN > 0:
                                                BTN = BN - (A*BT).sum()
                                                try:
                                                    JT[n,m] = (A * BT).sum() / BTN
                                                except TypeError:
                                                    JT[n,m] = 0
                                            else:
                                                JT[n,m] = 0
                                            #print A.sum(),B.sum()
                                            if BN > 0:
                                                BN -= (A*B).sum()
                                                try:
                                                    J[n,m] = (A * B).sum() / BN
                                                except TypeError:
                                                    J[n,m] = 0
                                                    
                                            else:
                                                J[n,m] = 0
                                            I[n,m] = len(nodelist)
                                            #x = zeros((10**3,))
                                            #for q in range(10**3):
                                            #    shuffle(B)
                                            #    x[q] = minimum(A,B).sum() / maximum(A,B).sum()
                                            #PJ[n,m] = 1. - sum(x <= J[n,m]) / 10.**3
                                                
                                        except AttributeError:
                                            col_to_delete.append(m)
                                        m += 1
                    except AttributeError:
                        row_to_delete.append(n)
                    
                    n+=1
    
    to_delete = intersect1d(row_to_delete,col_to_delete)    
    #print 'deleted row / cols: ' + str(len(to_delete))
    
    J = delete(J,to_delete,0)
    J = delete(J,to_delete,1)

    JT = delete(JT,to_delete,0)
    JT = delete(JT,to_delete,1)

    I = delete(I,to_delete,0)
    I = delete(I,to_delete,1)

    C = delete(C,to_delete,0)
    C = delete(C,to_delete,1)

    CT = delete(CT,to_delete,0)
    CT = delete(CT,to_delete,1)

    
    #PJ = delete(PJ,to_delete,0)
    #PJ = delete(PJ,to_delete,1)

    #PC = delete(PC,to_delete,0)
    #PC = delete(PC,to_delete,1)
    
    
    nodes = delete(nodes,to_delete)
    links = delete(links,to_delete)

    filename = ''

    if len(args) > 0:
        for i in args:
            filename += i
    else:
        filename += 'total'
    
    try:
        os.chdir('similarity_results')
    except OSError:
        os.mkdir('similarity_results')
        os.chdir('similarity_results')
        
    
    savetxt(filename + '.wmatrix',C,fmt = '%.4f')
    savetxt(filename + '.wTmatrix',CT,fmt = '%.4f')
    savetxt(filename + '.jTmatrix',JT,fmt = '%.4f')
    savetxt(filename + '.jmatrix',J,fmt = '%.4f')
    savetxt(filename + '.intersection',I,fmt = '%.4f')
    #savetxt(filename + '.jpvalues',PJ,fmt = '%.4f')
    #savetxt(filename + '.wpvalues',PC,fmt = '%.4f')
    savetxt(filename + '.nodes',nodes,fmt = '%.4f')
    savetxt(filename + '.links',links,fmt = '%.4f')
    
    output = open(filename + '.labels','wb')
    output.writelines(labels)
    output.flush()    
    
    #print 'shape of similarity matrix:'
    #print J.shape
                                    
    os.chdir('..')
                    
    
    
if __name__ == "__main__":
    main()