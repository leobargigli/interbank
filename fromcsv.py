# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:06:33 2013

@author: leonardo
"""

import os
import sys
from optparse import OptionParser
from numpy import *

USAGE = "%prog [-d(omestic)] [-f(oreign)] FILENAME"
USE_DESCRIPTION = "The FILENAME file must be a csv file."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity = None )

parser.add_option( "-d" , "--d",
    action = "store_const",
    const = 1, 
    dest ="domestic",
    help = "Only domestic exposures are included."
)

parser.add_option( "-f" , "--f",
    action = "store_const",
    const = 2, 
    dest ="foreign",
    help = "Only foreign exposures are included."
)


def main():

    ( opts , args ) = parser.parse_args()
    filename = args[0]
    edgetype = dtype([('date','|S10'),
                      ('weight',float32),
                    ('source','|S10'),
                    ('dest','|S10' ),
                    ('natura_rapporto','|S10'),
                    ('cred_deb','|S10'),
                    ('infragruppo',int32),
                    ('location','|S10'),
                    ('statoctp','|S10'),
                    ('maturity','|S10')])

    try:
        edgelist = loadtxt(filename, delimiter = ',',dtype = edgetype)
    except IndexError:
        edgelist = loadtxt(filename, delimiter = ';',dtype = edgetype)
        
    dates = unique(edgelist['date'])
    savetxt('dates.csv',dates,fmt ='%s')
    

    for i in dates:
        
        reporterfile = i + '_reporters.csv'
        indices = where(edgelist['date'] == i)
        y_edgelist = edgelist[indices]
        reporters = unique(y_edgelist['source'])
        print i + '_totreporters:' + str(len(reporters))
        savetxt(reporterfile,reporters,fmt ='%s')
        
    
    
    if opts.domestic == 1:
        
        indices = where(edgelist['location'] == 'domest')
        edgelist = edgelist[indices]
        indices = where(edgelist['cred_deb'] == 'DEBITI')
        edgelist = edgelist[indices]
        opt = 'dom'
        

    elif opts.foreign == 2:
        indices = where(edgelist['location'] == 'estero')
        edgelist = edgelist[indices]
        indices = where(edgelist['cred_deb'] == 'IMPIEGHI')
        source = edgelist['dest'][indices].copy()
        dest = edgelist['source'][indices].copy()
        edgelist['source'][indices] = source
        edgelist['dest'][indices] = dest
        opt = 'for'
        
    else:
        
        indices_l = where(edgelist['location'] == 'domest')
        indices_d = where(edgelist['cred_deb'] == 'IMPIEGHI')
        indices = intersect1d(indices_l[0],indices_d[0])
        edgelist = delete(edgelist,indices)
        indices = where(edgelist['cred_deb'] == 'IMPIEGHI')
        source = edgelist['dest'][indices].copy()
        dest = edgelist['source'][indices].copy()
        edgelist['source'][indices] = source
        edgelist['dest'][indices] = dest
        opt = 'tot'
        
    for i in dates:
        
        indices = where(edgelist['date'] == i)
        y_edgelist = edgelist[indices]
        filename = i + '_' + opt 
        print i,len(y_edgelist),y_edgelist['weight'].sum()
        save(filename,y_edgelist)
                                
    
if __name__ == "__main__":
    main()

    
