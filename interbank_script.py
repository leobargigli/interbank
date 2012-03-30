from interbank_classes import *
import numpy as np
import os
import sys
from optparse import OptionParser
from plfit import plfit
from plplot import plplot



USAGE = "%prog --n(odes) NODEFILE FILENAME "
USE_DESCRIPTION = "The FILENAME file must be a weighted edgelist."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity=None )
parser.add_option( "--nodes", '--n',
    dest = "nodelist",
    help = "a list of nodes for the analysis", 
)


def main():
    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    filename = args[0]
    
    
    
    if opts.nodelist:
        nodelist = np.loadtxt(opts.nodelist, dtype = str,delimiter = ',')
        label = '_adj'
    else:
        nodelist = None
        label = '_unadj'
    
    if filename.find('dom') <> -1:
        label = ''
    
    Y = Year(filename)

    try:
        os.chdir(Y.filename + label +'_stats')
    except OSError:
        os.mkdir(Y.filename  + label + '_stats')
        os.chdir(Y.filename  + label + '_stats')
  
    Y.stats(nbunch = nodelist)
    G = Y.Net
    if len(args) > 1:
        img = args[-1]
    else:
        img = 'png'

    filename = Y.filename + label + '_stats.dat'
    output = open(filename, 'a')
    output.write('\n')
    distG = dists(G, nbunch = nodelist)
    fmts = {'out-weight': '%.2f', 'in-degree': '%.1i', 'in-weight': '%.2f', 'out-degree':'%.1i', 'cells': '%.2f'}
   
    for i in distG:
        
        x = distG[i]
        distfile = Y.filename + '_' + i  + label +'.distr'
        np.savetxt(distfile, x, fmt = fmts[i])
        if i <>'cells':
            knnk = k_vs_nnk(G, i, nbunch = nodelist)
            knnk = knnk.items()
            knnk = np.array(knnk)
            k = knnk[:, 0]
            nnk = knnk[:, 1]
            labels=[r'$k$',r'$k_{nn}$','Average neighbor %s vs node %s'%(i, i)]
            scatter(k, nnk, labels, Y.filename + label,fmt = img)

        indices = np.where(x == 0)
        x = np.delete(x, indices)
        x = list(x)
        try:
            alpha,xmin,ntail = plfit(x)
        except ValueError:
            pass
        try:
            h = plplot(x,xmin,alpha, ntail, i, Y.filename + label, format = img)
            output.write('estimated exponent of %s distribution: %.2f\n'%(i, alpha))
            output.write('# datapoints in the tail of %s distribution: %i\n'%(i , ntail))
            output.write('min %s value in the tail: %.2f\n'%(i, xmin))

        except TypeError:
            pass
    
    options = [('out',True,0),('out',False,50), ('in',True,0),('in',False,50)]
    for i,j,k in options:
        part_dict = participation_ratio(G, i, nbunch = nodelist, degree = j, quant = k)
        part_dict = part_dict.items()
        part_dict = np.array(part_dict)
        epart = part_dict[:, 0]
        part = part_dict[:, 1]
        labeldict = {True: 'inverse degree', False: 'weight'}
        ydict = {True: r'$1 / k$', False : ''} 
        labels = [ydict[j],r'Part. ratio',r'Participation ratio vs %s (%s)' % (labeldict[j],i)]
        scatter(epart, part, labels, Y.filename + label, diag = True,fmt = img)

            
        
    clust_vs_degree(G, Y.filename, format = img, nbunch = nodelist)
    clust_vs_degree(G, Y.filename, weight = 'weight', format = img, nbunch = nodelist)
    
    os.chdir('..')
        
    
  
    
    
def __validate_opts_and_args( opts , args ):
    """Makes sure that one file has been passed and that it exists.
    """
    if len(args)<1:
        parser.print_help()
        sys.exit( 1 )
    if not os.path.exists( args[0] ): 
        parser.print_help()
        sys.stderr.write( "FILE {0} DOES NOT EXISTS\n".format(args[0]) )
        sys.exit( 1 )

if __name__ == "__main__":
    main()


