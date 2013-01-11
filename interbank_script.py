from interbank_classes import *
import numpy as np
import os
import sys
from optparse import OptionParser
from plfit import plfit
from plplot import plplot
from basic_stats_funcs import *


USAGE = "%prog --n(odes) NODEFILE --m --r FILENAME "
USE_DESCRIPTION = "The FILENAME file must be a weighted edgelist."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity=None )
parser.add_option( "--nodes", '--n',
    dest = "nodelist",
    help = "a list of nodes for the analysis", 
)

parser.add_option( "--maturity", '--m',
    dest = "maturity",
    help = "must be either longterm or overnight", 
)

parser.add_option( "--rapporto", '--r',
    dest = "rapporto",
    help = "must be either SECURED or UNSECURED", 
)


def main():
    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    filename = args[0]

    if opts.nodelist:
        nodelist = np.loadtxt(opts.nodelist, 
                              dtype = str,
                              delimiter = ',')
        label = '_adj'
        
    else:
        nodelist = None
        label = '_unadj'
    
    if filename.find('dom') <> -1:
        label = ''

    if opts.rapporto is None:
        opts.rapporto = 'total'
    
    if opts.maturity is None:
        opts.maturity = 'total'

    
    Y = Year(filename,
             rapporto = opts.rapporto, 
             maturity = opts.maturity)

    # this is to clean from fake domestic links    
    if nodelist is not None:
        nodelist = np.intersect1d(nodelist,Y.nodes)
    
   
    try:
        os.chdir(Y.filename + opts.rapporto + opts.maturity + label + '_stats')
    except OSError:
        os.mkdir(Y.filename + opts.rapporto + opts.maturity + label + '_stats')
        os.chdir(Y.filename + opts.rapporto + opts.maturity + label + '_stats')
  
    G = Y.Net
    
    distG = dists(G, nbunch = nodelist)
    Y.stats(distG,nbunch = nodelist)

#    if len(args) > 1:
#        img = args[-1]
#    else:
#        img = 'png'

#    fmts = {
#    'gross out-weight': '%.2f', 
#    'net out-weight': '%.2f', 
#    'in-degree': '%.1i', 
#    'gross in-weight': '%.2f', 
#    'net in-weight': '%.2f', 
#    'out-degree':'%.1i', 
#    'gross cells': '%.2f',
#    'net cells': '%.2f'
#     }
#    #print distG
#    
#    filename = Y.filename + label + opts.rapporto + opts.maturity + '_stats.dat'
#    
#    output = open(filename, 'a')
#    output.write('\n')    
#    
#    for i in distG:
#        
#        x = distG[i]
#
#        
#        distfile = Y.filename + '_' + i + label + opts.rapporto + opts.maturity + '.distr'
#        np.savetxt(distfile, x, fmt = fmts[i])
#        if i <>'net cells' and i <> 'gross cells':
#            knnk = k_vs_nnk(G, i, nbunch = nodelist)
#            knnk = knnk.items()
#            knnk = np.array(knnk)
#            k = knnk[:, 0]
#            nnk = knnk[:, 1]
#            labels = [r'$k$',r'$k_{nn}$','Average_neighbor_%s_vs_node_%s'%(i, i)]
#            scatter(k, nnk, labels, opts.rapporto + opts.maturity + Y.filename + label,fmt = img)
#
#        indices = np.where(x == 0)
#        x = np.delete(x, indices)
#        x = list(x)
#        
#        try:
#            alpha,xmin,ntail = plfit(x)
#        except ValueError:
#            pass
#        try:
#            h = plplot(x,xmin,alpha, ntail, i, opts.rapporto + opts.maturity + Y.filename + label, format = img)
#            output.write('estimated exponent of %s distribution: %.2f\n'%(i, alpha))
#            output.write('# datapoints in the tail of %s distribution: %i\n'%(i , ntail))
#            output.write('min %s value in the tail: %.2f\n'%(i, xmin))
#
#        except TypeError:
#            pass
#        except UnboundLocalError:
#            pass
#        except ValueError:
#            pass
#    
#    options = [
#        ('out',True,0),
#        ('out',False,50), 
#        ('in',True,0),
#        ('in',False,50)]
#        
#    for i,j,k in options:
#        try:
#            part_dict = exp_part_ratio(G, i, nbunch = nodelist, degree = j, quant = k)
#            part_dict = part_dict.items()
#            part_dict = np.array(part_dict)
#            epart = part_dict[:, 0]
#            part = part_dict[:, 1]
#            labeldict = {True: 'inverse degree', False: 'weight'}
#            ydict = {True: r'$1 / k$', False : ''} 
#            labels = [ydict[j],r'Part. ratio',r'Participation_ratio_vs_%s%s' % (i,labeldict[j])]
#            scatter(epart, part, labels, opts.rapporto + opts.maturity + Y.filename + label, diag = True,fmt = img)
#
#        except IndexError:
#            pass
#        except ValueError:
#            pass
#            
#        
#    clust_vs_degree(G, opts.rapporto + opts.maturity + Y.filename, format = img, nbunch = nodelist)
#    clust_vs_degree(G, opts.rapporto + opts.maturity + Y.filename, format = img, nbunch = nodelist, directed = True)
    
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


