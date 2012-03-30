
import numpy as np
import os
import sys
from optparse import OptionParser
from plfit import plfit
from plplot import plplot
from interbank_classes import *



USAGE = "%prog --bip --imgfmt FILENAME METHOD"
USE_DESCRIPTION = "The FILENAME file must be a weighted edgelist."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity = None )
parser.add_option( '--fmt',
    dest = "fmt",
    help = "fileformat for graphs", 
)

parser.add_option( '--bip',
    dest = "bip",
    action = 'store_const',
    const = 1,
    help = "fileformat for graphs", 
)






def main():
    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    filename = args[0]


    Y = Year(filename)
    SV = SVnet(Y)
    

    method = args[1]
    
    fmt = opts.fmt
    if fmt is None:
        fmt = 'png'
    
    label = ''

    try:
        os.chdir(Y.filename + label +'_comm_detect')
    except OSError:
        os.mkdir(Y.filename  + label + '_comm_detect')
        os.chdir(Y.filename  + label + '_comm_detect')
 
    

    if method == 'infomap':
        
        SV.to_pajek()
        SV.Infomap()
        cluster = SV.from_pajek()
    
    elif method == 'spectral':
        
        W = Y.Adj
        q, pdiff,p, svs = n_of_modules_markov(W)
        cluster = spectral_partition(W,q)
        plot_svs(svs,pdiff,q,filename, fmt = fmt)
        
    
    cover,M = SV.makecover(cluster)
    Y.saveDigraph()
    SV.edgelist()
    SV.make_Digraph()
    community_stats(filename)
    
    
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


