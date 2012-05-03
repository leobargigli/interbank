
import numpy as np
import os
import sys
from optparse import OptionParser
from plfit import plfit
from plplot import plplot
from interbank_classes import *



USAGE = "%prog --fmt --sl FILENAME METHOD"
USE_DESCRIPTION = "The FILENAME file must be a weighted edgelist."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity = None )
parser.add_option( '--fmt',
    dest = "fmt",
    help = "fileformat for graphs", 
)

parser.add_option( '--sl',
     dest = "sl",
     action = 'store_const',
     const = True,
     help = "selfloops included", 
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

    selfloops = opts.sl

    if selfloops is None:
        selfloops = False
    print selfloops

    selfloopdict = { 1: '_withselfloops', 0: ''}
    
    try:
        os.mkdir(Y.filename  + label + selfloopdict[selfloops]+'_comm_detect')
    except OSError:
        pass

    if method == 'infomap':
        
        SV.to_pajek()
        SV.Infomap(selfloops = selfloops)
        cluster = SV.from_pajek()
        os.chdir(Y.filename + label + selfloopdict[selfloops] +'_comm_detect')
    
    elif method == 'spectral':
        
        W = Y.Adj
        if selfloops is False:
            W = W.todense()
            indices = np.diag_indices_from(W)
            W[indices] = 0.
            W = csc_matrix(W)
            
        q, pdiff,p, svs = n_of_modules_markov(W)
        cluster = spectral_partition(W,q)
        plot_svs(svs,pdiff,q,filename, fmt = fmt)
        os.chdir(Y.filename + label +selfloopdict[selfloops] +'_comm_detect')
    
    cover,M = SV.makecover(cluster)
    Y.saveDigraph()
    SV.edgelist()
    SV.make_Digraph()
    community_stats(filename,method)
    
    
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


