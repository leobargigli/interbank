import os
import sys
from optparse import OptionParser
from numpy import loadtxt

USAGE = "%prog [-d(omestic)] FILENAME"
USE_DESCRIPTION = "The FILENAME file must be a csv file."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity=None )
parser.add_option( "-d" , "--d",
    action="store_const",
    const= 1, 
    dest="domestic",
    help="Only domestic exposures are included."
)

def main():

    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    filename = args[0]
    optdict = {1: 'dom', None: 'tot'}
    tablename = os.path.splitext(filename)[0]

    optdict = {1: 'dom', None: 'tot'}

    selfloops = {0:'',1:'_sl'}
    for j in range(2):
        countfile = tablename + '_'+ optdict [opts.domestic] +selfloops[j]+'.count'
        sumfile = tablename + '_'+ optdict [opts.domestic] +selfloops[j]+'.sum'
        where = {1: 'WHERE location = "domest" AND infragruppo = "%i" '% (j), 
             None: 'WHERE NATURA_RAPPORTO LIKE "%%DEBITI%%" AND infragruppo = "%i" '% (j)} 
        command  = 'python querycsv.py -i %s -o %s \'SELECT DATA_CONTABILE,NATURA_RAPPORTO,maturity,sum(importo) from %s %s GROUP BY DATA_CONTABILE,NATURA_RAPPORTO,maturity\'' %(filename,sumfile, tablename, where[opts.domestic])
        os.system(command)
        command  = 'python querycsv.py -i %s -o %s \'SELECT DATA_CONTABILE,NATURA_RAPPORTO,maturity,count(importo) from %s %s GROUP BY DATA_CONTABILE,NATURA_RAPPORTO,maturity\'' %(filename,countfile, tablename, where[opts.domestic])
        os.system(command)

        
        #~ output = open(tablefile, 'r')
        #~ lines = output.readlines()
        #~ lines.pop(0)
    #~ 
        #~ if opts.domestic is None:
#~ 
            #~ command  = 'python querycsv.py -i %s -o %s \'SELECT ctp_controp,CTP_CAPOGRU,sum(importo) from %s WHERE DATA_CONTABILE = "%s" AND NATURA_RAPPORTO LIKE "%%IMPIEGHI%%" GROUP BY ctp_controp,CTP_CAPOGRU\'' %(filename,edgefile, tablename, i)
            #~ os.system(command)
            #~ output = open(edgefile, 'r')
            #~ lines2 = output.readlines()
            #~ lines2.pop(0)
            #~ for  j in lines2:          
                #~ lines.append(j)
#~ 
        #~ for j in xrange(len(lines)):
            #~ lines[j] = lines[j].replace('"', '')
        #~ 
        #~ output = open(edgefile, 'wb')
        #~ output.writelines(lines)
       #~ output.flush()
        
        
    

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
    
