import os
import sys
from optparse import OptionParser
from numpy import loadtxt,setdiff1d

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
    __validate_opts_and_args( opts , args )
    filename = args[0]
    tablename = os.path.splitext(filename)[0]
    datefile = 'date.csv'

    if opts.domestic is None and opts.foreign is not None:
        opts.domestic = 0
        
    if opts.foreign is None and opts.domestic is not None:
        opts.foreign = 0
        
    if opts.foreign is None and opts.domestic is None:
        opts.domestic = 1
        opts.foreign = 2
        
    optdict = {1: 'dom', 2: 'for' , 3: 'tot'}
    
    opt = opts.domestic + opts.foreign # possible values: 1 (domestic) ; 2 (foreign); 3 (total) 
    
    command  = 'python querycsv.py -i %s -o %s \"SELECT DATA_CONTABILE from %s GROUP BY DATA_CONTABILE\"' %(filename,datefile, tablename)
    os.system(command)
    dates = open(datefile, 'r')
    dates = dates.readlines()
    dates.pop(0)
    dates = [j.replace('"','') for j in dates]
    dates = [j.replace('\r\n','') for j in dates]

    for i in dates:
        
        reporterfile = i + '_reporters.csv'
        command  = 'python querycsv.py -i %s -o %s \"SELECT CTP_CAPOGRU from %s WHERE DATA_CONTABILE = \'%s\' GROUP BY CTP_CAPOGRU\"' % (filename , reporterfile , tablename, i)
        os.system(command)
        output = open(reporterfile, 'r')
        lines = output.readlines()
        lines.pop(0)

        lines = [j.replace('"', '') for j in lines]

        output = open(reporterfile, 'wb')
        output.writelines(lines)
        output.flush()
        
        edgefile = i + '_'+ optdict[opt] + '.edgelist'
        where = {
                1: 'WHERE DATA_CONTABILE = \'%s\' AND location = \'domest\' AND (NATURA_RAPPORTO = \'DEBITI UNSECURED\' OR NATURA_RAPPORTO = \'DEBITI SECURED\')'% (i), 
                3: 'WHERE DATA_CONTABILE = \'%s\' AND (NATURA_RAPPORTO = \'DEBITI UNSECURED\' OR NATURA_RAPPORTO = \'DEBITI SECURED\')'% (i),
                2: 'WHERE DATA_CONTABILE = \'%s\' AND location = \'estero\' AND (NATURA_RAPPORTO = \'DEBITI UNSECURED\' OR NATURA_RAPPORTO = \'DEBITI SECURED\')'% (i)
                 } 
        command  = 'python querycsv.py -i %s -o %s \"SELECT CTP_CAPOGRU,ctp_controp,sum(importo),location from %s %s GROUP BY CTP_CAPOGRU,ctp_controp\"' % (filename,edgefile, tablename, where[opt])
        os.system(command)
        output = open(edgefile, 'r')
        lines = output.readlines()
        lines.pop(0)
    
        if opt <> 1:

            command  = 'python querycsv.py -i %s -o %s \"SELECT ctp_controp,CTP_CAPOGRU,sum(importo),location from %s WHERE DATA_CONTABILE = "%s" AND location = \'estero\' AND (NATURA_RAPPORTO = \'IMPIEGHI UNSECURED\' OR NATURA_RAPPORTO = \'IMPIEGHI SECURED\') GROUP BY ctp_controp,CTP_CAPOGRU\"' %(filename,edgefile, tablename, i)
            os.system(command)
            output = open(edgefile, 'r')
            lines2 = output.readlines()
            lines2.pop(0)
            for  j in lines2:          
                lines.append(j)

        lines = [j.replace('"', '') for j in lines]
        
        output = open(edgefile, 'wb')
        output.writelines(lines)
        output.flush()

        
    

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
  #  if opts.domestic and opts.foreign:
  #      sys.stderr.write( "*** choose either the --d or the --f option ***\n\n")
  #      parser.print_help()
  #      sys.exit( 1 )

if __name__ == "__main__":
    main()
    
