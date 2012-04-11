import os
import sys
from optparse import OptionParser
from numpy import loadtxt

USAGE = "%prog [-d(omestic)] FILENAME"
USE_DESCRIPTION = "The FILENAME file must be a csv file."
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity=None )

parser.add_option( "-d" , "--d",
    action = "store_const",
    const = 1, 
    dest ="domestic",
    help = "Only domestic exposures are included."
)

parser.add_option( "-f" , "--f",
    action = "store_const",
    const = 1, 
    dest ="foreign",
    help = "Only domestic exposures are included."
)


def main():

    ( opts , args ) = parser.parse_args()
    __validate_opts_and_args( opts , args )
    filename = args[0]
    optdict = {1: 'dom', 2: 'for' , None: 'tot'}
    tablename = os.path.splitext(filename)[0]
    datefile = 'date.csv'
    opt = opts.domestic + opts.foreign
    
    command  = 'python querycsv.py -i %s -o %s \'SELECT DATA_CONTABILE from %s GROUP BY DATA_CONTABILE\'' %(filename,datefile, tablename)
    os.system(command)
    dates = open(datefile, 'r')
    dates = dates.readlines()
    dates.pop(0)
    for i in xrange(len(dates)):
        dates[i] = dates[i].replace('"', '')
        dates[i] = dates[i].replace('\r\n', '')

    for i in dates:
        
        reporterfile = i + '_reporters.csv'
        command  = 'python querycsv.py -i %s -o %s \'SELECT CTP_CAPOGRU from %s WHERE DATA_CONTABILE = "%s" GROUP BY CTP_CAPOGRU\'' % (filename , reporterfile , tablename, i)
        os.system(command)
        output = open(reporterfile, 'r')
        lines = output.readlines()
        lines.pop(0)
        for j in xrange(len(lines)):
            lines[j] = lines[j].replace('"', '')
        output = open(reporterfile, 'wb')
        output.writelines(lines)
        output.flush()
        
        edgefile = i + '_'+ optdict [opts.domestic] + '.edgelist'
        output = open(edgefile, 'wb')
        where = {1: 'WHERE DATA_CONTABILE = "%s" AND location = "domest"'% (i), 
                 None: 'WHERE DATA_CONTABILE = "%s" AND NATURA_RAPPORTO LIKE "%%DEBITI%%"'% (i)} 
        command  = 'python querycsv.py -i %s -o %s \'SELECT CTP_CAPOGRU,ctp_controp,sum(importo) from %s %s GROUP BY CTP_CAPOGRU,ctp_controp\'' %(filename,edgefile, tablename, where[opts.domestic])
        os.system(command)
        output = open(edgefile, 'r')
        lines = output.readlines()
        lines.pop(0)
    
        if opts.domestic is None:

            command  = 'python querycsv.py -i %s -o %s \'SELECT ctp_controp,CTP_CAPOGRU,sum(importo) from %s WHERE DATA_CONTABILE = "%s" AND NATURA_RAPPORTO LIKE "%%IMPIEGHI%%" GROUP BY ctp_controp,CTP_CAPOGRU\'' %(filename,edgefile, tablename, i)
            os.system(command)
            output = open(edgefile, 'r')
            lines2 = output.readlines()
            lines2.pop(0)
            for  j in lines2:          
                lines.append(j)

        for j in xrange(len(lines)):
            lines[j] = lines[j].replace('"', '')
        
        output = open(edgefile, 'wb')
        output.writelines(lines)
        output.flush()
        
#        try:
#                edgelist = loadtxt(edgefile, delimiter =',')
#        except IOError:
#                volume = 0.
#        try:
#            volume = edgelist[:, 2].sum()
#        except IndexError:
#            if len(edgelist) > 0:
#                volume = edgelist[-1]
#            else:
#                volume = 0.
#        print volume
        
    

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
    if opts.domestic and opts.foreign:
        sys.stderr.write( "choose either the --d or the --f option".format(args[0]) )
        parter.print_help()

if __name__ == "__main__":
    main()
    
