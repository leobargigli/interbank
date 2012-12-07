# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:31:52 2012

@author: leonardo
"""

from numpy import *
import os
import sys
from optparse import OptionParser


USAGE = "%prog FILENAME"
USE_DESCRIPTION = "The FILENAME file must be a csv file."
parser = OptionParser(USAGE , description= " " + USE_DESCRIPTION) 
#parser.set_defaults( verbosity = None )


def main():
    
    args  = parser.parse_args()
    filename = args[1][0]
    datefile = 'date.csv'
    tablename = os.path.splitext(filename)[0]
    
    # here we select the dates in the dataset
    command  = 'python querycsv.py -i %s -o %s \"SELECT DATA_CONTABILE from %s GROUP BY DATA_CONTABILE\"' %(filename,datefile, tablename)
    os.system(command)
    dates = open(datefile, 'r')
    dates = dates.readlines()
    dates.pop(0)
    dates = [j.replace('"', '') for j in dates]
    dates = [j.replace('\r\n', '') for j in dates]

    ######

    for i in dates:
        
        # here we find the reporters
        
        reporterfile = i + '_reporters.csv'    
        command  = 'python querycsv.py -i %s -o %s \"SELECT CTP_CAPOGRU from %s WHERE DATA_CONTABILE = \'%s\' GROUP BY CTP_CAPOGRU\"' % (filename , reporterfile , tablename, i)
        os.system(command)
        output = open(reporterfile, 'r')
        reporters = output.readlines()
        reporters.pop(0)

        ####        
        
        # here we extract the counterparts of foreign credit relations (i.e. foreign debtors)
        
        debtfile = i + '_foreign.debt'
        
        command  = 'python querycsv.py -i %s -o %s \"SELECT ctp_controp from %s WHERE DATA_CONTABILE = "%s" AND location = \'estero\' AND NATURA_RAPPORTO = \'IMPIEGHI UNSECURED\' OR NATURA_RAPPORTO = \'IMPIEGHI SECURED\' GROUP BY ctp_controp\"' %(filename,debtfile, tablename, i)
        
        os.system(command)
        output = open(debtfile, 'r')
        foreign_debtors = output.readlines()
        foreign_debtors.pop(0)

        ####

        # here we find non reporting debtors 
        
        non_reporters = setdiff1d(foreign_debtors,reporters) # non reporting nodes are a subset of foreign counterparts
        savetxt(debtfile,non_reporters,fmt = '%30s')        
        
        
if __name__ == "__main__":
    main()

    