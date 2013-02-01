import os
import sys
from optparse import OptionParser
from numpy import intersect1d,minimum,maximum,delete,savetxt,zeros,loadtxt
import interbank_classes as bdi

USAGE = "%prog YEAR"
USE_DESCRIPTION = "YYYYMMDD must be provided as an argument"
parser = OptionParser( USAGE , description= " " + USE_DESCRIPTION ) 
parser.set_defaults( verbosity=None )
parser.add_option( "-d" , "--d",
    action="store_const",
    const= 1, 
    dest="domestic",
    help="Only domestic exposures are included."
)

def main():

    args = parser.parse_args()[1]

    rapporto = [
    'SECURED',
    'UNSECURED']
        
    maturity = [
    'overnight',
    'longterm',
    'nonsignif']
    
    location = [
    'dom',
    'tot_unadj']
    
    dates = open('date.csv', 'r')
    dates = dates.readlines()
    dates.pop(0)
    dates = [j.replace('"','') for j in dates]
    dates = [j.replace('\r\n','') for j in dates]

    for i in args:

        select_dates = list()
        try:
            j = dates.index(i)
            select_dates.append(dates[j])
        except ValueError:
            pass
        if len(select_dates) > 0:
            dates = select_dates
        
    for k in dates:
        for h in location:
            lines = list()
            lines.append('label \t gross volume \t net volume \r\n')
            filename = k + '_' + h.split('_')[0] + '.npy'
            for i in rapporto:
                for j in maturity:
                    G = bdi.Year(filename, rapporto = i, maturity = j)                    
                    gross_vol =  G.Adj.sum()
                    net_vol = gross_vol - sum(G.Adj.diagonal())
                    line = i + '_' + j + '\t' + str(gross_vol) + '\t' + str(net_vol) + '\r\n'
                    lines.append(line)
                    
            output = open(filename.split('.')[0] + '.volumes','wb')
            output.writelines(lines)
            output.flush()
            
    

if __name__ == "__main__":
    main()
    
