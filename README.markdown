# Istruzioni

1. Copiare i dati da analizzare nella cartella in cui sono stati salvati i file .py. I file di dati devono essere in formato csv.
 
    **N.B.** I nomi dei campi devono essere gli stessi del sample gia' inviato.

2. Convertire il file .csv in una edgelist.
    
    Da terminale digitare > `python fromcsv.py NOMEFILE.csv`
    
    **N.B.**:    
    * Usare l'opzione -d o --d per includere solo le transazioni contrassegnate come "domestiche".    
    * La procedura crea tanti file .edgelist quanti sono i diversi valori nel campo DATA_CONTABILE.
    Se NON viene utilizzata l'opzione -d, la procedura restituisce un file denominato DATACONT_tot.edgelist
    In caso contrario restituisce un file denominato DATACONT_dom.edgelist
    
3. Calcolare le statistiche per la rete.

Da terminale digitare:

 `python interbank_script.py DATACONT_tot(dom).edgelist`

## Dipendenze

Oltre a python, si richiede l'installazione dei seguenti moduli: 

* numpy & scipy http://numpy.scipy.org/
* matplotlib http://matplotlib.sourceforge.net/
* networkx http://networkx.lanl.gov/



