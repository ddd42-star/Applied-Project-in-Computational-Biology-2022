# DNA Multiplication
Generalmente per poter sequenziare il genoma di un batterio/virus si deve moltiplicare il genoma.(aggiungere dettagli qui). Pero con le tecniche sviluppate negi anni non e' piu necessario amplificare il DNA. Grazie a cio' basta solo replicare il DNA e poi verra' analizzato. Si vuole allora scoprire in che modo dopo ogni replicazione, il DNA cambia, quali sono le Rimozioni,le Inserzioni e le Sostituzioni.\
Iniziando in una provetta con un determinato ambiente si moltiplica il DNA. Poi dopo un certo tempo t, viene presa una parte di soluziona di una provetta e viene trasferita in un'altra provetta. Il genoma si moltiplica ancora e dopo un certo periodo, viene presa ancora una parte di soluzione e viene trasferita in un'a;tra provetta. IN questo modo, viene amplificato il genoma a diversi tempi.\
# DNA Sequencing tools
Ci sono diverse tecniche di sequencing: Illumina e Nanopore
- Illumina interpreta sequenze corte (300 bp)
- Nanopore interpreta sequenze piu lunghe (300 kbp)
Per questo progetto il genoma e' stato moltiplicato e poi letto con nanopore (citare da chi vengono i dati).
# NANOPORE
![nanoporeImage](https://user-images.githubusercontent.com/80390025/179186649-401cc01d-233e-4672-a170-da7aadf058e1.png)

- Alle estremita' del canale sono presenti due elettrodi, in modo di far passare la corrente al suo interno.
- Si fa passare il DNA/RNA in un canale dei nanopori. Durante il passagio della molecola la corrente ionica cambia (aggiungere dettagli, spiegare come funziona).
- Questi cambiamenti di corrente vengono misurati e producono differenti segnali per i nucleotidi A,C,G,T.
![Color-online-Principle-of-nanopore-DNA-sequencing-technology-a-A-representative](https://user-images.githubusercontent.com/80390025/179193920-84cc087c-9acb-4839-97de-b4390db5dd18.png)
- La macchina raccoglie i differenti signali ed essi vengono trasformati in nucleotidi grazie a Guppy (spiegare come funziona)
- Alla fine si avranno per ogni sequenza di genoma, diverse letture.
# Ricostruzione del genoma
A partire dalle varie letture il genoma viene ricostruito. La ricostruzione viene fatta a partire da un algoritmo sviluppato da gruppo Neher (https://github.com/neherlab/genome-assembly).

# Minimap2
Dopo aver ottenuto il genoma di riferimento si vuole osservare quali errori sono presenti nelle varie reads. Si usa la libreria minimap2 (https://github.com/lh3/minimap2).
- minimap2 -a ref.fa query.fq > alignment.sam
- poi si usa Samtools(referenza) per ordinare le reads che leggono il genoma di riferimento -> samtools -@ {CPU} file.sam (input) -> file.bam (output)(http://www.htslib.org/doc/samtools-sort.html)
- samtools index -@ {CPU} file.bam -> it compress the sorted file to quick access (http://www.htslib.org/doc/samtools-index.html)
- Infine si puo' usare l-interfaccia IGV per avere una visione ordinata delle read ordinate e della sequenza di riferimento.
Ci sono due tipi di letture, quelle chiamate Forwards e quelle chiamate Reverse.
Ci sono due tipi di letture poiche' normalmente la macchina non sa se sta leggendo una sequenza che va da 3->5 o da 5->3. Percio' la macchina interpreta le sequenze Forward quelle seguono per esempio una direzione (3->5 o 5->3). Invece le reads che non hanno una corrispondenza con la refernza in Forward ma hanno una corrispondeza con la reverse strand vengono denominate reverse. Forward hanno i nucletotidi 'A','C','G','T' invece reverse hanno i nucleotdi 'a','c','g','t'.
# Cigar String
Invece di memorizare tutte i nucleotidi che combaciano o  non con la referenza viene usata la string cigar. Essa permette di rispamiare memoria poiche' salva i nucleotidi cosi':
- se c'e' un match con M
- se c'e' una rimozione con D
- se c'e' un inserimento con I
- se c'e' un parziale match con Hard clipping H
- soft clipping
# Pysam
In questo progetto per allineare
