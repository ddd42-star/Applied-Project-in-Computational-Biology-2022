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
In questo progetto per accedere alle informazioni presenti nel file bam usiamo la libreria pysam. Essa permette di accedere a tutte le informazioni necessarie salvate nel file come ad esempio i Nucleotidi delle letture che sono allineate con il genoma di riferimento, la qualita' con cui vengono assengate ad una determinata lettura, il numero di letture che sono allineate ad una determinata posizione del genoma e cosi via.\

## Task 1- Error vs quality
Si vuole conoscere in quale proporzione i match tra letture e genoma di riferimento sono veramente match e quali invece sono errori. Come anticipato, la macchina da' uno score ad un determinato nucleotide. Esso serve a stabilire quanto la macchina e' sicura che il nucleotide indicato e' effetivamente quello letto. Uno score basso che va da 1-10 sottolinea come la macchina sia molto insicura, invece uno score alto sottolinea che la macchina e' molto sicura della previsione. Per poter stabilire la distribuzione tra match e qualita' dobbiamo trovare per ogni lettura il nucleotide query, confrontarlo con il nucleotide di riferimento e raccogliere i numeri degli score lungo tutto il genoma.\
Prima di iniziare pero' abbiamo ricercato la distribuzione delle qualita' dei Nucleotidi in tutte le letture, per capire in che modo i Nucleoti A,C,G,T venivano predetti. Nel file bam le qualita' dei Nucleotidi e' salvata con i segni char ASCII.
"In FASTQ files, quality scores are encoded into a compact form, which uses only 1 byte per quality value. In this encoding, the quality score is represented as the character with an ASCII code equal to its value + 33. The following table demonstrates the relationship between the encoding character, its ASCII code, and the quality score represented." (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm). Percio' per traformare i caratteri trovati in cifre bisogna sottrarre ad ogni simbolo il valore di 33.\
Ad esempio: ord('!') - 33 = 0\
Sapendo come decodificare i byte trovati, abbiamo guardato in ogni lettura il numero di qualita' per i nucleotidi A,C,G e T. Per ogni nucleotide abbiamo raccolto tutte le cifre delle qualita' lette e le abbiamo visualizzate grazie ad un istogramma.
![Quality_distribution_All_reads](https://user-images.githubusercontent.com/80390025/179225174-3b1fa30c-0a02-4bc2-96ca-3672bb8e2eaf.png)
Come si puo vedere la distribuzione ha un andamento binomiale. C'e' un primo picco con un valore di qualita' pari a 10 ed un picco con una valore di qualita' pari a 30. Il secondo picco e' piu' alto percio significa che avremo molti nucleotidi indovinati correttamente (qualita' della mappa alta) ed invece avremo meno nucleotidi che saranno indovinati correttamente a qualita' piu' bassa.\
Secondariamente vogliamo capire in che modo di quanti dati di lettura disponiamo. Avere tanti dati e' molto utile in quanto ti permette di avere risultati piu' accurati ed effettuare un'analisi piu' approfondita. Per trovare questo pero' avevamo bisogno di uno strumento che contasse per ogni posizione sul genoma le letture presenti. Abbiamo percio' fatto affidamento a *pileup*, una funzione di pysam che fa l'operazione di mettere insieme tutte le reads che mappano in una data posizione, e ritorna un iteratore su pileup columns. Grazie a pileup possiamo ricavare ad ogni posizione sul genoma:(https://pysam.readthedocs.io/en/latest/api.html?highlight=Column#pysam.PileupColumn.get_query_names)
- get_mapping_qualities
- get_num_aligned
- get_query_names
- get_query_positions
- get_query_qualities
- get_query_sequences
Oltre a cio' possiamo conoscere molte altre cose, tra cui il numero di letture ad ogni posizione nel genoma. utilizzare pileup e' molto utile poiche' si tratta di un IteratorColoumn. Un iterator e' un oggetto che permette di salvare oggetti molto grossi utilizzando pochissima memoria. infatti invece di salvare l'oggetto per intero viene salvato solamente la referenza di una posizione con quella successiva. In questo modo eseguendo un iterazione possiamo avere accesso a tutti i dati contenuti nel file per tutte le letture nello stesso tempo.\
Iterando lungo tutte le posizioni, salviamo tutte le quantita' di letture presenti ad ogni poszione e otteniamo il sequente grafico:
![Numer_of_reads](https://user-images.githubusercontent.com/80390025/179231010-023526f9-292f-4664-a363-0b709778f687.png)
L'andamento di questo grafico e' possibile descriverlo con una 'u'. Il genoma che stiamo analizzando appartiene ad un batterio. i batteri presentano un genoma chiuso e rotondo. Dato che la replicazione inizia ad un estreimita' nelle due direzioni e si conclude all'estremita' opposto osserveremo come ben visualizzato nel grafico che all'inizio abbiamo molte letture invece verso la fine, ossia il centro nel grafico, le letture saranno inferiori.(inserire disegno dna batterio)\
Per ultimo abbiamo controllato lungo tutto il genoma il valore medio delle qualita'. Abbiamo iterato lungo tutto il genoma sempre utilizzando pileup e abbiamo ottenuto il seguente grafico:
![Average_quality](https://user-images.githubusercontent.com/80390025/179232196-30ad0612-5228-4950-b91b-35f0af972449.png)
Il grafico soprea raffigurato mostra che lungo tutto il genoma il valore medio delle qualita' assegnate in ogni posizione varia da 18 a 22 circa.

Per poter trovare la distribuzione degli errori per ogni qualita' abbiamo costruito una matrice che conteneva\

[(N_i di riferimento, N_j della query) : {valore della qualita':numero di volte che si osserva la qualita'}]

N_i rappresenta il **nucleotide di riferimento** sul genoma invece N_j rappresenta **nucleotide che si osserva sulle reads**. Per creare questa matrice abbiamo utilizzato i *defaultdict* (https://docs.python.org/3/library/collections.html#collections.defaultdict). Questo tipo di oggetto fa parte della classe dei containers e ha la caratteristica di cacolare il numero di volte che una chiave viene salvata. Per noi e' ideale poiche' dobbiamo salvare per ogni coppia nucleotide di riferimento e nucleotide query il numero di qualita' osservate lungo tutto il genoma. Alla fine quindi avremo una matrice composta da 16 coppie di nucleotidi con i rispettivi dizionari contentneti le qualita' e il numero di volte che una qualita' e' presente. Questi sono i grafici che abbiamo ottenuto:
![plot_A_reference_log_quality](https://user-images.githubusercontent.com/80390025/179235788-cdd613f5-41b2-4797-a76a-ac13d6a46a5f.png)
![plot_A_reference_normed_quality](https://user-images.githubusercontent.com/80390025/179235793-4ccf86a1-0d53-47b9-9543-f23ad7a92005.png)
![plot_C_reference_log_quality](https://user-images.githubusercontent.com/80390025/179235795-832afabf-b2e0-4514-93d7-367fa803e38f.png)
![plot_C_reference_normed_quality](https://user-images.githubusercontent.com/80390025/179235796-a08c84c1-a8e9-4c54-816d-504a8b3c9915.png)
![plot_G_reference_log_quality](https://user-images.githubusercontent.com/80390025/179235798-adbc22d9-543a-460d-b343-050b477af5a0.png)
![plot_G_reference_normed_quality](https://user-images.githubusercontent.com/80390025/179235800-afbbd8e4-b5d1-4ee9-8a5f-a46e2430f0a9.png)
![plot_T_reference_log_quality](https://user-images.githubusercontent.com/80390025/179235803-c45b9f83-dff4-4ee4-a4c5-b7aecb9bbaab.png)
![plot_T_reference_normed_quality](https://user-images.githubusercontent.com/80390025/179235804-38db5ad3-367d-4a11-a5a2-8747b38dac43.png)

Per ogni set di Nucleotide di riferimento troviamo tutte le coppie con l'asse logaritmico e l'asse normalizzato. Questi grafici ci confermano che quando la macchina indica che un nucleotide fa match con la referenza si tratta il piu' delle volte del giusto nucleotide. Invece a qualita' basse e' molto probabile di ottenere un errore ma in proporzione questo errore e' relativamente basso, come indicato in questo grafico.
![plot_A_reference_quality](https://user-images.githubusercontent.com/80390025/179237101-80791128-6578-4a7b-98bd-39d420a7930f.png)
![plot_C_reference_quality](https://user-images.githubusercontent.com/80390025/179237196-374cab05-8564-4a0f-8c1a-404e9af7b459.png)
![plot_G_reference_quality](https://user-images.githubusercontent.com/80390025/179237229-7fd92372-75fa-4ff7-93d7-d13e31e9b1bd.png)
![plot_T_reference_quality](https://user-images.githubusercontent.com/80390025/179237261-cb2bbd49-2a92-40d3-9475-3d739fba4ec4.png)

Interessente da osservare e' il fatto che nei casi AG e GA ci sono piu' errori, come confermato anche dagli autori della tecnologia nanopore in questo articolo ([journal.pone.0257521.pdf](https://github.com/ddd42-star/Applied-Project-in-Computational-Biology-2022/files/9121381/journal.pone.0257521.pdf))

Infine per poter utilizzare la matrice creata e non doverla rifare ogni volta l'abbiamo salvato con la libreria pickle, che permette di trasformare un oggetto in python in una stringa byte e all'occorrenza utilizzarla.(https://docs.python.org/3/library/pickle.html)

## Matrice delle probabilita' per tre diversi intervalli
utilizzando la matrice calcolata in precedenza siamo interessati alle probabilita' dei singoli errori per ogni nucleotide letto versus i nucletotidi veri.
La matrice ha questa forma:\
|   | A | C | G | T |
|---|---|---|---|---|
| A | ? | ? | ? | ? |
| C | ? | ? | ? | ? |
| G | ? | ? | ? | ? |
| T | ? | ? | ? | ? |

Sulle diagonali troviamo le percentuali degli errori che sono letti correttamente, poi si trovano gli errori presenti nelle letture. Questa matrice l'abbiamo fatto per tre diversi intervalli:
- per la qualita' minore di 20
- per la qualita' maggiore di 20
- per la qualita' totale
![Matrix_error](https://user-images.githubusercontent.com/80390025/179247479-60dd72e8-7041-44c2-8c81-6450e9946b0c.png)
Le matrici mostrano che gli errori maggiori si trovano con le coppie in cui sono presenti AG o GA e infatti hanno un errore che e' quasi il doppio rispetto a tutti gli altri.

## Frequency of error at each level of quality
Poi ci siamo domandati quale fosse la frequenza di errore per ogni livello di qualita' trovato. Avevamo ipotizzato di avere a qualita' basse un livello di errore elevato e invece a qualita' basse un livello di errore piu' piccolo.

![Frequency of error Forward read](https://user-images.githubusercontent.com/80390025/179258931-de684e07-89bf-46c8-ac4c-c5e204fd24e7.png)
Per le letture Forward troviamo che come ipotizzato gli errori piu' frequenti si trovano a qualita' piu' basse. La linea grigia rappresenta l'errore lineare (da chiarire!) e come possiamo vedere tutte le frequenze sono al di sotto, anche per quelle piu' basse. Cio' significa che si' le frequenze misurate hanno errori alti pero' sono al di sotto dell'errore medio. A qualita' piu' alte invece l'errore e' molto piu' basse, ogni 100 nucleotidi sono 1 e' errato.
![Frequency of error reverse read](https://user-images.githubusercontent.com/80390025/179260326-041695b0-826c-434f-9e8e-4142a7798cab.png)
Similmente anche per le letture reverse troviamo lo stesso comportamento, percio' siamo sicuri poiche' e' stato confermato che a qualita' piu' basse e' piu' frequente che la macchina sbagli invece a qualita' piu' alte la frequenza e' molto minore.

## Error with sequences
Dopo aver analizzato coppie di nucleotidi ci siamo domandati in che modo sequence corte di tre o cinque nucleotidi si comportano. Otteniamo alcune sequenze piu' facili da sbagliare oppure otteniamo sequenze che non presentano nessun vantaggio rispetto ad altre. Per rsilvere questa domanda abbiamo fatto una nuova matrice, in cui al posto del nucleotide di referenza, abbiamo selezionato tutti i 64 possibili tripletti. Poi abbiamo calcolcato il numero di errore per ogni tripletto, sommando le qualita' e dividendo per il titale numero di possibilita'. Dopo aver ottenuto questa percentuale, abbiamo diviso ogni tripletto per la reciproca base che corrispondeva alla referenza:\
P_e('triplet')/P(single nucleotide)\
Abbiamo calcolato il logaritmo di questa divisione in modo da ottenere la distribuzione degli errori.
![logRatioTriplet](https://user-images.githubusercontent.com/80390025/179264109-46334017-f254-46ee-b582-6c235ce9117f.png)
Dove troviamo zero ci troviamo nella zona dove o possiamo fare giusto o dove possiamo sbagliare. Invece nella coda a sinistra siamo sicuri di sbagliare molto poco e nella zona destra siamo sicuri di sbagliare molto. Non ci sono pero' sequenze che sbagliano con molta maggiore evidenza rispetto ad altre. Abbiamo allora calcolato la matrice con sequenza di riferimento con cinque nucleotidi e poi abbiamo anche in questo caso calcolato il rapporto logaritmico con la percentuale di errore del singolo nucleotide.\
P_e('quintuple')/P(single nucleotide)\
![logRatioQuintuple](https://user-images.githubusercontent.com/80390025/179264933-c5184deb-b84f-41bb-bf23-eb22a0cdea7f.png)
In questo caso abbiamo moltee piu' sequenze e quindi la distribuzione e' molto piu' facile da osservare. I casi importanti sono alle due estremita'. All'estremita' vicino a -1 troviamo la sequenza GGATC e questa e' quasi sempre riconosciuta correttamente. Invece all'estremita' di destra troviamo due sequenze nei pressi di 0.75, ossia la sequenza CCTGG e la sequenza TCGGG. Queste due sequenze quasi sicuramente verrano erroneamente identificate dalla macchina come ci viene anche confermato dall'articole degli autori di questa tecnologia ([NanoporeArticle.pdf](https://github.com/ddd42-star/Applied-Project-in-Computational-Biology-2022/files/9122424/NanoporeArticle.pdf)).

## Contesto  nelle sequenze
l'ultima domanda che ci siamo posti in questa prima parte riguardava i contesti. In che modo la presenza di uno o piu' nucleotidi prima o dopo influisce sulla corretta identificazione del nucleotide di riferimento nella query? Per rispondere a questa domanda abbiamo utilizzato un concetto fisico conosciuto nominato come entropia. (fare disegno dei due sistemi)\
Questa idea nasce dalla inferenza nella teoria dell'informazione e vuole valutare in che modo la presenza o meno dei nucleotidi inferisca sul riconoscimento di altri nucleotidi da parte dell atecnologia nanopore.\
(Scrivere derivazione formule trovate)\
La prima domanda che abbiamo risposto riguardava l'inferenza di un nucleotide i + l rispetto ad un nucleotide a posizione i. Abbiamo ipotizzato che il primo o il secondo nucleotide avesse un effetto sul nucleotide di riferimento ma poi piu' si andava distanti piu' l'influenza scemeva via via.
![Entropy_of_signal](https://user-images.githubusercontent.com/80390025/179267409-eb7c653e-d4c0-4dd0-b18f-93bb554ad8c7.png)
E quello che abbiamo trovato, rappresentato in questo grafico conferma la nostra ipotesi. Come possiamo vedere il nucleotide a distanza +1 influisce sull'informazione del nucleotide di riferimento. Invece il nucleotide a distanza -1 influisce in modo minore pero' e' anche esso importante. Successivamente piu' la distanza cresce piu' l'influenza dei nucleotidi dopo diminuisce fino a saturarsi dopo una distanza di 5. E questo risultato e' anche confermato dal fatto che al massimo in un nanoporo ci stanno 5 nucleotidi percio' l'influenza oltre questa distanza e' pari a nulla.

La seconda domanda che abbiamo risposto riguarda in che modo l'informazione del contesto influenza delle sequenze di nucleotidi? Ora non siamo piu' interessati al singolo nucleotide bensi' ad una sequenza di 3,5 o 7 nucleotidi. Per questa domanda crediamo che diminuisca con l'aumento della sequenza, ossia con sequenze piu' basse sia piu' piu' grande ma poi diminuisca fino a raggiungere livelli vicini allo 0.
![Information_of_the error](https://user-images.githubusercontent.com/80390025/179268985-f2f4071e-0afb-461a-9df9-cd2cf1854ecf.png)
Il risultato ottenuto non ha soddisfatto le nostre aspettative. Infatti la curva dell'entropia non diminuisce in modo esponenziale. Fino a che siamo di fronte a sequenze di lunghezza di 3,5 o 7 siamo sicuri di ottenere risultati giusti, pero' dopo siamo di fronte a sequenze sempre piu' grandi che raggiungono la sequenza totale del genoma. Poi si ottiene un ulteriore diminuzione che e' riconducibile all'overfitting, ovvero quando non si hanno abbasntanza letture e percio' i risultati ottenuti non sono piu' attendibili.\
Similmente otteniamo lo stesso risultato se partiamo da un nucleotide e guardiamo l'influenza su sequenze successive (Forward) o precedenti (Reverse). Fino a sequenze di lunghezza 3 o 5 ci possiamo fidare dei risultati e poi le sequenze diventano troppo grandi. Il pericolo e' che magari certe sequenze non le vediamo mai e a loro viene assegnata probabilita' di sbagliare di 0 quando in verita' non lo sono.
(chiedere dettagli per spiegare meglio questa parte)

la terza domanda che abbiamo analizzato riguarda calcolare l;entropia dell'errore a partire da sequenze che partono da un nucleotide e poi vanno in una direzione (verso destra o verso sinistra) (inserire disegno --------------------> oppure <-------------------)
Inserire grafici\
In questo caso si vede meglio la curva esponenziale che avevamo ipotizzato di ottenere. Abbiamo distinto tra Reverse e Forward e abbiamo notato come le due curve siano differenti.E' piu' informativo il nucleotide dopo il nucleotide di riferimento rispetto a quello prima, come confermato dal grafico ottenuto nella prima parte di questo ultimo esercizio.







