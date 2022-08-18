import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import Counter, defaultdict
from Bio import SeqIO
import re







# Find length of the insertion in all the genome

def find_length_insertion(file):
    
    bamfile = pysam.AlignmentFile(file,'rb')
    
    insertion_for = defaultdict(list)
    insertion_rev = defaultdict(list)

    for col in bamfile.fetch():
        
        if col.is_secondary or col.is_supplementary:
            continue
        
        insertions = insertion_for if col.is_forward else insertion_rev
        for tipo,lunghezza in col.cigartuples:
            insertions[tipo].append(lunghezza)
    return insertion_for,insertion_rev





# Find Insertion quality
def quality_insertion(file):
    
    bamfile = pysam.AlignmentFile(file,'rb')
    #qualita = defaultdict(list)
    insertion_for = defaultdict(list)
    insertion_rev = defaultdict(list)
    #qual = {}
    for col in bamfile.fetch():
        n = 0
        if col.is_secondary or col.is_supplementary:
            continue
        insertions_qual = insertion_for if col.is_forward else insertion_rev
        #if col.is_forward
        for tipo, lunghezza in col.cigartuples:
            insertions_qual[lunghezza].append(np.mean(col.query_qualities[n:n+lunghezza]))    
    return insertion_for,insertion_rev


forward_qual, reverse_qual = quality_insertion('reads.sorted.bam')

# Count of a certain quality insertion
def intervall (dizionario,min,max):

    x = dizionario.keys()
    w = []
    for l in x:
        qs = np.array(dizionario[l])
        qmin , qmax = min,max
        w.append(np.sum((qs <= qmax) & (qs >= qmin)))
    
    return x,w



# Search in all the genome for all the homopolymer

def find_homopolymer_nucleotide(filesequence):
    # Read true sequence
    true_sequence = []
    for seq in SeqIO.read(filesequence, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)


    matrice = defaultdict(list)
    altra = defaultdict(list)
    n = true_string[0]
    L = 1
    p = 0
    for x in range(1,len(true_string)-1):

        if p == len(true_string) - 1:
            continue
        nucleotide = true_string[x]

        if n == nucleotide:
            L += 1
            continue
        else:
            #save old status
            if p != 0:
                if L > 1:
                    for space in range(p,p+L):
                        matrice[space].append((n,L))
                    #matrice[(n,L)].append((p,p+L-1))
                    #matrice[(p,p+L-1)].append((n,L))
                    altra[(L*n,L)].append((p,p+L-1))
                else:
                    #matrice[(n,L)].append((p,p))
                    matrice[p].append((n,L))
                    altra[(n,L)].append((p,p))
            #reset status

            n = nucleotide
            L = 1
            p = x

    return matrice,altra
# Function 2.0
def find_homopolymer_nucleotide_new(filesequence):
    # Read true sequence
    true_sequence = []
    for seq in SeqIO.read(filesequence, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)


    #matrice = defaultdict(list)
    matrice  = {}
    altra = defaultdict(list)
    n = true_string[0]
    L = 1
    p = 0
    for x in range(1,len(true_string)-1):

        if p == len(true_string) - 1:
            continue
        nucleotide = true_string[x]

        if n == nucleotide:
            L += 1
            continue
        else:
            #save old status
            if p != 0:
                if L > 1:
                    space = np.arange(p,p+L)
                    matrice[space[0]] = (n,L)#.append((n,L))
                    matrice[space[-1]] = (n,L)#.append((n,L))
                    #for space in range(p,p+L):
                    #    matrice[space].append((n,L))
                    #matrice[(n,L)].append((p,p+L-1))
                    #matrice[(p,p+L-1)].append((n,L))
                    altra[(n,L)].append((p,p+L-1))
                else:
                    #matrice[(n,L)].append((p,p))
                    matrice[p] = (n,L)#.append((n,L))
                    altra[(n,L)].append((p,p))
            #reset status

            n = nucleotide
            L = 1
            p = x
    #matrix = {}
    #for x in matrice.keys():
    #    matrix[x] = matrice[x]

    return matrice,altra


# Search in all the genome the position i and position i+1 where an insertion happen
def insertion_position_and_sequence(file):

    # Reading file
    bamfile = pysam.AlignmentFile(file,'rb')
    # List for Forward and Reverse read    
    insertion_for = []
    insertion_rev = []
    tot_for = np.zeros(bamfile.get_reference_length(bamfile.references[0]))
    tot_rev = np.zeros(bamfile.get_reference_length(bamfile.references[0]))

    # Iterate over all the reads
    for read in bamfile.fetch():
        
        if read.is_secondary or read.is_supplementary:
            continue
        
        insertions = insertion_for if read.is_forward else insertion_rev
        total = tot_for if read.is_forward else tot_rev
        # Finding the start position of the reference and the start position of the query
        ref_beg = read.pos
        position_reference = read.pos
        position_query = 0 #read.query_alignment_start
        # read.query_alignment_start give the first position that is not a soft clip.
        # If there is a soft clip at first, we have to skip the it
        #if read.cigartuples[0][0] == 4:
        #    cigarstring = read.cigartuples[1:]
        #else:
        #    cigarstring = read.cigartuples
        #print(read.cigartuples[1:])
        for tipo,lunghezza in read.cigartuples:
            if tipo == 0: # match
                position_reference += lunghezza
                position_query += lunghezza
            elif tipo == 1: # Insertion
                ins_pos = (position_reference-1,position_reference)
                ins_sequence = read.query_sequence[position_query:position_query+lunghezza]
                insertions.append([ins_pos, ins_sequence])
                position_reference += 0
                position_query += lunghezza
            elif tipo == 2: # Deletion
                position_reference += lunghezza
                position_query += 0
            elif tipo == 4: # Soft clip
                position_reference += 0
                position_query += lunghezza
            else: # Hard clip
                position_reference += 0
                position_query += 0

        ref_end = position_reference
        total[ref_beg+1:ref_end-1] += 1

    return insertion_for, insertion_rev,tot_for,tot_rev


# Count the length of the inserted homopolymer starting from the nucleotide i
def search_homoploymer_forward_start(forward, homopolymer_sequence):

    dizionario = defaultdict(list)
    for x in forward:
        # ex. x = [(18,19), 'TT']
        start_pos = x[0][0] # position nucleotides i, in ex. is 18
        nucleotides = homopolymer_sequence[start_pos][0][0] # nucleotides that there is in start_pos, in ex. {18: [('C,1)]}
        inserted_sequence = x[1] # sequence that has been inserted, in ex. 'TT'
        if nucleotides == inserted_sequence[0]: # if the nucleotides in position i is the same of the first of the insertion
                                                # we find the maximum length of the homopolymer inserted using regex
            object_to_search = '^[' + nucleotides + ']+'
            seq_match = re.search(object_to_search,inserted_sequence).group()# sequence of the max homopolymer
            dizionario[homopolymer_sequence[start_pos][0]].append(len(seq_match))# we save the length of the homopolymer
        else:
            continue

    return dizionario
# Count the length of the inserted homopolymer starting from the nucleotide i+1
def search_homoploymer_forward_end(forward, homopolymer_sequence):

    dizionario = defaultdict(list)
    for y in forward:
        end_pos = y[0][1]
        nucleotides_end = homopolymer_sequence[end_pos][0][0]
        inserted_sequence_end = y[1]
        if nucleotides_end == inserted_sequence_end[-1]:

            object_to_search_at_the_end = '[' + nucleotides_end + ']+$'
            seq_end_match = re.search(object_to_search_at_the_end, inserted_sequence_end).group()
            dizionario[homopolymer_sequence[end_pos][0]].append(len(seq_end_match))
        else:
            continue
    return dizionario

# search insertion 2.0
def search_one_after_another(forward, homopolymer):

    #Search for homopolymer for i
    dizionario = defaultdict(list)
    #for x in forward:
        # ex. x = [(18,19), 'TT']
        #start_pos = x[0][0]# position nucleotides i, in ex. is 18
        #if start_pos in homopolymer:

        #    nucleotides = homopolymer[start_pos][0][0] # nucleotides that there is in start_pos, in ex. {18: [('C,1)]}
        #    inserted_sequence = x[1] # sequence that has been inserted, in ex. 'TT'
        #    if nucleotides == inserted_sequence[0]: # if the nucleotides in position i is the same of the first of the insertion
                                                    # we find the maximum length of the homopolymer inserted using regex
        #        object_to_search = '^[' + nucleotides + ']+'
        #        seq_match = re.search(object_to_search,inserted_sequence).group()# sequence of the max homopolymer
        #        dizionario[homopolymer[start_pos]].append(len(seq_match))# we save the length of the homopolymer
        #    else:
        #        continue
        #else:
        #    continue
    #Search for homopolymer for i+1        
    for y in forward:
        end_pos = y[0][1]
        if end_pos in  homopolymer:

            nucleotides_end = homopolymer[end_pos][0][0]
            inserted_sequence_end = y[1]
            if nucleotides_end == inserted_sequence_end[-1]:

                object_to_search_at_the_end = '[' + nucleotides_end + ']+$'
                seq_end_match = re.search(object_to_search_at_the_end, inserted_sequence_end).group()
                dizionario[homopolymer[end_pos]].append(len(seq_end_match))
            else:
                continue

        else:
            continue

    return dizionario

# Find length of the homopolymer 2.0
def find_length_homopolymer_new(forward_seq,reverse_seq,homopolymer):


    dizionario_forward= search_one_after_another(forward_seq,homopolymer)

    dizionario_reverse = search_one_after_another(reverse_seq,homopolymer)

    return dizionario_forward,dizionario_reverse

# Count the legnth of the homopolymer simultaneosly for forward and reverse reads
def find_length_homopolymer(forward_seq,reverse_seq,homopolymer):


    dizionario_lunghezze_start = search_homoploymer_forward_start(forward_seq,homopolymer)

    reverse_lunghezze_start = search_homoploymer_forward_start(reverse_seq, homopolymer)

    dizionario_lunghezze_end = search_homoploymer_forward_end(forward_seq,homopolymer)

    reverse_lunghezze_end = search_homoploymer_forward_end(reverse_seq,homopolymer)

    return dizionario_lunghezze_start, dizionario_lunghezze_end, reverse_lunghezze_start, reverse_lunghezze_end



# Function for doing plots

#INSERTION
def plot_insertion_length(dict_forward, dict_reverse):

    # M = 0
    # I = 1
    # D = 2
    # N = 3 (ref_skip)
    # S = 4 (soft clip)
    # H = 5 (hard clip)

    plt.figure(figsize=(10,8))
    plt.hist(dict_reverse[1], bins=np.arange(800),histtype='step',label='Insertion Reverse')
    plt.hist(dict_forward[1], bins=np.arange(800),histtype='step',label='Insertion Forward')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel('length of the insertion sequence')
    plt.legend()
    return plt.show()

# DELETION
def plot_deletion_length(dict_forward, dict_reverse):

    # M = 0
    # I = 1
    # D = 2
    # N = 3 (ref_skip)
    # S = 4 (soft clip)
    # H = 5 (hard clip)

    plt.figure(figsize=(10,8))
    plt.hist(dict_reverse[2], bins=np.arange(800),histtype='step',label='Deletion Reverse')
    plt.hist(dict_forward[2], bins=np.arange(800),histtype='step',label='Deletion Forward')
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel('length of the deletion sequence')
    plt.legend()
    return plt.show()

# HARD AND SOFT CLIP
def plot_Soft_and_Hard_clip_length(dict_forward, dict_reverse):

    # Plot Hard Clip
    # M = 0
    # I = 1
    # D = 2
    # N = 3 (ref_skip)
    # S = 4 (soft clip)
    # H = 5 (hard clip)

    kwargs = {
        "bins" : np.logspace(0,5,100),
        "histtype" : "step",
    #     "density" : True,
    }

    plt.figure(figsize=(10,8))
    plt.hist(dict_reverse[5],label='Hard Clip Reverse', **kwargs)
    plt.hist(dict_forward[5],label='Hard Clip Forward', **kwargs)
    plt.hist(dict_reverse[4],label='Soft Clip Reverse', **kwargs)
    plt.hist(dict_forward[4],label='Soft Clip Forward', **kwargs)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel('length of the Hard and Soft clip sequence')
    plt.legend()
    return plt.show()



# PLOT OF THE QUALITY INTERVALL
def plot_of_the_quality_of_insertion(quality_forward, quality_reverse):

    plt.figure(figsize=(10,8))
    bin = np.logspace(0,4,100)
    bad_x, bad_w = intervall(quality_forward,0,10)
    rbad_x,rbad_w = intervall(quality_reverse,0,10)
    plt.hist(bad_x,weights=bad_w, bins=bin, histtype='step', label='bad quality forward insertion')
    plt.hist(rbad_x,weights=rbad_w,bins=bin,histtype='step',label='bad quality reverse insertion')
    good_x, good_w = intervall(quality_forward,30,40)
    rgood_x,rgood_w = intervall(quality_reverse,30,40)
    plt.hist(rgood_x,weights=rgood_w,bins=bin,histtype='step', label='good quality reverse insert')
    plt.hist(good_x,weights=good_w, bins=bin, histtype='step', label='good quality forward insert')
    plt.xlabel('length of insertion')
    plt.ylabel('mean quality for every length')
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    return plt.show()


# PLOT LENGHT INSERTION IN HOMOPOLYMER
def plot_homopolymer_length(nucleotide, dizionario):

    plt.figure(figsize=(10,8))
    for x in range(1,9):
        if dizionario[(nucleotide,x)]:
            plt.hist(dizionario[(nucleotide, x)],density=True,bins = np.arange(100),label='(' + nucleotide + ', '+ str(x) +')', histtype='step')
        else:
            continue
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Length of the homopolymer inserted')
    plt.ylabel('Count of the homopolyer')
    plt.legend()
    return plt.show()


# PLOT COMPARISON START AND END
def plot_comparison_between_start_end(dict_start, dict_end, nucleotide, lunghezza):

    name = '(' + nucleotide + ',' + str(lunghezza) + ') '
    plt.figure(figsize=(10,8))
    plt.hist(dict_start[(nucleotide, lunghezza)],density=True, bins=np.arange(100), label= name +'position i', histtype='step')
    plt.hist(dict_end[(nucleotide,lunghezza)], density=True, bins=np.arange(100), label= name + 'position i+1', histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Length of the homopolymer inserted')
    plt.ylabel('Count of the homopolyer')
    plt.legend()
    return plt.show()


# COUNT DIFFERENCES BETWEEN i and i+1
def delta_function(lista):
    lungo = len(lista)
    r = 0
    ritorno = np.zeros(len(lista))
    while r < lungo-1:
        tmp = lista[r] - lista[r+1]
        ritorno[r] = tmp
        r += 1
    # add last value
    ritorno[len(lista)-1] = lista[-1]

    return ritorno


# Creating unique dictonary for insertion and deletion
def deletionsAnsInsertion_dict(matrice_of_deletion, normal,dict_for,tot_for):

    # Count number of deletion
    number_of_deletion = {}

    for x in normal.keys():
        tmp = []
        for y in normal[x]:
            nn = matrice_of_deletion[''][y[0]:y[1]+1]
            
            tmp.append(nn)
        number_of_deletion[x] = tmp

    # Calculate differencies
    delta_deletions = defaultdict(lambda:defaultdict(int))
    for x in number_of_deletion.keys():
        tmp = []
        intervallo = np.arange(-1,-x[1]-1,-1)
        if x[1] == 1:
            for m in number_of_deletion[x]:
                
                #delta_deletions[x][intervallo[0]].append(m[0])
                delta_deletions[x][intervallo[0]] += m[0] 
        else:
            for g in number_of_deletion[x]:
                sorted_list = sorted(g,reverse=True)
                #print(g)
                differenze = delta_function(sorted_list)
                #print(differenze)
                for n,y in enumerate(intervallo):
                    #delta_deletions[x][y].append(differenze[n])
                    delta_deletions[x][y] += differenze[n]
    # Making sure that the dictionary of insertion and deletion has common keys
    inserzioni = {}

    for x in delta_deletions.keys():
        if x in dict_for:
            tmp = dict(sorted(Counter(dict_for[x]).items()))
        else:
            for l in range(1,x[1]+1):
                tmp[l] = 0
        inserzioni[x] = tmp


    # Adding toghether insertion and deletion data
    totale_dizio = {}
    for x in delta_deletions.keys():
        totale_dizio[x] = {**dict(sorted(delta_deletions[x].items())), **inserzioni[x]}

    # Findind delta 0

    # First summing all the count for insertion and delition for each nucleotide
    totale_read = defaultdict(int)

    for x in normal.keys():
        for y in normal[x]:
            totale_read[x] += tot_for[y[0]]

    # Creating new dicotionry includin also delta 0
    final_dizio = {}
    for x in totale_dizio.keys():
        n = 0
        for y in totale_dizio[x].keys():
            somma = totale_dizio[x][y]
            n += somma
        totale_dizio[x][0] = totale_read[x] - n
        final_dizio[x] = dict(sorted(totale_dizio[x].items()))

    return final_dizio