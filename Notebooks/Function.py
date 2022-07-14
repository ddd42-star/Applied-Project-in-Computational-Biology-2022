from ast import main
from hashlib import new
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import Counter, defaultdict
from Bio import SeqIO
import pickle
from scipy.stats import binom_test


def histogramACGT(file):
    
    # read bam file
    bamfile = pysam.AlignmentFile(file,'rb')
    dizio = defaultdict(list)
    for ind in bamfile:
        
        if ind.is_unmapped or ind.is_secondary or ind.is_supplementary:
            continue
            
        if ind.is_reverse:
            continue
            
        # quality number at each position
        ind_qual = []
        if ind.qual is None:
            continue
        else:    
            for x in ind.qual:
                    ind_qual.append(ord(x) - 33)
        # list contaning the nucleotide
        nucleotide = [x for x in ind.query]
        # tuple of (Nucleotide, quality)
        list_of_quality = list(zip(nucleotide,ind_qual))

        # save for each nucleotide the quality score
        for x,y in list_of_quality:
            dizio[x].append(y)
    # find uniques
    nucleotides = unique_dict(dizio)
    
    return nucleotides#dict_A,dict_C,dict_T,dict_G

def unique_dict(dizionario):
    
    new_dict = {}
    for x in ["A","C","G","T"]:
        new_dict[x] = Counter(dizionario[x])
    
    return new_dict

def average_number_of_reads(filename):

    bamfile = pysam.AlignmentFile(filename,'rb')

    number_of_reads = []
    for col in bamfile.pileup(min_base_quality=0):
        number_of_reads.append(col.get_num_aligned())
    return number_of_reads

def average_quality(filename):
    
    bamfile = pysam.AlignmentFile(filename,'rb')

    avg_quality = []
    for col in bamfile.pileup(min_base_quality=0):
        avg_quality.append(np.mean(col.get_query_qualities()))

    return avg_quality

def plot_error(base,dizionario, filename, log=False, densita = False,reverse = None):
    """
    @log is True if one wants the y axis log scaled or False not log scaled
    @ densita True if one wants the distribution normalized or False not normalized
    @reverse True if one wants Forward or False reverse
    This function plot the distribution of quality from a dictionary containing:
    ("Nucleotide of reference", "Nucleotide of query", Forward(True)/Reverse(False) : quality count [])
    
    ------------------------------------------------------
    It Return the plot if the dictonary with the option of log scale, normalized, Forward/Reverse
    """
    leggenda = []
    plt.figure(figsize=(10, 8), dpi=80)
    for b in ['A','C','G','T']:
        plt.hist(dizionario[(base,b, reverse)].keys(),
                 weights=dizionario[(base,b,reverse)].values(),
                 bins=np.arange(100),
                 histtype='step',
                 density=densita,
                 label=base+b)
            
    plt.legend()
    plt.xlabel("Quality score")
    plt.ylabel("# of single nucleotides reads")
    #plt.savefig(base+'_reference_log')
    if log == True:
        plt.yscale('symlog')
    plt.savefig(filename)    
    return plt.show()
    
def create_int_defaultdict():
    return defaultdict(int)

    
def read_quality_seq_single(file,true_file):
    
    bamfile = pysam.AlignmentFile(file,'rb')
    true_sequence = []
    # saving in a list the true sequence
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    
    # Create the matrix of quality for forward and reverse
    #matrice = defaultdict(lambda : defaultdict(int))
    matrice = defaultdict(create_int_defaultdict)
    # dictionary for the reverse strand
    nucleotidi_reverse = {'a': 'T', 'c': 'G', 'g' : 'C', 't': 'A','':0}
    nucleotide_coniugato = {'A': 'T', 'C': 'G', 'G' : 'C', 'T': 'A'}
    
    # iterate over all the genome
    for col in bamfile.pileup(min_base_quality=0):
        # quality at n position
        q = col.get_query_qualities()
        # nucleotide at n position
        b = col.get_query_sequences()
        # make a list of tuple with (quality value, nucleotide)
        assembly = list(zip(q,b))
        # iterate over the list and save in the dictonary
        for q, unique in assembly:
            if unique == '':
                continue
            #Reverse
            if unique in ['a','g','t','c']:
                reverse_nucleotide = nucleotide_coniugato[true_sequence[col.pos]]
                matrice[(reverse_nucleotide,nucleotidi_reverse[unique], False)][q] += 1     
            
            #Forward
            else:
                true_nucleotide = true_sequence[col.pos]
                matrice[(true_nucleotide,unique, True)][q] += 1
    
    new_matrix = {}
    
    for x in matrice.keys():
        new_matrix[x] = matrice[x]
            
    return new_matrix


def percentage_error(base,dizionario,reverse):
    
    lista_percetage, lista_error = [], []
    Nset = set(['A','C','G','T'])
    for y in range(1,101,1):
        somma_tot = 0
        somma = 0
        for x in ["A","C","G","T"]:
            somma_tot += dizionario[(base,x, reverse)][y]
            
        compl_nucleotide = Nset - set([base])
        for l in compl_nucleotide:
            somma += dizionario[(base,l,reverse)][y]
        if somma_tot < 1000:
            percentage = np.NaN
            error = np.NaN
        else:
            percentage = somma/somma_tot
            error = ((percentage * (1- percentage)) / somma_tot)   
        lista_percetage.append(percentage)
        lista_error.append(np.sqrt(error))
        

    return lista_percetage, lista_error

def plot_error_rate_pro_quality(dizionario, errorbar):
    
    
    plt.figure(figsize=(10,8))
    i = np.arange(40)
    for x in ["A","C","G","T"]:
        #plt.plot(np.arange(1,101,1),dizionario[x],'.--',label=x)
        plt.errorbar(np.arange(1,101,1), dizionario[x],errorbar[x], marker='.', linestyle='',label=x)
    plt.plot(i, np.power(10,(-i/10)),"--", label='10^(-q/10)', color='gray')
    plt.xlabel("Single quality score")
    plt.ylabel("Fraction of error")
    plt.yscale('symlog', linthresh=1e-5)
    plt.legend()
    return plt.show()
def sum_quality (base,query, dizionario,reverse, scale =None):
    
    somma = 0
    # quality < 20
    if scale == True:
        for s in {x : y for x,y in dizionario[(base,query,reverse)].items() if x < 20}.values():
            somma += s
    elif scale == False:# quality >=20
        for s in {x : y for x,y in dizionario[(base,query,reverse)].items() if x >= 20}.values():
            somma += s
    else:
        for s in dizionario[(base,query,reverse)].values():
            somma += s
    return somma
def sum_quality_tot (base, dizionario, reverse, scale = None):
    n = 0
    # quality < 20
    if scale == True:
        for x in ["A", "C", "G", "T"]:
            for y in {x : y for x,y in dizionario[(base,x,reverse)].items() if x < 20}.values():
                n += y
    elif scale == False:# quality >= 20
        for x in ["A", "C", "G", "T"]:
            for y in {x : y for x,y in dizionario[(base,x,reverse)].items() if x >= 20}.values():
                n += y
    else:
        for x in ["A", "C", "G", "T"]:
            for y in dizionario[(base,x,reverse)].values():
                n += y
    return n

def matrix_of_error_single(dizionario, reverse, scale):
    
    #initialized a 4x4 matrix
    matrice = np.zeros((4,4))
    #quality < 20
    if scale == True:
        for n,m in enumerate(["A", "C","G","T"]):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                #if tot == 0:
                #    continue
                matrice[n][y] = somma/tot  
    elif scale == False:# quality >= 20
        for n,m in enumerate(["A", "C","G","T"]):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                #if tot == 0:
                #    continue
                matrice[n][y] = somma/tot 
    else:# quality >= 20
        for n,m in enumerate(["A", "C","G","T"]):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                #if tot == 0:
                #    continue
                matrice[n][y] = somma/tot 
    return matrice

def plot_matrix_error(lista):
    
    fig, axs = plt.subplots(2,3, figsize=(12,10))
    norm = mpl.colors.LogNorm(vmax=1, vmin=1e-4)
    i,j = 0,0
    index = ["A","C","G","T"]
    titolo = ['Forward_less_20','Forward_plus_20', 'Forward_all',
              'Reverse_less_20','Reverse_plus_20', 'Reverse_all']
    
    lungo = 0
    for x in lista:
        m = axs[i,j].matshow(x,norm=norm,cmap="jet", alpha=0.5)
        for l in range(len(index)):
                for k in range(len(index)):
                    axs[i,j].text(k,l,x[l,k].round(decimals=4), ha="center", va="center", color="black")
        axs[i,j].set_title(titolo[lungo])
        lungo += 1
        j += 1
        if j == 3:
            j = 0
            i += 1
    plt.colorbar(m, label='Frequency', ax=axs)
    plt.setp(axs, xticks=range(len(index)), xticklabels=index,
        yticks=range(len(index)), yticklabels=index)
    

    

    
    return plt.show()

def matrix_of_error_triplet(dizionario, reverse, scale):
    
    #initialized a 64x4 matrix
    matrice = np.zeros((64,4))
    possibilities = []
    for x in ["A","T","C","G"]:
        for y in ["A","T","C","G"]:
            for z in ["A","T","C","G"]:
                possibilities.append(y+x+z)
    #quality < 20
    if scale == True:
          for n,m in enumerate(possibilities):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                if tot == 0:
                    continue
                matrice[n][y] = somma/tot  
    elif scale == False:# quality >= 20
        for n,m in enumerate(possibilities):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                if tot == 0:
                    continue
                matrice[n][y] = somma/tot 
    else:# quality >= 20
        for n,m in enumerate(possibilities):
            tot = sum_quality_tot(m, dizionario, reverse, scale)
            for y,x in enumerate(["A", "C","G", "T"]):
                somma = sum_quality(m,x, dizionario, reverse, scale)
                if tot == 0:
                    continue
                matrice[n][y] = somma/tot 
    return matrice


def read_quality_seq_triplet(file,true_file):
    
    bamfile = pysam.AlignmentFile(file,'rb')
    true_sequence = []
    # saving in a list the true sequence
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    
    true_string = ''.join(true_sequence)

    matrice = defaultdict(lambda : defaultdict(int))
    for col in bamfile.pileup(min_base_quality=0):
        # quality at n position
        q = col.get_query_qualities()
        # nucleotide at n position
        b = col.get_query_sequences()
        # make a list of tuple with (quality value, nucleotide)
        assembly = list(zip(q,b))
        
        for quality, unique in assembly:
            if unique == '':
                continue
            if unique in ["A", "C", "G", "T"]:
                if col.pos == len(true_string) - 1:
                    triplet = true_string[col.pos-1]+true_string[col.pos]+true_string[0]
                else:
                    triplet = true_string[col.pos-1]+true_string[col.pos]+true_string[col.pos+1] 
                
                matrice[(triplet,unique,True)][quality] += 1
    
    
    new_matrix = {}
    
    for x in matrice.keys():
        new_matrix[x] = matrice[x]
            
    return new_matrix 


def plot_of_error_matrix_triplet(matrice,possibilities):

    norm = mpl.colors.LogNorm(vmax=1, vmin=1e-4)
    plt.figure(figsize=(10,8))
    plt.matshow(matrice,norm=norm,cmap="jet")
    plt.yticks(range(64),possibilities)
    plt.xticks(range(4),['A','C','G','T'])
    plt.colorbar()

    return plt.show()

def somma_of_quality_triplet(base,dizionario,reverse):
    somma = 0
    Nset = set(['A','C','G','T'])
    
    compl_nucleotide = Nset - set([base[1]])
    
    for l in compl_nucleotide:
        value = sum(dizionario[(base,l,reverse)].values())
        somma += value
    return somma

def somma_tot(base, dizionario, reverse):
    
    n = 0
    
    for x in ["A", "C", "G", "T"]:
        for g in dizionario[(base,x,reverse)].values():
            n += g
    return n


def percentage_error_triplet(dizionario,reverse,possibilities):
    
    lista_percetage = {}
    for m in possibilities:
        tot = somma_tot(m, dizionario, reverse)
        somma = somma_of_quality_triplet(m, dizionario, reverse)
        lista_percetage[m] = somma/tot  
    return lista_percetage

def percentage_error_single(dizionario,reverse,possibilities):
    
    lista_percetage = {}
    
    for m in possibilities:
        tot = somma_tot(m, dizionario, reverse)
        somma = somma_of_quality_single(m, dizionario, reverse)
        lista_percetage[m] = somma/tot  
    return lista_percetage

def somma_of_quality_single(base,dizionario,reverse):
    somma = 0
    Nset = set(['A','C','G','T'])
    
    compl_nucleotide = Nset - set([base])
    
    for l in compl_nucleotide:
        value = sum(dizionario[(base,l,reverse)].values())
        somma += value
    return somma

def read_quality_seq_quintuple(file,true_file):
    
    bamfile = pysam.AlignmentFile(file,'rb')
    true_sequence = []
    # saving in a list the true sequence
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    
    true_string = ''.join(true_sequence)

    matrice = defaultdict(create_int_defaultdict)
    for col in bamfile.pileup(min_base_quality=0):
        # quality at n position
        q = col.get_query_qualities()
        # nucleotide at n position
        b = col.get_query_sequences()
        # make a list of tuple with (quality value, nucleotide)
        assembly = list(zip(q,b))
        
        for quality, unique in assembly:
            if unique == '':
                continue
            if unique in ["A", "C", "G", "T"]:
                if col.pos == 0:
                    quintuple = true_string[-2:]+true_string[:3]
                elif col.pos == len(true_string)-1:
                    quintuple = true_string[-3:]+true_string[:2]
                elif col.pos == 1:
                    quintuple = true_string[-1]+true_string[:4]
                elif col.pos == (len(true_string) - 2):
                    quintuple = true_string[-4:]+true_string[0]
                else:
                    quintuple = true_string[col.pos-2:col.pos+3] 
                
                
                matrice[(quintuple,unique,True)][quality] += 1
    
    
    new_matrix = {}
    
    for x in matrice.keys():
        new_matrix[x] = matrice[x]
            
    return new_matrix


def somma_of_quality_quintuple(base,dizionario,reverse):
    
    somma = 0
    Nset = set(['A','C','G','T'])
    compl_set = Nset - set([base[2]])
    
    for l in compl_set:
        values = dizionario[(base,l,reverse)].values()
        somma += sum(values)
        
    return somma

def percentage_error_quintuple(dizionario,reverse,possibilities):
    
    lista_percetage = {}
    
    for m in possibilities:
        tot = somma_tot(m, dizionario, reverse)
        somma = somma_of_quality_quintuple(m, dizionario, reverse)
        lista_percetage[m] = somma/tot  
    return lista_percetage
 
def dictonary_Contex(file,true_file, distance):
    #parse bam file containing reads
    bamfile = pysam.AlignmentFile(file,'rb')
    #parse reference genome file
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)
    L_genoma = len(true_string)
    
    dizionario = defaultdict(int)
    percentage_dict = {}
    lunghi_dizionari = {}
    
    for col in bamfile.pileup(min_base_quality=0):
        # nucleotide at n position
        nucleotides_read = col.get_query_sequences()
        
        for x in nucleotides_read:
            if x in ['A','C','G','T']:

                nucleotide_i = true_string[col.pos]
                il = (col.pos + distance) % L_genoma
                nucleotide_il = true_string[il]

                if nucleotide_i == x:
                    dizionario[(nucleotide_i+nucleotide_il,True)] += 1
                else:
                    dizionario[(nucleotide_i+nucleotide_il,False)] += 1
                        
    # Find Percentage of P(e,a)
    sum_tot = 0
    for x in dizionario.values():
        sum_tot += x

    #update values
    for y in dizionario.keys():
        percentage_dict[y] = dizionario[y] / sum_tot
        
    return percentage_dict


def entropy(dizionario):
    
    #assemble dictonary with total probabilities of contest
    def dizionario_prob_A(alpha):
        p = 0
        for t in [True, False]:
            k = (alpha, t)
            if k in dizionario:
                p += dizionario[k]
        return p
        
    # Find entropy    
    Entropy = 0 
    
    for x in dizionario.keys():
        Entropy += (dizionario[x] * np.log2(dizionario[x]/dizionario_prob_A(x[0])))
    
    result = -Entropy
    return result

def Contex(filename,true_name,distances):

    results = []
    list_of_dict = {}
    for x in distances:
        list_of_dict[x] = dictonary_Contex(filename,true_name,x)
        print(x)

    #for x in list_of_dict:
    #    results.append(entropy(x))
    return list_of_dict

def Find_Entropy(dizionari, distanze):
    risultato = []
    for x in distanze:
        risultato.append(entropy(dizionari[x]))
    return risultato
def dictonary_Contex_sequence(file,true_file, distance):
    #parse bam file containing reads
    bamfile = pysam.AlignmentFile(file,'rb')
    #parse reference genome file
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)
    L_genoma = len(true_string)
    
    dizionario_forward = defaultdict(int)
    dizionario_reverse = defaultdict(int)
    dizionario_mix = defaultdict(int)
    percentage_dict = {}
    lunghi_dizionari = {}
    
    for col in bamfile.pileup(min_base_quality=0):
        # nucleotide at n position
        nucleotides_read = col.get_query_sequences()
        read_i = true_string[col.pos]
        
        I = range(col.pos - distance + 1, col.pos + distance)
        I = [i % L_genoma for i in I]
        tag = "".join([true_string[i] for i in I])
        tag_complementary = reverse_Complementary(tag)
        read_complementary = reverse_Complementary(read_i)
        
        
        for x in nucleotides_read:
            
            if x == '':
                continue
            nucleotide_complementary = reverse_Complementary(x)    
            if x in ['A','C','G','T']:
                dizionario_forward[(tag,read_i == x)] += 1
                dizionario_mix[(tag,read_i == x)] += 1
            else:
                
                dizionario_reverse[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
                dizionario_mix[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
            
    forward_percentage = find_percentage(dizionario_forward)
    reverse_percentage = find_percentage(dizionario_reverse)
    mix_percentage = find_percentage(dizionario_mix)
    
    
    return forward_percentage,reverse_percentage,mix_percentage
 
def find_percentage(dizionario):
    
    
    sum_tot = 0
    new_dizionario = {}
    
    for x in dizionario.values():
        sum_tot += x
        
    #update values
    for y in dizionario.keys():
        new_dizionario[y] = dizionario[y] / sum_tot
    
    return new_dizionario

def Contex_multiple_sequence(filename,true_name,distances):

    results_forward,results_reverse,results_mix = [], [],[]
    
    for x in distances:
        forward,reverse,mix = dictonary_Contex_sequence(filename,true_name,x)
        print(x)
        results_forward.append(entropy(forward))
        results_reverse.append(entropy(reverse))
        results_mix.append(entropy(mix))
        
    return results_forward,results_reverse,results_mix


def dictonary_Contex_long_sequence(file,true_file, distance):
    #parse bam file containing reads
    bamfile = pysam.AlignmentFile(file,'rb')
    #parse reference genome file
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)
    L_genoma = len(true_string)
    
    dizionario_forward = defaultdict(int)
    dizionario_reverse = defaultdict(int)
    dizionario_mix = defaultdict(int)
    percentage_dict = {}
    lunghi_dizionari = {}
    
    for col in bamfile.pileup(min_base_quality=0):
        # nucleotide at n position
        nucleotides_read = col.get_query_sequences()
        read_i = true_string[col.pos]
        
        I = range(col.pos,col.pos + distance)
        I = [i % L_genoma for i in I]
        tag = "".join([true_string[i] for i in I])
        I_r = range(col.pos - distance + 1,col.pos + 1)
        I_r = [x % L_genoma for x in I_r]
        tag_r = ''.join([true_string[j] for j in I_r])
        tag_complementary = reverse_Complementary(tag_r)
        read_complementary = reverse_Complementary(read_i)
        
        
        for x in nucleotides_read:
            
            if x == '':
                continue
            nucleotide_complementary = reverse_Complementary(x)    
            if x in ['A','C','G','T']:
                dizionario_forward[(tag,read_i == x)] += 1
                dizionario_mix[(tag,read_i == x)] += 1
            else:
                
                dizionario_reverse[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
                dizionario_mix[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
            
    forward_percentage = find_percentage(dizionario_forward)
    reverse_percentage = find_percentage(dizionario_reverse)
    mix_percentage = find_percentage(dizionario_mix)
    
    
    return forward_percentage,reverse_percentage,mix_percentage
def dictonary_Contex_reverse_long_sequence(file,true_file, distance):
    #parse bam file containing reads
    bamfile = pysam.AlignmentFile(file,'rb')
    #parse reference genome file
    true_sequence = []
    for seq in SeqIO.read(true_file, 'fasta'):
        true_sequence.append(seq)
    true_string = ''.join(true_sequence)
    L_genoma = len(true_string)
    
    dizionario_forward = defaultdict(int)
    dizionario_reverse = defaultdict(int)
    dizionario_mix = defaultdict(int)
    percentage_dict = {}
    lunghi_dizionari = {}
    
    for col in bamfile.pileup(min_base_quality=0):
        # nucleotide at n position
        nucleotides_read = col.get_query_sequences()
        read_i = true_string[col.pos]
        
        I = range(col.pos - distance + 1,col.pos + 1)
        I = [i % L_genoma for i in I]
        tag = "".join([true_string[i] for i in I])
        I_r = range(col.pos,col.pos + distance)
        I_r = [x % L_genoma for x in I_r]
        tag_r = ''.join([true_string[j] for j in I_r])
        tag_complementary = reverse_Complementary(tag_r)
        read_complementary = reverse_Complementary(read_i)
        
        
        for x in nucleotides_read:
            
            if x == '':
                continue
            nucleotide_complementary = reverse_Complementary(x)    
            if x in ['A','C','G','T']:
                dizionario_forward[(tag,read_i == x)] += 1
                dizionario_mix[(tag,read_i == x)] += 1
            else:
                
                dizionario_reverse[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
                dizionario_mix[(tag_complementary,read_complementary == nucleotide_complementary)] += 1
            
    forward_percentage = find_percentage(dizionario_forward)
    reverse_percentage = find_percentage(dizionario_reverse)
    mix_percentage = find_percentage(dizionario_mix)
    
    
    return forward_percentage,reverse_percentage,mix_percentage

def reverse_Complementary(string):
    
    
    reverse = {'A':'T','C':'G','G':'C','T':'A','a':'T','c':'G','g':'C','t':'A'}
    l = []
    for x in string[::-1]:
        l.append(reverse[x])
    string_result = ''.join(l)
    return string_result


def Contex_long_Sequence(filename,true_name,distances):

    results_forward,results_reverse,results_mix = [], [],[]
    
    for x in distances:
        forward,reverse,mix = dictonary_Contex(filename,true_name,x)
        print(x)
        results_forward.append(entropy(forward))
        results_reverse.append(entropy(reverse))
        results_mix.append(entropy(mix))
        
    return results_forward,results_reverse,results_mix

def Contex_reverse_long_Sequence(filename,true_name,distances):

    results_forward,results_reverse,results_mix = [], [],[]
    
    for x in distances:
        forward,reverse,mix = dictonary_Contex_reverse(filename,true_name,x)
        print(x)
        results_forward.append(entropy(forward))
        results_reverse.append(entropy(reverse))
        results_mix.append(entropy(mix))
        
    return results_forward,results_reverse,results_mix
