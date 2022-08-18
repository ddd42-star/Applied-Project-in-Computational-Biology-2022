# %%
import pysam
import numpy as np
from collections import Counter, defaultdict
import pickle
# %%
bamfile = pysam.AlignmentFile('reads.sorted.bam','rb')

# %%
# We create the matrix where at each position of the reference genome, there is the content of the reads
n = 0
matrice = defaultdict(list)

for colonna in bamfile.pileup(min_base_quality=0,flag_filter=2320): # Flag filter for onlyForward read. To see Reverse read see https://davetang.org/muse/2014/03/06/understanding-bam-flags/
        
    seq = colonna.get_query_sequences()
    #seq_for = [x for x in seq if x in ['A','C','G','T','']]
    tmp = Counter(seq)
    for l in ['A','C','G','T','']:
        if l in tmp:
            matrice[l].append(tmp[l])
        else:
            matrice[l].append(0)
# %%
# we save the matrix so we can use it again in other Notebooks
with open('matrice_deletion.pkl','wb') as fh:
    pickle.dump(matrice,fh)

# %%
# We check our matrix with the program of Prof. Richard Neher
ll = np.load('pileup/allele_counts/arr_0.npy')
# %%
