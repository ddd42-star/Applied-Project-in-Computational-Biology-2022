# %%
import homopolymerLength as hL
import pickle
import pandas as pd
import seaborn as sea
# %%
# First we wanto to find how the insertion along all the reads are distributed.


dizionario_for,dizionario_rev = hL.find_length_insertion('reads.sorted.bam')


# Now we are interested in looking at the distribution of the length of INSERTION, DELETION, HARD AND SOFT CLIP 
# in all the reads along the entire genome


# INSERTION
hL.plot_insertion_length(dizionario_for, dizionario_rev)


# DELETION
hL.plot_deletion_length(dizionario_for, dizionario_rev)

# HARD and SOFT CLIP
hL.plot_Soft_and_Hard_clip_length(dizionario_for, dizionario_rev)

# %%
# Secondly we look the quality of the insertion. As we did for the matches, we want to know how good or bad the insertion are guessed
# from the machine.


forward_qual, reverse_qual = hL.quality_insertion('reads.sorted.bam')



# Plot of this quality intervall

hL.plot_of_the_quality_of_insertion(forward_qual, reverse_qual)

# %%
# Finally we look into the stratification of the insertion in all the reads, to see at the distribution in all the reads
# First we search in all the genome each homopolymer combination for A,C,G and T. Then we save in a list the information
#  of the position i and position i + 1 where on the genome happen the insertion and we also save the sequence of the insertion.
# To prove that the inserted sequence is the actual one, we checked with the tool IGV and see if in each site saved, there is 
# the right insertion. Now we can find the length of the insertion in the homopolymer, stratifing by length of the homopolymer
# found in the genome. Finally, we plot the count of homopolymer inserted along thei length, stratifing by the different type of 
# homopolymer find in the genome.
homopolymer_nucletides_sequence,normal = hL.find_homopolymer_nucleotide_new('assembled_genome/assembly.fna')

forward_seq,reverse_seq,tot_for,tot_rev = hL.insertion_position_and_sequence('reads.sorted.bam')

# check with IGV the insertion

dict_for,dict_rev = hL.find_length_homopolymer_new(forward_seq,reverse_seq,homopolymer_nucletides_sequence)

# %%
# Plot of Distribution length of the homopolymer
# Number of homopolymer for each class of homopolymer length
dis_A,dis_C,dis_G,dis_T = {},{},{},{}
for x in normal.keys():
    if 'A' in x[0]:
        dis_A[x] = len(normal[x])
    elif 'C' in x[0]:
        dis_C[x] = len(normal[x])
    elif 'G' in x[0]:
        dis_G[x] = len(normal[x])
    else:
        dis_T[x] = len(normal[x])
# Plot the Distribution
hL.plt.figure(figsize=(10,8))
hL.plt.hist(hL.np.arange(1,10),weights= dis_A.values(),histtype='step', label='A', bins=hL.np.arange(11)+0.5)
hL.plt.hist(hL.np.arange(1,10),weights= dis_C.values(),histtype='step', label='C',bins=hL.np.arange(11)+0.5)
hL.plt.hist(hL.np.arange(1,10),weights= dis_G.values(),histtype='step', label='G',bins=hL.np.arange(11)+ 0.5)
hL.plt.hist(hL.np.arange(1,10),weights= dis_T.values(), histtype='step', label='T',bins=hL.np.arange(11) + 0.5)
hL.plt.yscale('log')
hL.plt.xlabel('length of the homopolymer')
hL.plt.ylabel('Number of of homopolymer for each nucleotide')
hL.plt.legend()
hL.plt.show()
#hL.plt.savefig('Length of homopolymer.png')

# %%
# Find real length of the homopolymer
# Now Considering deletion
with open('matrice_deletion.pkl','rb') as fh:
    matrice_deletion = pickle.load(fh)

# Find unique dictonary with real length distribution
final_dictonary = hL.deletionsAnsInsertion_dict(matrice_deletion,normal,dict_for,tot_for)

# Plot delta distribution
# create dataframe
m = pd.DataFrame(final_dictonary)

# create new indexes
nuova_forma = m.unstack().reset_index()

# create new dataframe in order to plot the dictonary
df = {
    "nucl" : [],
    "L" : [],
    "Delta" : []
}
for idx, row in nuova_forma.iterrows():
    N = row[0]
    if hL.np.isnan(N):
        continue
    N = int(N)
    nucl, L, delta = row.iloc[0:3]
    df['nucl'] += [nucl] * N
    df['L'] += [L] * N
    df['Delta'] += [delta] * N

df = pd.DataFrame(df)

# The last steps takes lot of times, so we calculated once and saved the data in a pickle object
with open('DataFrame.pkl','wb') as fh:
    pickle.dump(df,fh)

# Open the pickle object
with open('DataFrame.pkl','rb') as f:
    dati_finali = pickle.load(f)

# plot the new data
hL.plt.figure(figsize=(10,8))
sea.boxenplot(x='L',y='Delta',hue=dati_finali['nucl'],data=dati_finali) # here we show only the plot with the differencies cases for ACGT
#sea.boxenplot(x='L',y='Delta',data=dati_finali) # here we show the general plot
hL.plt.ylim(bottom=-20,top=20)
hL.plt.show()
#hL.plt.savefig('DeltaHomoPolymerNoDivision.png')