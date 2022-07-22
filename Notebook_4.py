# %%
import homopolymerLength as hL



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
homopolymer_nucletides_sequence,altro = hL.find_homopolymer_nucleotide('assembled_genome/assembly.fna')

forward_seq,reverse_seq = hL.insertion_position_and_sequence('reads.sorted.bam')

# check with IGV the insertion

dict_start_for,dict_end_for, dict_start_rev, dict_end_rev = hL.find_length_homopolymer(forward_seq,reverse_seq,homopolymer_nucletides_sequence)

# Plot (you can choose which nucleotides and which dictonary to use)

hL.plot_homopolymer_length('A', dict_start_for)

hL.plot_comparison_between_start_end(dict_start_for,dict_end_for,'A', 1)



# %%
