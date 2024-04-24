##########   LU2SV483 - Project 2024   ##########
import sys
from Project_Library import basic_seq_analysis as SeqAnalysis

##########          STEP 1           ##########
print('--------------STEP 1--------------')

# QUESTION 1
print('\n----------Question 1----------\n')
# SeqTest : str ; sequence for the question 1
SeqTest = 'GCGCATGCTTACATAGGCCTAACAAACGGCTTATGACTAG'

print('The sequence to analyse:\n',SeqTest)
print(f'The length of the sequence is {SeqAnalysis.len_seq(SeqTest)} nucleotides')
print('The GC content is ',SeqAnalysis.GC_content(SeqTest))

# revcompl : str, reverse complement sequence
revcompl = SeqAnalysis.rev_compl(SeqTest)
print('The reverse complement is:\n',revcompl)
# atgs_seq : List[int] ; list of positions of ATG in the sequence
# stops_seq : List[int] ; list of positions of STOP in the sequence
# atgs_rev : List[int] ; list of positions of ATG in the reverse complement sequence
# stops_rev : List[int] ; list of positions of STOP in the reverse complement sequence
atgs_seq, stops_seq = SeqAnalysis.ORF_position(SeqTest)
atgs_rev, stops_rev = SeqAnalysis.ORF_position(revcompl)

print('ATG positions for the sequence:',atgs_seq)
print('STOP positions for the sequence:',stops_seq)
print('ATG positions for the reverse complement:',atgs_rev)
print('STOP positions for the reverse complement:',stops_rev)

# QUESTION 2
print('\n----------Question 2----------\n')
# seqal : List[str] ; sequences for the question 2
seqal = ['AGTTAGAG--TAGGGCCAGCCAGCATAGCA-----GGTT',
         'AGTTCCAGTATAGGCCCA---AGCAAAGCAGTACCGGTT']
# match : int ; number of matches
# subst : int ; number of substitutuions
# indel : int ; number of indels
match, subst, indel = SeqAnalysis.compare_seqs(seqal[0],seqal[1])
# edit_dist : int ; edit distance (the sum of indels and substitutions)
edit_dist = SeqAnalysis.edit_distance(seqal[0],seqal[1])
# al_score : float ; alignment score (weighted sum of matches, substitutuions and indels)
al_score = SeqAnalysis.align_score(seqal[0],seqal[1],+2,-1,-2)
# id_ratio : float ; identity ratio (the proportion of matches)
id_ratio = SeqAnalysis.identity_ratio(seqal[0],seqal[1])
print('Comparison of two sequences:\n',seqal[0],'\n', seqal[1])
print(f'Matches: {match}  Substitutions: {subst}  Indels: {indel}')
print('Alignment score:', al_score)
print('Identity ratio:', round(id_ratio,2))
