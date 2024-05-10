from Project_Library import basic_seq_analysis as SeqAnalysis

'''Question 1'''

print('The sequence to analyze:\n',SeqAnalysis.TEST_SEQ, sep='')
print(f'\nThe length of the sequence is {SeqAnalysis.len_seq(SeqAnalysis.TEST_SEQ)} nucleotides')
print('The GC content is ',SeqAnalysis.GC_content(SeqAnalysis.TEST_SEQ))

# revcompl : str, reverse complement sequence
revcompl = SeqAnalysis.rev_compl(SeqAnalysis.TEST_SEQ)
print('\nThe reverse complement is:\n',revcompl, sep='')

# atgs_seq : List[int] ; list of positions of ATG in the sequence
# stops_seq : List[int] ; list of positions of STOP in the sequence
# atgs_rev : List[int] ; list of positions of ATG in the reverse complement sequence
# stops_rev : List[int] ; list of positions of STOP in the reverse complement sequence
atgs_seq, stops_seq = SeqAnalysis.ORF_position(SeqAnalysis.TEST_SEQ)
atgs_rev, stops_rev = SeqAnalysis.ORF_position(revcompl)

print('\nATG positions for the sequence:', *atgs_seq)
print('STOP positions for the sequence:', *stops_seq)
print('ATG positions for the reverse complement:', *atgs_rev)
print('STOP positions for the reverse complement:', *stops_rev)

'''Question 2'''

# match : int ; number of matches
# subst : int ; number of substitutions
# indel : int ; number of indels
match, subst, indel = SeqAnalysis.compare_seqs(SeqAnalysis.TEST_ALIGNED_SEQ[0],SeqAnalysis.TEST_ALIGNED_SEQ[1])
# edit_dist : int ; edit distance (the sum of indels and substitutions)
edit_dist = SeqAnalysis.edit_distance(SeqAnalysis.TEST_ALIGNED_SEQ[0],SeqAnalysis.TEST_ALIGNED_SEQ[1])

# Wm, Ws, Wi : float ; weights for a match, substitution and indel respectively 
# W : str ; entered weight
Wm, Ws, Wi = [float(W) for W in input('\nEnter weights for a match, substitution and indel (separated by space)\n').split()]

# al_score : float ; alignment score (weighted sum of matches, substitutions and indels)
al_score = SeqAnalysis.align_score(SeqAnalysis.TEST_ALIGNED_SEQ[0],SeqAnalysis.TEST_ALIGNED_SEQ[1],Wm,Ws,Wi)
# id_ratio : float ; identity ratio (the proportion of matches)
id_ratio = SeqAnalysis.identity_ratio(SeqAnalysis.TEST_ALIGNED_SEQ[0],SeqAnalysis.TEST_ALIGNED_SEQ[1])

print('\nComparison of two sequences:\n',SeqAnalysis.TEST_ALIGNED_SEQ[0],'\n', SeqAnalysis.TEST_ALIGNED_SEQ[1], sep='')
print(f'\nMatches: {match}  Substitutions: {subst}  Indels: {indel}')
print('Edit distance:', edit_dist)
print('Alignment score is', al_score, f'(for Wm={Wm} Ws={Ws} Wi={Wi})')
print('Identity ratio:', round(id_ratio,2))
