import sys
sys.path.append('Project_library')

from matplotlib import pyplot as plt
from copy import deepcopy
import embl_analysis as EMBL
import numpy as np
import msa

def write_msa_as_fasta(path_to_unaligned_seqs, path_to_substitutuion_matrix, path_to_output):
    '''
    str * str * str -> None
    Assumtion: unaligned sequences are in fasta format file, substitution matrix is in mat format file
    Writes in an output fasta format file a reslut of Multiple Sequence Alignement of unaligned sequences
    using a given  substitution matrix
    '''
    # ProtSeqs : Dict{str, str} ; unaligned sequences {UniProtID : sequence}
    ProtSeqs = EMBL.read_fasta(path_to_unaligned_seqs)
    # SubMat : Dict{(str,str):float} ; substitution matrix
    SubMat = read_score_matrix(path_to_substitutuion_matrix)
    # MSA_score : float ; multiple alignement score
    # MSA_seqs : Dict{str:str} ; aligned sequences {UniProtID : sequence}
    MSA_score, MSA_seqs = msa.star_align(ProtSeqs, SubMat)
    print(len(ProtSeqs), 'sequences, score =', MSA_score)

    # stream_w : stream ; output file stream
    with open(path_to_output, 'w') as stream_w:
        # Note : ';' is used as a comment for fasta format
        stream_w.write(f'; {len(ProtSeqs)} sequences, score = {MSA_score}\n')
        # uniprot : str ; UniProtID
        for uniprot in MSA_seqs:
            stream_w.write(f'>{uniprot}\n{MSA_seqs[uniprot]}\n')

def msa_alignment_score(MSA_dict, w_match, w_mismatch, w_indel):
    '''
    Dict{str:str} * float * float * float -> float, List[float]
    Assumtion: MSA_dict contains same sized aligned sequences as values
    Returns an overall alignment score and a list of scores by position
    Note : a score for each positition is a sum of all pairwise scores for this position,
    when a gap matches another gap - the score is zero
    '''
    # MSA : List[str] ; list of aligned sequences
    MSA = list(MSA_dict.values())
    # overall_score : float
    overall_score = 0.0
    # scores_by_pos : List[float] ; list of scores by position
    scores_by_pos = []
    # pos : int ; position (index) considered across the sequences
    for pos in range(len(MSA[0])):
        # score : float ; a score to be calculated for one position
        score = 0.0

        # A nested 'for' loop is used to compare all possible paires of sequences
        # i : int ; index of the first sequence to consider
        for i in range(len(MSA)):
            # j : int ; index of the second sequence to consider
            for j in range(i+1,len(MSA)):
                # seq1, seq2 : str ; first and second sequence
                seq1 = MSA[i]
                seq2 = MSA[j]
                # When a gap matches another gap - the score is zero
                if not (seq1[pos] == '-' and seq2[pos] == '-' ):
                    if seq1[pos] == seq2[pos]:
                        score += w_match
                    elif seq1[pos] == '-' or seq2[pos] == '-':
                        score += w_indel
                    else:
                        score += w_mismatch
        scores_by_pos.append(score)
        overall_score += score

    return overall_score, scores_by_pos

def read_score_matrix(path):
    '''
    str -> Dict{(str,str):float}
    Assumtion: path leads to a mat format file with a substitution matrix,
    the first not empty and not commentary line is concidered as a header
    Note: '*' is considered as gap and substituted by '-' in an output dictionary
    Returns a dictionary {(aminoacid1, aminoacid2):score}
    '''
    # Score_dict : Dict{(str,str):float} ; dictionary {(aminoacid1, aminoacid2):score}
    Scores_dict = {}
    # stream_r : stream ; input file stream
    with open(path,'r') as stream_r:
        # line : str ; a line read from input file
        for line in stream_r:
            # The first not empty line which is not commentary is concidered as a header
            if line.strip()[0] != '#' and line.strip() != '' :
                # header : List[str] ; a header of the table
                header = line.split()
                break
        # line : str ; a line read from input file after a header
        for line in stream_r:
            if line.strip() != '':
                # aa_i : str ; aminoacid corresponding to a current row
                aa_i = line.split()[0]
                if aa_i == '*':
                    aa_i = '-'
                # j : int ; index of a current column
                j=0
                for score in line.split()[1:]:
                    # aa_j : str ; aminoacid corresponding to a current column
                    aa_j = header[j]
                    if aa_j == '*':
                        aa_j = '-'
                    Scores_dict[(aa_i,aa_j)]= float(score)
                    j+=1
    return Scores_dict

def dissimilarity_score(MSA_dict,Score_dict):
    '''
    Dict{str:str} * Dict{(str,str):float} -> float, List[float]
    Assumtion: MSA_dict contains same sized aligned porotein sequences as values,
    Score_dict is of a {(aminoacid1, aminoacid2):score} format
    Returns an overall alignment score and a list of scores by position
    '''
    # MSA : List[str] ; list of aligned sequences
    MSA = list(MSA_dict.values())
    # overall_score : float
    overall_score = 0
    # scores_by_pos : List[float] ; list of scores by position
    scores_by_pos = []
    # pos : int ; position (index) considered across the sequences
    for pos in range(len(MSA[0])):
        # score : float ; a score to be calculated for one position
        score = 0

        # A nested 'for' loop is used to compare all possible paires of sequences
        # i : int ; index of the first sequence to consider
        for i in range(len(MSA)):
            # j : int ; index of the second sequence to consider
            for j in range(i+1,len(MSA)):
                # seq1, seq2 : str ; first and second sequence
                seq1 = MSA[i]
                seq2 = MSA[j]
                score += Score_dict[(seq1[pos],seq2[pos])]
        scores_by_pos.append(score)
        overall_score += score

    return overall_score, scores_by_pos

def compute_dissimilarity_matrix(MSA_dict):
    '''
    Dict{str:str} -> List[List[float]], List[str]
    Assumtion: MSA_dict contains same sized aligned porotein sequences as values
    Note: positions with indels don't count for a total number of positions
    Returns a lower triangular half of a dissimilarity matrix (diagonal included)
    and a list of sequence names (header of the matrix)
    '''
    # seq_names : List[str] ; a header of dissimilarity table with sequence names
    seq_names = []
    # MSA : List[str] ; list of aligned sequences
    MSA = []

    # seq_name : str ; sequence name - a key of MSA_dict
    for seq_name in MSA_dict :
        seq_names.append(seq_name)
        MSA.append(MSA_dict[seq_name])
    # length : int ; the size of aligned sequences (must be the same)
    length = len(MSA[0])
    # We don't count positions with indels for a total number of positions, to calculate a propostion after
    # DistMat : List[List[float]]
    DistMat=[]
    # i : int ; index of a 1st sequence
    for i in range(len(MSA)):
        # seq1 : str ; the 1st sequence
        seq1 = MSA[i]
        DistMat.append([])
        # j : int ; index of a 2nd sequence
        for j in range(i+1):
            # seq2 : str ; the 2nd sequence
            seq2 = MSA[j]
            DistMat[-1].append(0)
            # nb_pos_wo_gaps : int ; a counter of positions without gaps
            nb_pos_wo_gaps = 0
            # nb_pos_w_diff_aa : int ; a counter of positions with different aminoacids
            nb_pos_w_diff_aa = 0
            # pos : int ; a position in a sequence
            for pos in range(length):
                if seq1[pos] != '-' and seq2[pos] != '-':
                    nb_pos_wo_gaps += 1
                    if seq1[pos] != seq2[pos] :
                        nb_pos_w_diff_aa += 1
            DistMat[-1][-1]+=nb_pos_w_diff_aa/nb_pos_wo_gaps

    return DistMat, seq_names

def triangular_to_rectangular(triang_Mat):
    '''
    List[List[int]] -> List[List[int]]
    Assumption: triang_Mat is a a lower triangular half of a matrix
    Returns a complete rectangular matrix
    Note: if rectangular matrix is given as an argument, the same matrix will be returned
    '''
    # Mat : List[List[int]] ; a triangular matrix to be completed
    Mat = deepcopy(triang_Mat)
    # i : int ; a row index
    for i in range(len(Mat)):
        # j : int ; a column index
        for j in range(len(Mat[i]),len(Mat)):
            Mat[i].append(Mat[j][i])
    return Mat

def plot_alignment_score(score_list, title, path_to_output):
    '''
    List[float] * str * str -> None
    Creates and saves a matplotlib figure of alignment score by position
    (I've chosen a heatmap as a nice representation)
    '''
    # fig : matplotlib.figure.Figure object  
    # ax : matplotlib.axes._axes.Axes object
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_yticks([])   # to get rid of useless y-axis ticks
    plt.imshow([score_list], cmap = 'grey',aspect='auto')
    plt.colorbar()
    plt.title(title)
    plt.savefig(path_to_output)

def plot_dissimilarity_matrix_as_heatmap(triang_Mat, seq_names,path_to_output):
    '''
    List[List[float]] * str * str -> None
    Assumptions: triang_Mat is a lower triangular half of dissimilarity matrix (diagonal included)
    Creates and saves a matplotlib figure of a dissimilarity matrix as a heatmap
    '''
    # fig : matplotlib.figure.Figure object  
    # ax : matplotlib.axes._axes.Axes object
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.imshow(triangular_to_rectangular(triang_Mat), cmap = 'YlOrRd', aspect='auto')
    plt.colorbar()
    plt.title('Dissimilarity Martix')
    ax.set_xticks(np.arange(len(seq_names)), labels=seq_names) ; ax.tick_params(axis='x', labelsize=5)
    ax.set_yticks(np.arange(len(seq_names)), labels=seq_names) ; ax.tick_params(axis='y', labelsize=5)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")  
    plt.savefig(path_to_output, dpi = 500)

def reorder_dissimilarity_matrix(triang_Mat, old_seq_names, new_seq_names):
    '''
    List[List[int]] * List[str] * List[str] -> List[List[int]]
    Assumtion: triang_Mat is a lower triangular half of original dissimilarity matrix (diagonal included),
    old_seq_names and new_seq_names contain the same elements, just in a different order
    Returns a lower triangular half of new dissimilarity matrix (diagonal included), where elements are
    reordered according to the new header (new_seq_names)
    '''
    # new_Mat : List[List[int]] ; a triangular matrix to be reordered
    new_Mat = deepcopy(triang_Mat)

    # No necessity to do this indexing stuff if you use .index() method
    # old_to_new_index : Dict{int:int} ; a dictionnary of correspondance of indices {old index : new index}
    old_to_new_index = {}
    # new_index : int
    for new_index in range(len(new_seq_names)):
        # old_index : int
        for old_index in range(len(old_seq_names)):
            if new_seq_names[new_index] == old_seq_names[old_index]:
                old_to_new_index[old_index] = new_index

    for i in range(len(triang_Mat)):
        for j in range(len(triang_Mat[i])):
            # Note: since it's a lower triangular matrix, column index can never exceed a row index!
            #       and they're interchangeble since this matrix is symmetrical
            # new_i : int ; new row index
            new_i = max(old_to_new_index[i],old_to_new_index[j])
            # new_j : int ; new column index
            new_j = min(old_to_new_index[i],old_to_new_index[j])

            # the more succinct way, by using built-in method .index()
            #new_i = max(new_seq_names.index(old_seq_names[i]),new_seq_names.index(old_seq_names[j]))
            #new_j = min(new_seq_names.index(old_seq_names[i]),new_seq_names.index(old_seq_names[j]))

            new_Mat[new_i][new_j] = triang_Mat[i][j]
    return new_Mat


if __name__ == '__main__':

    # Test_seqs : Dict[str:str] ; test sequences {id:sequence}
    Test_seqs = EMBL.read_fasta('Files_for_Project/SeqTest.fasta')
    # Score_dict : Dict{(str,str):float} ; {(aminoacid1, aminoacid2):score}
    Score_dict = read_score_matrix('Files_for_Project/blosum62.mat')
    # profs_MSA_score : float ; MSA score obtained by prof's module
    # MSA_seqs : Dict[str] ; aligned sequences {id:sequence} 
    profs_MSA_score, MSA_seqs = msa.star_align(Test_seqs, Score_dict)
    # my_MSA_score : float ; MSA score obtained by my module
    # scores_by_pos : List[float] ; list of scores by position
    my_MSA_score, scores_by_pos = dissimilarity_score(MSA_seqs,Score_dict)
    assert profs_MSA_score == my_MSA_score

    # DisMat : List[List[float]] ; dissimilarity matrix
    # seq_names : List[str] ; sequences' names (header of dissimilarity matrix)
    DisMat , seq_names = compute_dissimilarity_matrix(MSA_seqs)

    # new_seq_names : List[str] ; reordered sequences' names
    new_seq_names = seq_names[::-1] 
    # DisMat : List[List[float]] ; reordered dissimilarity matrix
    newDisMat = reorder_dissimilarity_matrix(DisMat, seq_names, new_seq_names)

    print('Test sequences:')
    print(*['>'+name+'\n'+ Test_seqs[name] for name in Test_seqs.keys()], sep='\n')
    print('Aligned:')
    print(*['>'+name+'\n'+ MSA_seqs[name] for name in MSA_seqs.keys()], sep='\n')
    print('Overall score:',my_MSA_score)
    print('Scores by position:\n', scores_by_pos)
    print('Dissimilarity matrix:')
    print(seq_names,*triangular_to_rectangular([[round(elt,2) for elt in row] for row in DisMat]), sep='\n')
    print('Reordered matrix:')
    print(new_seq_names,*triangular_to_rectangular([[round(elt,2) for elt in row] for row in newDisMat]), sep='\n')
    
    