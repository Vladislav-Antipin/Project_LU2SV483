from Project_Library import msa   # msa.pyc file was compiled with python3.10, so use this version to run the code
from Project_Library import basic_seq_analysis as SeqAnalysis
from Project_Library import embl_analysis as EMBL
from Project_Library import msa_analysis as MSA
from Project_Library import upgma as UPGMA

from matplotlib import pyplot as plt

from Bio import Phylo
from io import StringIO


# path_to_score_matrix : str ; path to score (aka substitution) matrix
path_to_score_matrix = './Files_for_Project/blosum62.mat'
# path_to_embl_summary : str ; path to an embl summary (generated in step 2)
path_to_embl_summary = "./Output_Files/Bartonella_embl_summary_usingLISTOFLISTS.txt"

# path_to_prots : str ; path to fasta file with protein sequences to be analyzed
path_to_prots = './Output_Files/selected_proteins.fasta'
# path_to_msa_result : str ; path to msa result (to be saved or used directly)
path_to_msa_result = './Output_Files/aligned_selected_proteins.fasta'

# path_to_alignment_score_plot : str ; path to a plot of MSA score by position to be saved
path_to_alignment_score_plot = './Output_Files/Alignement_score_plot_by_position.png'
# path_to_unordered_heatmap : str ; path to a heatmap of unordered sequences to be saved
path_to_unordered_heatmap = './Output_Files/Heatmap_for_unordered_sequences.png'
# path_to_ordered_heatmap : str ; path to a heatmap of ordered sequences to be saved
path_to_ordered_heatmap = './Output_Files/Heatmap_for_ordered_sequences.png'
# path_to_tree : str ; path to a tree to be saved
path_to_tree = './Output_Files/Tree.png'

# Use once to save the alignment result, otherwise you're gonna wait for eternity...
#MSA.write_msa_as_fasta(path_to_prots,path_to_score_matrix, path_to_msa_result)

# MSA_dict : Dict[str:str] ;  dictionnary of MSA result {sequence ID : aligned sequence}
MSA_dict = EMBL.read_fasta(path_to_msa_result)
# Score_dict : Dict{(str,str):float} ; dictionary {(aminoacid1, aminoacid2):score}
Score_dict = MSA.read_score_matrix(path_to_score_matrix)

# overall_score : float ; overall MSA score
# scores_by_pos : List[float] ; list of scores by position
overall_score, scores_by_pos = MSA.dissimilarity_score(MSA_dict,Score_dict)
print(f'BLOSUM62 Alignment Score:{overall_score}')   


MSA.plot_alignment_score(scores_by_pos, 'BLOSUM62 Alignment Score', path_to_alignment_score_plot)

# DistMat : List[List[float]] ; lower triangular half of a dissimilarity matrix (diagonal included)
# seq_names : List[str] ; list of sequence names (header of the matrix)
DistMat, seq_names = MSA.compute_dissimilarity_matrix(MSA_dict)

# TODO: if I include species, they name is too long to be well represented on a tree/heatmap
# Modifies seq_names so that they contain species

# stream_r : stream ; input file stream 
with open(path_to_embl_summary, 'r') as stream_r:
        # line : str ; a read line
        for line in stream_r:
                # i : int ; index of a sequence name
                for i in range(len(seq_names)):
                        if seq_names[i] in line:
                                seq_names[i] += ' B. ' + line.split()[2]

MSA.plot_dissimilarity_matrix_as_heatmap(DistMat, seq_names, path_to_unordered_heatmap)

# tree : Dict{ str : Tuple(int, Tuple(str, float), Tuple(str, float) ) } ; tree as a 
#        dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))}
# root_name : str ; name of the root node
tree, root_name = UPGMA.get_tree_upgma(DistMat, seq_names)

# ordered_seq_names : List[str] ; list of sequence names ordered by their position in the paranthesized form of tree
# newik_format : str ; tree represented in a newik (paranthesized) format
ordered_seq_names = UPGMA.get_ordered_names_from_tree(tree, root_name)
newik_format = UPGMA.convert_tree_from_dict_to_newik(tree, root_name)

# OrderedDistMat : List[List[int]]; a lower triangular half of new dissimilarity matrix (diagonal included),
#                  where elements are reordered according to the new header
OrderedDistMat = MSA.reorder_dissimilarity_matrix(DistMat, seq_names, ordered_seq_names)
MSA.plot_dissimilarity_matrix_as_heatmap(OrderedDistMat, ordered_seq_names, path_to_ordered_heatmap)

# Draw a tree with BioPython Module
# handle : StringIO object
handle = StringIO(newik_format)
# tree_to_build : Bio.Phylo.Newick.Tree object
tree_to_build = Phylo.read(handle, "newick")
Phylo.draw(tree_to_build, do_show=False)
plt.savefig(path_to_tree)

###### DEBUGGING #####

if True:
    DEFT_SEQ = {
            'SeqA': 'MSAGAEKQDRWIALHG',
            'SeqB': 'MSNGAEQERWIALVDHG',
            'SeqC': 'MSNAKGERLALVCHG',
            'SeqD': 'MSAVHKQDRLALVDIG' }

    MSA_score, MSA_seqs = msa.star_align(DEFT_SEQ, Score_dict)
    my_score_list = MSA.dissimilarity_score(MSA_seqs ,Score_dict)[1]

    print(my_score_list)
    
