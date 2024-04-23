import sys
sys.path.append('../Project_Library')
import msa   # msa.pyc file was compiled with python3.10, so use this version to run the code
import basic_seq_analysis as SeqAnalysis
import embl_analysis as EMBL
import msa_analysis as MSA
import upgma as UPGMA

from matplotlib import pyplot as plt

from Bio import Phylo
from io import StringIO

path_to_prots = '../Output_Files/selected_proteins.fasta'
path_to_msa_result = '../Output_Files/aligned_selected_proteins.fasta'
path_to_score_matrix = '../Files_for_Project/blosum62.mat'
path_to_alignment_score_plot = '../Output_Files/Alignement_score_plot_by_position.png'
path_to_unordered_heatmap = '../Output_Files/Heatmap_for_unordered_sequences.png'
path_to_tree = '../Output_Files/Tree.png'
path_to_ordered_heatmap = '../Output_Files/Heatmap_for_ordered_sequences.png'

# Use once to save the alignment result, otherwise you're gonna wait for eternity...
#MSA.write_msa_as_fasta(path_to_prots,path_to_score_matrix, path_to_msa_result)

MSA_dict = EMBL.read_fasta(path_to_msa_result)
Score_dict = MSA.read_score_matrix(path_to_score_matrix)

overall_score, scores_by_pos = MSA.dissimilarity_score(MSA_dict,Score_dict)
print(f'BLOSUM62 Alignment Score:{overall_score}')   # I DONT HAVE THE SAME SCORE as the one obtained by msa.py


MSA.plot_alignment_score(scores_by_pos, 'BLOSUM62 Alignment Score', path_to_alignment_score_plot)

DistMat, seq_names = MSA.compute_dissimilarity_matrix(MSA_dict)

MSA.plot_dissimilarity_matrix_as_heatmap(DistMat, seq_names, path_to_unordered_heatmap)

tree, root_name = UPGMA.get_tree_upgma(DistMat, seq_names)

ordered_seq_names = UPGMA.get_ordered_names_from_tree(tree, root_name)
newik_format = UPGMA.convert_tree_from_dict_to_newik(tree, root_name)

OrderedDistMat = MSA.reorder_dissimilarity_matrix(DistMat, seq_names, ordered_seq_names)
MSA.plot_dissimilarity_matrix_as_heatmap(OrderedDistMat, ordered_seq_names, path_to_ordered_heatmap)

# Build a tree with BioPython Module
handle = StringIO(newik_format)
tree_to_build = Phylo.read(handle, "newick")
Phylo.draw(tree_to_build, do_show=False)
plt.savefig(path_to_tree)

###### DEBUGGING #####

if False:
    DEFT_SEQ = {
            'SeqA': 'MSAGAEKQDRWIALHG',
            'SeqB': 'MSNGAEQERWIALVDHG',
            'SeqC': 'MSNAKGERLALVCHG',
            'SeqD': 'MSAVHKQDRLALVDIG' }

    ProtSeqs = EMBL.read_fasta(path_to_prots)

    MSA_score, MSA_seqs = msa.star_align(ProtSeqs, Score_dict)
    my_score = MSA.dissimilarity_score(MSA_seqs ,Score_dict)[0]

    print(f'Good score:{MSA_score}')
    print(f'My score:{my_score}')

    print(MSA_seqs)
    print('=====================================')
    print(MSA_dict)
