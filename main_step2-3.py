from Project_Library import embl_analysis as EMBL
from Project_Library import basic_seq_analysis as SeqAnalysis

'''          STEP 2           '''

# path_to_embl : str ; path to EMBL database file
path_to_embl = './Files_for_Project/Bartonella.dat'
# path_to_summary_output_from_dict : str
path_to_summary_output_from_dict = './Output_Files/Bartonella_embl_summary_usingDICTIONARY.txt'
# path_to_summary_output_from_list_of_lists : str
path_to_summary_output_from_listoflists = './Output_Files/Bartonella_embl_summary_usingLISTOFLISTS.txt'
# path_to_regex : str ; path to the file with regex - criteria to select CDSs
path_to_regex = './Files_for_Project/criteria_for_gene_selection.txt'
# path_to_fasta : str ; path to the Fasta file with genomic DNA
path_to_fasta = './Files_for_Project/Bartonella.fasta'

'''Question 1'''

# embl_sum : dict[str:List[str]] ; EMBL file summary {UniProtID : [Sequence ID, Species, Direction, Begin, End, Function]}
embl_sum_dict = EMBL.get_EMBL_summary_as_dict(path_to_embl)
EMBL.write_EMBL_summary_file_from_dict(embl_sum_dict, path_to_summary_output_from_dict )

# embl_sum : dict[str:List[str]] ; EMBL file summary {UniProtID : [Sequence ID, Species, Direction, Begin, End, Function]}
embl_sum_listoflists = EMBL.get_EMBL_summary_as_listoflists(path_to_embl)
EMBL.write_EMBL_summary_file_from_listoflists(embl_sum_listoflists, path_to_summary_output_from_listoflists )

'''Moral of the story : we didn't get exactly the same results using a list of lists
and a dictionary using UniProtID as a key - which means that
UniProtID is not a unique identifier for CDS and we better use a list of lists'''

print('EMBL summary saved to', path_to_summary_output_from_listoflists )

'''Question 2'''

# CDS_positions : Dict[str:List[str]] ; dictionary for each selected CDS {UniProtID : [Sequence ID, Direction, Begin, End]}
CDS_positions = EMBL.get_and_select_CDS_positions(path_to_embl, path_to_regex)

'''Question 3'''

# genome : Dict[str:str] ; a dictionary of genomic DNA {sequence ID : sequence}
genome = EMBL.read_fasta(path_to_fasta)

'''Question 4'''

# CDSs : Dict[str:str] ; a dictionary {UniProtID : CDS} for extracted CDSs
CDSs = EMBL.extract_CDSs_seq_from_fasta(CDS_positions,genome)


'''          STEP 3           '''

'''Question 1 (modified step 2, used dictionaries)'''

''' Question 1 (see SeqAnalysis.find_all_ORFs() function and it's asserts)'''

'''Question 3'''

# path_to_gencode : str ; path to file with the genetic code
path_to_gencode = './Files_for_Project/CodeGenetique.tab'
# path_to_protseqs : str ; path to save protein sequences
path_to_protseqs = './Output_Files/selected_proteins.fasta'

# gencode : Dict[str:str] ; dictionary for genetic code {Codon:Aminoacid}
gencode = SeqAnalysis.read_genetic_code(path_to_gencode)

# ProtSeqs : Dict[str:str] ; dictionary {UniProtID : Protein sequence}
ProtSeqs = {}
# uniprot : str ; UniProtID, a key of CDSs dictionary
for uniprot in CDSs :
    ProtSeqs[uniprot]=SeqAnalysis.translate(CDSs[uniprot],gencode)

with open(path_to_protseqs, 'w') as stream_w:
    # uniprot : str ; UniProtID, a key of ProtSeqs dictionary
    for uniprot in ProtSeqs:
        stream_w.write(f'>{uniprot}\n{ProtSeqs[uniprot][:-1]}\n')

print('Nucleotide sequences, selected with criteria from', path_to_regex,'\nare translated and saved to', path_to_protseqs)