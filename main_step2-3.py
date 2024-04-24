import sys
from Project_Library import embl_analysis as EMBL
from Project_Library import basic_seq_analysis as SeqAnalysis

##########          STEP 2           ##########
# path_to_embl : str ; path to EMBL database file
path_to_embl = '../Files_for_Project/Bartonella.dat'
# path_to_summary_output_from_dict : str
path_to_summary_output_from_dict = '../Output_Files/Bartonella_embl_summary_usingDICTIONARY.txt'
# path_to_summary_output_from_list_of_lists : str
path_to_summary_output_from_listoflists = '../Output_Files/Bartonella_embl_summary_usingLISTOFLISTS.txt'
# path_to_regex : str ; path to the file with regex - critirea to select CDSs
path_to_regex = '../Files_for_Project/criteria_for_gene_selection.txt'
# path_to_fasta : str ; path to the Fasta file with genomic DNA
path_to_fasta = '../Files_for_Project/Bartonella.fasta'

# embl_sum : dict[str:List[str]] ; EMBL file summary {UniProtID : [Sequence ID, Species, Direction, Begin, End, Function]}
embl_sum_dict = EMBL.get_EMBL_summary_as_dict(path_to_embl)
EMBL.write_EMBL_summary_file_ftom_dict(embl_sum_dict, path_to_summary_output_from_dict )

# embl_sum : dict[str:List[str]] ; EMBL file summary {UniProtID : [Sequence ID, Species, Direction, Begin, End, Function]}
embl_sum_listoflists = EMBL.get_EMBL_summary_as_listoflists(path_to_embl)
EMBL.write_EMBL_summary_file_ftom_listoflists(embl_sum_listoflists, path_to_summary_output_from_listoflists )

# Moral of the story : we didn't get exactly the same results using a list of lists
# and a dictionnary using UniProtID as a key - which means that
# UniProtID is not a unique identifier for CDS and we better use a list of lists

# CDS_positions : dict[str:List[str]] ; dictionary for each selected CDS {UniProtID : [Sequence ID, Direction, Begin, End]}
CDS_positions = EMBL.get_and_select_CDS_positions(path_to_embl, path_to_regex)

# genome : dict[str:str] ; a dictionary of genomic DNA {sequence ID : sequence}
genome = EMBL.read_fasta(path_to_fasta)

# CDSs : dict[str:str] ; a dictionary {UniProtID : CDS}
CDSs = EMBL.extract_CDSs_seq_from_fasta(CDS_positions,genome)


##########          STEP 3           ##########
# path_to_gencode : str ; path to txt file with the genetic code
path_to_gencode = '../Files_for_Project/CodeGenetique.tab'
# gencode : dict[str:str] ; dictionary for genetic code {Codon:Aminoacid}
gencode = {}

# streamr : stream ; input file stream
with open(path_to_gencode) as streamr:
    # line : str ; a read line of the file
    for line in streamr:
        if len(line.split())==2:
            # codon, aa : str ; codon and aminoacid read from genetic code file
            codon, aa = line.split()
            gencode[codon]=aa
# ProtSeqs : dict[str:str] ; dictionary {UniProtID : Protein sequence}
ProtSeqs = {}
# uniprot : str ; UniProtID, a key of CDSs dictionary
for uniprot in CDSs :
    ProtSeqs[uniprot]=SeqAnalysis.translate(CDSs[uniprot],gencode)

# uniprot : str ; UniProtID, a key of ProtSeqs dictionary
for uniprot in ProtSeqs:
    print(f'>{uniprot}\n{ProtSeqs[uniprot]}')

with open('../Output_Files/selected_proteins.fasta', 'w') as stream_w:
    for uniprot in ProtSeqs:
        stream_w.write(f'>{uniprot}\n{ProtSeqs[uniprot][:-1]}\n')
