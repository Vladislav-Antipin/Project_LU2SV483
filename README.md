# Bioinformatics Project LU2SV483

## Description

My first uni bioinformatics project with some basic sequence analysis and comparison tools. Some of them are used to compare hemin-binding proteins hbpC, hbpD, hbpE between 5 members of Bartonella spp. 

For this, corresponding genes were found in Bartonella's WGS sequence set from EMBL database, their coding sequences are localized in the genome and  translated into protein. Those protein sequences are aligned using MSA with star heuristic. Aligned protein sequences are used to build a neighbor joining tree with UPGMA algorithm.

## Usage

This repository already contains all outputs in Output_Files folder, but you can play with the code or change inputs in Files_for_Project folder. 

To create a virtual environment and install  there all required modules, use this command in the terminal:

Unix/MacOS:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Windows:
```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```
But note: star alignment of protein sequences takes some time so it was made once in step 4 and then commented out. Uncomment this line of code in main_step4.py if you want to analyze other protein sequences.
```python 
MSA.write_msa_as_fasta(path_to_prots,path_to_score_matrix, path_to_msa_result)
```
 
## Details

### [Step 1](main_step1.py)

For this step I've written some [basic functions for sequence analysis](Project_Library/basic_seq_analysis.py), here is the example output for the test nucleotide sequences:

```
The sequence to analyze:
GCGCATGCTTACATAGGCCTAACAAACGGCTTATGACTAG

The length of the sequence is 40 nucleotides
The GC content is  0.475

The reverse complement is:
CTAGTCATAAGCCGTTTGTTAGGCCTATGTAAGCATGCGC

ATG positions for the sequence: 4 32
STOP positions for the sequence: 13 19 33 37
ATG positions for the reverse complement: 26 34
STOP positions for the reverse complement: 1 7 19 29

Enter weights for a match, substitution and indel (separated by space)
2 -1 -2

Comparison of two sequences:
AGTTAGAG--TAGGGCCAGCCAGCATAGCA-----GGTT
AGTTCCAGTATAGGCCCA---AGCAAAGCAGTACCGGTT

Matches: 25  Substitutions: 4  Indels: 10
Edit distance: 14
Alignment score is 26.0 (for Wm=2.0 Ws=-1.0 Wi=-2.0)
Identity ratio: 0.64
```

### [Steps 2 and 3](main_step2-3.py)

*see* [*EMBL analysis module*](Project_Library/embl_analysis.py)

[EMBL database file](Files_for_Project/Bartonella.dat) *(European Molecular Biology Laboratory)* with  WGS *(Whole Genome Sequencing)* annotations was used to create a [summary table](Output_Files/Bartonella_embl_summary_usingLISTOFLISTS.txt) with essential information (EMBL ID, species, direction +/-, Uniprot ID and protein function) for each CDS *(CoDing Sequence)* from the file.

Note: I didn't get exactly the same result using a list of lists and a dictionary using Uniprot ID as a key. This means that Uniprot ID is not a unique identifier for CDS and we better use a list of lists or use another ID as a key.

The same EMBL file was employed for the selection of proteins and species of interest. [Criteria for this selection](Files_for_Project/criteria_for_gene_selection.txt) were formulated as regular expressions:
```
Bartonella (henselae|quintana|taylorii|vinsonii|washoeensis)
hbp(C|D|E)
OMP_b-brl
```
The short summary about [selected proteins](Output_Files/Bartonella_embl_summary_selective.txt):

| Species | Uniprot ID | Function |
|----------|----------|----------|
| Bartonella washoeensis  | J1JPX0   | OMP_b-brl domain |
| Bartonella washoeensis | J1JLY6  | OMP_b-brl domain |
| Bartonella washoeensis  | J0QJM5   | OMP_b-brl domain |
| Bartonella vinsonii | J0R567   | OMP_b-brl domain |
| Bartonella vinsonii | J0QWD4    | OMP_b-brl domain  |
| Bartonella vinsonii  | J1JRZ2  | OMP_b-brl domain  |
| Bartonella taylorii | J0RK79  | OMP_b-brl domain  |
| Bartonella taylorii | J0RIX3   | OMP_b-brl domain  |
| Bartonella taylorii | J0RAR5   | OMP_b-brl domain |
| Bartonella quintana   | Q8KP14   | hbpC |
| Bartonella quintana | Q8KP12   | hbpD |
| Bartonella quintana  | Q8KP11   | hbpE |
| Bartonella henselae   | X5M5M3   | hbpC |
| Bartonella henselae  | X5MEW9  | hbpD  |
| Bartonella henselae   | X5MHR9  | hbpE |

Then, thanks to found CDS boarders and directions, corresponding nucleotide sequences were localized in [the genome](Files_for_Project/Bartonella.fasta).

Obtained CDSs were translated into [protein sequences](Output_Files/selected_proteins.fasta) using a specific [genetic code](Files_for_Project/CodeGenetique.tab). 

In the scope of those steps I've also written a function [find_all_ORFs()](Project_Library/basic_seq_analysis.py) that finds all ORFs *(Open Reading Frames)* for a given sequence. 

### [Step 4](main_step4.py)

*see* [*MSA analysis*](Project_Library/msa_analysis.py) *and* [*UPGMA*](Project_Library/upgma.py) *modules*

First, obtained protein sequences were [aligned](Output_Files/aligned_selected_proteins.fasta) using [MSA](Project_Library/msa.py) *(Multiple Sequence Alignment)* with star heuristic.

With the help of [BLOSUM62](Files_for_Project/blosum62.mat) substitution matrix, the overall MSA score was computed (54255) as well as the score for each position of the alignment, which I've represented as a barcode:

![](Output_Files/Alignment_score_plot_by_position.png)

Based on aligned sequences, a dissimilarity matrix was computed and represented as a heatmap (sequences in a random order):

![](Output_Files/Heatmap_for_unordered_sequences.png)

The dissimilarity matrix was then put into the [UPGMA algorithm](Project_Library/upgma.py) *(Unweighted Pair Group Method with Arithmetic Mean)*. The resulting tree was at first stocked as a dictionary {parent : (weight, (child_i_name, dist_i), (child_j_name, dist_j))}.

This dictionary was unpacked using a stack data structure as an intermediate, and then converted to newick parenthesized format as well as to a list of ordered according to it sequences.

The neighbor joining (NJ) tree was drawn based on the newik representation with [Biopython](requirements.txt) module. 

![](Output_Files/Tree.png)

The dissimilarity matrix was rearranged so that the order of sequences matches the one in the tree:

![](Output_Files/Heatmap_for_ordered_sequences.png)

### General Conclusion

The comparison of hemin-binding proteins hbpC, hbpD, hbpE between 5 members of Bartonella spp. allows us to suggest that those proteins have diverged before the the speciation of *B. henselae, B. quintana, B. taylorii, B. vinsonii,* and *B. washoeensis*. 


