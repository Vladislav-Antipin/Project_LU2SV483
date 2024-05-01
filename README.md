# Project_LU2SV483

## Description
My first uni bioinformatics project with some basic sequence treatment and comparison tools. Some of them are used to compare hemin-binding proteins hbpC, hbpD, hbpE between 5 members of Bartonella spp. 

For this, corresponding genes were found in Bartonella's WGS sequence set from EMBL database, their coding sequences are localized in the genome and  translated into protein. Those protein sequences are aligned using MSA with star heuristic. Aligned protein sequences are used to build a tree with UPGMA algorithm.

## Usage
This repository already contains all outputs in Output_Files folder, but you can play with the code or change inputs in Files_for_Project folder. 

But note: star alignment of protein sequences takes some time so it was made once in step 4 and then commented out. Uncomment MSA.write_msa_as_fasta(path_to_prots,path_to_score_matrix, path_to_msa_result) in main_step4.py if you want to analyze other protein sequences.

## General Conclusion
The comparison of hemin-binding proteins hbpC, hbpD, hbpE between 5 members of Bartonella spp. allows us to suggest that those proteins have diverged before the the speciation of *B. henselae, B. quintana, B. taylorii, B. vinsonii,* and *B. washoeensis*. 