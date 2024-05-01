# Project_LU2SV483

## Description
My first uni bioinformatics project with some basic sequence treatment and comparison tools. Some of them are used to compare hemin-binding proteins hbpC, hbpD, hbpE between 5 members of Bartonella spp. 

## Usage
This repository already contains all outputs in Output_Files folder, but you can play with the code or change inputs in Files_for_Project folder. 

But note: star alignment of protein sequences takes some time so it was made once in step 4 and then commented out. Uncomment MSA.write_msa_as_fasta(path_to_prots,path_to_score_matrix, path_to_msa_result) in main_step4.py if you want to analyze other protein sequences.